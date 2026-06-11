#!/bin/sh
# run_all.sh INPUT_TSV OUT_DIR [PREC]
#
# High-precision solve + certified prove for a list of hard cases.
#
#   INPUT_TSV   rows: v <TAB> clers   (e.g. data/molasses_rejects_v4_50.tsv)
#   OUT_DIR     created; everything lands here
#   PREC        solver precision in decimal digits (default 1000)
#
# Env knobs:
#   JOBS    parallel chunk workers (default 4)
#   CHUNKS  chunk count (default 2*JOBS)
#   NICE    nice level; 0 disables (default 0)
#
# Stages:
#   1. clers decode             CLERS -> netcode (bulk, order-preserving)
#   2. solve_chunk.wls          lmSolve at PREC; OBJ per solved case
#   3. euclid_prover --batch    certified ACCEPT/REJECT on the OBJs as given
#   4. verdicts.tsv             one row per case, CLERS-keyed

set -eu

if [ $# -lt 2 ]; then
    echo "usage: $0 INPUT_TSV OUT_DIR [PREC]" >&2
    exit 2
fi

INPUT=$1
OUT=$2
PREC=${3:-1000}
JOBS=${JOBS:-4}
CHUNKS=${CHUNKS:-$(( JOBS * 2 ))}
NICE=${NICE:-0}

[ -f "$INPUT" ] || { echo "no such input: $INPUT" >&2; exit 2; }

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CLERS_BIN="$ROOT/submodules/clers/bin/clers"
PROVER="$ROOT/prover/euclid_prover"

# Mathematica front end: WS env override, else first of
# wolframscript / wolfram / math on PATH.
if [ -n "${WS:-}" ]; then
    : # caller-supplied, e.g. WS="wolframscript -file"
elif command -v wolframscript >/dev/null 2>&1; then
    WS="wolframscript -file"
elif command -v wolfram >/dev/null 2>&1; then
    WS="wolfram -script"
elif command -v math >/dev/null 2>&1; then
    WS="math -script"
else
    echo "no Mathematica front end (wolframscript/wolfram/math) on PATH" >&2
    exit 1
fi

NICE_PREFIX=""
[ "$NICE" != 0 ] && NICE_PREFIX="nice -n $NICE"

make -C "$ROOT/prover" >&2
make -C "$ROOT/submodules/clers" >&2

mkdir -p "$OUT/chunks" "$OUT/objs" "$OUT/solvelogs" \
         "$OUT/solvemanifest" "$OUT/reports"

{
    echo "host=$(hostname)"
    echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "input=$INPUT"
    echo "cases=$(grep -c . "$INPUT" | tr -d ' ')"
    echo "prec=$PREC"
    echo "jobs=$JOBS chunks=$CHUNKS nice=$NICE"
    echo "frontend=$WS"
    echo "repo_commit=$(git -C "$ROOT" rev-parse HEAD 2>/dev/null || echo unknown)"
} > "$OUT/manifest.txt"

# Stage 1: decode. Bulk, order-preserving; abort on any count mismatch.
cut -f2 "$INPUT" | "$CLERS_BIN" decode > "$OUT/netcodes.txt"
n_in=$(grep -c . "$INPUT")
n_out=$(grep -c . "$OUT/netcodes.txt")
if [ "$n_in" -ne "$n_out" ]; then
    echo "decode count mismatch: $n_in cases, $n_out netcodes" >&2
    exit 1
fi
paste "$INPUT" "$OUT/netcodes.txt" > "$OUT/cases.tsv"

# Stage 2: chunk (awk round-robin) and solve.
awk -v n="$CHUNKS" -v root="$OUT/chunks" '
    NF { print > sprintf("%s/chunk_%04d", root, (NR - 1) % n) }
' "$OUT/cases.tsv"

if command -v parallel >/dev/null 2>&1; then
    find "$OUT/chunks" -type f -name 'chunk_*' | sort | \
    parallel -j "$JOBS" --will-cite "
        $NICE_PREFIX $WS '$ROOT/solver/solve_chunk.wls' '$ROOT/solver' {} $PREC '$OUT/objs' \
            > '$OUT/solvemanifest/'{/}.tsv \
            2> '$OUT/solvelogs/'{/}.err
    " || echo "warning: parallel reported failed solver chunk(s)" >&2
else
    for chunk in "$OUT/chunks"/chunk_*; do
        [ -f "$chunk" ] || continue
        base=$(basename "$chunk")
        $NICE_PREFIX $WS "$ROOT/solver/solve_chunk.wls" "$ROOT/solver" \
            "$chunk" "$PREC" "$OUT/objs" \
            > "$OUT/solvemanifest/$base.tsv" \
            2> "$OUT/solvelogs/$base.err" \
            || echo "warning: solver chunk $base failed" >&2
    done
fi

cat "$OUT/solvemanifest"/chunk_*.tsv > "$OUT/solve_results.tsv"

# Completeness gate: every input case must have produced exactly one
# solver result row; otherwise abort rather than publish partial verdicts.
n_rows=$(grep -c . "$OUT/solve_results.tsv" || true)
if [ "$n_rows" -ne "$n_in" ]; then
    echo "solver result rows ($n_rows) != input cases ($n_in); aborting" >&2
    exit 1
fi

# Stage 3: prove the OBJs exactly as written.
awk -F '\t' '$3 == "SOLVED" { print $8 }' "$OUT/solve_results.tsv" \
    > "$OUT/objfiles.txt"
if [ -s "$OUT/objfiles.txt" ]; then
    $NICE_PREFIX "$PROVER" --batch --outdir "$OUT/reports" \
        < "$OUT/objfiles.txt" > "$OUT/prover_raw.tsv" || true
else
    : > "$OUT/prover_raw.tsv"
fi

# Stage 4: join, CLERS-keyed. Prover rows reference the OBJ path whose
# basename is <clers>.obj.
awk -F '\t' '
    FNR == NR {
        n = split($3, parts, "/");
        base = parts[n]; sub(/\.obj$/, "", base);
        verdict[base] = $2; message[base] = $5;
        next
    }
    {
        clers = $1; stage = $3;
        pv = (clers in verdict) ? verdict[clers] : "-";
        pm = (clers in message) ? message[clers] : "";
        print $2 "\t" clers "\t" stage "\t" pv "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" pm
    }
' "$OUT/prover_raw.tsv" "$OUT/solve_results.tsv" > "$OUT/verdicts.tsv"

echo "== summary"
echo "cases:        $n_in"
echo "solved:       $(awk -F'\t' '$3 == "SOLVED"' "$OUT/solve_results.tsv" | grep -c . || true)"
echo "solver fail:  $(awk -F'\t' '$3 != "SOLVED"' "$OUT/solve_results.tsv" | grep -c . || true)"
echo "prover accept: $(awk -F'\t' '$4 == "accept"' "$OUT/verdicts.tsv" | grep -c . || true)"
echo "prover reject: $(awk -F'\t' '$4 == "reject"' "$OUT/verdicts.tsv" | grep -c . || true)"
echo "prover fail:   $(awk -F'\t' '$4 == "fail"' "$OUT/verdicts.tsv" | grep -c . || true)"
echo "verdicts: $OUT/verdicts.tsv"
