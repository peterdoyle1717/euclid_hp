# euclid_hp

High-precision rescue stage for the neoplatonic existence proof.

The double-precision pipeline
([neoplatonic-proof-pipeline](https://github.com/peterdoyle1717/neoplatonic-proof-pipeline))
REJECTs 532 of the 8,239,684 prime 6-nets at v=4..50. Some of those are
true floppers — zero link turning, so no existence certificate can
exist. Others fail only because coordinates realized from a
double-precision solve are not accurate enough for the prover's
existence inequality. This repo rescues the second class: solve at high
precision in Mathematica, write coordinates that are correctly rounded
doubles, and let the unchanged interval prover decide.

The prover is the same single-file `euclid_prover.c` as the main
pipeline (vendored here). It takes the OBJ exactly as given and never
refines it; solving and proving stay separate.

## Stages

1. `clers decode` — CLERS → netcode (pinned submodule).
2. `solver/solve_chunk.wls` — `lmSolve[netcode, prec]` from
   `solver/lm_solve.wl`: all-bends Levenberg–Marquardt on edge bends,
   least-squares wish start, dent gate on every accepted step, target
   cone angle 60° directly, everything at `prec` decimal digits
   (spec in `solver/lm_solve.md`). Writes one OBJ per solved case,
   plus per-case link-turning diagnostics.
3. `prover/euclid_prover --batch` — certified interval-arithmetic
   ACCEPT/REJECT (existence, embedding, undentedness) on the OBJs.
4. `verdicts.tsv` — one CLERS-keyed row per case:
   `v clers solver_stage prover_verdict iters resid minTurn minTurnAt message`.

## Requirements

- Mathematica (`wolframscript`, `wolfram`, or `math` on PATH)
- C compiler; GNU parallel optional (serial fallback included)

## Run

```sh
git clone --recurse-submodules https://github.com/peterdoyle1717/euclid_hp.git
cd euclid_hp
./scripts/run_all.sh data/molasses_rejects_v4_50.tsv runs/rejects 1000
```

Env knobs: `JOBS` (parallel kernels, default 4), `CHUNKS` (default
2·JOBS), `NICE` (default 0). One Mathematica kernel runs per chunk, many
cases per kernel.

`data/molasses_rejects_v4_50.tsv` is the committed list of the 532
REJECT cases from the official v=4..50 run (`v <TAB> clers`).

## License

MIT.
