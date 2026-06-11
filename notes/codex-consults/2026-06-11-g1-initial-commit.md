# G1 design audit — euclid_hp initial commit

- Date: 2026-06-11
- Codex session: `019eb8f4-c3b7-7de2-99e9-10b791af6cad` (codex exec, gpt-5.5)
- Bundle: the full request below; codex additionally inspected the
  working tree itself (verified md5s, row counts, shell syntax, builds).

## Request (sent via stdin)

What this repo is: high-precision rescue stage — clers decode →
`lmSolve[netcode, prec]` (all-bends LM, wish start, dent gate, 60°
direct, arbitrary precision) → 50-digit OBJ → vendored `euclid_prover`
(interval arithmetic, never refines input) → CLERS-keyed verdicts.
Default input: the 532 molasses REJECTs at v=4..50.

Provenance of every staged convention (nothing re-derived here):
`prover/euclid_prover.c` byte-identical to
neoplatonic-proof-pipeline@81b034d (md5 ef04a05dd6fd968499d02b9a6dece0e8),
validated 2026-06-11 by fresh-clone reproduce runs on doob: v=4..20
compare=PASS (14/14), v=4..30 compare=PASS (68/68) vs the committed
molasses reference. `solver/lm_solve.wl`/`.md` byte-identical to
claude/sandbox/mma_play (md5 df4df2d75b32389f32282a7ebcb89fe8),
validated May 2026 (bit-equivalent behavior to prior production solver
at prec=200/1000 on v=8, v=26 cases; recorded in claude/MIGRATION.md).
New code: scripts/run_all.sh (orchestration) and solver/solve_chunk.wls
(chunk driver; turning formula lifted from claude/sandbox/prove_one.wls).

Questions: (1) is the vendoring provenance sufficient to treat staged
conventions as inherited-and-validated? (2) defects in the new code
paths? (3) anything that must change before commit?

## Codex findings (verbatim)

1. **BLOCKER** scripts/run_all.sh:102: serial fallback masks solver
   worker failure with `|| true`, then concatenates whatever manifests
   exist. A crashed/aborted chunk can silently drop cases and still
   publish `verdicts.tsv` with fewer rows than input. Add a post-solve
   count check against `n_in`, and preferably remove `|| true` or
   record chunk failure explicitly.
2. **Minor robustness** solver/solve_chunk.wls:80: non-`SOLVER_TOL`
   paths assume `lmIters` and `lmFinalResid` exist. `PARSE_FAIL`,
   `TOPOLOGY_FAIL`, `WISH_FAIL` do not guarantee those keys. Official
   decoded inputs likely avoid this, but failure reporting should use
   `Lookup`/defaults.
3. Vendoring provenance is sufficient for this initial commit as
   stated. I verified md5s: `prover/euclid_prover.c` =
   `ef04a05dd6fd968499d02b9a6dece0e8`; `solver/lm_solve.wl` =
   `df4df2d75b32389f32282a7ebcb89fe8`.

Other checks: input has 532 rows, no duplicate CLERS keys, shell syntax
passes, prover and clers build. Basename-keyed join is acceptable for
this alphabet and unique input.

Verdict line: `BLOCKED: serial fallback can silently publish incomplete
results after solver worker failure.`

## Resolution (this commit)

- Finding 1: fixed. Both parallel and serial paths now report failed
  chunks to stderr, and a completeness gate aborts before the prover
  stage unless solver result rows == input cases
  (scripts/run_all.sh, "Completeness gate" block).
- Finding 2: fixed. Failure rows use `Lookup[res, ..., "-"]` and only
  format the residual when present (solver/solve_chunk.wls).
- Re-validated: smoke (prec=200, 2 cases) re-run after the fixes;
  identical verdicts (CCAE accept; v12 reject via sigma certificate),
  completeness gate 2/2.

## Round 2 (commit-gate review, codex session 019eb8f8-9e96-7642-9ead-4d6521ad4932)

Finding: `solve_chunk.wls` adds its own link-turning/ring traversal even
though the staged prover carries the ccw-outside cycle and turning
convention; minTurn lands in verdicts.tsv, so convention drift could
make diagnostics misleading. (Codex's own numerical cross-check printed
identical turning totals for the two traversals on its test case.)

Disposition (advisory, accepted-in-part): keep the WL helper, document
it. It is verbatim from the validated May sweep driver (prove_one.wls);
it is diagnostic-only — verdicts come exclusively from the prover, whose
per-case report files carry the certified sin(T/2) turning; and it
provides the one quantity the double prover cannot: minTurn at the
~10^-(prec/2) flopper floor, far below double range. Sign anchored
empirically: CCAE gives minTurn=+5.7319 and the prover ACCEPTs the same
OBJ. Documentation added at the helper definition in solve_chunk.wls.
