# `lm_solve.wl` — Specification

A Mathematica package that finds a Euclidean polyhedral realization of a
triangulated combinatorial net by Levenberg–Marquardt iteration on edge
**bends**, with a dent gate to keep iterates in the undented cone.
**Intended** port of `euclid_clean.c` from the neoplatonic-proof-pipeline.
The file's header comment says so, but the C source was not diffed against
the Mathematica file when this spec was written, so an implementer should
treat the C file as ground truth on any discrepancy.

The driver `lmSolve[netcode, prec, maxIter]` runs the whole pipeline at an
arbitrary fixed numerical precision.

## 1. Inputs, outputs, top-level driver

### Entry point
```
lmSolve[netcode_String, prec_Integer, maxIter_Integer:200] -> Association
```

- `netcode` — `"a1,b1,c1;a2,b2,c2;..."`, semicolon-separated 1-based
  triangle vertex triples. The first triangle is the **base face**.
- `prec` — required numeric precision in decimal digits. There is no
  machine-precision entry point.
- `maxIter` — outer LM iteration cap. Default 200.

### Returned Association keys
- `status` — one of `"PARSE_FAIL"`, `"TOPOLOGY_FAIL"`, `"WISH_FAIL"`,
  `"REALIZE_FAIL"`, `"SOLVER_TOL"`,
  `"LAMBDA_SATURATED_POSITIVE_RESIDUAL"`, `"STALLED_POSITIVE_RESIDUAL"`.
- `failureReason` — coarse tag: `"-"`, `"LM_LAMBDA"`, `"LM_STALLED"`,
  `"LM_RETRIES"`, `"LM_MAX_ITER"`, `"LM_OTHER"`.
- `NV`, `NE`, `prec`, `alphaTargetDeg=60`.
- `lmIters`, `lmFinalResid`, `lmFinalLambda`, `lmMsg`.
- `bend` — final per-edge bend list, length NE, at the requested precision.
- `bendInit` — wish-start bend list.
- `topology` — the Association from `buildTopology`.
- `V` — NV-by-3 list of vertex coordinates from realize, if available.

### Precision discipline (REQUIRED)
The entire driver runs inside

```
Block[{$MinPrecision = prec + 30}, ...]
```

with a cushion of 30 digits. **Do not set `$MaxPrecision`.** Inside the
block, `bend0 = wishStart[top, prec+30]`, `alpha = N[Pi/3, prec+30]`, and
`tol = 10^(-prec)`.

## 2. Constants

```
$lmTargetDeg     = 60         (* target dihedral cone angle, degrees *)
$lmFinalTol      = 1/10^12    (* unused in current code; tol = 10^(-prec) *)
$lmLambdaInit    = 1
$lmMaxIter       = 200        (* default outer cap *)
$lmLambdaDown    = 3/10       (* exact rationals — propagate precision *)
$lmLambdaUp      = 10
$lmLambdaMax     = 10^12
$lmTinyRel       = 1/10^12
$lmMaxLmRetries  = 20
```

All as exact rationals (or symbolic `Pi`) so that high-precision arithmetic
is preserved.

## 3. Pipeline summary

```
netcode --parseNetcode-->  parsed
parsed  --buildTopology--> top
top     --setupAllBends--> layout
top     --wishStart-----> bend0          (HP, length NE)
                           +
                           alpha = N[Pi/3, prec+30]
                           tol   = 10^(-prec)
       --solveLMInner---> sol            (final bend, residual, status)
sol.bend --realize-----> V               (NV by 3 vertex positions)
```

## 4. Netcode parsing — `parseNetcode[s]`

Strip, split on `";"`, drop empty tokens, split each on `","` and
`ToExpression`. Validate every face is a length-3 list of positive
integers; return `$Failed` otherwise. Output:

```
<|"faces" -> List[List[Int,Int,Int]],
  "NV"    -> Max @@ Flatten[faces],
  "NF"    -> Length[faces]|>
```

## 5. Topology — `buildTopology[parsed]`

Builds canonical edges and a per-vertex "flower" (cyclically ordered ring
of incident faces and the corresponding "other" vertex).

### Edges
- Canonicalized as `{min,max}`. `edgeIdx[{lo,hi}] = e` records the
  assignment in **first-encountered face order**: scan faces in order,
  for each face register its three undirected edges if not seen yet,
  assigning the next integer index. `edges` is `Transpose[{ea, eb}]`,
  a length-NE list of `{lo,hi}`.
- `dirFace[{u,v}] = fi` for **each oriented edge** of face `fi`, i.e.
  `{a,b}`, `{b,c}`, `{c,a}`. If any oriented edge is reused, or any face
  is degenerate (`a==b` etc.), abort with `$Failed`. This enforces a
  2-manifold with consistent orientation: each undirected edge appears in
  exactly two faces, in opposite orientations.

### Vertex degrees
`vertDeg[[v]]` = count of edges incident to `v`, from `ea`/`eb`.

### Flowers (per vertex `v`)
The flower of `v` is the cyclic sequence of faces incident to `v`, walked
via dual-graph BFS through the opposite-edge link.

Algorithm:
1. Find `start` = first face containing `v`. If none, fail.
2. Set `cur = start`, `seenCount = 0`, collect `flT` (third vertex per
   step), `flE` (edge-index per step).
3. Loop: read `{a,b,c} = faces[cur]`. Set `third = c` if `v==a`, `a` if
   `v==b`, else `b` (i.e. third = the vertex that follows `v` in the
   face's cyclic order `a->b->c->a`). Append `third` to `flT`; append
   `edgeIdx[{min(v,third),max(v,third)}]` to `flE`. Increment
   `seenCount`. Lookup `nxt = dirFace[{v, third}]` (the face across edge
   `(v,third)` on the side reached by orientation `v->third`). If
   missing, fail. If `nxt == start`, break. Else `cur = nxt`.
4. Sanity: `seenCount` must equal `vertDeg[[v]]`; otherwise fail.
5. Guard: if `seenCount > NV + 6`, abort.

Returns
`<|"NV","NF","NE","faces","edges","edgeIdx","dirFace","vertDeg","flowerLen","flowerE","flowerThird"|>`.

## 6. Layout — `setupAllBends[top]`

Trivial **all-bends** layout: every edge is a free variable, every vertex
contributes a residual.

```
<|"base"   -> top["faces"][[1]],
  "varOfE" -> Range[NE],   (* variable column for each edge *)
  "resVs"  -> Range[NV],   (* residual rows are 3 per vertex *)
  "gateVs" -> Range[NV]|>  (* dent gate checks every vertex *)
```

`varOfE[e]` returns the column index of edge `e`'s bend in the LM step
vector; here it's identity, but the indirection lets future layouts pin
edges by setting `varOfE[e] = 0`.

## 7. Wish init — `wishStartExact`, `wishStart`

Closed-form (exact rational) seed for the bend vector. Idea: per-edge
"inv-avg" target `g_e = 1/d_i + 1/d_j` (where `d_v` = vertex degree, `i,j`
are the edge endpoints), KKT-projected onto the linear subspace `B x = c`
where `B` is the (NV × NE) vertex-edge incidence (all 1's, undirected),
`c = 1` (per-vertex sum-to-1 in revolutions).

KKT solve (exact arithmetic over rationals):

```
BBt    = B . Transpose[B]
rhs    = B . g - 2 c
lambda = LinearSolve[BBt, rhs]
x      = (g - Transpose[B] . lambda) / 2
```

Verify `Max[Abs[B . x - c]] == 0` exactly; fail otherwise. Then
`wishStart[top, prec] = N[2 Pi x, prec]` — convert revolutions to radians
at the requested precision. Output length NE.

## 8. Quaternions

Quaternions are 4-element lists `{w, x, y, z}`.

```
qMul[a,b] = {a1 b1 - a2 b2 - a3 b3 - a4 b4,
             a1 b2 + a2 b1 + a3 b4 - a4 b3,
             a1 b3 - a2 b4 + a3 b1 + a4 b2,
             a1 b4 + a2 b3 - a3 b2 + a4 b1}

qStep[alpha,beta] = { cos(alpha/2) cos(beta/2),
                     -cos(alpha/2) sin(beta/2),
                     -sin(alpha/2) sin(beta/2),
                      sin(alpha/2) cos(beta/2)}

qStepDBeta[alpha,beta] = d/dbeta qStep
                       = {-cos(alpha/2) sin(beta/2)/2,
                          -cos(alpha/2) cos(beta/2)/2,
                          -sin(alpha/2) cos(beta/2)/2,
                          -sin(alpha/2) sin(beta/2)/2}
```

Note the explicit `/2` (exact rational) inside `qStepDBeta` so that
precision propagates.

## 9. Holonomy residual — `holonomyResidualQuat`

For each residual vertex `v` (in `layout["resVs"]`), compute the product
of `qStep[alpha, bend[e]]` over the flower edges of `v` in flower order.
The residual contribution is the **vector part** (components 2,3,4) of the
product quaternion — three real numbers per vertex. Total residual length
`M = 3 * Length[resVs]`.

```
Q_v = Product_{s=1..flowerLen[v]}  qStep[alpha, bend[flowerE[v,s]]]
r_v = (Q_v[2], Q_v[3], Q_v[4])
```

The scalar part (`Q_v[1]`) is discarded; the constraint is that the
rotation around vertex `v` is the identity (modulo sign). At
`alpha = pi/3` the realized cone angle around each vertex is the full
turn `2 pi`.

## 10. Dense analytic Jacobian — `quatAnalyticJacobianDense`

J is `(3 NR) x NVAR` where NR = `Length[resVs]`, NVAR = `Max[varOfE]`.

For each residual vertex `v` with flower length `k`:

1. Compute per-step quaternions
   `qs[t] = qStep[alpha, bend[flowerE[v,t]]]`, `t=1..k`.
2. Compute per-step derivatives
   `dqs[t] = qStepDBeta[alpha, bend[flowerE[v,t]]]`.
3. Build prefix and suffix products with the identity quaternion
   `id = {1,0,0,0}` at the requested precision:
   - `P[1] = id`, `P[t+1] = qMul[P[t], qs[t]]`.
   - `S[k+1] = id`, `S[t] = qMul[qs[t], S[t+1]]` (built right-to-left).
4. For each flower step `t`: edge `e = flowerE[v,t]`, column
   `col = varOfE[e]`. If `col <= 0`, skip (variable is pinned).
   - `dQ = qMul[qMul[P[t], dqs[t]], S[t+1]]`.
   - Add `dQ[2..4]` to rows `3(v_index-1)+1 .. 3(v_index-1)+3`, column
     `col`.

J is initialized as a dense matrix of zeros at the input precision.

## 11. Dent gate

```
vertexTurn[v, bend, top]       = Sum over e in flowerE[v] of bend[e]
hasDentFirst[bend, top, layout] = first v in gateVs with vertexTurn[v]<0,
                                  else 0
```

A vertex with negative total bend is "dented." The LM step is rejected
whenever the trial bend would dent any gate vertex.

## 12. LM inner — `solveLMInner`

Levenberg–Marquardt with **Marquardt scaling** (D = diag(JᵀJ)),
retry-on-reject, dent gate, tiny-step stall detector, and no D-floor.

State: `bend`, residual `r`, residual norm `norm`, damping `lam`, outer
iter `it`, success flag.

Initialization:
```
r    = holonomyResidualQuat[alpha, bend, top, layout]
norm = Sqrt[r.r]
lam  = lambdaInit            (* default 1 *)
it   = 0
```

Outer loop while `it < maxIter`:

1. If `norm <= tol`, set `success = True`, `msg = "tol"`, break.
2. Compute `J`, `A = Jᵀ J`, `g = Jᵀ r`, `D = Diagonal[A]`.
3. **No D-floor.** If `Min[D] < $lmTinyRel * Max[D]`, `Abort[]` with a
   diagnostic print. (The C version floors; this port aborts loudly
   instead so missing guards are visible.)
4. Inner retry loop, up to `$lmMaxLmRetries` (20) times:
   - Solve `(A + lam * DiagonalMatrix[D]) delta = -g`. (Marquardt
     scaling: damping is added relative to each variable's diagonal.)
     `Quiet[Check[LinearSolve[...], $Failed]]`. If failed: `lam *= 10`.
     If `lam > 1e12`, set `msg = "lambda_saturated_linsolve"`, go to done.
     Else continue.
   - Apply step: `bendTrial[e] = bend[e] + delta[varOfE[e]]` for each
     edge with `varOfE[e] > 0`.
   - Compute `rTrial`, `nTrial = Sqrt[rTrial.rTrial]`,
     `dentV = hasDentFirst[bendTrial,...]`.
   - **Accept condition:** `nTrial < norm` AND `dentV == 0`. If accepted:
     - Tiny-step stall test: `dnorm = Sqrt[delta.delta]`,
       `xnorm = Sqrt[sum of bend[e]^2 over active edges]`. If
       `dnorm < (1+xnorm)/1e15` AND `nTrial >= norm`, commit the trial
       state, set `msg = "stalled"`, go to done.
     - Else commit: `bend = bendTrial`, `r = rTrial`, `norm = nTrial`,
       `lam *= 3/10`, mark accepted, break inner.
   - **Reject:** `lam *= 10`. If `lam > 1e12`, set
     `msg = "lambda_saturated"`, go to done. Continue inner.
5. If inner loop exits without accepting, set
   `msg = "lm_retries_exhausted"`, go to done.
6. `it += 1`.

After outer loop: if exited via `it >= maxIter`, set
`success = (norm <= tol)`, `msg = "max_iter_tol"` or `"max_iter"`.

Return
`<|"success", "iters", "finalResid"=norm, "finalLambda"=lam, "msg", "bend"|>`.

Note: lambda decrease on accept is **not** clipped to a minimum in this
port (the C comment mentions `>=1e-30` clipping, but the code does not
enforce it).

## 13. Realize — `reconstruct`, `realize`

Standard-gauge BFS placement of NV vertices from face 1.

### Base gauge (unit edges)
Let `{b0,b1,b2} = faces[1]`. Set:
```
V[b0] = (0, 0,  1/2)
V[b1] = (0, 0, -1/2)
V[b2] = (Sqrt[3]/2, 0, 0)
```
This pins translation, rotation, and reflection. The base triangle is
equilateral with edge length 1.

### BFS
- `placedF[1] = True`; queue starts at `{1}`.
- Dequeue `fi`, scan its three oriented edges `(a,b)` in face order. The
  face across that edge is `otherFi = dirFace[{b,a}]` (the reverse
  orientation). If missing, skip (boundary; not expected for a closed
  mesh). If already placed, skip.
- In `otherFi`, find the third vertex `c` (not `a`, not `b`).
- In `fi`, find `p` (the third vertex of `fi`, not `a`, not `b`) — the
  "previous" vertex used to orient the fold direction.
- `e = edgeIdx[{min(a,b),max(a,b)}]`; `theta = bend[e]`.
- `V[c] = placeThird[V, a, b, p, theta]`.
- Mark `c` placed, mark `otherFi` placed, enqueue `otherFi`.

### `placeThird[V, a, b, p, theta]`
Place `c` on the unit-edge equilateral triangle on edge `(a,b)`, folded
by signed dihedral angle `theta` away from `p`'s plane.

```
m     = (V[a] + V[b]) / 2
eHat  = V[b] - V[a]                 (* length 1 in normal use; not renormalized *)
pPerp = V[p] - m                    (* in-plane vector perpendicular to ab on p's side *)
uHat  = pPerp / Norm[pPerp]
vHat  = Cross[uHat, eHat]
vHat  = vHat / Norm[vHat]
cPerp = -Cos[theta] uHat + Sin[theta] vHat
V[c]  = m + (Sqrt[3]/2) cPerp
```

Intuition: rotate `-uHat` (the in-plane unit perpendicular pointing away
from `p`) by `theta` toward `vHat = uHat x eHat`. `theta = 0` plants `c`
flat across edge `ab` from `p` (zero bend, coplanar). Positive `theta`
folds `c` toward the `vHat` side. The factor `Sqrt[3]/2` is the apothem
of an equilateral triangle of side 1.

Note `eHat` is not renormalized; the code relies on edges being unit
length to numerical precision in the realized state. (Edge lengths can be
verified post hoc with `edgeLengths`.)

If any vertex remains unplaced, return
`<|"error" -> "unplaced_vertices", "V", "vPlaced"|>`; otherwise return
the `V` list.

`realize[bend, top]` wraps `reconstruct` and tags failure as
`<|"error" -> "reconstruct_failed"|>` or success as `<|"V" -> V|>`.

## 14. Driver assembly — `lmSolve`

```
Block[{$MinPrecision = prec + 30},
  parsed = parseNetcode[netcode]            -- "PARSE_FAIL" on null
  top    = buildTopology[parsed]            -- "TOPOLOGY_FAIL"
  layout = setupAllBends[top]
  bend0  = wishStart[top, prec + 30]        -- "WISH_FAIL"
  alpha  = N[Pi/3, prec + 30]
  tol    = 10^(-prec)
  sol    = solveLMInner[alpha, bend0, tol, maxIter, top, layout, 1]
  rz     = realize[sol.bend, top]           -- "REALIZE_FAIL" if no V
  status = "SOLVER_TOL" iff sol.success, else lambda-saturation or
           "STALLED_POSITIVE_RESIDUAL"
  return assoc of all of the above plus V from rz
]
```

## 15. `edgeLengths`, `writeObjHP`

Utilities:

```
edgeLengths[V, top] = Table[Norm[V[edges[e,2]] - V[edges[e,1]]], e=1..NE]

writeObjHP[path, V, faces, digits:50]:
  For each vertex v: write "v {x} {y} {z}\n"
                     via ToString[CForm[N[..., digits]]]
  For each face f:   write "f {f[1]} {f[2]} {f[3]}\n"
```

`CForm` is used so the OBJ is parseable by the C downstream tools (e.g.
`euclid_prover`).

## 16. Invariants and conventions an implementer must preserve

- **No machine-precision entry point.** Every numeric quantity threaded
  through the solver carries `prec + 30` digits or more. The
  `Block[{$MinPrecision = prec + 30}, ...]` enforces a precision floor
  for the whole driver. Do not introduce a `$MaxPrecision` setting; let
  Mathematica's adaptive precision raise where needed for things like
  `Sqrt[3]/4`.
- **`alpha = N[Pi/3, prec+30]`, `tol = 10^(-prec)`.** Tolerance is at the
  requested precision, not half.
- **Rational constants.** All LM tuning constants (lambda factors,
  tinyrel, finaltol-if-used) are exact rationals.
- **Triangle orientation is an input invariant.** The topology builder
  uses `dirFace` keyed by oriented edges and assumes the netcode already
  has consistent winding. Do not validate or "fix" orientation by
  computing signed volume or flipping faces.
- **Dent gate uses bends, not coordinates.** A vertex is dented iff its
  bend sum is negative. The C realize step never re-checks dentedness
  from geometry.
- **No D-floor in the solver.** If
  `min(diag(JᵀJ)) < 1e-12 * max(diag)`, abort with diagnostic. Add the
  floor back only when an observed case demands it.
- **Base face pinning** (`V[b0]`, `V[b1]`, `V[b2]` as in §13) is the
  canonical gauge; all reported residuals/edge lengths/dihedrals are in
  this frame.
- **All-bends layout.** Every edge bend is a free variable, every vertex
  contributes a 3-vector residual, every vertex is dent-gated.
  `varOfE[e] = 0` would pin edge `e`, but that hook is unused in this
  entry point.

## 17. Behavioral sanity check (unverified expectations)

The solver was not run while this spec was written — the numbers below
are *expectations* an implementer can use as a smell test, not
measurements. Verify against the C reference (or against `lmSolve` on
the original `lm_solve.wl`) before trusting them.

```
res = lmSolve["1,2,3;1,3,4;1,4,2;2,4,3", 200]   (* CCAE / tetrahedron, prec=200 *)
```

Expected qualitatively:
- `res["status"]` == `"SOLVER_TOL"`.
- `res["NV"] == 4`, `res["NE"] == 6`.
- `res["lmFinalResid"]` below `tol = 10^(-200)`.
- `Max[Abs[edgeLengths[res["V"], res["topology"]] - 1]]` small relative
  to the realized scale; the exact magnitude depends on how much
  precision the realize step retains and was *not* measured here.
- `lmIters` small — the tetrahedron is the easy case, with no flat
  vertices.

The iteration-count regime (linear vs near-quadratic) for nontrivial
cases is an open question in this session; do not assume one from this
spec. For a nontrivial sanity case, pick a `v=26` net and confirm
`status == "SOLVER_TOL"` with `lmFinalResid < tol` and
`min link turn > 0` — but those values, too, should be measured rather
than predicted.

Other claims in this spec (data-structure shapes, the dent-gate rule,
the gauge formulas, the constants, `tol = 10^(-prec)`, the
`$MinPrecision = prec + 30` block with no `$MaxPrecision`) are read
directly from the source file and stand without separate measurement.
