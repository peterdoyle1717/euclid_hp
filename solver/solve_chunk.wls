(* solve_chunk.wls SOLVER_DIR CHUNK_TSV PREC OBJ_DIR

   Chunk worker: one Mathematica kernel, many cases.

   CHUNK_TSV rows:  v <TAB> clers <TAB> netcode
   For each row, run lmSolve[netcode, prec, 20000] from lm_solve.wl.
   On SOLVER_TOL: compute per-vertex link turning from the realized
   coordinates, write OBJ_DIR/<clers>.obj (50-digit CForm), and emit

     clers v SOLVED iters resid minTurn minTurnAt objpath wall

   On any other status:

     clers v SOLVER_FAIL_<status> iters resid - - - wall

   One TSV row per case on stdout; diagnostics go to stderr. *)

argv = $ScriptCommandLine;
If[Length[argv] < 5,
  Print["usage: solve_chunk.wls SOLVER_DIR CHUNK_TSV PREC OBJ_DIR"];
  Exit[2]];
solverDir = argv[[2]];
chunkPath = argv[[3]];
prec      = ToExpression[argv[[4]]];
objDir    = argv[[5]];

SetOptions[$Output, PageWidth -> Infinity];
Get[FileNameJoin[{solverDir, "lm_solve.wl"}]];

fmtNum[x_] := If[x == 0, "0", ToString[CForm[N[x, 6]]]];

(* Per-vertex link turning from realized coordinates V: at each vertex,
   sum the spherical-link geodesic turning via the ArcTan[den, num]
   formula over the cyclically ordered unit directions to neighbors.

   Deliberately kept as a separate WL implementation (verbatim from the
   May 2026 sweep driver, claude/sandbox/prove_one.wls) rather than a
   port of the prover's cycle code: this diagnostic runs at solver
   precision and resolves flopper floors near 10^-(prec/2), far below
   anything the double-precision prover can represent. It gates
   nothing — verdicts come from the prover alone, and the prover's
   report files carry the certified sin(T/2) turning value. Sign
   anchor: the tetrahedron (CCAE) gives minTurn = +5.7319 here, and
   the prover ACCEPTs the same OBJ. *)
linkTurnings[V_List, faces_List] := Module[
  {nv = Length[V], em = Association[], nbrs, cyclicNbrs, turning},
  Do[Module[{a = f[[1]], b = f[[2]], c = f[[3]]},
     em[{a, b}] = c; em[{b, c}] = a; em[{c, a}] = b], {f, faces}];
  nbrs = Table[{}, {nv}];
  Do[Module[{a = f[[1]], b = f[[2]], c = f[[3]]},
     nbrs[[a]] = Union[nbrs[[a]], {b, c}];
     nbrs[[b]] = Union[nbrs[[b]], {a, c}];
     nbrs[[c]] = Union[nbrs[[c]], {a, b}]], {f, faces}];
  cyclicNbrs[v_] := Module[{ring = {nbrs[[v]][[1]]}, cur, nxt},
    cur = ring[[1]];
    While[True,
      nxt = Lookup[em, Key[{v, cur}], 0];
      If[nxt === ring[[1]] || nxt === 0, Break[]];
      AppendTo[ring, nxt]; cur = nxt];
    ring];
  turning[v_] := Module[{ring, dirs, k, A, B, C, num, den, t = 0},
    ring = cyclicNbrs[v]; k = Length[ring];
    If[k < 3, Return[Indeterminate]];
    dirs = Table[
      Module[{d = V[[ring[[i]]]] - V[[v]]}, d/Sqrt[d . d]], {i, k}];
    Do[
      A = dirs[[Mod[i - 2, k] + 1]];
      B = dirs[[i]];
      C = dirs[[Mod[i, k] + 1]];
      num = B . Cross[A, C];
      den = (A . B) (B . C) - (A . C);
      t += ArcTan[den, num],
      {i, k}];
    t];
  Table[turning[v], {v, nv}]
];

rows = Select[StringSplit[ReadString[chunkPath], "\n"], # =!= "" &];

Do[
  Module[{fields, v, clers, netcode, t0, res, wall, V, faces,
          turnings, minTurn, minTurnAt, objpath},
    fields = StringSplit[row, "\t"];
    If[Length[fields] < 3,
      WriteString["stderr", "bad row: " <> row <> "\n"]; Continue[]];
    {v, clers, netcode} = fields[[1 ;; 3]];
    t0 = AbsoluteTime[];
    res = lmSolve[netcode, prec, 20000];
    wall = AbsoluteTime[] - t0;
    If[res["status"] =!= "SOLVER_TOL",
      (* PARSE_FAIL/TOPOLOGY_FAIL/WISH_FAIL rows lack lmIters/lmFinalResid *)
      Print[StringRiffle[{clers, v,
        "SOLVER_FAIL_" <> ToString[res["status"]],
        ToString[Lookup[res, "lmIters", "-"]],
        Module[{rr = Lookup[res, "lmFinalResid", "-"]},
          If[rr === "-", "-", fmtNum[rr]]],
        "-", "-", "-", ToString[N[wall, 4]]}, "\t"]];
      Continue[]];
    V = res["V"]; faces = res["topology"]["faces"];
    turnings = linkTurnings[V, faces];
    minTurn = Min[turnings];
    minTurnAt = Position[turnings, minTurn][[1, 1]];
    objpath = FileNameJoin[{objDir, clers <> ".obj"}];
    writeObjHP[objpath, V, faces, 50];
    Print[StringRiffle[{clers, v, "SOLVED",
      ToString[res["lmIters"]], fmtNum[res["lmFinalResid"]],
      fmtNum[minTurn], ToString[minTurnAt], objpath,
      ToString[N[wall, 4]]}, "\t"]]],
  {row, rows}];
