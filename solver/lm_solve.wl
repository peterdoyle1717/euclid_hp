(* lm_solve.wl
   Regenerated from lm_solve.md.  See that file for the spec; this is
   meant to be reproducible from it.

   Driver:  lmSolve[netcode, prec, maxIter:200]
*)

(* ------------------------------------------------------------------
 * Constants (exact rationals so HP propagates).
 * ------------------------------------------------------------------ *)

$lmTargetDeg    = 60;
$lmFinalTol     = 1/10^12;
$lmLambdaInit   = 1;
$lmMaxIter      = 200;
$lmLambdaDown   = 3/10;
$lmLambdaUp     = 10;
$lmLambdaMax    = 10^12;
$lmTinyRel      = 1/10^12;
$lmMaxLmRetries = 20;

(* ------------------------------------------------------------------
 * parseNetcode  "a,b,c;d,e,f;..."  ->  <|faces, NV, NF|>
 * ------------------------------------------------------------------ *)

parseNetcode[s_String] := Module[{toks, faces},
  toks  = Select[StringSplit[StringTrim[s], ";"], # =!= "" &];
  faces = (ToExpression /@ StringSplit[#, ","]) & /@ toks;
  If[!AllTrue[faces,
      Length[#] === 3 && AllTrue[#, IntegerQ[#] && # >= 1 &] &],
     Return[$Failed]];
  <|"faces" -> faces,
    "NV"    -> Max @@ Flatten[faces],
    "NF"    -> Length[faces]|>
];

(* ------------------------------------------------------------------
 * buildTopology  edges canonical (min,max); dirFace on oriented
 * edges; per-vertex flower walked via dual graph.
 * ------------------------------------------------------------------ *)

buildTopology[parsed_Association] := Module[
  {faces = parsed["faces"], NV = parsed["NV"], NF, NE,
   edgeIdx, dirFace, ne, ea, eb, addEdge, edges, vertDeg,
   flowerLen, flowerE, flowerThird,
   fi, a, b, c, i, v, start, cur, third, nxt, flE, flT, seenCount},

  NF = Length[faces];
  edgeIdx = <||>; dirFace = <||>;
  ne = 0; ea = {}; eb = {};
  addEdge[u_, w_] := Module[{lo = Min[u, w], hi = Max[u, w]},
    If[!KeyExistsQ[edgeIdx, {lo, hi}],
       ne += 1; AppendTo[ea, lo]; AppendTo[eb, hi];
       edgeIdx[{lo, hi}] = ne]];

  Do[
    {a, b, c} = faces[[fi]];
    If[a == b || b == c || a == c, Return[$Failed]];
    If[KeyExistsQ[dirFace, {a, b}], Return[$Failed]];
    If[KeyExistsQ[dirFace, {b, c}], Return[$Failed]];
    If[KeyExistsQ[dirFace, {c, a}], Return[$Failed]];
    dirFace[{a, b}] = fi; dirFace[{b, c}] = fi; dirFace[{c, a}] = fi;
    addEdge[a, b]; addEdge[b, c]; addEdge[c, a],
    {fi, NF}];

  NE      = ne;
  edges   = Transpose[{ea, eb}];
  vertDeg = ConstantArray[0, NV];
  Do[vertDeg[[ea[[i]]]] += 1; vertDeg[[eb[[i]]]] += 1, {i, NE}];

  flowerLen   = ConstantArray[0, NV];
  flowerE     = Table[{}, {NV}];
  flowerThird = Table[{}, {NV}];

  Do[
    start = -1;
    Do[If[MemberQ[faces[[fi]], v], start = fi; Break[]], {fi, NF}];
    If[start < 0, Return[$Failed]];

    cur = start; seenCount = 0; flE = {}; flT = {};
    While[True,
      If[seenCount > NV + 6, Return[$Failed]];
      {a, b, c} = faces[[cur]];
      third = Which[v == a, c, v == b, a, True, b];
      AppendTo[flT, third];
      AppendTo[flE, edgeIdx[{Min[v, third], Max[v, third]}]];
      seenCount += 1;
      nxt = Lookup[dirFace, Key[{v, third}], Missing[]];
      If[MissingQ[nxt], Return[$Failed]];
      If[nxt === start, Break[]];
      cur = nxt];

    flowerLen[[v]]   = seenCount;
    flowerE[[v]]     = flE;
    flowerThird[[v]] = flT;
    If[seenCount =!= vertDeg[[v]], Return[$Failed]],
    {v, NV}];

  <|"NV" -> NV, "NF" -> NF, "NE" -> NE,
    "faces" -> faces, "edges" -> edges,
    "edgeIdx" -> edgeIdx, "dirFace" -> dirFace,
    "vertDeg" -> vertDeg,
    "flowerLen" -> flowerLen,
    "flowerE" -> flowerE,
    "flowerThird" -> flowerThird|>
];

(* ------------------------------------------------------------------
 * Trivial all-bends layout.
 * ------------------------------------------------------------------ *)

setupAllBends[top_Association] := <|
  "base"   -> top["faces"][[1]],
  "varOfE" -> Range[top["NE"]],
  "resVs"  -> Range[top["NV"]],
  "gateVs" -> Range[top["NV"]]|>;

(* ------------------------------------------------------------------
 * Wish init.  Exact KKT projection of inv-avg seed onto B x = 1
 * (B = vertex-edge incidence, all 1's).
 * ------------------------------------------------------------------ *)

wishStartExact[top_Association] := Module[
  {NV = top["NV"], NE = top["NE"], edges = top["edges"],
   vertDeg = top["vertDeg"], B, g, c, BBt, rhs, sol, x, e},
  B = ConstantArray[0, {NV, NE}];
  Do[B[[edges[[e, 1]], e]] = 1; B[[edges[[e, 2]], e]] = 1, {e, NE}];
  g   = Table[1/vertDeg[[edges[[e, 1]]]]
            + 1/vertDeg[[edges[[e, 2]]]], {e, NE}];
  c   = ConstantArray[1, NV];
  BBt = B . Transpose[B];
  rhs = B . g - 2 c;
  sol = Quiet[Check[LinearSolve[BBt, rhs], $Failed]];
  If[sol === $Failed, Return[$Failed]];
  x = (g - Transpose[B] . sol)/2;
  If[Max[Abs[B . x - c]] =!= 0, Return[$Failed]];
  x
];

wishStart[top_Association, prec_Integer] := Module[
  {xExact = wishStartExact[top]},
  If[xExact === $Failed, Return[$Failed]];
  N[2 Pi xExact, prec]
];

(* ------------------------------------------------------------------
 * Quaternions  {w, x, y, z}.
 * ------------------------------------------------------------------ *)

qMul[a_List, b_List] := {
  a[[1]] b[[1]] - a[[2]] b[[2]] - a[[3]] b[[3]] - a[[4]] b[[4]],
  a[[1]] b[[2]] + a[[2]] b[[1]] + a[[3]] b[[4]] - a[[4]] b[[3]],
  a[[1]] b[[3]] - a[[2]] b[[4]] + a[[3]] b[[1]] + a[[4]] b[[2]],
  a[[1]] b[[4]] + a[[2]] b[[3]] - a[[3]] b[[2]] + a[[4]] b[[1]]};

qStep[alpha_, beta_] := Module[{ca, sa, cb, sb},
  ca = Cos[alpha/2]; sa = Sin[alpha/2];
  cb = Cos[beta/2];  sb = Sin[beta/2];
  {ca cb, -ca sb, -sa sb, sa cb}];

qStepDBeta[alpha_, beta_] := Module[{ca, sa, cb, sb},
  ca = Cos[alpha/2]; sa = Sin[alpha/2];
  cb = Cos[beta/2];  sb = Sin[beta/2];
  {-ca sb/2, -ca cb/2, -sa cb/2, -sa sb/2}];

(* ------------------------------------------------------------------
 * Residual:  per vertex v, vector part of product of qStep over
 * flower edges of v.
 * ------------------------------------------------------------------ *)

holonomyResidualQuat[alpha_, bend_List, top_, layout_] := Module[
  {resVs = layout["resVs"], flowerE = top["flowerE"],
   flowerLen = top["flowerLen"], r = {}, t, v, kk, s, e, Q},
  Do[
    v = resVs[[t]]; kk = flowerLen[[v]];
    Q = qStep[alpha, bend[[flowerE[[v, 1]]]]];
    Do[
      e = flowerE[[v, s]];
      Q = qMul[Q, qStep[alpha, bend[[e]]]],
      {s, 2, kk}];
    AppendTo[r, Q[[2]]]; AppendTo[r, Q[[3]]]; AppendTo[r, Q[[4]]],
    {t, Length[resVs]}];
  r
];

(* ------------------------------------------------------------------
 * Dense analytic Jacobian.  Prefix/suffix products around each
 * vertex's flower.  Rows = 3*(residual vertex index); columns =
 * varOfE[e].  varOfE[e] = 0 pins edge e.
 * ------------------------------------------------------------------ *)

quatAnalyticJacobianDense[alpha_, bend_List, top_, layout_] := Module[
  {resVs = layout["resVs"], varOfE = layout["varOfE"],
   flowerE = top["flowerE"], flowerLen = top["flowerLen"],
   nvar, NR, J, id, zero, i, v, k, t, e, qs, dqs, P, S, col, dQ, row},
  id   = SetPrecision[{1, 0, 0, 0}, Precision[alpha]];
  zero = SetPrecision[0, Precision[alpha]];
  nvar = Max[varOfE];
  NR   = Length[resVs];
  J    = ConstantArray[zero, {3 NR, nvar}];
  Do[
    v = resVs[[i]]; k = flowerLen[[v]];
    qs  = Table[qStep[alpha,     bend[[flowerE[[v, t]]]]], {t, k}];
    dqs = Table[qStepDBeta[alpha, bend[[flowerE[[v, t]]]]], {t, k}];
    P   = ConstantArray[id, k + 1];
    Do[P[[t + 1]] = qMul[P[[t]], qs[[t]]], {t, k}];
    S   = ConstantArray[id, k + 1];
    Do[S[[t]] = qMul[qs[[t]], S[[t + 1]]], {t, k, 1, -1}];
    Do[
      e   = flowerE[[v, t]];
      col = varOfE[[e]];
      If[col <= 0, Continue[]];
      dQ  = qMul[qMul[P[[t]], dqs[[t]]], S[[t + 1]]];
      row = 3 (i - 1);
      J[[row + 1, col]] += dQ[[2]];
      J[[row + 2, col]] += dQ[[3]];
      J[[row + 3, col]] += dQ[[4]],
      {t, k}],
    {i, NR}];
  J
];

(* ------------------------------------------------------------------
 * Dent gate:  vertex_turn(v) = sum of bend on flowerE[v].
 * ------------------------------------------------------------------ *)

vertexTurn[v_Integer, bend_List, top_Association] :=
  Total[bend[[#]] & /@ top["flowerE"][[v]]];

hasDentFirst[bend_List, top_, layout_] := Module[
  {gv = layout["gateVs"], t, v},
  Do[
    v = gv[[t]];
    If[vertexTurn[v, bend, top] < 0, Return[v]],
    {t, Length[gv]}];
  0
];

(* ------------------------------------------------------------------
 * LM inner.  Marquardt scaling, retry-on-reject, dent gate, tiny-
 * step stall detector, no D-floor (loud Abort if D would collapse).
 * ------------------------------------------------------------------ *)

solveLMInner[alphaIn_, bendIn_List, tol_, maxIter_, top_, layout_,
             lambdaInitIn_] := Module[
  {alpha = alphaIn, bend = bendIn,
   varOfE = layout["varOfE"], NE = top["NE"],
   r, norm, lam, it, rt, J, A, g, D, maxDiag, floorV,
   bendTrial, rTrial, nTrial, dentV, delta, dnorm, xnorm,
   accepted, success, msg, iters, finalLambda, e},

  r    = holonomyResidualQuat[alpha, bend, top, layout];
  norm = Sqrt[r . r];
  lam  = lambdaInitIn;
  it   = 0;
  success = False; msg = "uninit";

  While[it < maxIter,
    If[norm <= tol, success = True; msg = "tol"; Break[]];

    J = quatAnalyticJacobianDense[alpha, bend, top, layout];
    A = Transpose[J] . J;
    g = Transpose[J] . r;
    D = Diagonal[A];
    maxDiag = Max[D];
    floorV  = $lmTinyRel * maxDiag;
    If[Min[D] < floorV,
      Print["ABORT solveLMInner: D-floor would have engaged."];
      Print["  iter = ", it, "  min diag(JᵀJ) = ", Min[D]];
      Print["  max diag(JᵀJ) = ", maxDiag, "  floor = ", floorV];
      Abort[]];

    accepted = False;
    Do[
      delta = Quiet[Check[
        LinearSolve[A + DiagonalMatrix[lam D], -g], $Failed]];
      If[delta === $Failed,
        lam *= $lmLambdaUp;
        If[lam > $lmLambdaMax,
           msg = "lambda_saturated_linsolve"; Goto[doneInner]];
        Continue[]];

      bendTrial = bend;
      Do[If[varOfE[[e]] > 0,
            bendTrial[[e]] += delta[[varOfE[[e]]]]], {e, NE}];
      rTrial = holonomyResidualQuat[alpha, bendTrial, top, layout];
      nTrial = Sqrt[rTrial . rTrial];
      dentV  = hasDentFirst[bendTrial, top, layout];

      If[nTrial < norm && dentV == 0,
        dnorm = Sqrt[delta . delta];
        xnorm = Sqrt[Total[
          (If[varOfE[[#]] > 0, bend[[#]], 0])^2 & /@ Range[NE]]];
        If[dnorm < (1 + xnorm)/10^15 && nTrial >= norm,
          bend = bendTrial; r = rTrial; norm = nTrial;
          msg = "stalled"; Goto[doneInner]];
        bend = bendTrial; r = rTrial; norm = nTrial;
        lam *= $lmLambdaDown;
        accepted = True; Break[]];

      lam *= $lmLambdaUp;
      If[lam > $lmLambdaMax,
         msg = "lambda_saturated"; Goto[doneInner]],
      {rt, $lmMaxLmRetries}];

    If[!accepted, msg = "lm_retries_exhausted"; Goto[doneInner]];
    it += 1;
  ];

  If[it >= maxIter,
    success = norm <= tol;
    msg = If[success, "max_iter_tol", "max_iter"]];

  Label[doneInner];
  iters = it; finalLambda = lam;
  <|"success" -> success, "iters" -> iters,
    "finalResid" -> norm, "finalLambda" -> finalLambda,
    "msg" -> msg, "bend" -> bend|>
];

(* ------------------------------------------------------------------
 * placeThird.  Unit-edge equilateral triangle on (a,b) folded by
 * signed dihedral theta away from p.
 *
 *   m     = (Va+Vb)/2
 *   eHat  = Vb-Va
 *   uHat  = (Vp-m)/Norm
 *   vHat  = Cross[uHat,eHat]/Norm
 *   Vc    = m + (Sqrt[3]/2) (-Cos[theta] uHat + Sin[theta] vHat)
 * ------------------------------------------------------------------ *)

placeThird[V_List, a_Integer, b_Integer, p_Integer, theta_] := Module[
  {Va = V[[a]], Vb = V[[b]], Vp = V[[p]],
   m, eHat, pPerp, uHat, vHat, cPerp},
  m     = (Va + Vb)/2;
  eHat  = Vb - Va;
  pPerp = Vp - m;
  uHat  = pPerp/Norm[pPerp];
  vHat  = Cross[uHat, eHat];
  vHat  = vHat/Norm[vHat];
  cPerp = -Cos[theta] uHat + Sin[theta] vHat;
  m + (Sqrt[3]/2) cPerp
];

(* ------------------------------------------------------------------
 * reconstruct:  BFS from face 1 in the standard gauge.
 *   V[b0] = (0,0,1/2), V[b1] = (0,0,-1/2), V[b2] = (Sqrt[3]/2,0,0).
 * ------------------------------------------------------------------ *)

reconstruct[bend_List, top_Association] := Module[
  {NV = top["NV"], NF = top["NF"], faces = top["faces"],
   dirFace = top["dirFace"], edgeIdx = top["edgeIdx"],
   V, vPlaced, placedF, queue, qh, b0, b1, b2,
   fi, fa, fb, fc, vs, i, a, b, otherFi, oa, ob, oc, c, p, e},

  V       = ConstantArray[{0, 0, 0}, NV];
  vPlaced = ConstantArray[False, NV];
  {b0, b1, b2} = faces[[1]];
  V[[b0]] = {0, 0,  1/2};
  V[[b1]] = {0, 0, -1/2};
  V[[b2]] = {Sqrt[3]/2, 0, 0};
  vPlaced[[b0]] = True;
  vPlaced[[b1]] = True;
  vPlaced[[b2]] = True;

  placedF = ConstantArray[False, NF];
  placedF[[1]] = True;
  queue = {1}; qh = 1;

  While[qh <= Length[queue],
    fi = queue[[qh]]; qh += 1;
    {fa, fb, fc} = faces[[fi]];
    vs = {fa, fb, fc};
    Do[
      a = vs[[i]]; b = vs[[Mod[i, 3] + 1]];
      otherFi = Lookup[dirFace, Key[{b, a}], Missing[]];
      If[MissingQ[otherFi], Continue[]];
      If[otherFi == fi || placedF[[otherFi]], Continue[]];
      {oa, ob, oc} = faces[[otherFi]];
      c = Which[oa =!= a && oa =!= b, oa,
                ob =!= a && ob =!= b, ob,
                True, oc];
      p = Which[fa =!= a && fa =!= b, fa,
                fb =!= a && fb =!= b, fb,
                True, fc];
      e = edgeIdx[{Min[a, b], Max[a, b]}];
      V[[c]] = placeThird[V, a, b, p, bend[[e]]];
      vPlaced[[c]] = True;
      placedF[[otherFi]] = True;
      AppendTo[queue, otherFi],
      {i, 3}]];

  If[!AllTrue[vPlaced, # &],
     <|"error" -> "unplaced_vertices",
       "V" -> V, "vPlaced" -> vPlaced|>,
     V]
];

realize[bend_List, top_Association] := Module[{V = reconstruct[bend, top]},
  If[AssociationQ[V],
     Return[<|"error" -> "reconstruct_failed", "details" -> V|>]];
  <|"V" -> V|>
];

edgeLengths[V_List, top_Association] := Module[{edges = top["edges"]},
  Table[Norm[V[[edges[[e, 2]]]] - V[[edges[[e, 1]]]]], {e, top["NE"]}]
];

(* ------------------------------------------------------------------
 * Driver.  All HP work inside Block[{$MinPrecision = prec + 30}, ...].
 * No $MaxPrecision setting.  tol = 10^(-prec).
 * ------------------------------------------------------------------ *)

lmSolve[netcode_String, prec_Integer, maxIter_Integer:200] := Block[
  {$MinPrecision = prec + 30},
  Module[
   {parsed, top, layout, bend0, alpha, tol, sol, rz, status, reason},

   parsed = parseNetcode[netcode];
   If[parsed === $Failed, Return[<|"status" -> "PARSE_FAIL"|>]];

   top = buildTopology[parsed];
   If[top === $Failed, Return[<|"status" -> "TOPOLOGY_FAIL"|>]];

   layout = setupAllBends[top];

   bend0 = wishStart[top, prec + 30];
   If[bend0 === $Failed,
     Return[<|"status" -> "WISH_FAIL", "topology" -> top|>]];

   alpha = N[Pi/3, prec + 30];
   tol   = 10^(-prec);

   sol = solveLMInner[alpha, bend0, tol, maxIter, top, layout,
                      $lmLambdaInit];
   rz  = realize[sol["bend"], top];

   status = Which[
     !KeyExistsQ[rz, "V"], "REALIZE_FAIL",
     sol["success"],       "SOLVER_TOL",
     MemberQ[{"lambda_saturated", "lambda_saturated_linsolve"},
             sol["msg"]],  "LAMBDA_SATURATED_POSITIVE_RESIDUAL",
     True,                 "STALLED_POSITIVE_RESIDUAL"];

   reason = Switch[sol["msg"],
     "tol",                       "-",
     "max_iter_tol",              "-",
     "lambda_saturated",          "LM_LAMBDA",
     "lambda_saturated_linsolve", "LM_LAMBDA",
     "stalled",                   "LM_STALLED",
     "lm_retries_exhausted",      "LM_RETRIES",
     "max_iter",                  "LM_MAX_ITER",
     _,                           "LM_OTHER"];

   <|"status" -> status, "failureReason" -> reason,
     "NV" -> top["NV"], "NE" -> top["NE"], "prec" -> prec,
     "alphaTargetDeg" -> $lmTargetDeg,
     "lmIters"      -> sol["iters"],
     "lmFinalResid" -> sol["finalResid"],
     "lmFinalLambda"-> sol["finalLambda"],
     "lmMsg"        -> sol["msg"],
     "bend"         -> sol["bend"],
     "bendInit"     -> bend0,
     "topology"     -> top|>
   ~Join~ (If[KeyExistsQ[rz, "V"], rz, <||>])
  ]
];

(* ------------------------------------------------------------------
 * writeObjHP — emit HP OBJ in C-parseable CForm.
 * ------------------------------------------------------------------ *)

writeObjHP[path_String, V_List, faces_List, digits_Integer:50] := Module[
  {fh, fmt, v, f},
  fh = OpenWrite[path];
  fmt[x_] := ToString[CForm[N[x, digits]]];
  Do[WriteString[fh,
       "v ", fmt[V[[v, 1]]], " ", fmt[V[[v, 2]]],
            " ", fmt[V[[v, 3]]], "\n"],
     {v, Length[V]}];
  Do[WriteString[fh,
       "f ", faces[[f, 1]], " ", faces[[f, 2]],
            " ", faces[[f, 3]], "\n"],
     {f, Length[faces]}];
  Close[fh];
  path
];
