(* Stage B trace driver for the 1-loop bubble pipeline.
   Runs MMA's AMFlow on bubble (same invariants as the C++ test),
   intercepting AMFlow's Layer 11-16 key functions and dumping every
   call's input/output as JSONL to $TRACE_OUT.

   Usage: TRACE_OUT=...jsonl math -script trace_bubble.wl

   The wrappers use the Block-trick to safely recurse into the original
   definitions without disturbing DownValues long-term.  To avoid the
   `ToExpression` context-shadowing pitfall, we reference every symbol
   by its literal context-qualified name (AMFlow`X, not "AMFlow`X").  *)

rootdir = ExpandFileName[FileNameJoin[{DirectoryName[$InputFileName], "..", ".."}]];
amfdir = FileNameJoin[{rootdir, "reference", "amflow-master"}];
Get[FileNameJoin[{amfdir, "AMFlow.m"}]];

SetReductionOptions["IBPReducer" -> "Kira"];
$KiraExecutable   = "/usr/local/bin/kira";
$FermatExecutable = "/usr/share/Ferl7/fer64";

SetDefaultOptions[];
SetGlobalOptions["SilentMode" -> True];
SetExpansionOptions["XOrder" -> 100, "ExtraXOrder" -> 20];

AMFlowInfo["Family"]       = bubblefam;
AMFlowInfo["Loop"]         = {l};
AMFlowInfo["Leg"]          = {p};
AMFlowInfo["Conservation"] = {};
AMFlowInfo["Replacement"]  = {p^2 -> s};
AMFlowInfo["Propagator"]   = {l^2 - msq, (l - p)^2 - msq};
AMFlowInfo["Numeric"]      = {s -> 5, msq -> 1};
AMFlowInfo["NThread"]      = 4;

(* --- Trace infrastructure --- *)

$stepCounter = 0;
nextStep[] := (++$stepCounter);
toStr[e_]  := ToString[e // InputForm, CharacterEncoding -> "ASCII"];

$traceRecords = {};
addRec[r_] := AppendTo[$traceRecords, r];

(* Manual helper wrappers for functions where the generic Block-trick
   proved fragile.  We copy the AMFlow.m bodies verbatim in structure,
   then add trace records around the exact semantics we care about. *)

Begin["AMFlow`Private`"];

Unprotect[AnalyzeTopology];
AnalyzeTopology[denominators_] := Module[{s = Global`nextStep[], u, f, massterm, components, info},
  Global`addRec[<|"step" -> s, "fn" -> "AnalyzeTopology", "phase" -> "in",
           "denominators" -> Global`toStr[denominators]|>];
  {u, f, massterm} = EvaluateUF[denominators];
  components = If[Head[#] === Times, List @@ #, {#}] & @ Factor[u];
  info = AnalyzeComponent[#, f, massterm] & /@ components;
  info = Select[info, #[[3]] > 0 &];
  info = Sort[info, #1[[3]] >= #2[[3]] &];
  Global`addRec[<|"step" -> s, "fn" -> "AnalyzeTopology", "phase" -> "out",
           "u" -> Global`toStr[u], "f" -> Global`toStr[f], "massterm" -> Global`toStr[massterm],
           "components" -> Global`toStr[components], "result" -> Global`toStr[info]|>];
  info
];

Unprotect[AMFCandidateComponent];
AMFCandidateComponent[info0_, mode_] := Module[{s = Global`nextStep[], vac, single, result},
  Global`addRec[<|"step" -> s, "fn" -> "AMFCandidateComponent", "phase" -> "in",
           "info0" -> Global`toStr[info0], "mode" -> Global`toStr[mode]|>];
  vac = VacuumQ[info0];
  single = SingleMassQ[info0];
  result = GetFirst @ Which[
    !vac, AllPossiblePosition[info0, mode],
    !single, AllPossiblePosition[info0, "Branch"],
    True, {}
  ];
  Global`addRec[<|"step" -> s, "fn" -> "AMFCandidateComponent", "phase" -> "out",
           "vacuumQ" -> Global`toStr[vac], "singleMassQ" -> Global`toStr[single],
           "result" -> Global`toStr[result]|>];
  result
];

Unprotect[FindAllRegion];
FindAllRegion[topposi_] := Module[
  {s = Global`nextStep[], branch0, branch, cutbranch, cutposi, trans, patt, regions,
   pres, regionsByScale, regionsAfterPres, result},
  Global`addRec[<|"step" -> s, "fn" -> "FindAllRegion", "phase" -> "in",
           "topposi" -> Global`toStr[topposi]|>];
  branch0 = ToSquareAll[ReducedPropagator][[1]] /. Thread[ReducedLeg -> 0];
  branch = DeleteDuplicates[branch0[[topposi]]];
  cutbranch = If[Head[Cut] === List, Pick[branch0, Cut, 1], {}];
  cutposi = Flatten@Table[Position[Expand[branch - cutbranch[[i]]], 0],
                          {i, Length@cutbranch}];
  trans = Select[BranchToLoop /@ Tuples[branch, Length@Loop], # =!= $Failed &];
  patt = Tuples[{0, 1}, Length@Loop];
  regions = Join @@ Table[
    {trans[[i]], patt[[j]], BranchScale[branch /. trans[[i]], patt[[j]]]},
    {i, Length@trans}, {j, Length@patt}
  ];
  regions = Select[regions, AllTrue[#[[-1, cutposi]], # === 0 &] &];
  regionsByScale = GatherBy[regions, Last];
  regions = regionsByScale[[All, All, ;; 2]];
  pres = PrescriptionOf /@ branch;
  regionsAfterPres = Select[#, PrescriptionOf /@ Expand[branch /. #[[1]]] === pres &] & /@ regions;
  Global`addRec[<|"step" -> s, "fn" -> "FindAllRegion", "phase" -> "mid",
           "branch0" -> Global`toStr[branch0], "branch" -> Global`toStr[branch],
           "cutbranch" -> Global`toStr[cutbranch], "cutposi" -> Global`toStr[cutposi],
           "trans" -> Global`toStr[trans], "patt" -> Global`toStr[patt],
           "regions_pre_gather" -> Global`toStr[regionsByScale],
           "pres" -> Global`toStr[pres],
           "regions_post_pres" -> Global`toStr[regionsAfterPres]|>];
  If[AnyTrue[regionsAfterPres, # === {} &],
    Print["FindAllRegion: some regions are prohibited by prescriptions."];
    Abort[]
  ];
  result = regionsAfterPres[[All, 1]];
  Global`addRec[<|"step" -> s, "fn" -> "FindAllRegion", "phase" -> "out",
           "result" -> Global`toStr[result]|>];
  result
];

Unprotect[ZeroRegionQ];
ZeroRegionQ[region_, integrals_] := Block[
  {trans, scale, fullde, factor, expde, s = Global`nextStep[], sortedExpde, groupIdx,
   groupedExpde, result},
  Global`addRec[<|"step" -> s, "fn" -> "ZeroRegionQ", "phase" -> "in",
           "region" -> Global`toStr[region], "integrals" -> Global`toStr[integrals]|>];
  {trans, scale} = region;
  fullde = Expand[ReducedPropagator[[GetTopPosition[integrals]]] /. trans /. RegionRule[scale]];
  factor = If[FreeQ[#, $Eta], 1, $Eta] & /@ fullde;
  fullde = Expand[fullde*factor^-1];
  expde = Coefficient[fullde, $Eta, 0];
  sortedExpde = SortBy[expde, (# /. Thread[Join[Loop, ReducedLeg] -> 0]) === 0 &];
  groupIdx = MaximalGroup[Coefficient[#, SPList[]] & /@ sortedExpde];
  groupedExpde = sortedExpde[[groupIdx]];
  result = ZeroSectorQ[groupedExpde];
  Global`addRec[<|"step" -> s, "fn" -> "ZeroRegionQ", "phase" -> "out",
           "fullde" -> Global`toStr[fullde], "factor" -> Global`toStr[factor],
           "expde_raw" -> Global`toStr[expde], "expde_sorted" -> Global`toStr[sortedExpde],
           "maximal_group" -> Global`toStr[groupIdx],
           "expde_grouped" -> Global`toStr[groupedExpde],
           "result" -> Global`toStr[result]|>];
  result
];

Protect[AnalyzeTopology, AMFCandidateComponent, FindAllRegion, ZeroRegionQ];
End[];

(* Generic wrap for 1..5 argument symbols.  Caller provides the literal
   symbol (not a string) so we avoid ToExpression/Symbol context
   ambiguities. *)
SetAttributes[wrap1, HoldFirst];
wrap1[sym_, name_String] := Module[{orig},
  orig = DownValues[sym];
  If[orig === {}, Print["wrap1: no DV for ", name]; Return[]];
  Unprotect[sym];
  DownValues[sym] = {};
  sym[a_] := Module[{s = nextStep[], result},
    addRec[<|"step" -> s, "fn" -> name, "phase" -> "in",
             "a" -> toStr[a]|>];
    result = Block[{sym}, DownValues[sym] = orig; sym[a]];
    addRec[<|"step" -> s, "fn" -> name, "phase" -> "out",
             "result" -> toStr[result]|>];
    result
  ];
];

SetAttributes[wrap2, HoldFirst];
wrap2[sym_, name_String] := Module[{orig},
  orig = DownValues[sym];
  If[orig === {}, Print["wrap2: no DV for ", name]; Return[]];
  Unprotect[sym];
  DownValues[sym] = {};
  sym[a_, b_] := Module[{s = nextStep[], result},
    addRec[<|"step" -> s, "fn" -> name, "phase" -> "in",
             "a" -> toStr[a], "b" -> toStr[b]|>];
    result = Block[{sym}, DownValues[sym] = orig; sym[a, b]];
    addRec[<|"step" -> s, "fn" -> name, "phase" -> "out",
             "result" -> toStr[result]|>];
    result
  ];
];

SetAttributes[wrap3, HoldFirst];
wrap3[sym_, name_String] := Module[{orig},
  orig = DownValues[sym];
  If[orig === {}, Print["wrap3: no DV for ", name]; Return[]];
  Unprotect[sym];
  DownValues[sym] = {};
  sym[a_, b_, c_] := Module[{s = nextStep[], result},
    addRec[<|"step" -> s, "fn" -> name, "phase" -> "in",
             "a" -> toStr[a], "b" -> toStr[b], "c" -> toStr[c]|>];
    result = Block[{sym}, DownValues[sym] = orig; sym[a, b, c]];
    addRec[<|"step" -> s, "fn" -> name, "phase" -> "out",
             "result" -> toStr[result]|>];
    result
  ];
];

SetAttributes[wrap4, HoldFirst];
wrap4[sym_, name_String] := Module[{orig},
  orig = DownValues[sym];
  If[orig === {}, Print["wrap4: no DV for ", name]; Return[]];
  Unprotect[sym];
  DownValues[sym] = {};
  sym[a_, b_, c_, d_] := Module[{s = nextStep[], result},
    addRec[<|"step" -> s, "fn" -> name, "phase" -> "in",
             "a" -> toStr[a], "b" -> toStr[b],
             "c" -> toStr[c], "d" -> toStr[d]|>];
    result = Block[{sym}, DownValues[sym] = orig; sym[a, b, c, d]];
    addRec[<|"step" -> s, "fn" -> name, "phase" -> "out",
             "result" -> toStr[result]|>];
    result
  ];
];

(* Wrap only high-level Layer-16 pipeline entry points by literal
   symbol reference.  Low-level helpers (RegionRule, ZeroRegionQ, etc.)
   are NOT wrapped to avoid disturbing internal usage patterns that
   Block[{sym}, ...] scoping tends to break. *)
wrap2[AMFlow`RegionPower,          "RegionPower"];
wrap3[AMFlow`BoundaryIntegrands,   "BoundaryIntegrands"];
wrap2[AMFlow`BoundaryIntegrals,    "BoundaryIntegrals"];
wrap4[AMFlow`ReduceBoundary,       "ReduceBoundary"];
wrap1[AMFlow`AMFSystemBoundaryCondition, "AMFSystemBoundaryCondition"];
wrap2[AMFlow`AMFSystemSolution,    "AMFSystemSolution"];

Print["wrap complete — trace count: ", Length[$traceRecords]];

(* --- Run SolveIntegrals for bubble at eps=1/100 --- *)

Print["==BEGIN=="];
precision = 30;
epsorder  = 4;
targets   = {j[bubblefam, 1, 1], j[bubblefam, 1, 0], j[bubblefam, 0, 1]};

t0 = AbsoluteTime[];
sol = SolveIntegrals[targets, precision, epsorder];
t1 = AbsoluteTime[];
Print["==DONE== time = ", t1 - t0, " s"];
Print["==RESULT=="];
Do[Print[sol[[i]]], {i, Length[sol]}];

Print["==TRACE_RECORDS_COUNT== ", Length[$traceRecords]];

outfile = Environment["TRACE_OUT"];
If[outfile === $Failed || outfile === None,
  outfile = FileNameJoin[{rootdir, "build", "trace", "mma_bubble.jsonl"}]];
If[!DirectoryQ[DirectoryName[outfile]], CreateDirectory[DirectoryName[outfile]]];
stream = OpenWrite[outfile];
Do[
  WriteLine[stream,
    ExportString[rec, "RawJSON", "Compact" -> True]],
  {rec, $traceRecords}];
Close[stream];
Print["==WROTE== ", outfile];

Quit[];
