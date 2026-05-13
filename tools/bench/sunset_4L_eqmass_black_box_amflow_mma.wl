(* ::Package:: *)

(*
  4-loop equal-mass sunset (banana) — Phase 3 oracle (loop axis = 4 +
  equal-mass uniformity).  5 internal lines all carrying the SAME
  mass^2 = msq, between two external vertices (p enters one, exits
  the other).  Distinct from:
    * banana_4L_mixed (4L, ALTERNATING A,B,A,B,A masses)
    * cutbanana_4L    (4L, massless, all 5 cuts on top sector)
    * sunset_2L       (2L, equal mass, lower loops)
  Numerics: psq = -3 (Euclidean, below threshold 5*msq = 5), msq = 1.
  BlackBoxDot = 5 to avoid the "inconsistent masters" abort pattern at
  low dot in 4-loop banana topology.
*)

current = If[$FrontEnd === Null, $InputFileName, NotebookFileName[]] // DirectoryName;
repo = FileNameJoin[{current, "..", ".."}];
Get[FileNameJoin[{repo, "reference", "amflow-master", "AMFlow.m"}]];

SetReductionOptions[
  "IBPReducer" -> "Kira",
  "BlackBoxRank" -> 0,
  "BlackBoxDot" -> 5,
  "DeleteBlackBoxDirectory" -> False
];
Kira`$KiraExecutable = "/usr/local/bin/kira";
Kira`$FermatExecutable = "/usr/share/Ferl7/fer64";

SetAMFOptions[
  "AMFMode" -> {"Prescription", "Mass", "Propagator"},
  "EndingScheme" -> {"Tradition", "SingleMass"},
  "WorkingPre" -> 120, "ChopPre" -> 20, "RationalizePre" -> 100,
  "XOrder" -> 240, "ExtraXOrder" -> 280,
  "LearnXOrder" -> -1, "TestXOrder" -> 5
];

AMFlowInfo["Family"] = sunset4eq;
AMFlowInfo["Loop"] = {l1, l2, l3, l4};
AMFlowInfo["Leg"] = {p1};
AMFlowInfo["Conservation"] = {};
AMFlowInfo["Replacement"] = {p1^2 -> psq};
AMFlowInfo["Propagator"] = {
  l1^2 - msq,
  l2^2 - msq,
  l3^2 - msq,
  l4^2 - msq,
  (l1 + l2 + l3 + l4 - p1)^2 - msq,
  (l1 + p1)^2,
  (l2 + p1)^2,
  (l3 + p1)^2,
  (l4 + p1)^2,
  (l1 - l2)^2,
  (l1 - l3)^2,
  (l1 - l4)^2,
  (l2 - l3)^2,
  (l2 - l4)^2
};
AMFlowInfo["Numeric"] = {psq -> -3, msq -> 1};
AMFlowInfo["NThread"] = 4;

targets = {j[sunset4eq, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]};
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "sunset_4L_eqmass_mma_cache"}]];

Print["==BENCH== mma 4-loop equal-mass sunset BlackBoxAMFlow"];
Print["==CACHE== ", cache];
Print["==EPS== ", epslist // InputForm];

elapsed = AbsoluteTiming[
  sol = BlackBoxAMFlow[targets, epslist, cache];
][[1]];

Print["==TIME== ", elapsed];
Print["==RESULT_INPUTFORM== ", sol // InputForm];

Do[
  key = Keys[sol][[i]];
  val = Values[sol][[i, 1]];
  Print["==VALUE== ", key // InputForm, " :: ",
    N[Re[val], 50] // InputForm, " :: ", N[Im[val], 50] // InputForm],
  {i, Length[Keys[sol]]}
];

Quit[];
