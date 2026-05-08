(* Stage A trace driver: run MMA's AMFlow on the 1x1 PowerLaw case,
   capturing all key Layer-6 intermediates via Sow/Reap + Block-trick.

   Usage:
     TRACE_OUT=/path/to/trace.jsonl math -script trace_minimal.wl

   Produces a JSON-Lines file (one JSON object per line) suitable for
   side-by-side diff against the C++ dump produced by DESOLVER_DEBUG_BC. *)

scriptDir = DirectoryName[$InputFileName];
repoRoot = ExpandFileName[FileNameJoin[{scriptDir, "..", ".."}]];
dsdir = FileNameJoin[{repoRoot, "reference", "amflow-master", "diffeq_solver"}];
Get[FileNameJoin[{dsdir, "DESolver.m"}]];

SetDefaultOptions[];
SetGlobalOptions["SilentMode" -> True];
(* Match C++ defaults (XOrder=100, ExtraXOrder=20) so the sparse matrix
   dimensions and fid cutoffs agree with the C++ side. *)
SetExpansionOptions["XOrder" -> 100, "ExtraXOrder" -> 20];

(* ---------------------------------------------------------------------------- *)
(*  Tracing infrastructure                                                     *)
(* ---------------------------------------------------------------------------- *)

$stepCounter = 0;
nextStep[] := (++$stepCounter);

(* Convert an expression to a canonical string suitable for diffing. *)
toStr[expr_] := ToString[expr // InputForm,
                         CharacterEncoding -> "ASCII"];

(* Convert a numeric (possibly with extended precision) to a shortened
   decimal string for easier reading. *)
toN[expr_, digits_:30] := ToString[N[expr, digits] // InputForm,
                                   CharacterEncoding -> "ASCII"];

(* Emit a trace record. *)
sowRec[data_] := Sow[data, "TR"];

(* ---------------------------------------------------------------------------- *)
(*  Wrap each target function                                                   *)
(* ---------------------------------------------------------------------------- *)

(* AMFlow (top-level entry) *)
origAMFlow = DownValues[AMFlow];
Unprotect[AMFlow]; DownValues[AMFlow] = {};
AMFlow[de_, bcs_] := Module[{s = nextStep[], result},
  sowRec[<|"step" -> s, "fn" -> "AMFlow", "phase" -> "in",
           "de" -> toStr[de], "bcs" -> toStr[bcs]|>];
  result = Block[{AMFlow}, DownValues[AMFlow] = origAMFlow;
                           AMFlow[de, bcs]];
  sowRec[<|"step" -> s, "fn" -> "AMFlow", "phase" -> "out",
           "result" -> toStr[result], "result_N" -> toN[result, 30]|>];
  result
];

(* CalcInf *)
origCalcInf = DownValues[CalcInf];
Unprotect[CalcInf]; DownValues[CalcInf] = {};
CalcInf[de_, bcs_] := Module[{s = nextStep[], result},
  sowRec[<|"step" -> s, "fn" -> "CalcInf", "phase" -> "in",
           "de" -> toStr[de], "bcs" -> toStr[bcs]|>];
  result = Block[{CalcInf}, DownValues[CalcInf] = origCalcInf;
                            CalcInf[de, bcs]];
  sowRec[<|"step" -> s, "fn" -> "CalcInf", "phase" -> "out",
           "result" -> toStr[result]|>];
  result
];

(* BuildTaylor *)
origBT = DownValues[DESolver`Private`BuildTaylor];
Unprotect[DESolver`Private`BuildTaylor];
DownValues[DESolver`Private`BuildTaylor] = {};
DESolver`Private`BuildTaylor[mat_, ini_] := Module[{s = nextStep[], result},
  sowRec[<|"step" -> s, "fn" -> "BuildTaylor", "phase" -> "in",
           "mat" -> toStr[mat], "ini" -> toStr[ini]|>];
  result = Block[{DESolver`Private`BuildTaylor},
    DownValues[DESolver`Private`BuildTaylor] = origBT;
    DESolver`Private`BuildTaylor[mat, ini]];
  sowRec[<|"step" -> s, "fn" -> "BuildTaylor", "phase" -> "out",
           "result" -> toStr[result]|>];
  result
];

(* PoincareRank *)
origPR = DownValues[DESolver`Private`PoincareRank];
Unprotect[DESolver`Private`PoincareRank];
DownValues[DESolver`Private`PoincareRank] = {};
DESolver`Private`PoincareRank[mat_] := Module[{s = nextStep[], result},
  sowRec[<|"step" -> s, "fn" -> "PoincareRank", "phase" -> "in",
           "mat" -> toStr[mat]|>];
  result = Block[{DESolver`Private`PoincareRank},
    DownValues[DESolver`Private`PoincareRank] = origPR;
    DESolver`Private`PoincareRank[mat]];
  sowRec[<|"step" -> s, "fn" -> "PoincareRank", "phase" -> "out",
           "result" -> toStr[result]|>];
  result
];

(* NHEquations *)
origNHE = DownValues[NHEquations];
Unprotect[NHEquations]; DownValues[NHEquations] = {};
NHEquations[mat_, mode_] := Module[{s = nextStep[], result},
  sowRec[<|"step" -> s, "fn" -> "NHEquations", "phase" -> "in",
           "mat" -> toStr[mat], "mode" -> toStr[mode]|>];
  result = Block[{NHEquations}, DownValues[NHEquations] = origNHE;
                                NHEquations[mat, mode]];
  sowRec[<|"step" -> s, "fn" -> "NHEquations", "phase" -> "out",
           "result" -> toStr[result]|>];
  result
];

(* NHEquationsNum *)
origNHEN = DownValues[DESolver`Private`NHEquationsNum];
Unprotect[DESolver`Private`NHEquationsNum];
DownValues[DESolver`Private`NHEquationsNum] = {};
DESolver`Private`NHEquationsNum[nheqs_] := Module[{s = nextStep[], result},
  sowRec[<|"step" -> s, "fn" -> "NHEquationsNum", "phase" -> "in",
           "nheqs" -> toStr[nheqs]|>];
  result = Block[{DESolver`Private`NHEquationsNum},
    DownValues[DESolver`Private`NHEquationsNum] = origNHEN;
    DESolver`Private`NHEquationsNum[nheqs]];
  sowRec[<|"step" -> s, "fn" -> "NHEquationsNum", "phase" -> "out",
           "result" -> toStr[result]|>];
  result
];

(* ReadBCS *)
origRB = DownValues[DESolver`Private`ReadBCS];
Unprotect[DESolver`Private`ReadBCS];
DownValues[DESolver`Private`ReadBCS] = {};
DESolver`Private`ReadBCS[bcs_, region_] := Module[{s = nextStep[], result},
  sowRec[<|"step" -> s, "fn" -> "ReadBCS", "phase" -> "in",
           "bcs" -> toStr[bcs], "region" -> toStr[region]|>];
  result = Block[{DESolver`Private`ReadBCS},
    DownValues[DESolver`Private`ReadBCS] = origRB;
    DESolver`Private`ReadBCS[bcs, region]];
  sowRec[<|"step" -> s, "fn" -> "ReadBCS", "phase" -> "out",
           "result" -> toStr[result]|>];
  result
];

(* ReverseBCS / UnionBCS *)
origRevB = DownValues[DESolver`Private`ReverseBCS];
Unprotect[DESolver`Private`ReverseBCS];
DownValues[DESolver`Private`ReverseBCS] = {};
DESolver`Private`ReverseBCS[bcs_] := Module[{s = nextStep[], result},
  sowRec[<|"step" -> s, "fn" -> "ReverseBCS", "phase" -> "in", "bcs" -> toStr[bcs]|>];
  result = Block[{DESolver`Private`ReverseBCS},
    DownValues[DESolver`Private`ReverseBCS] = origRevB;
    DESolver`Private`ReverseBCS[bcs]];
  sowRec[<|"step" -> s, "fn" -> "ReverseBCS", "phase" -> "out", "result" -> toStr[result]|>];
  result
];

origUnB = DownValues[DESolver`Private`UnionBCS];
Unprotect[DESolver`Private`UnionBCS];
DownValues[DESolver`Private`UnionBCS] = {};
DESolver`Private`UnionBCS[bcs_] := Module[{s = nextStep[], result},
  sowRec[<|"step" -> s, "fn" -> "UnionBCS", "phase" -> "in", "bcs" -> toStr[bcs]|>];
  result = Block[{DESolver`Private`UnionBCS},
    DownValues[DESolver`Private`UnionBCS] = origUnB;
    DESolver`Private`UnionBCS[bcs]];
  sowRec[<|"step" -> s, "fn" -> "UnionBCS", "phase" -> "out", "result" -> toStr[result]|>];
  result
];

(* CalcTaylor *)
origCT = DownValues[DESolver`Private`CalcTaylor];
Unprotect[DESolver`Private`CalcTaylor];
DownValues[DESolver`Private`CalcTaylor] = {};
DESolver`Private`CalcTaylor[mat_, bc_] := Module[{s = nextStep[], result},
  sowRec[<|"step" -> s, "fn" -> "CalcTaylor", "phase" -> "in",
           "mat" -> toStr[mat], "bc" -> toStr[bc]|>];
  result = Block[{DESolver`Private`CalcTaylor},
    DownValues[DESolver`Private`CalcTaylor] = origCT;
    DESolver`Private`CalcTaylor[mat, bc]];
  sowRec[<|"step" -> s, "fn" -> "CalcTaylor", "phase" -> "out",
           "result" -> toStr[result]|>];
  result
];

(* ConstructMatrix *)
origCM = DownValues[DESolver`Private`ConstructMatrix];
Unprotect[DESolver`Private`ConstructMatrix];
DownValues[DESolver`Private`ConstructMatrix] = {};
DESolver`Private`ConstructMatrix[dx_, ax_, totalorder_] := Module[{s = nextStep[], result},
  sowRec[<|"step" -> s, "fn" -> "ConstructMatrix", "phase" -> "in",
           "dx" -> toStr[dx], "ax" -> toStr[ax], "totalorder" -> totalorder|>];
  result = Block[{DESolver`Private`ConstructMatrix},
    DownValues[DESolver`Private`ConstructMatrix] = origCM;
    DESolver`Private`ConstructMatrix[dx, ax, totalorder]];
  sowRec[<|"step" -> s, "fn" -> "ConstructMatrix", "phase" -> "out",
           "result" -> toStr[result]|>];
  result
];

(* SparseGaussian *)
origSG = DownValues[DESolver`Private`SparseGaussian];
Unprotect[DESolver`Private`SparseGaussian];
DownValues[DESolver`Private`SparseGaussian] = {};
DESolver`Private`SparseGaussian[sp_, nh_] := Module[{s = nextStep[], result},
  sowRec[<|"step" -> s, "fn" -> "SparseGaussian", "phase" -> "in",
           "sp" -> toStr[sp], "nh" -> toStr[nh]|>];
  result = Block[{DESolver`Private`SparseGaussian},
    DownValues[DESolver`Private`SparseGaussian] = origSG;
    DESolver`Private`SparseGaussian[sp, nh]];
  sowRec[<|"step" -> s, "fn" -> "SparseGaussian", "phase" -> "out",
           "result" -> toStr[result]|>];
  result
];

(* ---------------------------------------------------------------------------- *)
(*  Run the case                                                                *)
(* ---------------------------------------------------------------------------- *)

de  = {{ (99/100)/(eta - 3) }};
bcs = {{ 99/100 -> 100 }};

{finalResult, sown} = Reap[AMFlow[de, bcs], "TR"];
traceList = If[sown === {}, {}, sown[[1]]];

Print["==N_TRACE_RECORDS== ", Length[traceList]];
Print["==FINAL_N== ", toN[finalResult, 30]];

(* Serialize to JSON Lines. *)
outfile = Environment["TRACE_OUT"];
If[outfile === $Failed || outfile === None,
  outfile = FileNameJoin[{repoRoot, "build", "trace", "trace_minimal.jsonl"}]];
If[!DirectoryQ[DirectoryName[outfile]], CreateDirectory[DirectoryName[outfile]]];

stream = OpenWrite[outfile];
Do[WriteLine[stream, ExportString[rec, "RawJSON",
                                  "Compact" -> True,
                                  "ConversionRules" -> {}]],
   {rec, traceList}];
Close[stream];

Print["==WROTE== ", outfile];

Quit[];
