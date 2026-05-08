If[$FrontEnd === Null, $InputFileName, NotebookFileName[]] // DirectoryName // SetDirectory;

Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "reference", "amflow-master",
  "diffeq_solver", "DESolver.m"}]];

cacheDir = FileNameJoin[{DirectoryName[$InputFileName], "cache", "sunrise_amflow", "1"}];

masters = Get[FileNameJoin[{cacheDir, "masters"}]];
diffeq = Get[FileNameJoin[{cacheDir, "diffeq"}]];
boundary = Get[FileNameJoin[{cacheDir, "boundary"}]];
direction = Get[FileNameJoin[{cacheDir, "direction"}]];
epslist = Get[FileNameJoin[{cacheDir, "epslist"}]];
boundaryMI = Get[FileNameJoin[{cacheDir, "boundarymi"}]];
{regions, powers, pattern} = Get[FileNameJoin[{cacheDir, "bpattern"}]];

pattern = Select[pattern, AnyTrue[Flatten[#[[1]] - boundary[[All, 1]]], IntegerQ] &];

evaluate[tab_, list_] := Map[# . list &, tab /. epsRule, {2}];
tobcs[beh_, tab_] := MapThread[Thread[#1 - Range[0, Length @ #2 - 1] -> #2] &, {beh /. epsRule, tab}];

n = 1;
epsRule = {eps -> epslist[[n]]};
bmi = Map[N[#, $MinPrecision] &, Values[boundaryMI][[All, All, n]], {2}];
de = Together[diffeq /. epsRule];
bc = MapThread[tobcs[#1, evaluate[#2, #3]] &, {boundary[[All, 1]], boundary[[All, -1]], bmi}];
bc = Join @@@ Transpose[bc];
bc = MapThread[Join[#1, #2] &, {bc, Transpose[Thread /@ Thread[(pattern /. epsRule) -> 0]]}];

SetDefaultOptions[];
SetGlobalOptions["WorkingPre" -> 232, "ChopPre" -> 20, "RationalizePre" -> 100];
SetExpansionOptions["XOrder" -> 464, "ExtraXOrder" -> 20, "LearnXOrder" -> -1, "TestXOrder" -> 5];
SetRunningOptions["RunDirection" -> direction];

LoadSystem[$InternalSystem, de, bc, Infinity];
run = RunEta[GetPoles[DE[$InternalSystem]]];

Print["RUN=", InputForm[N[run, 30]]];

allrule = CalcInf[de, bc];
Print["CALCINF_FIRST5=", InputForm@Table[N[allrule[[i, 1, 2, 1, 1 ;; 5]], 30], {i, Length[allrule]}]];

InfToRegular[$InternalSystem, run[[1]]];
Print["BC_AFTER_INF=", InputForm[N[BC[$InternalSystem], 30]]];

RegularRun[$InternalSystem, run];
Print["BC_AFTER_RUN=", InputForm[N[BC[$InternalSystem], 30]]];

SolveAsyExp[$InternalSystem];
Print["ZERO_SOL=", InputForm[N[PickZeroRuleS /@ AsyExp[$InternalSystem], 30]]];

ClearSystem[$InternalSystem];
Quit[];
