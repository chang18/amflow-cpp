If[$FrontEnd === Null, $InputFileName, NotebookFileName[]] // DirectoryName // SetDirectory;

Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "reference", "amflow-master",
  "AMFlow.m"}]];

cacheDir = FileNameJoin[{DirectoryName[$InputFileName], "cache", "sunrise_amflow", "1"}];
epslist = Get[FileNameJoin[{cacheDir, "epslist"}]];
sol = Get[FileNameJoin[{cacheDir, "solution"}]];

vals = j[sunrise, 1, 1, 1, 0, 0] /. sol;

SetAMFOptions["WorkingPre" -> 232, "ChopPre" -> 20];

expr = FitEps[epslist, vals, -4];
series = Normal @ Series[expr, {$Eps, 0, -1}];

Print["FIT_SERIES=", InputForm[N[series, 40]]];
Print["COEFF_M2=", InputForm[N[SeriesCoefficient[expr, {$Eps, 0, -2}], 40]]];
Print["COEFF_M1=", InputForm[N[SeriesCoefficient[expr, {$Eps, 0, -1}], 40]]];

Quit[];
