If[$FrontEnd===Null,$InputFileName,NotebookFileName[]]//DirectoryName//SetDirectory;
Get["/home/chang/workspace/agent/trials/amflowRefactor/reference/amflow-master/diffeq_solver/DESolver.m"];
SetDefaultOptions[];
SetExpansionOptions["ExtraXOrder" -> 480];


If[!FileExistsQ["diffeq"], Print["error: diffeq not found."]; Abort[]];
diffeq = Get["diffeq"];
{regions, powers, pattern} = Get["bpattern"];


epsrule = {eps -> 10^-4};
deinf = Together[-Power[eta,-2](diffeq/.epsrule/.eta -> 1/eta)];
npattern = -pattern/.epsrule;
orders = DetermineBoundaryOrder[deinf, #]&/@npattern;


border = Table[
k = Select[Range[Length[pattern]], IntegerQ[pattern[[#, 1]]-powers[[i, 1]]]&][[1]];
orders[[k]]-(pattern[[k]]-powers[[i]]),  {i, Length[powers]}];
border = Map[If[#<0,-1,#]&, border, {2}];


Put[border, "border"];
Quit[];