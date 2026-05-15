# upstream `AMFlow.m` — symbol map

(Line numbers track <https://gitlab.com/multiloop-pku/amflow> at the
time of the v1.0.0 release.)

1489 lines. Outer structure:

- 7-167: `BeginPackage` + usage declarations.
- 169: `Begin["`Private`"]`.
- 172-401: common package state, options, generic utilities, integral notation.
- 398-444: scalar-product / propagator algebra.
- 446-525: UF polynomials and topology / family analysis.
- 527-621: eta scheme (AMFMode candidate enumeration).
- 623-812: boundary / region machinery.
- 815-824: single-mass vacuum lookup table.
- 827-1204: AMF system setup + iterative setup + solve.
- 1207-1306: black-box wrappers (BlackBoxReduce / BlackBoxDiffeq / BlackBoxAMFlowSingle / BlackBoxAMFlow).
- 1309-1361: `SolveIntegrals` (the Laurent-fitted top-level entry).
- 1364-1476: `SolveIntegralsGaugeLink` (Layer 17b, **out of scope**).

## Globals & options

`AMFlow.m:181-183`: `$Eps -> Global`eps`, `$Eta -> Global`eta`, `$JHead -> Global`j`.

`AMFlow.m:188`: `$AMFlowInfoKeywords = {Family, Loop, Leg, Conservation, Replacement, Propagator, Numeric, NThread, Prescription, Cut}`.

`AMFlow.m:191-200`: shorthand accessors `Family / Loop / Leg / Conservation / Replacement / Propagator / Numeric / NThread / Prescription / Cut := AMFlowInfo[...]`.

`AMFlow.m:203-212`: derived helpers
- `ReducedLeg` — legs not eliminated by Conservation.
- `ReducedReplacement` — replacement rules expressed only in `ReducedLeg . ReducedLeg` scalar products.
- `ReducedPropagator := Propagator/.Conservation`.

`AMFlow.m:215-222`: AMFlowInfo state ops `DefinedKeywords / ExtractGlobalVariables / DefineGlobalVariables / UnsetGlobalVariables / SaveGlobalVariables / LoadGlobalVariables / PrintGlobalVariables`. Save/Load use `$Global[n]` keyed snapshots.

`AMFlow.m:258` `Options[SetAMFOptions]` defaults:
- `AMFMode -> {"Prescription", "Mass", "Propagator"}`
- `EndingScheme -> {"Tradition", "Cutkosky", "SingleMass"}`
- `D0 -> 4`
- `WorkingPre -> 100`, `ChopPre -> 20`, `RationalizePre -> 100`
- `XOrder -> 100`, `ExtraXOrder -> 20`, `LearnXOrder -> -1`, `TestXOrder -> 5`

`AMFlow.m:277` `Options[SetReductionOptions]` defaults:
- `IBPReducer -> "FiniteFlow+LiteRed"` (the project uses `Kira` instead — see `ibp_interface/Kira/interface.m`)
- `BlackBoxRank -> 3`, `BlackBoxDot -> 0`
- `ComplexMode -> True`, `DeleteBlackBoxDirectory -> False`

## Integral notation (361-395)

- `ToJ[ind]` → `j[Family, ind...]`.
- `FromJ[jint]` → indices.
- `JSector` (sector mask 0/1), `JProp` (#props with index>0), `JDot` (sum of dot indices−1), `JRank` (sum of −indices for index<0).
- `IntegralWeight[jint] = {-JProp, -JDot, -JRank, Min[indices], jint}`; sort with `SortIntegrals`.
- `GetTopSector[jints]` — sector pattern with `_` at every position that has any positive index across the list.
- `GetTopPosition` — list of those positions.
- `IntegralQ`, `CheckIntegrals` — shape sanity gates.
- `GetTopSectorList[jlist]` — distinct top sectors among a list (each maximal among the others).
- `SplitTarget[jlist]` (386) — partitions integrals by which top sector contains them. Used by `BlackBoxAMFlow` to chunk targets per top sector.

## SP/D algebra (402-443)

- `SPList[]` — all loop-loop and loop-leg scalar products.
- `DListSymbol[]` — `Denominators[1..n]`.
- `CheckCompleteness[denoms]` — rank of SP coefficient matrix must equal `Length[SPList[]]`.
- `SPListToDListSymbol[denoms]` — invert the SP matrix to express SP in terms of `Denominators[i]`.
- `ToCompleteExplicit[denoms]` (425) — pad partial denominators with squared loop+leg combinations until SP-complete; uses `MaximalGroup`.
- `ToSquare / ToSquareAll / SquaredDenominators` — write a propagator as `momentum^2 + mass`.

## UF polynomials (453-525)

- `FeynmanPara[n]` (450) → `FeynmanParameters[1..n]`.
- `EvaluateABC[denominators]` (453) — A (loop quadratic), B (loop linear), C (constant), each as `Σ x_i (...)`.
- `EvaluateUF[denominators]` (466) — `u = Det[A]`, `f0 = - mass . FeynmanPara`, `f = Expand[Together[u(-Cx + Bx . Inverse[Ax] . Bx)] / .ReducedReplacement] - u*f0`. Returns `{u, f, f0}`. `u===0` → `{0,0,0}`.
- `ZeroSectorQ[denominators]` (476) — Lee–Pomeransky scaleless test on `g = u + f + u*massterm` after substituting `Numeric`.
- `AnalyzeComponent[u0, f, massterm]` (488) — for one factor of u: returns `{u0, vacQ, loopnum, var, mass}`.
- `AnalyzeTopology[denominators]` (499) — factor `u` and analyze each monomial component, sort by loopnum.
- `FactorizeFamily[denominators, patts]` (510) — for single-mass children: per component, redefine loops so each component owns its own loop subset, return `{loop, prop, mas}` triples.

## Eta scheme (527-622)

- `PrescriptionOfLoop / PrescriptionOf` (532, 535) — read `Prescription` per loop momentum then propagate to a propagator (mixed → `$Failed`).
- `AnalyzeTopSector[topposi]` (545) — `AnalyzeTopology[ReducedPropagator[[topposi]]]` then attach cut + prescription per component.

`info0` shape after AnalyzeTopSector: `{u0, vacQ, loopnum, var, mass, cut, pres}`.

- `VacuumQ` (556) — `vacQ && all cut == 0`.
- `SingleMassQ` (557) — `VacuumQ && exactly one mass == 1, rest == 0`.
- `PhaseVolumeQ` (558) — `all cut == 1 && Length[var] == loopnum + 1`.
- `EndingQ` (559) — `SingleMassQ || PhaseVolumeQ`.

`AllPossiblePosition[info0, mode]`:
- `Prescription` (562) — uncut props that have `pres = -1`.
- `Mass` (566) — group props by independent mass.
- `Propagator` (575) — uncut props sorted by branch length.
- `Branch` (581) — uncut props grouped by being in the same branch.
- `Loop` (589) — uncut props grouped by being in the same loop.
- `All` (598) — flatten `Propagator`.

`AMFCandidateComponent[info0, mode]` (604) — picks first candidate, falling back to Branch if not single-mass.
`AMFCandidate[info, mode]` (610), `AMFPosition[topposi, mode]` (617) — final eta-position selection across components/modes.
`AMFEtaC[posi]` (621) — coefficient vector with -1 at each chosen position, 0 elsewhere.

## Boundary / region (624-812)

- `BranchToLoop[candidate]` (628) — solve for loop = inverse(coeff) . loop; `$Failed` if singular.
- `BranchScale[branch, scale]` (631) — per branch, 0 if free of all "large" loops (scale==1), else 1.
- `FindAllRegion[topposi]` (634) — enumerate all `(trans, patt)` with branch transformations and binary scales, dedup by branch scale, drop scale-incompatible cuts and prescription-violating regions.
- `RegionRule[patt]` (651) — `Loop -> eta^(patt/2) Loop`.
- `ZeroRegionQ[region, integrals]` (654) — leading subsector of the region is scaleless.
- `ExplicitRetionRule[region]` (666) — combine trans + RegionRule for printing.
- `RegionPower[integrals, region]` (669) — leading-power formula `(2-eps) Total[scale] - Sum[powers indexed by props with eta factor]` for each integral.
- `BoundaryPattern[powers]` (680) — group rows by integer-shift equivalence on their first column, return one representative per group with the per-row maxima above the representative. **This is the actual leading-power per region used to set BC for the ODE.**

`BoundaryIntegrands[integrals, border, region]` (704) — **the sub-leading integrand engine**:
1. Apply trans + RegionRule to ReducedPropagator → `fullde`.
2. Strip `eta` factors → `factor`. `expde = Coefficient[fullde, $Eta, 0]`.
3. Complete the leading subsector with `ToCompleteExplicit[expde[topposi]]`.
4. `sptoD = SPListToDListSymbol[completede]`. Re-express `fullde / expde` in terms of D-symbols.
5. `structure = Coefficient[fullde, $Eta, {0, -1/2, -1}]` per propagator.
6. Build `integrand = ∏ (x[i,1] + x[i,2] eta^(-1/2) + x[i,3] eta^(-1))^(-p[i])`. Substitute zero-coefficient rules and equality rules.
7. `integrand = Together[SeriesCoefficient[..., {eta, ∞, k}]]` for `k = 0..Max[border]`.
8. Re-substitute the `x[i,j]` to actual structure entries.
9. For each integral i, take `integrand[[1..border[k]+1]]` and substitute the indices `FromJ[integrals[k]]` for `p`.
10. Return `{SquaredDenominators[completede], SimplifyIntegrand[integrands]}`.

`ApartRationals[rationals]` (732) — repeated `Apart` over each D-symbol.
`LaportaIntegrals[term]` (743) — denominator must be a monomial in D-symbols; returns `{D-shift -> coef}`.
`BoundaryIntegrals[completede, integrands]` (753) — partition by `analyze[term]` (which D-symbols appear in the denominator), apart, group by D-signature → returns `{family, partitionedIntegrals}` per D-signature, where each `family = completede + shifts`.
`ReduceBoundary[integrals, border, region, dir]` (784) — for each new family produced by `BoundaryIntegrals`: snapshot AMFlowInfo, set new propagators (and projected Cut), call `BlackBoxReduce`, restore. Return `{powers, config, masters, table}` per family. The `table` rows already encode the boundary integrand reduction in masters.

## Single-mass vacuum table (819-824)

- `Vacuum[1, 1]`, `Vacuum[2, 3]`, `Vacuum[3, 4]`, `Vacuum[3, 5]`, `Vacuum[4, 5]` — closed-form Gamma expressions.
- `Vacuum[a___]` — Abort.

## AMF system file map (835-859)

`AMFSystemPath` keys: `Config`, `Preferred`, `GlobalPreferred`, `SubSystem`, `DiffeqSetup`, `Masters`, `Diffeq`, `BPattern`, `BOrder`, `BoundaryReduce`, `Boundary`, `Direction`, `BoundaryMI`, `Solution`, `EpsList`.

## Per-system setup (866-991)

Order of operations inside `AMFSystemSetup[sysid, preferred, etac]` (972):

1. `AMFSystemInitialize` (866) — set `Propagator = ReducedPropagator + eta * etac`, write Config + Preferred.
2. If `etac == 0...0`: ending — stop.
3. Otherwise:
   - `AMFSystemDifferentialEquation` (876) — `BlackBoxDiffeq[Preferred, {eta}, ...]`, write Masters + Diffeq.
   - `AMFSystemBoundaryOrder` (890) — `regions = FindAllRegion[GetTopPosition[Masters]]` filtered by `!ZeroRegionQ`, `powers = RegionPower[Masters, region]`, `pattern = BoundaryPattern[powers]`, write `BPattern`. Then template a DESolver script (`SetExpansionOptions[ExtraXOrder]`, `Get[diffeq]`, `DetermineBoundaryOrder[deinf, npattern[k]]`) and run via `RunCommand[$WolframPath, ...]`. Result `border` per region per integral, written to BOrder.
   - `AMFSystemBoundaryCondition` (937) — for each region, `ReduceBoundary[Masters, border[i], region[i], path]` and concat. Written to Boundary.
   - `AMFSystemDirection` (959) — pick `Im` if eta-touching propagators have `pres == -1`, else `NegIm`.

## Ending scheme dispatch (994-1061)

`AMFSystemEndingQ[preferred, scheme]`:
- `Tradition` (995) — `AMFPosition[GetTopPosition[preferred], $AMFMode] === {}`.
- `Cutkosky` (996) — NOT (`all components ending` AND `exactly one PhaseVolume component`).
- `SingleMass` (1002) — NOT (`all ending` AND `no PhaseVolume` AND `Loop nonempty` AND `all preferred indices >= 0`).

`AMFSystemSetupMaster[sysid, preferred, schemes]` (1010):
- If `AllTrue[endingQ]` for all schemes: regular `AMFSystemSetup` with `etac = 0`. `GlobalPreferred = {preferred, 1, preferred, Identity}`. Returns `{sysid}`.
- Otherwise: pick **first** scheme where `endingQ === False` and recurse into the scheme-specific overload.

Per-scheme overloads:

- `Tradition` (1019) — `etac = AMFEtaC[AMFPosition[posi, $AMFMode]]`. `GlobalPreferred = {preferred, 1, preferred, Identity}`. Returns `{sysid}`.

- `Cutkosky` (1026) — **clears Prescription + Cut globally** (the legendary `AMFlow.m:1031`), prefactor `2 (Pi^(2-eps) (2Pi)^(2eps-4))^L (-1)^(L+1)` per master with `L = cutcom[1, 3]` = loopnum of the unique phase-volume component, eta inserted into top, op = `Im`. `GlobalPreferred = {preferred, prefactor, preferred, Im}`. Returns `{sysid}`.

- `SingleMass` (1039) — `FactorizeFamily` to per-component children; for each component:
  1. Find `drop` = index of the unique massive prop in component (mass == -1 after `ToSquareAll`).
  2. New AMFlowInfo: `Loop \ {droploop}`, `Leg = {droploop}`, `Conservation = {}`, `Replacement = {droploop^2 -> -1}`, `Propagator = ToCompleteExplicit[Drop[prop, drop]]`, clear Prescription + Cut.
  3. `jpre = ToJ[PadRight[Drop[indices, drop], OriginalLen]]` for each preferred master in this component.
  4. Prefactor `(-1)^(-n) Gamma[2-eps-m] Gamma[-2+eps+m+n] / (Gamma[2-eps] Gamma[n])` where `n = mas[All, drop]`, `m = Total - mas[All, drop] - (loopnum - 1)*(2 - eps)`.
  5. Recursive `AMFSystemSetup[sysid+i-1, jpre, etac]` and `GlobalPreferred = {preferred, prefactor, jpre, Identity}`.
  6. Returns `{sysid, sysid+1, ..., sysid+nComponents-1}` — multi-root recursion.

  Mirrors `AMFlow.m:1056-1057`: the immediate child uses `$AMFMode` to compute `etac` (rather than re-entering scheme dispatch), but the next level inherits the full ending-scheme list via `AMFSystemSetupMaster`.

## Iterative system setup (1068-1093)

`AMFSystemsSetup[preferred]`:
- Run `AMFSystemSetupMaster` to set up the head systems.
- For each head system, read `Boundary`, and for each boundary entry call `AMFSystemsSetup` recursively against its `boundary[[i, 3]]` (the new family masters to solve). Each AMFlowInfo is staged via `DefineGlobalVariables[boundary[[i, 2]]]` per boundary entry.
- Write `SubSystem = list of subsystem id lists per boundary entry`.

`AMFSystemsTree[headsysid]` (1087) — recursive tree of subsystems.
`AMFSystemsEnding[headsysid]` (1093) — leaf systems.

## Solve (1100-1204)

`AMFSystemSolution[sysid, epslist]` (1100):
- Read Preferred + SubSystem; write EpsList.
- Trivial cases:
  - `subsystem = {}` and `GetTopPosition[preferred] = {}` → preferred all 1's (constants).
  - `subsystem = {}` and top sector nonempty → fail with hint to set Solution by hand.
- Else: solve each subsystem first (`AMFSystemsSolution`), collect into `BoundaryMI`. Then template a DESolver script:
  - Read Masters / Diffeq / Boundary / Direction / EpsList / BoundaryMI / BPattern.
  - For each eps point: `de = diffeq /. eps`, `bc = MapThread[tobcs[#1, evaluate[#2, #3]]&, {boundary[All,1], boundary[All,-1], bmi}]`, append pattern-zero rows, call `AMFlow[de, bc]` (the DESolver entry).
  - `ParallelTable` over the eps grid.
- Run via subprocess `$WolframPath -script ...`.

`AMFSystemCombineSolution[sysids]` (1188):
- For each sysid: `{global, coe, local, op} = GlobalPreferred`, `sol = Solution`, `epslist = EpsList`. Build `Thread[global -> coe(eps) * op[local /. sol]]` (across the eps grid).
- Combine across systems by **multiplying** values at the same global key (this is how SingleMass children combine into the parent).

`AMFSystemsSolution[sysids, epslist]` (1200) — drive each system, then combine.

## Black-box wrappers (1218-1306)

- `BlackBoxReduce[jints, jpreferred, dir]` (1218):
  - `top = GetTopSector[Join[jints, jpreferred]]`.
  - `rank = Max[$BlackBoxRank, JRank/@all]`, `dot = Max[$BlackBoxDot, JDot/@all]`.
  - `IBPSystem[top, rank, dot, jpreferred, $ComplexMode, dir]`, `AnalyticReduction[jints]`.
  - Special cases: empty jints → `{{}, {}}`; `Loop = {}` → `{{ToJ[{}]}, {ToJ[{}] -> {1}}}`.
- `BlackBoxDiffeq[jpreferred, vars, dir]` (1235):
  - `dot = Max[$BlackBoxDot, JDot/@jpre + 1]` (note the `+1` — needed for derivative of an integral with dots).
  - `IBPSystem`, `DifferentialEquation[vars]`.
- `BlackBoxAMFlowSingle[jints, epslist, dir]` (1250):
  - Pre-reduce `BlackBoxReduce[jints, {}, .../system 0]` → `{masters, rules}`.
  - Trivial `masters = {}` → `sol = {}`.
  - Else: `AMFSystemsSetup[masters]`, `AMFSystemsSolution[sysids, epslist]` → master grid, then apply the pre-reduce `rules` to project back to `jints`.
  - Padding: integrals not in `Keys[rules]` (zero coefficients) → fill with zeros across the eps grid.
- `BlackBoxAMFlow[jints, epslist, dir]` (1293):
  - `targets = SplitTarget[jints]`.
  - For each part: `BlackBoxAMFlowSingle[part, epslist, dir]`, concat.
  - All parts reuse the same `dir` because `BlackBoxAMFlowSingle`
    itself sets `$AMFDirectory = dir`; upstream `BlackBoxAMFlow` does
    **not** scope per-part directories.  This C++ port's
    `solve_integrals_runner` scopes per-part by `root + "/part_" + i`,
    so a single `dir` can carry multiple top-sector groups without
    Kira cache collisions.

## Top-level entries (1313-1361)

`GenerateNumericalConfig[goal, order]` (1313):
- `loop = Max[Length[Loop], 1]`.
- `number = Ceiling[5/2 * order + 2*loop]` eps samples (must be ≤ 100).
- `eps0 = Power[10, -1/2*loop - goal/(order+1)]` then `Rationalize[N[eps0, MachinePrecision], eps0/100]`.
- `epslist = eps0 + Range[number] * eps0/100`.
- `singlepre = Max[Ceiling[(number + 2*loop)*(1/2*loop + goal/(order+1))], 30]`.
- `workpre = 2 * singlepre`, `xorder = 4 * singlepre`.
- Returns `{epslist, workpre, xorder}`.

`FitEps[epslist, vlist, leading]` (1327) — `Chop[Fit[data, eps^(leading + Range[0, number-1]), eps], 10^-$ChopPre]`.

`SolveIntegrals[jints, goal, order]` (1333):
- If `eps` is set in `Numeric`: 1-sample mode. `epslist = {eps_user + (4-D0)/2}`; `workpre = 2*goal`, `xorder = 4*goal`. Run `BlackBoxAMFlow`, return numbers.
- Otherwise: `{epslist0, workpre, xorder} = GenerateNumericalConfig[goal, order]`, `epslist = epslist0 + (4-D0)/2` (the **internal** eps grid). Run `BlackBoxAMFlow`. `leading = -2 Length[Loop]`. For each integral, fit Laurent `FitEps[epslist0, values, leading]` against `epslist0` (the **user-facing** grid), then `Series[..., {eps, 0, order + leading}]` truncated to `goal` digits.

## Cross-reference into the C++ port

- `BlackBoxAMFlow split` ↔ `cli::run_solve_integrals_json` part loop;
  per-part workdir `root/part_<i>` is created inside
  `solve_integrals_runner` to avoid Kira-cache collisions when one
  invocation covers multiple top-sector groups.
- `BlackBoxAMFlowSingle` ↔ `pipeline::black_box_amflow_single`.
- Cutkosky cleared globals (`:1031`) ↔ `pipeline::without_cut_or_prescription`.
- SingleMass child scheme inheritance (`:1056-1057`) ↔ `make_single_mass_child_requests` (full ending-scheme list).
- BlackBoxReduce dot/rank floor (`:1226-1227`) ↔ `ibp::black_box_reduce` `JDot/JRank` floor.
- D0 internal eps shift (`:1342, :1351`) ↔ `solve_ordinary_integrals_laurent` shift logic.
