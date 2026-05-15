# upstream `diffeq_solver/DESolver.m` — symbol map

(Line numbers track <https://gitlab.com/multiloop-pku/amflow> at the
time of the v1.0.0 release.)

1149 lines. Outer structure:

- 7-57: `BeginPackage` + usage declarations.
- 59: `Begin["`Private`"]`.
- 62-129: option setters (Global / Expansion / Running) + defaults.
- 132-187: helpers (FitEps, EvaluateAsymptoticExpansion, PowerSeries / RationalExpansion).
- 189-304: equation analysis (AnalyzeBlock / AsymptoticBehavior / NHEquations / NHEquationsNum / PoincareRank).
- 308-446: diagonal normalization (Balance / FindProjector / ToFuchsian / NormalizeEigen / NormalizeDiagonal).
- 449-490: off-diagonal normalization (SolveOffDiagonal / ToFuchsianGlobal).
- 493-515: NormalizeMat (composition).
- 518-587: integration contour (GetPoles / RunUnit / RunSegment / RunEta).
- 590-680: sparse solve (Sparsify / ConstructMatrix / SparseGaussian).
- 683-731: boundary conditions (BuildTaylor / DetermineBoundaryOrder / ReverseBCS / UnionBCS / ReadBCS).
- 734-795: expand at infinity (CalcTaylor / CalcInf).
- 798-849: regular point (Calcx1x2 / CalcRun).
- 852-979: zero-boundary (Calcx00 / FindLogPower / LearnFromRuleS / ExtendExpansion / PlusExpansion / PSMapRuleS / ToPS / CalcZero / TimesAsyExp / PlusAsyExp).
- 982-1080: system handle (LoadSystem / ClearSystem / InfToRegular / RegularRun / RegularInterpolation / SolveAsyExp).
- 1083-1139: `AMFlow[de, bc]` — top-level orchestrator.

## Globals & options

- `eta` is the auxiliary squared mass.
- `Options[SetGlobalOptions]` (66): `WorkingPre -> 100`, `ChopPre -> 20`, `SilentMode -> False`, `RationalizePre -> 20`. Setter calls `SetWorkingPrecision[WorkingPre]`.
- `Options[SetExpansionOptions]` (82): `XOrder -> 100`, `ExtraXOrder -> 20`, `LearnXOrder -> -1`, `TestXOrder -> 5`.
- `Options[SetRunningOptions]` (97): `RunRadius -> 2`, `RunLength -> 200`, `RunCandidate -> 10`, `RunDirection -> "NegIm"`.
- `SetWorkingPrecision[p]` (112): clamps `$MinPrecision = $MaxPrecision = p`.
- `AMFN[a]` (115) `:= N[a, WorkingPre]`. `AMFChop[a]` (116) `:= Chop[a, 10^-ChopPre]`. `AMFRationalize` (118) `:= Rationalize[a, 10^-RationalizePre]`.

## Helpers (143-187)

- `MatrixDensity / MatrixDegree`.
- `FitEps[rule, leading]` (154) — local DESolver Laurent fit (different signature from AMFlow.m's `FitEps[epslist, vlist, leading]`).
- `EvaluateExpansion[exp, x0]` (161): `exp` is a 2-D table indexed `[logpower+1, etapower+1]`, evaluates at given `x0` as `Σ x0^etapower * Log[x0]^logpower`.
- `EvaluateAsymptoticExpansion[mu -> exp, x0]` (163): `x0^mu * EvaluateExpansion`.
- `PickElement / PickList / PickMat` (166-168): `list[[order+1]]` with default 0.
- `TimesPoly[f, poly]` (173): truncated × precise polynomial multiplication.
- `InversePoly[f, poly]` (176): truncated 1/poly given leading nonzero.
- `RationalExpansion[num, de, expansion] = InversePoly[TimesPoly[expansion, num], de]`.
- `MapRationalExpansion[nummap, demap, expansions]` (183): per-row `InnerRationalExpansion`.

## Block analysis (193-304)

- `AnalyzeBlockOld[mat]` (193) — banded-form analysis (right-most nonzero column running max).
- `AnalyzeBlock[mat]` (204): `detest = Boole[#=!=0]&` then `AnalyzeBlock0`. Computes nonzero-position closure per row, gathers + sorts by dependency. Throws `"AnalyzeBlock: bad blocks."` if blocks don't cover all rows.
- `SubBlockID[mat]` (240): for each block i, the list of strictly-lower blocks j that mat couples into.
- `UnionBehavior[behs]` (247): `{(eta-power, max log-power)}` per distinct eta-power.
- `UnionBehavior2[behs]` (248): same key, but `Total[log-powers] + 1` if multiple.
- `AsymptoticBehavior[mat]` (251): per block, derive essential behavior from leading Jordan structure of `mat[block, block]*eta /. eta -> 0`, union with subblock behavior, store recursively. **Demands input matrix is already in Fuchsian + Jordan-leading form**.
- `PoincareRank[mat]` (266): `Max[Exponent[mat // Denominator, eta, Min]] - 1`.
- `NHEquations[mat, mode]` (271): for each block, factor out an η^k front:
  - `Singular`: `factor = eta` (for `x dx D[f,x] == A f + B g`).
  - `Regular`: `factor = 1`.
  - `Taylor`: `factor = eta^(PoincareRank+1)`.
  Returns `Transpose@{dxlist, axlist, bxlist, blocks, sub}` where `dx = LCM[denominators of factor*mat[block,block]]`, `ax = factor*dx*mat[block,block] // Together`, `bx = factor*dx*mat[block,sub] // Together`. For non-Singular modes, `dx` is multiplied by `factor` after construction.
- `NHEquationsNum[nheqs]` (295): `ToNum` each polynomial (`CoefficientList[poly, eta] // N`).

## Normalization (308-515)

- `PBar[P] = I - P`. `Balance[P] = (I-P) + P/eta`. `InvBalance[P] = (I-P) + P*eta`.
- `ReduceL0[L0, r, lblock]` (317): off-diagonal cleanup helper inside FindProjector.
- `FindProjector[cp, cp1]` (345): for irreducible Poincaré rank > 0 case, builds projector Q from Jordan blocks of `cp = mat*eta^(p+1) /. eta -> 0` and `cp1 = D[mat,eta]*eta^(p+1) + (p+1)*eta^p*mat /. eta -> 0`. Errors out if any Jordan diagonal is nonzero (irreducible Poincaré rank).
- `ToFuchsian[mat]` (372): single-block Fuchsian reduction. Iterates `Balance[Q]` until PoincaréRank == 0.
- `NormalEigen[feps] = {feps - Floor[Re[feps]], Floor[Re[feps]]}` (389).
- `NormalizeEigenQ[mat]` (393): True if any leading-residue eigenvalue has `Floor[Re] != 0`.
- `ShearingTransformation[mat]` (397): build T = Jordan-eigenvectors scaled by `eta` (or `1/eta`) per eigenvalue floor, apply.
- `NormalizeEigen[mat]` (410): loop ShearingTransformation until eigenvalues in `[0, 1)` strip.
- `NormalizeDiagonal[mat]` (422): for each block:
  1. ToFuchsian.
  2. NormalizeEigen.
  3. Constant Jordan rotation on leading residue.
  Then update lower off-diagonal `B[rows,cols] = invT[rows,rows] . mat[rows,cols] . T[cols,cols]` — block-diagonal first-step transform.
- `SolveOffDiagonal[a0, b0, c0, p]` (454): solve Sylvester `b0 + c0.G - G.a0 + p*G == 0`.
- `ToFuchsianGlobal[mat]` (464): off-diagonal cleanup. For each (i, j) with i > j, in reverse j order, while PoincaréRank > 0 of off-diagonal block: `transTB[i, j, G]` adjusts `B[i,j]`, `T[*,j]`, `invT[i,*]`.
- `NormalizeMat[mat]` (497): `NormalizeDiagonal` → `ToFuchsianGlobal`, then compose `T2[rows,cols] = T[rows,rows] . T2[rows,cols]` and `invT2[rows,cols] = invT2[rows,cols] . invT[cols,cols]` for the lower-triangular slots. Returns `{T2, invT2, B}`.

## Contour generation (522-587)

- `AllFactors[exp]` (522): factor list of joined denominators.
- `Zeros[poly]` (523): `NSolve[poly == 0, eta, WorkingPrecision -> WorkingPre]`.
- `GetPoles[mat]` (524): rationalized roots of all denominator factors.
- `FirstStep[polelist]` (528): `RunRadius * Max[|polelist|]` (or 1 if only pole is 0).
- `LastStep[polelist]` (532): `Min[|nonzero polelist|] / RunRadius`.
- `RunUnit[polelist]` (536): generate path `[0, 1]` avoiding integer-real poles. Step size `min[point]/RunRadius`. Bails if `RunLength` exceeded.
- `RunSegment[polelist, ini, fin]` (555): `RunUnit` rescaled to `[ini, fin]`.
- `RunEtaDirection[polelist, 1]` (558): real-axis run from `FirstStep` to `LastStep`.
- `RunEtaDirection[polelist, direction_]` (563): rotated path along the unit-direction.
- `RunEtaDirection[polelist, mode_String]` (568): for `Re` / `Im` / `NegRe` / `NegIm`: try `RunCandidate * 2 + 1` directions slightly tilted, pick the shortest non-empty.
- `RunEta[polelist]` := `RunEtaDirection[polelist, RunDirection]`.

## Sparse solve (593-680)

- `Sparsify[matrix]` (593): `{dim, sp[matrix]}` where each row stores `{col -> value}` only for nonzero entries.
- `SparseTranslate[sp, n]` (601): shift column indices by `n`.
- `ConstructMatrix[dx, ax, totalorder]` (607): build the truncated linear system for `dx D[f,x] - ax f = g`. `f` has `(totalorder + 2) * Length[ax]` unknowns ordered by `(power, integral)`.
- `SplitSystem[sys, m]` (626): split rows by whether they reference variable index `m`.
- `PlusSparse / ScalarSparse` (634-635).
- `ForwardSparseGaussian[sparse, vset]` (638): Gauss-eliminate over variable list `vset` in order, returning `{rset, residueSystem}`.
- `SparseGaussian[{dim, sparselist}, nh]` (660): block-walking variant. Augments `nh` as the trailing column. Iterates `ForwardSparseGaussian` per block, threading residue forward.

## Boundary order & boundary input (687-731)

- `BuildTaylor[mat, ini]` (687): `mat -> mat'[i,j] = eta^-ini[i] * mat[i,j] * eta^ini[j]` for off-diagonal, and `mat'[i,i] = mat[i,i] - ini[i]/eta` for diagonal. This converts a singular system at `eta = 0` with prescribed leading powers `ini` into a Taylor system on the residual.
- `DetermineBlockBoundaryOrder[mat, power]` (692): per block, build Taylor matrix, solve with `XOrder=ExtraXOrder`. Find unsolved coefficients (free Taylor parameters), pick the order at which there are no unsolved variables, return `{order_per_integral}` (or -1 if none).
- `DetermineBoundaryOrder[mat, power]` (720): apply `DetermineBlockBoundaryOrder` per block and merge by integrating row index back.
- `ReverseBCS[bcs] = Thread[-Keys[bcs] -> Values[bcs]]` (723): flip eta-power signs.
- `UnionBCS[bcs]` (724): sum per distinct eta-power (after the integer-shift dedup done by ReadBCS).
- `ReadBCS[bcs, region]` (727): for an integer-shift equivalence class `region`, find min eta-power, return `{ini, normalized bcs starting at 0}`.

## CalcTaylor / CalcInf (739-795)

- `CalcTaylor[mat, bc]` (739): integrand of `CalcInf`. Builds NHEquations[mat, "Taylor"], constructs matrix to total order `XOrder + ExtraXOrder`, sparse-solves. Inserts boundary conditions for unsolved coefficients (free Taylor parameters, indexed by `(integral, order)`). Returns `{f0[i, 0..XOrder]}` per integral.
- `CalcInf[de, bcs]` (782):
  1. `deinf = -1/eta^2 * (de /. eta -> 1/eta)`.
  2. `bcinf = UnionBCS /@ ReverseBCS /@ bcs`.
  3. `regions = NormalEigen[Keys[bcinf]][[1]] // Flatten // DeleteDuplicates` — distinct eta-power equivalence classes.
  4. For each region: `{ini, bc} = Transpose[ReadBCS[#, region] & /@ bcinf]`, `mat = BuildTaylor[deinf, ini]`, `rule = Thread[ini -> CalcTaylor[mat, bc]]`, accumulate.
  Returns `allrule[[i]] = {ini -> Taylor coeffs, ...}` per integral — i.e. an `AsyExpansion` keyed by leading eta-power.

## Regular point (803-849)

- `ExpandNHEquationsNum[nheqn, x0]` (804): Taylor-shift each polynomial via `Binomial[i,j] x0^(i-j)` table.
- `Calcx1x2[nheqn, bc, x0]` (816): one-segment forward Taylor recurrence. Closed form for `f0[block[[i]], n+1]` from `f0[*, 0..n]`. Refuses if `dx[0] == 0`.
- `CalcRun[de, bc, run]` (842): multi-segment chain. Per step `i`: `Calcx1x2`, then evaluate the resulting Taylor at `run[[i+1]] - run[[i]]` to get bc for next step.

## Zero boundary (856-979)

- `Calcx00[nheq, nheqn, bc, x0, behavior]` (856) — **the largest function in DESolver.m, the one split into `amflow/ode/zero.cpp`'s 4-phase pipeline**.

  Per block k (with own behavior list):
  1. `a00 = (axlist[k] / dxlist[k] /. eta -> 0)` — leading residue.
  2. For each behavior (`spbeh`, `logk`):
     - Compute the source `nh = MapRationalExpansion[bxexpn, bxexpd, ...]` with f0 of subblocks.
     - For each `n = 0 .. XOrder`:
       - `posi = Position[Together[(n + spbeh) - Diagonal[a00]], 0]` — resonance positions.
       - If `posi != {}`: **essential** behavior. Init `kk[spbeh, 0, n] = I`, `bb[spbeh, 0, n] = 0`. Use **rec1** to fill `bb / kk[spbeh, p+1, n]` for p ascending. Record `pair[m]` for the resonance constraints.
       - Else: invertible. Use **rec2** with `inv = ((n + spbeh) I - a00)^-1 / dx[0]`. Fill `bb / kk[spbeh, p, n]` for `p = logk, logk-1, ..., 0` descending.
  3. Match equations: at `x0`, evaluate `Σ_m essentialset[m] kk[m] . variables[m] == bc[block] - Σ bb[m]`.
  4. Pair constraints: from `pair[m]`, equate `kk[m, logk+1, n] . variables[m] + bb[m, logk+1, n]` to 0 on the non-resonant rows.
  5. Solve `compensate = Length[allvar] - Length[pairtoequ]` matchequ rows + all pairtoequ. Drop residual matchequ rows.
  6. Substitute back to obtain `f0[block[[i]], spbeh, p, n]` for each behavior, log power, eta power.

  Returns: `Range[Length[bc]] /. Thread[Join@@blocks -> Join@@final]` — per integral, list of `{spbeh -> {p_table_of_eta_coefficients}}`.

  `rec1` (858): the standard rec for non-invertible block — resolves `func[spbeh, p, n]` at `(p, n)` using lower (p, n).
  `rec2` (863): inverts `(n + spbeh) I - a00` to solve for `func[spbeh, p, n]`.

- `FindLogPower[exp]` (920): largest log power `k` with `exp[k+1] != 0`.
- `LearnFromRuleS[rules]` (921): for each rule, drop log levels that are all-zero up to `Min[TestXOrder, LearnXOrder]+1`. Used in the learning phase.
- `LearnFromRuleSAll[ruleslist, blocks]` (922): per block, union learned behaviors.

- `ExtendExpansion[exp, n]` (928): pad to `n` log rows.
- `PlusExpansion[exps]` (929): sum a list of expansions, padding to common max log row count.
- `RescaleExpansion[exp, order]` (930): shift eta-rows by `order` (i.e. multiply expansion by `eta^order`, truncate).

- `PlusRuleS[rules0]` (935): for a list of rules over the same eta-power equivalence class, take the min eta-power, rescale and add.
- `UnionRuleS[rules]` (936): `PlusRuleS` per integer-shift equivalence class.
- `PSTimesRuleS[ps, rules]` (938): "PowerSeries × RuleS" — multiply a (poly, exp) pair into a rule list.
- `PSMapRuleS[psmap, ruleslist]` (939): given a matrix-of-PowerSeries map (e.g. `T`), apply to a list of asymp expansions (one per integral).

- `ToPS[poly]` (943): convert a Laurent poly in eta to `{coeffs, leading_eta_exponent}`.

- `CalcZero[de, bc, x0]` (954) — **top-level zero entry**:
  1. `{T, invT, B} = NormalizeMat[de]`.
  2. `bcT = invT(x0) . bc`.
  3. `nheq = NHEquations[B, "Singular"]`. `nheqn = NHEquationsNum[nheq]`.
  4. `behavior = AsymptoticBehavior[B]`.
  5. If `LearnXOrder >= 0`: redo step 6 with `XOrder = LearnXOrder` to learn behavior, then `LearnFromRuleSAll`.
  6. `PSMapRuleS[ToPS[T], Calcx00[nheq, nheqn, bcT, x0, behavior]]` — push back through `T`.

- `TimesAsyExp[rational, asyexp]` (969): multiply rational × asymp.
- `PlusAsyExp[asyexplist]` (978) `:= UnionRuleS[Join @@ asyexplist]`.

## System handle (985-1080) — the mutable side that this C++ port explicitly removes

- `LoadSystem[sysid, de, bc, point]` (989) — sets `DE[sysid] = de; BC[sysid] = bc; P[sysid] = point; AsyExp[sysid] = Null`.
- `ClearSystem[sysid]` (1010).
- `InfToRegular[sysid, x0]` (1025): require `P[sysid] == Infinity`, set `AsyExp = CalcInf[DE, BC]`, `BC = Σ_eta_power EvaluateAsymptoticExpansion[asyexp, AMFN[1/x0]]`, `P = x0`.
- `RegularRun[sysid, run]` (1034): `BC = CalcRun[DE, BC, AMFN[run]]`, `P = run[[-1]]`.
- `RegularInterpolation[sysid, samples]` (1042): for each requested sample, either single-Taylor reuse (`Calcx1x2`) if within `radius`, else `RunSegment + RegularRun + Calcx1x2` chain. Restores caller's `BC, P` afterward. Returns the value list.
- `SolveAsyExp[sysid]` (1078) := `AsyExp[sysid] = CalcZero[DE, BC, AMFN[P]]`.

## Top-level orchestrator (1091-1139)

- `PickZeroRuleS[rules]` (1091): from a per-integral `AsyExp`, pick the eta-power-0 rule. Validate that exactly one integer-eta-power region survives. If `key > 0` → return 0. If `key == 0` → return `(key /. rules)[[1, 1]]` (constant log-row, constant eta-row coefficient). If `key < 0` → drop the negative coefficients with `Print["dropped terms when pickzero"]`, return `[[1, -key+1]]`.
- `AMFlow[de, bc]` (1103) — **the entry called by AMFlow.m's AMFSystemSolution template**:
  1. `LoadSystem[$InternalSystem, de, bc, Infinity]`.
  2. `run = RunEta[GetPoles[DE]]`. If empty → abort.
  3. `InfToRegular[$InternalSystem, run[[1]]]`.
  4. `RegularRun[$InternalSystem, run]`.
  5. `SolveAsyExp[$InternalSystem]`.
  6. `sol = PickZeroRuleS /@ AsyExp[$InternalSystem]`.
  7. `ClearSystem[$InternalSystem]`. Return `sol`.

## Cross-reference

- `AMFlow[de, bc]` ↔ `ode::solve_ode(de, boundaries, options)` (`amflow/ode/solve_ode.hpp`).
- `GetPoles` ↔ `ode::singular_points`.
- `RunEta / RunSegment / RunEtaDirection` ↔ `ode::PathPlan / make_path_plan / first_step / last_step / run_unit / run_segment / run_eta_direction / run_eta`.
- `CalcInf` ↔ `ode::solve_at_infinity` (composes `make_infinity_variable_matrix / boundary_regions / read_boundary_region / build_taylor_matrix / calc_taylor`).
- `CalcRun` ↔ `ode::calc_run`.
- `Calcx1x2` ↔ `ode::calcx1x2`.
- `RegularInterpolation` ↔ `ode::regular_interpolation`.
- `CalcZero` ↔ `ode::solve_at_zero` (with the `normalize_mat_for_calc_zero` numeric-Jordan-rotation fallback for irrational eigenvalues).
- `Calcx00` ↔ `ode::build_zero_recurrence_block + solve_zero_boundary_block + solve_zero_blocks` (4-phase split).
- `NormalizeMat` ↔ `ode::normalize_mat` (rational path) + `normalize_mat_for_calc_zero` (algebraic path).
- `LearnFromRuleS / FindLogPower` ↔ `ode::learn_from_rule_s / learn_from_rule_s_all`.
- `PickZeroRuleS` ↔ `ode::pick_zero_solution`.
- `LoadSystem / ClearSystem / DE / BC / P / AsyExp` ↔ **intentionally absent** in this port — replaced by value-typed `SolveOdeRequest / SolveOdeResult` so there is no mutable system handle.
- `DetermineBoundaryOrder` ↔ `ode::determine_boundary_order` (the sub-leading boundary-order chain consumed by `CalcInf`'s region table).
- `BuildTaylor` ↔ `ode::build_taylor_matrix`.
- `CalcTaylor` ↔ `ode::calc_taylor`.
- `ConstructMatrix / SparseGaussian / ForwardSparseGaussian` ↔ `ode::SparseSystem / construct_matrix / forward_sparse_gaussian / sparse_gaussian`.

## Open / surprising things

- `Calcx00` Solve has `compensate = Length[allvar] - Length[pairtoequ]` matchequ rows kept and the rest dropped. The number of matchequ rows is `Length[block]`, and `pairtoequ` is the number of resonant rows. So we must drop `Length[block] - compensate = pairtoequ.length` rows to make the system square. **The C++ `solve_boundary_linear_system` uses Dixon over Q[i]; upstream uses `Solve[..., allvar]` symbolically.  There is no explicit "rationalize via decimal then exact rational solve" pattern in upstream — the C++ port introduces it because Acb interval solves are too lossy for the actual matching system.**
- `AMFlow[de, bc]` always uses `$InternalSystem` as the sysid. So the system handle is in fact a global singleton inside `AMFlow[]` calls; the parallel sampling in `AMFSystemSolution` uses `ParallelTable` with `DistributedContexts -> All`, which evaluates `AMFlow` per kernel — each kernel has its own `$InternalSystem`. **The value-style `solve_ode` call is correct: it uses no shared mutable state.**
- DESolver.m's `FitEps[rule, leading]` (154) is **different** from AMFlow.m's `FitEps[epslist, vlist, leading]` (1327). The DESolver one is for internal asymptotic series fits; the AMFlow one is the user-facing eps Laurent fit. This C++ port's `solve_ordinary_integrals_laurent` uses the latter pattern.
