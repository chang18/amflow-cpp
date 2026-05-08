# upstream `ibp_interface/Kira/interface.m` — symbol map

(Line numbers track <https://gitlab.com/multiloop-pku/amflow> at the
time of the v1.0.0 release.)

520 lines. Implements the four contract functions declared by AMFlow.m
(`SetReducerOptions`, `IBPSystem`, `AnalyticReduction`, `DifferentialEquation`)
against the Kira binary.

## Globals & WSL plumbing (10-58)

- `$KiraExecutable`, `$FermatExecutable` — executable paths (set by `install.m`).
- `wslQ / wsl / pathChange / fileNameJoin / kiraExecutable / fermatExecutable / runCommand` (18-36) — WSL-on-Windows adapter; `RunCommand` shells out via `wsl` if needed and exports `FERMATPATH`.
- Pulls `Family / Loop / IndepLeg / Propagator / Cut / NThread` from `AMFlowInfo` shorthand.
- `$AuxLeg` (61): a synthetic auxiliary leg `<Family>AuxLeg` symbol used for momentum conservation.
- `Leg / Conservation` (62-63): wrap `IndepLeg` + the synthetic auxiliary leg, and add `$AuxLeg -> -Total[IndepLeg]`. So Kira's input always has a closed conservation system.

## Numeric splits (43-79)

- `SPToSTU` (43): `complete = independent SP between IndepLeg^2 entries → invariant rules` derived from `AMFlowInfo["Replacement"]` after Conservation. Returns `{p1*p2 -> s, p2^2 -> 0, ...}` style rules in invariant variables.
- `IBPRule := If[!ComplexMode, Numeric, Select[Numeric, Im[#[[2]]] == 0&]]` (56) — **real-only numeric values go into Kira**.
- `CompensateRule := Complement[Numeric, IBPRule]` (57) — complex/imaginary numeric values; applied **after** Kira reduction, in `AnalyticReduction` (475) and `DifferentialEquation` (504).
- `STUToSP` (66): inverse of `SPToSTU`, maximal-rank subset.
- `Momentum / Mass` (75-76): `ToSquareAll[Propagator] /. IBPRule` — propagator written as `momentum^2 + mass`.
- `MassScale` (79): kinematic invariants + bare scalar masses appearing in `Momentum / Mass`.
- `ep = Symbol["Global`eps"]` (83) — local alias.

## Config files (90-169)

- `Config[dir]` (93): writes `config/integralfamilies.yaml` + `config/kinematics.yaml`.
  - `top_level_sectors` is a single integer `Σ 2^(i-1)` for positions where `TopSector[i] === _`.
  - `propagators: [["momentum", -mass], ...]` (sign flip!).
  - `cut_propagators` (cut indices).
  - `kinematics`: each invariant declared with mass dimension 2.
  - `scalarproduct_rules`: from `SPToSTU`, two cases: `Power[p, 2]` → `[[p, p], value]`, `Times[p1, p2]` → `[[p1, p2], value]`.
  - Optional `permutation_option: <perm>` if `$PermutationOption =!= None`.
- `Preferred[preferred, dir]` (163): writes the `preferred` file in Kira format, one integral per line, blank-line separated, `j[Family, ...]` rewritten to `Family[...]`.

## Job files (176-264)

- `ReductionJob[dir]` (176): writes `jobs.yaml` for one of `$ReductionMode ∈ {Masters, Kira, FireFly, Mixed, NoFactorScan}`.
  - **Masters** mode: `select_mandatory_recursively`, `run_initiate: masters`. No output (no `kira2math`).
  - All others: `select_mandatory_list: [[fam, target]]`, plus a per-mode reduction block, plus `kira2math` output emission.
  - Top sector: `Σ 2^(i-1)` over `_`-positions of `TopSector`.
  - `r = Length[TopSector] - Count[TopSector, 0] + IBPDot` (number of nonzero positions in top + dot).
  - `s = IBPRank`. `d = IBPDot`.
  - `integral_ordering = $IntegralOrder` (default 5).
- `Target[target, dir]` (256): writes the `target` file (`Family[...]` per line, blank-line separated).

## Run command (271-275)

- `RunKira[]` (271): builds `[kira, "-pNThread", "jobs.yaml", "-s<inv>=<value>" ...]`.
  - Each invariant in `IBPRule ∩ (MassScale ∪ {ep})` becomes a `-s` arg.
  - `eps` is special: `-sd=4-2*<value>` (substituting `d` instead of `eps`).
  - **All other complex `Numeric` values are NOT passed**; they go into `CompensateRule` and are substituted after Kira returns.

## Symbolic differentiation (282-405)

This block is "compute derivatives" code originally from LiteIBP.m by T. Peraro:

- `mp / mm / gD / g4 / LContract` (286-323): scalar-product algebra over abstract momenta `mm[p]`, with metrics `gD` (D-dim) and `g4` (4-dim).
- `MomDerivative[expr, mm[q][mu], dim, gmetric, finalcontraction]` (327): symbolic ∂/∂q^μ acting on a polynomial in `mp[*, *]`. Extension `MomDerivative[expr, mm[q][mu]] := MomDerivative[expr, mm[q][mu], MetricD, gD, LContract]` (334) — D-dim default.
- `LIBPLoopMomenta / LIBPIndepExtMomenta / SPRule / Props / LIBPIds / LIBPInvariants / LIBPDenoms / LIBPSps / LIBPSpsToJ / LIBPGetDerivatives / LIBPDerivivative / LIBPComputeDerivatives / LIBPDeriv / ComputeDerivative` (342-405): the LiteIBP path that:
  1. For a topology (Family), builds the SP↔invariant scalar product map from `SPToSTU`.
  2. Builds `LIBPDenoms[fam] = {j[fam, -e_i] -> Together[Props[i]]}` — i.e. expresses raising one denominator power in terms of scalar products.
  3. `LIBPGetDerivatives[topo]` solves the linear system that turns derivatives w.r.t. external momenta into derivatives w.r.t. the kinematic invariants.
  4. `LIBPDeriv[j[t, a__], s]` then differentiates one J integral w.r.t. invariant `s` symbolically, applying the chain rule on each propagator and reducing back to J integrals via `LIBPSpsToJ`.
  5. `ComputeDerivative[target, x] = Collect[Expand[LIBPDeriv[#, x]& /@ target] //. j[a, b__] j[a, c__] -> j[a, b+c]]`.

## SetReducerOptions (419-426)

`Options[SetReducerOptions]`:
- `IntegralOrder -> 5`
- `ReductionMode -> "Kira"` (one of Kira/FireFly/Mixed/NoFactorScan/Masters; the AMFlow.m caller forces Masters in IBPSystem then a real mode in AnalyticReduction).
- `PermutationOption -> None`.

## The four contract functions (429-505)

### `IBPSystem[top, rank, dot, preferred, complexmode, dir]` (429)

1. Forces `$ReductionMode = "Masters"`.
2. `CheckCompleteness[]` (defined back in AMFlow.m).
3. `$ReductionDirectory = dir`. Wipe and recreate `dir`.
4. Set `TopSector = top`, `IBPRank = rank`, `IBPDot = dot`, `ComplexMode = complexmode`.
5. `Config[dir]`, `ReductionJob[dir]`, `Preferred[preferred, dir]`.
6. Run Kira via `runCommand[RunKira[], ProcessDirectory -> dir, ProcessEnvironment -> ... + FERMATPATH]`.
7. Parse `dir/results/<Family>/masters` (one integral per line) → `j[Family, ...]` form.
8. Cache to `dir/results/masters_mma`.

### `AnalyticReduction[target]` (454)

1. `Target[target, dir]` writes the new target file.
2. `ReductionJob[dir]` rewrites jobs.yaml in the current `$ReductionMode`.
3. Wipe `firefly_saves`, `tmp` to avoid stale state.
4. Run Kira.
5. `masters = GetFile[results/masters_mma]`. If empty → `{{}, {}}`.
6. `rule = GetFile[results/<Family>/kira_target.m]` — the kira2math output. If missing → identity rules `masters -> e_i`.
7. Drop `rule[[2]] === 0` entries.
8. `{col, mat} = CoefficientArrays[Values[rule], masters]`. `col` must be all zero (otherwise wrong reduction). `mat` is the coefficient matrix.
9. **Substitute `d -> 4 - 2*ep` in the matrix** to convert Kira's `d` variable to the `eps`.
10. **Apply `CompensateRule`** to the matrix, threading the imaginary numeric values that bypassed Kira.
11. Append `IdentityMatrix[Length@masters]` so masters are also returned with identity rows. Returns `{masters, Thread[Join[Keys[rule], masters] -> mat]}`.

### `DifferentialEquation[vars]` (480)

1. Read `masters` from `results/masters_mma`.
2. `LIBPComputeDerivatives[Family]` (only once per family).
3. `der = ComputeDerivative[masters, var] /. IBPRule` per `var`.
4. `integrals = Cases[der, _j, Infinity]` — all J integrals appearing in derivatives.
5. `{masnew, rules} = AnalyticReduction[integrals]`. **Must equal `masters`** (else error).
6. Build `red[jint] = rules[jint]` lookup (default 0 row), substitute into `der` to get `diffeq[v][i][j]`.
7. **Sort masters**: by `SortIntegrals` (weight = `{-JProp, -JDot, -JRank, ...}`), then re-gather by `JSector` and reverse — so masters are ordered by sector (smallest first within sector) and the diffeq matrix becomes lower-triangular.
8. Permute rows/cols of each `diffeq[v]` accordingly.
9. **Apply `CompensateRule`** to the final matrix.
10. Returns `{sortmasters, vars, diffeq}`.

## Cross-reference

- `IBPSystem[top, rank, dot, preferred, complexmode, dir]` ↔ `ibp::black_box_reduce` first phase + `ibp::run_kira` (Masters mode), with `KiraConfig` carrying `dimension_base` for `d -> 4 - 2*eps`.
- `AnalyticReduction[target]` ↔ `ibp::black_box_reduce` second phase: writes target file, runs Kira in Kira/FireFly mode, reads masters/target, parses coefficients via `parse_kira_expression`. The the C++ parser already lifts `d -> 4 - 2*eps` through `dimension_base`. **CompensateRule** equivalent: passes the full numeric kinematics through `parse_kira_expression` against a context whose vars include the complex-valued ones too — so the imaginary part substitutes after the parse, equivalent to the upstream's "apply CompensateRule" step.
- `DifferentialEquation[vars]` ↔ `ibp::black_box_diffeq`, composing `ibp::libp_deriv` (the LiteIBP-style chain rule) with `ibp::black_box_reduce`.
- `LIBPDeriv` ↔ `ibp::libp_deriv` — currently restricted to differentiation w.r.t. **scalar invariants** (rejects loop / leg derivatives). the upstream's `LIBPDeriv` only supports `D[expr, s]` as a special case via `LIBPDeriv[expr_, s_] := D[expr, s]` (402); the differential-equation use case in DESolver.m only uses invariant-side derivatives + the AMFlow eta variable, so this matches the restriction.
- `SortIntegrals` master ordering ↔ the IBP master sorting (need to verify in `ibp::black_box_diffeq` / `ibp::black_box_reduce` that the same `JSector`-grouped reverse-weight order is used).

## Surprising things

- The auxiliary leg `<Family>AuxLeg` is a per-family fresh Symbol — so two families never share that symbol name. This C++ port must replicate this when writing `kinematics.yaml` (or use any other closure trick that gives Kira a complete momentum-conservation set).
- `d -> 4 - 2*ep` substitution is applied in `AnalyticReduction[target]` (475) after parsing the kira_target.m output. `IBPSystem` itself runs Masters mode and never sees `eps` at all — Kira's `d` is just whatever `dimension_base` it was given via `-sd=`. This C++ port's `KiraConfig::dimension_base` mirrors this exactly.
- `CompensateRule` is non-trivial: any `Numeric` entry whose value has a nonzero imaginary part is **not** passed to Kira, but is substituted in the coefficient table and DE matrix afterwards. This C++ port handles this implicitly because `parse_kira_expression` evaluates against the full numeric kinematics, including the complex ones.
- DifferentialEquation **re-runs `AnalyticReduction`** internally for the integrals that appear in the symbolic derivatives. Then it sorts and permutes the master basis. This C++ port's `ibp::black_box_diffeq` must reproduce that same sort to keep the master order consistent with subsequent boundary planning.
- `IBPSystem` and `AnalyticReduction` share state via the global `$ReductionDirectory`. This C++ port passes `dir` explicitly to every entry to remove that hidden coupling.
