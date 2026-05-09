# MMA Parity Audit — v1.0 retrospective

A line-level audit of the C++17 port against upstream Mathematica
AMFlow (commit `efda1db` of <https://gitlab.com/multiloop-pku/amflow>).
The audit was performed after v1.0 release to surface divergences that
the 12 oracle benchmarks listed in [`AUDIT.md`](../AUDIT.md) do not
cover.

## TL;DR

| Severity | Count | Action taken |
|---|---|---|
| 🟢 verified                  | 65 | — |
| 🟡 unverified (oracle gap)    | 21 | Documented in §3 below; tracked for future bench expansion. |
| 🔴 actual divergence          |  6 | 2 fixed in this release; 4 documented for follow-up (§2). |
| ⚪ intentionally not ported   | 17 | — |

Net assessment: **no oracle-validated path is wrong**.  All 🔴 items are
either now fixed or latent (no oracle exercises them).  The audit is
preventative, not corrective — its main output is the 🟡 list,
which scopes the next round of benchmark expansion.

---

## 1. Methodology

For each upstream file, an audit pass classified every public symbol
into one of four severities:

- **🟢 verified** — semantics match upstream, and an oracle benchmark
  in [`AUDIT.md`](../AUDIT.md) §4 covers the path.
- **🟡 unverified** — semantics look correct but no oracle exercises
  the path.  Latent risk — kept in this list until an oracle is added.
- **🔴 divergent** — actual semantic difference between MMA and C++.
  Each one is a mini bug report.
- **⚪ not ported** — explicitly out of scope (gauge link / HQET /
  SCET / Wilson / non-Kira IBP backends / WSL plumbing).

Three audit passes ran in parallel:

| Audit | Upstream file | C++ targets |
|---|---|---|
| ODE engine  | `diffeq_solver/DESolver.m` (1162 lines) | `src/ode/`, `include/amflow/ode/` |
| Pipeline    | `AMFlow.m` (1533 lines) | `src/qft/`, `src/ibp/` (excluding `kira_*`), `src/pipeline/`, `src/api/` |
| IBP backend | `ibp_interface/Kira/interface.m` (532 lines) | `src/ibp/kira_*.cpp`, `include/amflow/ibp/kira.hpp` |

Full per-pass reports: `/tmp/audit_desolver.md`, `/tmp/audit_amflow.md`,
`/tmp/audit_kira.md` (locally generated; not committed).

---

## 2. 🔴 Divergences

### D1. Default `RunLength = 200` (upstream: 1000) — **FIXED in v1.0.1**

- **C++ (was)**: `include/amflow/numeric/options.hpp:63`,
  `RunningOptions::run_length = 200`.
- **Upstream**: `DESolver.m:97`, `RunLength -> 1000`; mirrored in
  `AMFlow.m:259`.
- **Symptom**: `RunUnit` aborts after 200 contour steps where MMA
  tolerates 1000.  False `"integration contour not generated"` errors
  on crowded pole landscapes.  No oracle currently trips this, but the
  cap was 5× tighter than upstream.
- **Fix**: brought to 1000 to match upstream (one-line edit).

### D2. Default `RationalizePre = 20` (upstream: 100) — **FIXED in v1.0.1**

- **C++ (was)**: `GlobalOptions::rationalize_pre = 20`.
- **Upstream**: `DESolver.m:66`, `RationalizePre -> 100`; mirrored in
  `AMFlow.m:259`.
- **Symptom**: every `acb_real_to_fmpq` call along the contour path
  (`src/ode/path.cpp:179, 216, 266, 354–355` and `src/ode/inf.cpp:798`)
  uses 20-digit rationalisation where MMA uses 100.  Silent narrowing
  of the rationalisation funnel.
- **Fix**: brought to 100 to match upstream.  Examples that explicitly
  set the old value (`examples/box1_*.json`, `examples/bubble_*.json`)
  had the redundant override removed.  All 12 committed bench `*_cpp.json`
  already set 80 or 100 explicitly, so they are unaffected.

### D3. Boundary sub-system loses parent's `Cut` info — **deferred**

- **Upstream** `ReduceBoundary` (`AMFlow.m:790–803`): when the parent
  system has a non-empty `Cut`, projects the parent's cut propagators
  (after the region transform) onto the new boundary-sub-family's
  propagator basis and emits `cut_propagators:[...]` to Kira.  Aborts
  on `Count[cut, 1] != Count[Cut, 1]`.
- **C++** `AMFSystem::build_boundary` (`src/pipeline/amfsystem.cpp:1452`):
  hard-codes `cut={}` when constructing the boundary `FamilyConfig`.
  The detection assertion is also missing.
- **Latent** because every Cutkosky benchmark
  (`cutbubble_1L`, `cutsunrise_2L`, `cutbanana_3L`,
  `tt_cutkosky_probe`) clears `cut` at the parent level via the
  Cutkosky scheme, never exercising the Tradition-with-cut path.
- **Recommended fix**: project parent's `cut` through `region.transform`
  in `build_boundary`, pass the result to `qft::FamilyConfig::build`,
  raise on count mismatch.

### D4. `BoundaryIntegrands` Jacobian factor `|Det|^(4-2eps)` is silently dropped — **deferred (defended-by-precondition)**

- **Upstream** (`AMFlow.m:731-732`): multiplies integrands by
  `Abs[Det[Coefficient[#, Loop]&/@Values[trans]]]^(4-2*eps)` (the
  Jacobian of the loop redefinition).
- **C++** `qft::boundary_integrands` (`src/qft/boundary.cpp`): no
  Jacobian factor.
- **Currently moot** because `branch_momenta`
  (`src/qft/region.cpp:108–114`) asserts each propagator has unit
  leading loop coefficient → `|Det| = ±1` → factor = 1.  Effectively
  🟢-by-precondition but fragile if the precondition is ever relaxed.
- **Recommended fix**: defensive assert that the computed Jacobian
  equals `±1`, or compute it explicitly and raise the precondition
  to "Jacobian must be unit".

### D5. Kira `ComplexMode` / `CompensateRule` real-only filter unimplemented — **deferred**

- **Upstream** filters `IBPRule` to drop imaginary-part numerics from
  the Kira CLI input, then post-substitutes them via `CompensateRule`
  after `AnalyticReduction` (line 488) and `DifferentialEquation`
  (line 516).
- **C++** `KiraConfig::numeric_values` is a flat `map<string,string>`
  of pre-computed rational strings — no imaginary handling.  The 5th
  `IBPSystem` parameter `complexmode` has no C++ counterpart.
- **Latent**: all 12 oracles use purely-real `numeric_values`.
- **Recommended fix**: add an imaginary-extraction step in
  `kira_run.cpp` before assembling the `-s` arguments, plus a post-Kira
  substitution in `kira_parse.cpp`.  Non-trivial.

### D6. RHS-not-in-master is silently dropped — **deferred**

- **Upstream** (`Kira/interface.m:486`): `Abort[]`s when
  `CoefficientArrays` reports a non-master residue (i.e., a J-integral
  appears on the right-hand side that is not in the masters list).
- **C++** `reduce.cpp:344-360`: the diffeq matrix builder `continue`s
  past unknown rhs J integrals.
- **Diagnostic regression only** — healthy Kira runs are unaffected.
- **Recommended fix**: throw `std::runtime_error` on unknown rhs
  integral with the family / indices of the offender.

---

## 3. 🟡 Unverified branches (oracle gaps)

These branches look correct to inspection but no committed benchmark
exercises them.  Each is a candidate for a future oracle.

| Branch | Site | Why it's unverified |
|---|---|---|
| `analyze_block` uses `AnalyzeBlock0` topology | `src/ode/blocks.cpp` | Upstream default is `AnalyzeBlock1`; equivalent on nested/equal blocks (the IBP common case), could differ on overlapping non-nested closures. |
| `Calcx00` boundary linear-system selection | `src/ode/zero.cpp:1314-1703` | Heuristic match-row selection + Dixon over Q[i]; upstream uses one symbolic `Solve[…, allvar]`.  Equivalent on full-rank systems. |
| Acb-inverse failure fallback in `Calcx00` | `src/ode/zero.cpp:1199-1206` | C++ degrades to "all positions resonate"; upstream's symbolic inverse never fails. |
| `evaluate_taylor` strips arb radii | `src/ode/regular.cpp` | Documented in source; abandons rigorous error bars on the regular-running stage. |
| Jordan block ordering | `src/ode/jordan.cpp` | Groups by eigenvalue then descending depth; MMA sorts globally by descending size.  Same final algebra, different column permutation when ≥2 distinct eigenvalues. |
| `SolveIntegrals` single-eps fast path | `src/pipeline/solve_integrals.cpp:867-944` | `AMFlow.m:1364-1374` short-circuits when `eps` is in `Numeric`; C++ always Laurent-fits. |
| Per-system `AMFSystemDirection` | `src/ode/path.cpp:431` | Upstream `AMFlow.m:981-991` picks Im / NegIm per-system based on prescription; C++ uses global `RunDirection`.  `AMFSystemOptions::direction` field is unread. |
| `SingleMassQ` substitutes `Numeric` | `src/pipeline/amfsystem.cpp:795` | C++ enhancement; MMA tests literally.  Latent: if a user supplies `mass` symbolically and forgets `Numeric`, MMA refuses to ending and C++ accepts. |
| `factorize_family` mass=−1 detection substitutes `Numeric` | `src/pipeline/amfsystem.cpp:2204` | Same enhancement as previous row. |
| Cutkosky physical-mass safety check | `src/pipeline/amfsystem.cpp:2669-2727` | Upstream `AMFlow.m:1050` aborts on negative `cutcom[[1,5]]` masses; C++ skips.  Misuse with negative masses → garbage instead of abort. |
| Auto-applied `Vacuum[L,n]` table | `src/pipeline/amfsystem.cpp:631-678` | C++ enhancement: looks up the closed-form table at ending.  MMA aborts unless the user supplies a `Solution`.  Closed-forms verified entry-by-entry. |
| Ending-master Kira reduction loop | `src/pipeline/amfsystem.cpp:710-785` | C++-only path: lowers an arbitrary ending master through repeated `BlackBoxReduce` calls. |
| `zero_sector_q` substitutes generic primes | `src/qft/topology.cpp:42-48` | Both are generic-point substitutions sufficient to detect scalelessness; benign. |
| `region_power` skips `/.Numeric` | `src/qft/findregion.cpp:463-529` | Expression contains only `eps` + integer constants; benign by construction. |
| `factorize_family` no-redef fallback | `src/pipeline/factorize.cpp:268-285` | Dead code on tested inputs (square / well-conditioned systems). |
| `LIBPDeriv` multi-invariant case | `src/ibp/libp_deriv.cpp` | Single-invariant tested; multi-invariant unverified. |
| `MasterRank`/`MasterDot` filter | `src/ibp/reduce.cpp` | Default `Infinity` is no-op. |
| `kira_target.m` parser strictness | `src/ibp/kira_parse.cpp` | Custom tokenizer is not fuzzed. |
| Coefficient parser fragility | `src/ibp/kira_parse.cpp` | No decimal-point support; integer exponents only.  Kira's normal output never trips this. |
| Diffeq nested reduce raises rank/dot via second `apply_jdot_jrank_floor` | `src/ibp/reduce.cpp` | Conservative; never lower than upstream. |
| `r = nonzero(top) + IBPDot` arithmetic match | `src/ibp/kira_yaml.cpp` | Spot-checked but not exhaustively. |

---

## 4. ⚪ Not ported (intentional)

- `SolveIntegralsGaugeLinkSingle` / `SolveIntegralsGaugeLink` (gauge-link / HQET / SCET / Wilson lines).
- `ExpandGaugeX`, `GenerateSquare` (linear-propagator reshaping for gauge links).
- `IBPReducer = "FiniteFlow+LiteRed"` — Kira-only by design; no FiniteFlow / FIRE / LiteRed / Blade bridge.
- `WSL` Windows compatibility plumbing in `Kira/interface.m`.
- `BinaryFormat` Windows newline handling.
- `$PermutationOption` (Kira yaml field; never set in benches).
- `FireFly` / `Mixed` / `NoFactorScan` reduction modes (Kira-only `Masters` and `Kira` modes used).
- `install.m` discovery scripts.

---

## 5. Upstream drift since v1.0 port

The upstream master branch has progressed since the C++ port was
performed.  Diffs against our reference snapshot:

- **New `Trivial` ending scheme** (`AMFlow.m:259`, `1016-1095`):
  upstream auto-appends `"Trivial"` to the user's `EndingScheme` list
  (`AMFlow.m:1034`) so there is always a fallback.  Our port has only
  `Tradition`, `Cutkosky`, `SingleMass`.  Latent: a configuration that
  trips all three of ours would error where MMA would silently fall
  back.  **No oracle currently hits this.**
- **New `UseCache` and `SkipReduction` options** (`AMFlow.m:259`):
  caching/skip toggles for AMFSystem persistence.  Not ported.  Not
  blocking.
- **Source-comment line-number drift**: ~30 lines off in places (the
  upstream files have grown ~3 % since the port was performed).
  Citations like `// Mirrors AMFlow.m:1342/1351` now nominally point
  to `:1368/1377`.  All 39 in-source citations were collected; none
  affects correctness.  Recommend pinning citations to a specific
  upstream commit hash in a future cleanup pass.

These upstream-side changes are all **additions**; no upstream
behaviour we already ported has changed.

---

## 6. Follow-up items (in rough priority order)

1. **Port `Trivial` ending scheme** so the scheme dispatcher can never
   exhaust without a fallback (D-class — minor; latent).
2. **Project parent `cut` into boundary sub-system** (D3).
3. **Defensive Jacobian assert in `boundary_integrands`** (D4).
4. **`ComplexMode` / imaginary-numeric pipeline in Kira interface** (D5);
   bigger feature.
5. **Throw on RHS-not-in-master in `reduce.cpp`** (D6); one-line.
6. **Pin source-comment line citations** to an upstream commit hash
   and add a banner in `docs/REFERENCE_MAP.md` explaining the policy.
7. **Add oracles for the 🟡 unverified branches** that touch user-
   reachable behaviour (priority: per-system Direction, multi-invariant
   `LIBPDeriv`, Tradition-with-cut family, scaleless-via-`Numeric`).
8. **Stale comment cleanup**: `qft/topology.cpp:340` references
   `AnalyzeComponent` (the MMA name) for what is now `make_component`.

None of these are release-blockers for v1.0.  Items 1, 2, 3, 5, 8 are
small enough to bundle into a v1.1; items 4, 6, 7 are projects.
