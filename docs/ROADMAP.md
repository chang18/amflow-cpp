# Roadmap

This document records the post-v1.0 development plan agreed with the
maintainer.  Priority is **implementation completeness and bug
exposure first**, ahead of accessibility (Docker / Python bindings)
or new physics features.

The audit
([`docs/AUDIT_MMA_PARITY.md`](AUDIT_MMA_PARITY.md))
left two kinds of work behind:

1. **Known-incomplete implementation paths.**  Every "deferred" 🔴
   item or "abort-instead-of-compute" stub.  Each is either
   implemented properly, or — if the implementation cost is
   incommensurate with the use case — explicitly deferred with a
   *loud* entry-point rejection so users get a clear error rather
   than a silent wrong answer.  Nothing is papered over.
2. **Oracle-uncovered code paths.**  All 21 🟡 unverified branches.
   For each, either add an oracle benchmark that exercises it (turning
   🟡 → 🟢), or document it as a known-equivalent design choice
   (e.g. theoretical equivalence, intentional dead-code fallback).

Items are listed in execution order.

---

## Phase 1 — Close all known implementation gaps

### Phase 1A — small fixes (~½ day total)

| # | Item | Upstream reference | Status |
|---|---|---|---|
| 1 | Port `Trivial` ending scheme (auto-fallback when other schemes don't match) | `AMFlow.m:1016-1095` | ✅ done |
| 2 | `SolveIntegrals` single-eps fast path (skip Laurent fit when user supplies `eps -> value` in `Numeric`) | `AMFlow.m:1364-1374` | ✅ done |
| 3 | Cutkosky physical-mass safety check (abort on negative `cutcom[[1,5]]` masses) | `AMFlow.m:1050` | ✅ done |
| 4 | Per-system `AMFSystemDirection` (compute per-system Im / NegIm from prescriptions of η-touching loops) | `AMFlow.m:981-991` | ✅ done |
| 5 | Pin in-source upstream-line citations to commit `efda1db` via REFERENCE_MAP banner | — | ✅ done |

### Phase 1B — D3 proper Tradition-with-cut projection — ✅ done

Replaced the Phase 1A detect-and-throw guard at
`src/pipeline/amfsystem.cpp:1452` with the real projection:

- For each parent propagator with `cut[k] == 1`, transform via the
  bare `region.transform.map` (matches upstream `region[[1]]` — no
  half-eta scaling) using `qft::apply_region_rule`.
- Project the result back to `fc_with_eta_->ctx` via
  `project_mfrac_by_name` (the bare transform is `__amf_*`-free).
- For each `fam.prop[i]`, test
  `is_zero(fc_with_eta_->apply_replacement(fam.prop[i] − cutde[j]))`
  against each transformed parent cut prop; record a 1 in `sub_cut[i]`
  on first match.
- Raise on `Count(sub_cut, 1) != Count(parent.cut, 1)` (mirror of
  `AMFlow.m:801` — "eta may have been inserted to cut denominators").
- Pass `sub_cut` to `qft::FamilyConfig::build`.

**Acceptance gate**: new oracle
[`tools/bench/tradcut_phase_2L_eps001_*`](../tools/bench) (2-loop
Tradition-with-cut probe, family from upstream
`examples/automatic_phasespace/run.wl`).  C++ matches MMA at
relative error **2.68 × 10⁻³⁰** on `j[phase, 1, 0, 1, 0, 1, 0, 0]`
(`s = 100, msq = 1, eps = 1/100`).  This is the only oracle
covering Tradition-with-cut.

### Phase 1C — D5 ComplexMode / imaginary-numeric pipeline — **deferred indefinitely**

**Status (post-investigation, 2026-05-09)**: investigated and deferred;
the JSON entry-point now actively rejects complex-form numeric values
with a clear "not implemented" error rather than silently mishandling
them.

**Why deferred**: the upstream `IBPRule` / `CompensateRule` split
(`Kira/interface.m:50-57`, `488`, `519`) cannot be ported by the
originally-imagined small patch.  The reason is structural — the C++
algebra layer carries Mfracs over **Q** (FLINT `fmpz_mpoly_q_t`), so a
complex kinematic invariant like `s = 1 + 2*I` cannot be substituted
into an Mfrac as a value.  Two faithful-to-upstream paths were
identified, both substantial:

- **Direction A** — promote `algebra::Mpoly` / `algebra::Mfrac` to
  Q[i] (Gaussian rationals).  FLINT does *not* provide a
  multivariate-polynomial type over Z[i], so the entire algebra-layer
  arithmetic (add, mul, GCD, content/primitive, substitute, etc.)
  would need a new Q[i]-coefficient implementation.  Cascades into
  `numeric::FmpqPoly`/`RationalFunction`/`RationalMatrix`, since the
  fundamental-matrix path consumes Mfracs and produces Q-rational
  output.  Estimated **~2000 LOC + 1-2 weeks**.

- **Direction B** — keep Mfracs over Q; defer the substitution of
  complex invariants to the rational-function-construction step
  (`mfrac_to_eta_rational`, `boundary_acbmat`).  Add a parallel
  pipeline using `numeric::AcbPoly` / `AcbRationalFunction` /
  `AcbRationalMatrix` (Arb provides `acb_poly_t` natively).  Mathematically
  equivalent to Direction A by the universal property of polynomial
  rings (substitution at the end commutes with rational-function
  arithmetic).  Estimated **~1100 LOC + ~7 days**.

**Workaround for end users**: real-only kinematics is the only
supported regime.  The JSON dispatcher now rejects the object form
`{"re":..,"im":..}` in `amf_options.blackbox.numeric_values` with an
error that points to this audit divergence (D5).

**If revisited**: build Direction B unless complex Replacement rules
(not just complex Numeric) become a project goal.

---

## Phase 2 — 21 🟡 → 🟢 (oracle expansion / bug exposure) (~7 days)

Each 🟡 item from the audit gets either a new oracle benchmark or an
explicit "documented as theoretical equivalence" record.

### 2A — high-risk (algorithm differs; could harbour bugs) — ~4 days

Each item: 1 day to construct the case, generate MMA reference, commit
the triplet.

| 🟡 item | Case to construct | Status |
|---|---|---|
| `analyze_block` topology choice | non-nested overlapping-closure block matrix | ✅ done (2026-05-10) — `test_ode_blocks.cpp` `AnalyzeBlock.NonNestedOverlapping_*` (Y-shape / mutual-plus-dependents / diamond), basis-invariant correctness locked + AnalyzeBlock0 output hand-traced |
| `Calcx00` heuristic match-row vs symbolic Solve | rank-deficient boundary system | ✅ done (2026-05-10) for the **full-rank** path — 6 existing `Calcx00_*` unit tests in `test_ode_zero.cpp` lock equivalence with hand-derived exact solutions; 12 oracle benches confirm at rel ~10⁻³⁰.  The genuinely rank-deficient case (acb-inverse fallback path at `zero.cpp:1199-1206`) remains 🟡 — Phase 3 synthetic-bench candidate |
| Jordan block ordering | ODE block with ≥2 distinct eigenvalues | ✅ done (2026-05-10) — `test_ode_jordan.cpp` `*Distinct*` cases (2-/3-distinct-eigenvalue, mixed block sizes, similarity transform); end-to-end "no impact on final integral" independently locked by 12 oracle benches |
| `LIBPDeriv` multi-invariant | family with two kinematic invariants and a derivative bench | ✅ done (2026-05-09) — `test_ibp_libp_deriv.cpp` `*TwoMassBubble*` covers free-symbol multi-invariant chain rule; Replacement-defined invariant scope documented in `include/amflow/ibp/libp_deriv.hpp` |

### 2B — medium-risk — ~2 days

| 🟡 item | Case | Status |
|---|---|---|
| `SingleMassQ` literal-vs-Numeric | mass left symbolic without `Numeric` | ✅ done (2026-05-10) — qft-layer `single_mass_q` stays upstream-literal; pipeline wrapper `single_mass_q_numeric` is the deliberate Numeric-substitution enhancement, locked by `test_qft_amfmode.cpp` `SingleMassQ_*` (2) + `test_amflow_amfsystem.cpp` `SingleMassEnhancement_*` (2) |
| `factorize_family` mass=−1 detection | likewise | ✅ done (2026-05-11) — `find_mass_minus_one` applies Numeric before literal `-1` test (same enhancement pattern as SingleMassQ); end-to-end locked by `test_amflow_amfsystem.cpp` `FactorizeFamilyMassMinusOne_NumericResolvesSymbolic` |
| `MasterRank`/`MasterDot` non-default filter | explicit filter values | ✅ done (2026-05-10) — intentionally not exposed; C++ pins to upstream default `Infinity` (no filter); contract documented in `include/amflow/ibp/reduce.hpp` and audit row |
| `r = nonzero(top) + IBPDot` arithmetic | non-trivial `IBPDot` value | ✅ done (2026-05-10) — `test_ibp_kira.cpp` `WriteJobs_R_*` (zero/non-trivial IBPDot, mixed pattern, non-binary defensive); C++ formula tightened to `count_if(!= 0)` to match upstream literal |

### 2C — low-risk (theoretical equivalence; documentation only) — ✅ done (2026-05-11)

| 🟡 item | Closure |
|---|---|
| `evaluate_taylor` strips arb radii | ✅ `test_ode_regular.cpp` `EvaluateTaylor_StripsArbRadii_LocksMidpointOnlyContract` — feeds non-zero arb radius into a coefficient and asserts output radius is zero |
| `zero_sector_q` uses generic primes | ✅ Documented theoretical equivalence (generic-point evaluation suffices for scalelessness); locked by 5 existing `ZeroSectorQ_*` tests covering 1-loop and 2-loop scaleless / non-scaleless cases |
| `region_power` skips `/.Numeric` | ✅ `test_qft_region.cpp` `RegionPower_OutputIsEpsOnlyPlusIntegers_AuditRow200Equivalence` — multi-invariant box family, asserts every output term has zero exponent on every non-eps variable |
| `factorize_family` no-redef fallback | ✅ Defensive fallback documented in source comment + audit row; trigger condition is rank-deficient redef matrix; unreached on tested inputs but mathematically safe (leaves family unchanged) |
| Coefficient parser fragility | ✅ `test_ibp_kira.cpp` adds 4 negative-acceptance tests (`ParseExpression_DecimalLiteralThrows`, `_NonIntegerExponentThrows`, `_GarbageCharacterThrows`, `_UnclosedParenThrows`) on top of 7 existing positive-acceptance tests |

---

## Phase 3 — Tail audit closures + Active oracle expansion (ongoing, ~10+ days cumulative)

### 3.A–3.E — Final 5 🟡 closures (all done 2026-05-11)

The 5 🟡 audit rows that survived through Phase 2 are now all 🟢:

| 🟡 item | Closure approach |
|---|---|
| 3.A Diffeq nested reduce conservative rank/dot floor | Audit equivalence paragraph: inner `apply_jdot_jrank_floor` is applied to already-floored opts, so effective `(rank, dot) >=` upstream's analogous step.  Larger IBP system is a superset of smaller one — final reduction rules match. |
| 3.B `kira_target.m` parser strictness | 5 new negative-acceptance tests (`ReadTargetTable_*`): missing-file/no-brace → empty (silent OK), malformed LHS / RHS-without-J → throw, literal `0` RHS → filtered |
| 3.C Auto-applied `Vacuum[L,n]` table | All 5 (L, n) entries already had per-entry numeric reference tests in `test_qft_vacuum.cpp`; coverage probe + unknown-throws also locked.  Audit row promoted with reference. |
| 3.D Ending-master Kira reduction loop | New end-to-end test `EndingTadpole_J2_RecursiveKiraLoweringMatchesReference`: provides J[tad, 2] preferred without explicit_boundary, forces Kira-lowering loop to derive J[2] from J[1] via IBP identity; matches MMA Laurent at two eps points (Kira-gated, ~4.5s) |
| 3.E `Calcx00` acb-inverse fallback (rank-deficient) | Documented as **conservative numerical safety net**: after `normalize_mat`, `a00` is Jordan-form so the fallback is unreachable in well-behaved precision; the "all-resonate" treatment is mathematically conservative-but-correct (downstream pair constraints resolve spurious log terms to zero).  Synthetic trigger requires unrealistic chop_pre or bypassing normalize_mat — both outside production usage. |

After 3.A–E, the audit table contains **86 🟢 / 0 🟡 / 6 🔴 (5 fixed + 1 deferred D5) / 17 ⚪**.

### 3.F — L=4 banana oracle (2026-05-11) + D7 fix (2026-05-12) — **completed, matches MMA at rel 7.2 × 10⁻³¹**

First L=4 diversity-axis bench: 4-loop equal-mass banana sunrise
(family `banana4`, 4 loops + 1 external leg, 5 massive propagators
+ 9 ISPs, target `J[1,1,1,1,1,0,…,0]` at psq=-3, msq=1, eps=1/1000).

- **MMA reference** completed cleanly in 385 s.  10 sampled
  integral values committed in `tools/bench/banana_4loop_eps001_mma_reference.json`.
- **C++ failed** with `error: ibp::diffeq: master count differs
  between Masters and Reduce calls`.  Root cause: dual-Kira-call
  structural divergence from upstream's single-call pattern.
  Surfaced as **new audit divergence D7** in `docs/AUDIT_MMA_PARITY.md`.
- **Status**: bench triplet committed regardless (MMA reference is
  a useful invariant for verifying any future D7 fix; the cpp.json
  + mma.wl artefacts let the bench be re-run end-to-end once D7 is
  resolved).  Bench is in `run_perf_audit.sh` rotation with a
  marker comment.

D7 fix landed 2026-05-12 (`src/ibp/reduce.cpp` `reduce` and
`diffeq`):

- Both functions now do a Masters-mode preheat call (sector-wide
  enumeration via `select_mandatory_recursively`) followed by a
  Reduce-mode call (`select_mandatory_list` for the specific
  targets), at the SAME `(rank, dot)`.  Each call runs in its own
  subdirectory (`<work_dir>/masters_preheat/` and
  `<work_dir>/target_reduce/`) because our Kira 2.x release
  refuses to share a `$ReductionDirectory` between the two calls
  (the Masters-mode `run_initiate: masters` doesn't register the
  `-s` numeric substitutions, so the subsequent Reduce-mode call
  Aborts with `Kira::update_auxiliary_file: Last Kira run set 0
  variables to numeric values, this time you request N`).
  Upstream MMA shares one dir and relies on Kira's tmp/-based
  incremental state; we accept the small overhead of re-running
  the IBP setup in the second dir for the architectural-parity
  guarantee.
- `ibp::diffeq` no longer calls `reduce()` for the inner step
  (which would have re-floored `(rank, dot)` over the derivative
  integrals — the original L=4-banana bug).  Instead it inlines
  the Reduce-mode Kira invocation at `opts_eff`'s `(rank, dot)`,
  mirroring upstream's `AnalyticReduction` which inherits
  `IBPRank`/`IBPDot` globals from `IBPSystem`.
- Both functions add a SubsetQ guard on the Reduce-mode master
  file against the Masters-mode sector enumeration, mirroring
  upstream's `If[!SubsetQ[masters, str], Abort["inconsistent
  masters from Kira"]]`.

Locked by L=4 banana oracle (rel 7.2 × 10⁻³¹ / 1.6 × 10⁻³⁰, 399 s).
All 545 pre-existing tests remain green.

### 3.G — Multi-invariant oracle (1-loop electroweak box, 2026-05-12) — **completed, matches MMA at rel ≤ 3.8 × 10⁻³⁰**

First oracle for the **Invariants ≥3** diversity axis: 1-loop
electroweak box `ewbox`, 4 propagators with alternating W/Z masses,
4 massless external legs, 4 distinct kinematic invariants
{`s, t, mWsq, mZsq`} at `s = 7, t = -3, mWsq = 1, mZsq = 4/3`,
`eps = 1/1000`.

- Targets: corner `j[ewbox, 1, 1, 1, 1]` + 6 sub-masters (W/Z
  tadpoles, WW/ZZ bubbles, WZW/WZZ triangles).  All 7 sampled
  values pass at rel ≤ 4 × 10⁻³⁰ (well within project tolerance).
- MMA reference: 62 s.  C++ wallclock: 39 s (faster at 1-loop scale).
- **No bug surfaced.**  The multi-invariant code paths
  (`kira_yaml.cpp::mass_scale_names` → `kira_run.cpp` `-s` emission
  → `libp_deriv` multi-invariant chain rule) all behaved correctly
  end-to-end.  The WW bubble being above threshold (`s = 7 > 4mW² = 4`)
  gives physical imaginary parts that lock the cross-mass-scale
  branch-cut handling too.

Triplet committed: `tools/bench/ewbox_1loop_eps001_black_box_amflow_{cpp.json,mma.wl}`
+ `tools/bench/ewbox_1loop_eps001_mma_reference.json`.  Added to
`run_perf_audit.sh` rotation.

### 3.H — Multi-cut Cutkosky oracle (4-loop massless cutbanana, 2026-05-12) — **completed, matches MMA at rel ≤ 2.2 × 10⁻³⁰**

First oracle for the **multi-cut Cutkosky** diversity axis: 4-loop
massless equal-mass cutbanana `cutbanana4` with all 5 internal
lines simultaneously on-shell, corresponding to a 5-particle
Cutkosky cut at `s = 4, eps = 1/100`.  Extends the cut series
(cutbubble 2-cut → cutsunrise 3-cut → cutbanana_3L 4-cut) to 5
simultaneous cuts in a 4-loop family.

- Targets: corner `j[cutbanana4, 1,1,1,1,1, 0,...]` + 3 dotted
  masters (5th index 2, 3, 4).  All 4 sampled values pass at rel
  ≤ 2.2 × 10⁻³⁰.  All imag parts zero as expected for pure-real
  Cutkosky phase-space integrals.
- MMA reference: 149 s.  C++ wallclock: 141 s (within 6%).
- **No bug surfaced.**  MMA log shows 5-level subsystem recursion
  (4-loop top down through 3/2/1-loop intermediate subsystems to
  the vacuum ending) and the C++ `AMFSystemSolution` chain mirrors
  it exactly — Cutkosky ending scheme + `cut_propagators` wiring +
  multi-system recursion all behaved correctly end-to-end.

Triplet committed: `tools/bench/cutbanana_4L_eps001_black_box_amflow_{cpp.json,mma.wl}`
+ `tools/bench/cutbanana_4L_eps001_mma_reference.json`.  Added to
`run_perf_audit.sh` rotation.

### Oracle diversity expansion (ongoing, no fixed end)

Beyond the (now-empty) 🟡 list, **diversify** the oracle set so that
"we-didn't-think-of-this" bugs have a chance to surface.

| Diversity axis | Current coverage | Target additions |
|---|---|---|
| Loop number  | L=1, 2, 3, **4 (banana, 2026-05-11)** | further L=4 (e.g. 4-loop sunrise) |
| Invariants   | 1 (`s`) or 2 (`s, t`), **4 (ewbox `s, t, mWsq, mZsq`, 2026-05-12)** | further multi-invariant at higher L |
| Mass config  | all-massless / all-equal-mass      | mixed (some massive, some massless) |
| Cuts         | single phase-space (cutbubble 2-cut, cutsunrise 3-cut, cutbanana_3L 4-cut), **5-cut (cutbanana_4L, 2026-05-12)** | further multi-cut topologies |
| Sector size  | ≤ 9 propagators                    | ≥ 10 |
| `eps` extremes | `1/1000`, `1/100`                | `1/1`, `1/10000` (boundary cases) |

Each new oracle is ~1 day of work (construct family → MMA reference →
C++ run → commit triplet).  5–10 new oracles = 1–2 weeks.  If any new
oracle exposes a bug, fix it inline and the bug becomes a
parity-validated 🟢 going forward.

This phase has no fixed end — it is the project's living quality
loop.

---

## What is **not** in this roadmap

The following items are **explicitly deferred** until Phase 1+2 are
complete and external user feedback exists:

- Docker image with pre-installed FLINT 3.4 + Kira + Fermat.
- Python bindings (`pybind11` wrapper around `api::run_json`).
- Conda-forge packaging.
- Multi-threading over ε samples.
- New physics features beyond Phase 3 oracle diversity.
- Continuous benchmarking in CI.

These are accessibility / outreach work; they amplify whatever
correctness baseline the project has.  Doing them before correctness
is fully closed would risk shipping a polished interface to
under-validated code.

---

## Update policy

- Mark each completed item with a date and the commit hash.
- When a 🟡 → 🟢 conversion exposes a real bug, link the bug-fix
  commit + the oracle that caught it.
- When a Phase 3 new-oracle case exposes a bug, do the same.
- Do not silently delete items that turn out to be infeasible —
  document the reason.

Phase 1+2 are bounded work (~16 days total).  Phase 3 is a continuous
practice.  The project's parity-validation surface only grows.
