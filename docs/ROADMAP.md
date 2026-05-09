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
| `Calcx00` heuristic match-row vs symbolic Solve | rank-deficient boundary system | pending |
| Jordan block ordering | ODE block with ≥2 distinct eigenvalues | ✅ done (2026-05-10) — `test_ode_jordan.cpp` `*Distinct*` cases (2-/3-distinct-eigenvalue, mixed block sizes, similarity transform); end-to-end "no impact on final integral" independently locked by 12 oracle benches |
| `LIBPDeriv` multi-invariant | family with two kinematic invariants and a derivative bench | ✅ done (2026-05-09) — `test_ibp_libp_deriv.cpp` `*TwoMassBubble*` covers free-symbol multi-invariant chain rule; Replacement-defined invariant scope documented in `include/amflow/ibp/libp_deriv.hpp` |

### 2B — medium-risk — ~2 days

| 🟡 item | Case |
|---|---|
| `SingleMassQ` literal-vs-Numeric | mass left symbolic without `Numeric` |
| `factorize_family` mass=−1 detection | likewise |
| `MasterRank`/`MasterDot` non-default filter | explicit filter values |
| `r = nonzero(top) + IBPDot` arithmetic | non-trivial `IBPDot` value |

### 2C — low-risk (theoretical equivalence; documentation only) — ~1 day

For each item below: write a unit test or a paragraph in
`AUDIT_MMA_PARITY.md` documenting why the C++ choice is provably
equivalent to the upstream choice.

- `evaluate_taylor` strips arb radii (intentional, documented)
- `zero_sector_q` uses generic primes (any generic point works)
- `region_power` skips `/.Numeric` (no invariants in expression)
- `factorize_family` no-redef fallback (dead-code; add unit test)
- Coefficient parser fragility (add unit test for the supported grammar)

---

## Phase 3 — Active oracle expansion (ongoing, ~10+ days cumulative)

Beyond the 🟡 list, **diversify** the oracle set so that
"we-didn't-think-of-this" bugs have a chance to surface.

| Diversity axis | Current coverage | Target additions |
|---|---|---|
| Loop number  | L=1, 2, 3                          | L=4 (e.g. 4-loop banana, 4-loop sunrise) |
| Invariants   | 1 (`s`) or 2 (`s, t`)              | ≥3 (electroweak: `s, t, m_W, m_Z`) |
| Mass config  | all-massless / all-equal-mass      | mixed (some massive, some massless) |
| Cuts         | single phase-space (`cutbubble`/`cutsunrise`/`cutbanana`) | multi-cut Cutkosky |
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
