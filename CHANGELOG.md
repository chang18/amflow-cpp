# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## Unreleased — post-v1.0 audit-driven patches

Patch-level commits on top of v1.0.0 (no version bump yet).
Together they close all silent-wrong-result paths surfaced by the
post-release line-level MMA parity audit
([`docs/AUDIT_MMA_PARITY.md`](docs/AUDIT_MMA_PARITY.md)).

### Fixed (alignment with upstream)
- `RunningOptions::run_length` default raised from 200 to **1000** to
  match upstream `RunLength = 1000` (`AMFlow.m:259`,
  `DESolver.m:97`).  Prevents premature `RunUnit` aborts on crowded
  pole landscapes.  (Audit D1.)
- `GlobalOptions::rationalize_pre` default raised from 20 to **100**
  to match upstream `RationalizePre = 100`.  Removes the silent
  precision-narrowing on rationalization steps along the contour.
  (Audit D2.)
- `ibp::black_box_reduce` now throws `std::runtime_error` if Kira's
  reduction returns an RHS J-integral that is not in the master
  list, mirroring upstream `Kira/interface.m:486` `Abort`.  Prior
  behaviour silently dropped the offending row.  (Audit D6.)
- **D3 — Tradition-with-cut boundary projection.**  v1.0 silently
  dropped the parent system's `Cut` when constructing the boundary
  sub-family.  The Phase 1A patch converted that to a loud abort.
  **Phase 1B implements the full projection** mirroring upstream
  `ReduceBoundary` (`AMFlow.m:790-803`): each parent cut prop is
  passed through the bare `region.transform.map`; each fam.prop is
  matched against the transformed cut props modulo
  `reduced_replacement`; the new `sub_cut` is built and asserted to
  preserve `Count[cut, 1]`; the result is passed to the
  `qft::FamilyConfig::build` of the sub-family.  Backed by a new
  oracle benchmark (`tradcut_phase_2L_eps001_*`) — 2-loop
  Tradition-with-cut from upstream `examples/automatic_phasespace`,
  matches MMA at relative error 2.68 × 10⁻³⁰.  (Audit D3.)
- `branch_to_loop` (`src/qft/region.cpp`) gains a defensive assert
  that `det(A) = ±1` (constant), where `A` is the loop-redefinition
  matrix.  This catches any future relaxation of the `branch_momenta`
  unit-leading-loop-coefficient precondition before the missing
  `|Det|^(4-2eps)` Jacobian factor (upstream `AMFlow.m:731-732`)
  silently produces wrong boundary integrands.  (Audit D4.)

### Added
- [`docs/AUDIT_MMA_PARITY.md`](docs/AUDIT_MMA_PARITY.md): line-level
  upstream-parity audit report.  65 🟢 verified / 21 🟡 unverified /
  6 🔴 (5 fully fixed, 1 deferred to v1.1) / 17 ⚪ not ported.
- [`docs/PERFORMANCE.md`](docs/PERFORMANCE.md): wall-clock baseline
  vs MMA on the 12 oracle benchmarks.  Median speedup ~2× on
  Kira-light benches; converges to ~1× on Kira-heavy benches.
- [`docs/ROADMAP.md`](docs/ROADMAP.md): post-v1.0 development plan
  (Phase 1 implementation completeness, Phase 2 oracle expansion,
  Phase 3 ongoing diversification).
- [`tools/bench/run_perf_audit.sh`](tools/bench/run_perf_audit.sh):
  shell driver to reproduce the perf run.
- New oracle benchmark
  [`tools/bench/tradcut_phase_2L_eps001_*`](tools/bench): the first
  Tradition-with-cut parity test (D3 acceptance gate).
- New ending scheme: `EndingScheme::Trivial`, auto-appended to the
  user's `ending_schemes` list as a final fallback (mirror of
  upstream `AMFlow.m:1034`).
- New `solve_integrals` single-eps fast path: when
  `numeric_values["eps"]` is supplied, skip the Laurent fit and
  return the integral evaluated at that eps (mirror of
  `AMFlow.m:1364-1374`).
- New per-system path direction: `AMFSystem::setup()` computes the
  η-touching loops' prescription consensus and overrides the global
  `run_direction` for that system's ODE solve (mirror of
  `AMFlow.m:981-991`).
- Cutkosky setup now validates that all phase-volume component
  masses are non-negative after `Numeric` substitution; raises
  otherwise (mirror of `AMFlow.m:1050`).
- `apply_blackbox_options` (`src/api/run_json.cpp`) now actively
  rejects the complex-numeric object form
  `{"re":..,"im":..}` in `amf_options.blackbox.numeric_values`
  with a clear "not implemented" error pointing to audit divergence
  D5.  Prior behaviour would have failed deeper in the parser with an
  unclear message; the explicit rejection makes the feature gap
  visible at the input boundary.

### Known limitations (deferred indefinitely)
- **D5 — Kira `ComplexMode` / imaginary-numeric pipeline is not
  implemented.**  Users who need to evaluate at numeric kinematics
  with non-zero imaginary parts cannot do so via this port; the
  upstream filters such values through `IBPRule` / `CompensateRule`
  but the C++ port treats `numeric_values` as a flat real-valued
  map.  All 12 oracle benchmarks use purely-real numerics, so this
  surface is untested.  **Investigated post-v1.0 (2026-05-09) and
  deferred indefinitely** — see [`docs/AUDIT_MMA_PARITY.md`](docs/AUDIT_MMA_PARITY.md)
  D5 and [`docs/ROADMAP.md`](docs/ROADMAP.md) §"Phase 1C" for the
  architectural trade-off (Q[i] algebra extension vs. parallel
  acb-rational pipeline).
- ~~**`Trivial` ending scheme is not ported.**~~  Now ported (Phase 1A).
  Upstream auto-appends `"Trivial"` to the user's `EndingScheme` list
  as a fallback when none of `Tradition` / `Cutkosky` / `SingleMass` apply
  (`AMFlow.m:1016-1095`).  No oracle currently configures this
  fallback path.  Planned for v1.1.

## [1.0.0] — 2026-05-08

First public release.  C++17 reimplementation of the
auxiliary-mass-flow algorithm originally introduced in the Mathematica
package AMFlow by Liu and Ma
(<https://gitlab.com/multiloop-pku/amflow>; Comput. Phys. Commun. 283
(2023) 108565).  All algorithmic credit belongs to the upstream
authors; this project contributes only the C++ implementation and a
numerical-parity test/benchmark harness.

### Added
- Domain-oriented architecture under the `amflow::` namespace:
  `numeric`, `algebra`, `ode`, `qft`, `ibp`, `pipeline`, `api`, `cli`.
- Public headers under `include/amflow/<domain>/` with strict ABI/PIMPL
  boundaries between domains.
- ODE solver — reimplementation of the upstream `DESolver.m` covering
  the workflows exercised by the oracle benchmarks; matches the
  reference at relative precision ~1e-30 across all 12 oracle cases.
- AMFlow + Kira pipeline implementing the core algorithms of upstream
  `AMFlow.m`: `amflow`, `black_box_amflow`, and `solve_integrals` modes.
- `amflow_cli` — JSON-driven driver around `amflow::api::run_json`,
  exposing all three top-level modes.
- 500+ GoogleTest cases, including 12 sampled-parity oracle benchmarks
  under `tools/bench/` (`eps = 1/1000`).
- Mathematica reference drivers under `tools/math_ref/` and
  `tools/bench/` (require a local clone of upstream AMFlow; see
  `reference/README.md`), plus committed JSON reference outputs for
  regression tests.
- CMake install / export rules; downstream projects can use
  `find_package(AMFlowCpp)` and link `AMFlowCpp::amflow`.

### Notes
- The upstream Mathematica AMFlow source is **not** vendored.  See
  `reference/README.md` for how to obtain it locally to regenerate
  benchmark reference data.

### Out of scope
- `SolveIntegralsGaugeLink`, HQET / SCET / Wilson-line workflows — see
  `AUDIT.md`.

[1.0.0]: https://github.com/chang18/amflow-cpp/releases/tag/v1.0.0
