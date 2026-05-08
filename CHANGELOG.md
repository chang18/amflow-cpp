# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

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
