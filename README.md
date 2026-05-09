# AMFlow.cpp

[![CI](https://github.com/chang18/amflow-cpp/actions/workflows/ci.yml/badge.svg)](https://github.com/chang18/amflow-cpp/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)](CHANGELOG.md)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://en.cppreference.com/w/cpp/17)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20087172.svg)](https://doi.org/10.5281/zenodo.20087172)

A C++17 reimplementation of the auxiliary-mass-flow algorithm for
multi-loop Feynman integrals.

## Scientific attribution

This project re-implements the algorithms presented in the Mathematica
package **AMFlow** by Xiao Liu and Yan-Qing Ma:

- Upstream Mathematica package: <https://gitlab.com/multiloop-pku/amflow>
- Algorithm paper: X. Liu and Y.-Q. Ma, *AMFlow: A Mathematica package
  for Feynman integrals computation via auxiliary mass flow*, Comput.
  Phys. Commun. **283** (2023) 108565,
  [doi:10.1016/j.cpc.2022.108565](https://doi.org/10.1016/j.cpc.2022.108565).

This repository contributes **only** the C++17 implementation and a
numerical-parity test/benchmark harness against the Mathematica
reference; **all algorithmic credit belongs to the AMFlow authors**.
We claim no algorithmic novelty.

If you use this software in published research, please cite **both** the
upstream paper above (algorithm) and this repository (implementation).
The repository is archived on Zenodo with DOI
[10.5281/zenodo.20087172](https://doi.org/10.5281/zenodo.20087172)
(concept DOI; always resolves to the latest release).  Machine-readable
metadata is in [`CITATION.cff`](CITATION.cff); see [`NOTICE`](NOTICE)
for full attribution.

> **Note** — the upstream Mathematica source is **not** vendored.
> Documentation in this repo refers to symbols and line numbers in
> upstream files (`AMFlow.m`, `diffeq_solver/DESolver.m`,
> `ibp_interface/Kira/interface.m`); to follow those references and
> to regenerate Mathematica reference data, clone the upstream into
> `reference/amflow-master/` (gitignored).  See
> [`reference/README.md`](reference/README.md).

## Project provenance

This C++17 implementation was developed primarily by **AI coding
agents** under the direction of the maintainer, with every behavioural
change gated by numerical-parity verification against the upstream
Mathematica reference.  The 12 oracle benchmarks under `tools/bench/`
and the 508-case GoogleTest suite are the contract that the AI-led
implementation has to honour for any change to land.

Maintainer / contact: **3250800970@qq.com** (please open a GitHub
issue first when possible; use email for inquiries that don't fit a
public issue).

## Layout

Domain-oriented architecture:

```
amflow::numeric -> amflow::algebra -> {amflow::ode,
                                       amflow::qft -> {amflow::ibp,
                                                       amflow::pipeline -> amflow::api -> amflow::cli}}
```

Public headers under `include/amflow/<domain>/`; implementation under
`src/<domain>/`.

## Working rules

- If a workflow succeeds in upstream `AMFlow.m` and the corresponding
  branch is already ported here, a C++ mismatch is a parity bug.
- Do not invent extra correctness criteria beyond Mathematica parity.
- Sampled single-point parity at `eps = 1/1000` is the primary
  numerical contract; see `tools/bench/`.
- All Mathematica-computed benchmark values are kept committed as
  regression standards (raw caches are regenerable; see
  `tools/bench/README.md`).

## Status

| Area | State |
|---|---|
| Core ODE solver (port of upstream `DESolver.m`) | Implemented; 12/12 oracle cases match MMA reference at `rel ~1e-30` |
| AMFlow + Kira pipeline (upstream `AMFlow.m` core algorithms) | Implemented for the covered workflows |
| Top-level entries (`amflow`, `black_box_amflow`, `solve_integrals`) | Implemented; exposed via `amflow_cli` JSON modes |
| `SolveIntegralsGaugeLink`, HQET / SCET / Wilson lines | Out of scope |

See [`AUDIT.md`](AUDIT.md) for the validated families and benchmark
inventory.

## Quick start

On Ubuntu, the packaged dependencies needed for the default build are:

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake pkg-config libflint-dev \
  libgtest-dev nlohmann-json3-dev
```

```bash
cmake -S . -B build -DAMFLOW_BUILD_DRIVER=ON
cmake --build build -j32
ctest --test-dir build --output-on-failure -j 4

# Raw ODE example (no IBP)
./build/src/cli/amflow_cli examples/power_law.json

# Full AMFlow examples (require Kira + Fermat at runtime)
./build/src/cli/amflow_cli examples/box1_black_box_amflow.json
./build/src/cli/amflow_cli examples/bubble_solve_integrals.json
./build/src/cli/amflow_cli examples/box1_solve_integrals.json
```

To use as a CMake dependency once installed:

```cmake
find_package(AMFlowCpp 1.0 REQUIRED)
target_link_libraries(my_target PRIVATE AMFlowCpp::amflow)
```

## Mathematica reference generation

The repository includes Wolfram drivers for regenerating committed
reference data under `tests/data/math_ref/` and sampled benchmark data
under `tools/bench/`.  These drivers expect a local clone of the
upstream Mathematica AMFlow at `reference/amflow-master/`
(gitignored — see [`reference/README.md`](reference/README.md) for the
two acceptable layouts).

```bash
# Either clone upstream directly into reference/amflow-master/ ...
git clone https://gitlab.com/multiloop-pku/amflow.git reference/amflow-master

# ... or symlink to an existing clone:
# ln -s /path/to/your/clone reference/amflow-master

cd tools/math_ref
AMF_REF_CFG=$PWD/cfg_bubble_1L.wl math -script run_amflow_kira.wl
```

For current parity benchmarks at `eps = 1/1000` see
[`tools/bench/`](tools/bench/) and [`AUDIT.md`](AUDIT.md).

## Documentation map

**For users** — read in this order:

1. [`docs/USER_GUIDE.md`](docs/USER_GUIDE.md) — step-by-step walk-through of one no-IBP example and one full Feynman-integral example.
2. [`docs/JSON_SCHEMA.md`](docs/JSON_SCHEMA.md) — every field of every input/output, for all three CLI modes.
3. [`docs/FAQ.md`](docs/FAQ.md) — common build / runtime / parity questions.

**For contributors** — additionally:

| Document | What it covers |
|---|---|
| [`CONTRIBUTING.md`](CONTRIBUTING.md) | TL;DR contribution rules + maintainer contact |
| [`docs/CONTRIBUTING.md`](docs/CONTRIBUTING.md) | Detailed dev workflow: tests, benchmarks, references, debugging, env-var trace flags |
| [`AGENTS.md`](AGENTS.md) | AI-coding-agent onboarding (the AI-led development contract) |
| [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) | Domain DAG and end-to-end data flow |
| [`docs/INVARIANTS.md`](docs/INVARIANTS.md) | Project-wide rules that must stay true (precision, RAII, ODE-boundary traps) |
| [`docs/REFERENCE_MAP.md`](docs/REFERENCE_MAP.md) | Upstream Mathematica symbol → C++ symbol mapping |
| [`AUDIT.md`](AUDIT.md) | Current parity status, validated surface, benchmark inventory |
| [`notes/mma_*_map.md`](notes/) | Per-file upstream Mathematica source → C++ port maps |
| [`reference/README.md`](reference/README.md) | How to clone the upstream MMA AMFlow locally for reference data regeneration |
| [`tools/bench/README.md`](tools/bench/README.md) | Sampled-parity benchmark conventions and regeneration |

## Dependencies

| Component | Why |
|---|---|
| C++17 compiler | Core build |
| CMake >= 3.16 | Build system |
| FLINT / Arb | Rational, polynomial, and arbitrary-precision arithmetic |
| GoogleTest | Test suite |
| nlohmann/json | Driver JSON I/O |
| Kira + Fermat | IBP reduction backend (runtime; tests gracefully skip when absent) |
| Mathematica / Wolfram Engine + upstream AMFlow clone | Optional; only needed to regenerate reference data |

## Contributing

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for the GitHub-facing entry
point and [`docs/CONTRIBUTING.md`](docs/CONTRIBUTING.md) for the full
project workflow.  Behavioural changes must be backed by upstream
AMFlow parity evidence and focused C++ regression tests.

For questions, bug reports, or collaboration inquiries, open a GitHub
issue at <https://github.com/chang18/amflow-cpp/issues>, or contact
the maintainer at <3250800970@qq.com>.

## License

MIT — see [`LICENSE`](LICENSE) and [`NOTICE`](NOTICE).  The upstream
Mathematica AMFlow is also MIT-licensed; this project does not bundle
any upstream source.
