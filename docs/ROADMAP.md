# Roadmap

This document records the post-v1.0 development plan agreed with the
maintainer.  Priority is **implementation completeness and bug
exposure first**, ahead of accessibility (Docker / Python bindings)
or new physics features.

The audit
([`docs/AUDIT_MMA_PARITY.md`](AUDIT_MMA_PARITY.md))
left two kinds of work behind:

1. **Known-incomplete implementation paths.**  Every "deferred" 🔴
   item or "abort-instead-of-compute" stub.  These must be
   implemented properly, not papered over.
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

### Phase 1B — D3 proper Tradition-with-cut projection (~3 days)

Replace the current detect-and-throw guard at
`src/pipeline/amfsystem.cpp:1452` with the real projection:

- For each parent propagator with `cut[k] = 1`, transform via
  `region.transform`.
- For each new sub-family propagator, symbolic-compare modulo
  `reduced_replacement` against the transformed cut props.
- Set the new `cut[]` array; assert
  `Count[cut, 1] == Count[parent_cut, 1]` (mirroring `AMFlow.m:801`).
- Pass to `qft::FamilyConfig::build`.

**Mandatory acceptance gate**: a new oracle benchmark with a
Tradition-scheme cut family (Mathematica reference + C++ output
matches at `rel ~ 1e-30`).  Without an oracle, the implementation
cannot be trusted.

### Phase 1C — D5 ComplexMode / imaginary-numeric pipeline (~5 days)

Implement the upstream `IBPRule` / `CompensateRule` split for complex
numeric kinematics:

- Split `KiraConfig::numeric_values` into real-vs-complex during
  `kira_run`; pass only real parts as `-s` arguments to Kira.
- After Kira returns, apply the complex-numeric values via a
  `CompensateRule`-style substitution in the parsed coefficient table
  and DE matrix.
- Wire through `ibp::black_box_reduce` and `ibp::black_box_diffeq`.

**Mandatory acceptance gate**: a new oracle benchmark with a
non-real-numeric kinematic invariant (e.g., a Cutkosky probe with a
genuinely complex kinematic point).

---

## Phase 2 — 21 🟡 → 🟢 (oracle expansion / bug exposure) (~7 days)

Each 🟡 item from the audit gets either a new oracle benchmark or an
explicit "documented as theoretical equivalence" record.

### 2A — high-risk (algorithm differs; could harbour bugs) — ~4 days

Each item: 1 day to construct the case, generate MMA reference, commit
the triplet.

| 🟡 item | Case to construct |
|---|---|
| `analyze_block` topology choice | non-nested overlapping-closure block matrix |
| `Calcx00` heuristic match-row vs symbolic Solve | rank-deficient boundary system |
| Jordan block ordering | ODE block with ≥2 distinct eigenvalues |
| `LIBPDeriv` multi-invariant | family with two kinematic invariants and a derivative bench |

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
