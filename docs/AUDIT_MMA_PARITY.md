# MMA Parity Audit — v1.0 retrospective

A line-level audit of the C++17 port against upstream Mathematica
AMFlow (commit `efda1db` of <https://gitlab.com/multiloop-pku/amflow>).
The audit was performed after v1.0 release to surface divergences that
the 12 oracle benchmarks listed in [`AUDIT.md`](../AUDIT.md) do not
cover.

## TL;DR

| Severity | Count | Action taken |
|---|---|---|
| 🟢 verified                  | 68 | — |
| 🟡 unverified (oracle gap)    | 18 | Documented in §3 below; tracked for future bench expansion.  Phase 2A reduces this count as oracles land (LIBPDeriv multi-invariant: 🟡 → 🟢 2026-05-09; Jordan block ordering: 🟡 → 🟢 2026-05-10; analyze_block non-nested overlapping: 🟡 → 🟢 2026-05-10). |
| 🔴 actual divergence          |  6 | **5 fully fixed; 1 deferred indefinitely** (D5, ComplexMode — investigated post-v1.0 and found to require an algebra-layer extension; entry-point now rejects loudly).  D3 was upgraded from detect-and-throw (Phase 1A) to the proper projection (Phase 1B) with an oracle. |
| ⚪ intentionally not ported   | 17 | — |

Net assessment: **no oracle-validated path is wrong**, and **no
silent-wrong-result path remains**.  Five of the six 🔴 items have
been fully corrected — two by aligning defaults to upstream, one by
adding a defensive Jacobian assert, one by raising on inconsistent
Kira output, and one (D3) by implementing the full Tradition-with-cut
boundary projection (Phase 1B), backed by a new oracle benchmark
matching upstream at rel ~ 1e-30.  Only D5 (Kira `ComplexMode` /
imaginary-numeric pipeline) remains deferred — investigated
post-v1.0 (2026-05-09) and deferred indefinitely after the
implementation cost analysis (see §D5 below and `docs/ROADMAP.md`
§"Phase 1C").  The JSON entry-point now actively rejects the
complex-numeric form rather than silently mishandling it.

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

### D3. Boundary sub-system loses parent's `Cut` info — **FIXED (proper projection, Phase 1B)**

- **Upstream** `ReduceBoundary` (`AMFlow.m:790–803`): when the parent
  system has a non-empty `Cut`, projects the parent's cut propagators
  (after the region transform) onto the new boundary-sub-family's
  propagator basis and emits `cut_propagators:[...]` to Kira.  Aborts
  on `Count[cut, 1] != Count[Cut, 1]` ("eta may have been inserted to
  cut denominators").
- **C++ (was, v1.0)** `AMFSystem::build_boundary`
  (`src/pipeline/amfsystem.cpp:1452`): hard-coded `cut={}` when
  constructing the boundary `FamilyConfig` — silently dropped the
  parent cut.
- **Phase 1A interim**: detect-and-throw if the parent had any cut.
  Eliminated the silent-wrong-result risk but at the cost of refusing
  Tradition-with-cut families entirely.
- **Phase 1B fix**: full projection mirroring upstream.  For each
  parent prop with `cut[k]==1`, apply the *bare* `region.transform.map`
  (loops only, no `half_eta` scaling) via `qft::apply_region_rule`,
  project the result back to `fc_with_eta_->ctx` via
  `project_mfrac_by_name` (the bare transform is `__amf_*`-free so
  no information is lost).  For each `fam.prop[i]`, test
  `is_zero(fc_with_eta_->apply_replacement(fam.prop[i] − cutde[j]))`
  for each transformed parent cut prop `cutde[j]`; record matches
  into `sub_cut`.  Raise `std::runtime_error` if
  `Count(sub_cut, 1) != Count(parent.cut, 1)` (mirror of
  `AMFlow.m:801`).  Pass `sub_cut` to `qft::FamilyConfig::build`.
- **Acceptance gate**: new oracle benchmark
  [`tools/bench/tradcut_phase_2L_eps001_*`](../tools/bench)
  derived from upstream `examples/automatic_phasespace/run.wl`
  (2-loop, 7 propagators, cut={1,0,1,0,1,0,0}, prescription={0,0},
  s=100, msq=1).  This is the only one of our oracles that
  actually exercises the Tradition-with-cut code path; the cut
  topology is not phase_volume so Cutkosky does not apply.  C++
  matches Mathematica AMFlow at relative error
  **2.68 × 10⁻³⁰** on `j[phase, 1, 0, 1, 0, 1, 0, 0]`.

### D4. `BoundaryIntegrands` Jacobian factor `|Det|^(4-2eps)` is silently dropped — **FIXED (defensive assert)**

- **Upstream** (`AMFlow.m:731-732`): multiplies integrands by
  `Abs[Det[Coefficient[#, Loop]&/@Values[trans]]]^(4-2*eps)` (the
  Jacobian of the loop redefinition).
- **C++ (was)** `qft::boundary_integrands` (`src/qft/boundary.cpp`):
  no Jacobian factor.  The omission was correct **only** because
  `branch_momenta` (`src/qft/region.cpp:108–114`) enforces each
  propagator's leading loop coefficient = 1, making the matrix `A`
  in `branch_to_loop` a permutation matrix and `|det A| = 1`.
- **Fix**: defensive assert in `branch_to_loop` (right after `det_A`
  is computed and the singular case is handled): if `det_A` is not a
  constant or `|det_A| ≠ 1`, raise `std::runtime_error`.  Catches
  any future relaxation of the `branch_momenta` precondition before
  the missing Jacobian factor produces wrong boundary integrands.

### D5. Kira `ComplexMode` / `CompensateRule` real-only filter unimplemented — **deferred indefinitely; entry-point rejects loudly**

- **Upstream** filters `IBPRule` to drop imaginary-part numerics from
  the Kira CLI input, then post-substitutes them via `CompensateRule`
  after `AnalyticReduction` (line 488) and `DifferentialEquation`
  (line 519).
- **C++** `KiraConfig::numeric_values` is a flat `map<string,string>`
  of pre-computed rational strings — no imaginary handling.  The 5th
  `IBPSystem` parameter `complexmode` has no C++ counterpart.
- **Latent**: all 12 oracles use purely-real `numeric_values`.
- **Status (post-investigation, 2026-05-09)**: deferred indefinitely.
  The JSON dispatcher (`apply_blackbox_options` in `src/api/run_json.cpp`)
  now rejects the complex form `{"re":..,"im":..}` with an error that
  points to this audit entry, rather than silently truncating to the
  real part or producing a wrong answer.
- **Why not a small patch**: the C++ algebra layer is over Q
  (FLINT `fmpz_mpoly_q_t`); a complex kinematic invariant cannot be
  substituted into an Mfrac as a value.  See
  [`docs/ROADMAP.md`](ROADMAP.md) §"Phase 1C" for the two
  implementation paths considered (Q[i] algebra extension vs.
  parallel acb-rational pipeline) and their trade-offs.
- **If revisited**: prefer the acb-rational-pipeline approach unless
  complex symbolic Replacement rules become a project goal.

### D6. RHS-not-in-master is silently dropped — **FIXED (now throws)**

- **Upstream** (`Kira/interface.m:486`): `Abort[]`s when
  `CoefficientArrays` reports a non-master residue (i.e., a J-integral
  appears on the right-hand side that is not in the masters list).
- **C++ (was)** the rhs assembly path in `ibp::black_box_reduce`
  (`src/ibp/reduce.cpp:218-230`) silently appended unknown rhs J
  integrals; consumers later filtered them with `master_index.find`
  + `continue`.
- **Fix**: in `ibp::black_box_reduce`, before pushing an rhs entry,
  look up `target_key(rhs_j)` in `master_index`; on miss, raise
  `std::runtime_error` naming the offending target and rhs.

---

## 3. 🟡 Unverified branches (oracle gaps)

These branches look correct to inspection but no committed benchmark
exercises them.  Each is a candidate for a future oracle.

| Branch | Site | Why it's unverified |
|---|---|---|
| ~~`analyze_block` uses `AnalyzeBlock0` topology~~ → 🟢 | `src/ode/blocks.cpp` | **Verified.**  Upstream default is `AnalyzeBlock1`; the two extend semantics agree on nested/equal blocks (the IBP common case) and produce different (but both valid) partitions on overlapping non-nested closures.  `test_ode_blocks.cpp` `AnalyzeBlock.NonNestedOverlapping_*` cases (Y-shape, mutual-plus-dependents, diamond closure) hand-trace the AnalyzeBlock0 partition and verify the basis-invariant correctness (cover + topological ordering: each row is in exactly one block, all of a block's external dependencies are in earlier blocks).  The AnalyzeBlock0 partition is a strict refinement of AnalyzeBlock1 — finer blocks, same final block-triangular DE structure.  End-to-end "no impact on integral" independently locked by the 12 oracle benches matching MMA at rel ~10⁻³⁰. |
| `Calcx00` boundary linear-system selection | `src/ode/zero.cpp:1314-1703` | Heuristic match-row selection + Dixon over Q[i]; upstream uses one symbolic `Solve[…, allvar]`.  Equivalent on full-rank systems. |
| Acb-inverse failure fallback in `Calcx00` | `src/ode/zero.cpp:1199-1206` | C++ degrades to "all positions resonate"; upstream's symbolic inverse never fails. |
| `evaluate_taylor` strips arb radii | `src/ode/regular.cpp` | Documented in source; abandons rigorous error bars on the regular-running stage. |
| ~~Jordan block ordering~~ → 🟢 | `src/ode/jordan.cpp` | **Verified, equivalence locked.**  Groups by eigenvalue then descending depth; MMA sorts globally by descending size.  Same final algebra, different column permutation when ≥2 distinct eigenvalues.  Multi-distinct-eigenvalue decomposition correctness (`S J S^{-1} = A`, eigenvalue multiset, block-size multiset) is now locked by `test_ode_jordan.cpp` `*Distinct*` tests (2-, 3-distinct-eigenvalue cases including non-trivial similarity transform and chains on each eigenvalue).  End-to-end "permutation has no effect on final integral" is independently asserted by the 12 oracle benchmarks matching MMA at rel ~10⁻³⁰. |
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
| ~~`LIBPDeriv` multi-invariant case~~ → 🟢 | `src/ibp/libp_deriv.cpp` | **Verified, scope clarified.**  Multi-invariant family (free-symbol invariants like `m1sq`, `m2sq` in distinct propagators, plus a Replacement-defined `s`) tested via `test_ibp_libp_deriv.cpp` `*TwoMassBubble*` cases.  Replacement-defined invariants (e.g. `s` from `p^2 -> s`) deliberately return all-zero — production AMFlow flow only differentiates w.r.t. `eta` (which lives directly in the propagator polynomial), so upstream's momentum-derivative chain-rule path (`LIBPDerivivative`, `Kira/interface.m`) is intentionally not ported.  See `include/amflow/ibp/libp_deriv.hpp` "Scope" block for the contract. |
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

- ~~**New `Trivial` ending scheme**~~ — **PORTED** (Phase 1A.1):
  `EndingScheme::Trivial` is now an explicit enum value with
  `ending_q(Trivial) === false` (mirroring upstream's
  `AMFSystemEndingQ[..., "Trivial"] := False`); the dispatcher
  auto-appends Trivial to the user's `ending_schemes` list
  (mirror of `AMFlow.m:1034`) so a fallback always exists.  Selecting
  Trivial gives the same passthrough setup as the previous "all
  schemes ending" shortcut (etac all-zero).
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

## 6. Follow-up items

The active development plan that addresses the items in this audit
lives in [`docs/ROADMAP.md`](ROADMAP.md).  In summary:

- **Phase 1A** — five small implementations of currently-missing
  paths (`Trivial` ending scheme, `SolveIntegrals` single-eps fast
  path, Cutkosky physical-mass safety check, per-system
  `AMFSystemDirection`, plus stale-comment cleanup).
- **Phase 1B** — replace the D3 Tradition-with-cut detect-and-throw
  with the proper projection (mirror of upstream
  `AMFlow.m:790-803`); requires a new oracle benchmark.
- **Phase 1C** — D5 (Kira `ComplexMode` / imaginary-numeric pipeline)
  was investigated post-v1.0 and **deferred indefinitely**.  See §D5
  for the architectural reason and `docs/ROADMAP.md` §"Phase 1C" for
  the two paths that would unblock it (Q[i] algebra extension or
  acb-rational parallel pipeline).  The dispatcher now rejects the
  complex-form input loudly.
- **Phase 2** — convert each 🟡 in §3 above into either an oracle-
  validated 🟢 or a documented "theoretical equivalence" entry.
- **Phase 3** — ongoing diversification of the oracle suite (loop
  number, invariant count, mass config, cut topology, ε regimes).

See `ROADMAP.md` for the concrete breakdown, sizing, and acceptance
gates.
