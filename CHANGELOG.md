# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- **Multi-mass stress batch (2026-05-15)**: two additional oracle
  triplets probing axes that surfaced a latent bug.  Oracle total:
  40 → 42.
  - `tools/bench/vtx2_2L_3mass_eps001_*` — 2L 3-leg vertex with three
    distinct internal masses (mAsq=1, mBsq=4, mCsq=9).  **Surfaced
    audit divergence D13** in C++ `build_boundary`: the projection
    from the boundary's Laporta coefficient to the sub-system's
    reduced context aborted with `residual var 'mAsq' not in dst`
    because `to_complete_explicit` rank-filtered out one of the
    mass-bearing propagators.  Fixed in the same commit; corner
    matches MMA at rel 9.78e-31 (Re).  MMA wallclock 187 s.
  - `tools/bench/doublebox_2mass_single_2L_eps001_*` — 2L doublebox
    with SPARSE 2-mass placement (one mass per loop on corner-only
    propagators).  Different from doublebox2m (interleaved),
    doublebox_blockmass (block), and doublebox_diagmass (rung).
    Does NOT trigger D12 precision sensitivity nor D13 rank-filter
    issue.  Corner matches MMA at rel 3.51e-31 (Re).  MMA wallclock
    222 s.

- **Mass-and-topology diversity batch** (2026-05-15): three additional
  novel oracles probing previously-uncovered axes.  All three pass at
  default precision (no D12-style sensitivity).  Oracle total: 37 → 40.
  - `tools/bench/pentagon_1L_3mass_eps001_*` — 1L pentagon with three
    distinct internal masses (mAsq=1, mBsq=4, mCsq=9 on props 0, 2, 4).
    First 3-mass 1L 5-leg oracle.  Corner matches MMA at rel 6.2e-31.
    MMA wallclock 111 s.
  - `tools/bench/doublebox_blockmass_2L_eps001_*` — 2L doublebox with
    BLOCK 2-mass placement (complement to the interleaved D12 case;
    same family skeleton, mAsq on l1's two rails + mBsq on l2's two
    rails).  Confirms the audit's argument that D12 precision
    sensitivity is interleaved-specific — block converges at default
    `(working_pre, x_order, extra_x_order) = (120, 240, 280)`.
    Corner matches MMA at rel 8.6e-31.  MMA wallclock 424 s.
  - `tools/bench/xbox_crossmass_2L_eps001_*` — 2L non-planar xbox with
    mass on the (l1+l2) cross-rung propagator.  First xbox-with-mass
    oracle; exercises non-planar + mass code paths simultaneously.
    Corner matches MMA at rel 2.4e-30.  MMA wallclock 151 s.
- **Topology-diverse stress batch** (2026-05-15): three new oracle
  triplets exercising axes the previous suite did not cover.  All three
  match MMA at rel ~ 10⁻³⁰ on first run; oracle total rises 34 → 37.
  - `tools/bench/hexagon_1L_eps001_*` — 1L 6-leg hexagon, all-massless,
    pair-only Replacement.  **First 6-leg oracle** (previous max was
    pentabox 2L 5-leg).  Corner `j[hexagon, 1,1,1,1,1,1]` matches MMA
    at rel 1.65 × 10⁻³⁰ (Re).  MMA wallclock 205 s; 31 sampled values.
  - `tools/bench/doublebox_diagmass_2L_eps001_*` — 2L doublebox 4-leg
    with mass on the inter-loop `(l1-l2)` rung propagator.  Novel
    mass-placement axis (existing 2L mass benches put mass within one
    loop or block-symmetric; none on the shared rung).  Corner matches
    MMA at rel 5.10 × 10⁻³² (Re); does NOT exhibit D12-style precision
    sensitivity — converges at default `(working_pre, x_order,
    extra_x_order) = (120, 240, 280)`.  MMA wallclock 174 s.
  - `tools/bench/mercedes_3L_*` — 3L 2-leg Mercedes self-energy
    (triangle of 3 outer rails + 3 inner spokes from a central
    vertex).  Topologically distinct from existing 3L oracles
    (banana_3loop, dotted_3L_banana, bn3mix are parallel-banana
    2-vertex graphs).  Corner matches MMA at rel 4.06 × 10⁻³¹ (Re).
    MMA wallclock 306 s.
- **Doublebox 2L interleaved 2-mass oracle**
  (`tools/bench/doublebox2m_eps001_*`): two-loop doublebox 4-leg with
  cross-loop interleaved 2-mass placement (`mAsq` on `l1` prop 0 AND
  `l2` prop 3; `mBsq` on `l1` prop 1 AND `l2` prop 5), eps = 1/1000.
  Target is the corner `j[doublebox2m, 1,1,1,1,1,1,1,0,0]`.  Matches
  MMA at rel ~ 3.3 × 10⁻³¹ (Re) / Im at noise floor on both sides.
  Regression guard for D12 — see audit §D12.  C++ config requires
  `working_pre=200, x_order=400, extra_x_order=480` (about doubled
  vs the other doublebox benches) because the top-sector diffeq
  matrix has a complex-conjugate pole pair near the NegIm contour.
- **Pentabox 2L 5-leg oracle**
  (`tools/bench/pentabox_2L_eps001_*`): two-loop pentabox with five
  external legs, all-massless internals, eps = 1/1000.  Target is
  the corner `j[pentabox2L, 1,1,1,1,1,1,1,1,0,0,0]`.  Matches MMA at
  rel 9.35 × 10⁻¹² (Re) / 4.67 × 10⁻¹⁰ (Im); precision is
  bounded by the ~120-digit working precision minus the cancellation
  horizon inherent to the 76-master sub-system (intermediate values
  reach 10²⁰⁰⁺).  First 2L 5-leg massless oracle.  MMA reference
  JSON includes 172 sampled values (corner + 171 dotted masters).
- **Pentagon 1L massive oracle**
  (`tools/bench/pentagon_1L_W_mass_*`): one-loop pentagon with one
  massive propagator (`l^2 - msq`, `msq = 1`), eps = 1/1000.  Target
  is the corner `j[pentagon, 1, 1, 1, 1, 1]`.  Matches MMA at
  rel 9.15 × 10⁻³¹ (Re) / 1.27 × 10⁻³⁰ (Im) — full
  working-precision agreement.  First 1L 5-leg single-mass oracle;
  extends the all-massless `pentagon_1L_eps001` along the
  mass-configuration diversity axis.

D11 fix verified by an eight-variant pattern-diverse stress sweep
on 2026-05-14 (pentabox / doublebox / mercedes 3L / banana 4L across
ε, dotted-master, cross-threshold, and 1L/2L/3L/4L axes); all match
MMA at rel ~ 10⁻³⁰ except the inherently cancellation-bound
pentabox-extreme rows.  See audit §D11.

### Fixed
- **Audit divergence D13** — `build_boundary` projection failed when
  `to_complete_explicit` rank-filtered linearly-dependent mass-bearing
  propagators of the boundary sub-family.  The dropped propagator's
  mass scale was absent from the sub-family's reduced context but
  still present in the boundary's Laporta coefficient `lt.coeff`,
  causing `project_mfrac_by_name` to throw `residual var 'mAsq' not
  in dst`.  Fix in `src/pipeline/amfsystem.cpp::build_boundary`:
  apply numeric substitution (`mfrac_substitute(lt_coef_in_fc,
  numeric_q, sub_red_keep_names)`) before the projection, mirroring
  MMA `ReduceBoundary`'s `/. Numeric` at `AMFlow.m:817`.  No
  regression: ctest 547/547 pass, all 40 pre-existing oracles still
  match.  Discovered + resolved 2026-05-15 by the new
  `vtx2_2L_3mass` oracle.  Audit table is now
  **86 🟢 / 0 🟡 / 13 🔴 (12 fixed + 1 D5 out of scope) / 17 ⚪**.
- **Audit divergence D8** (`canonical_boundary_permutation` +
  `canonical_taylor_permutation` re-route SparseGaussian's free
  column): both functions in `src/ode/inf.cpp` previously sorted
  block rows by `int_offsets` before calling `BuildTaylor` /
  `ConstructMatrix` / `SparseGaussian`.  Upstream MMA's `DESolver`
  uses NO row permutation in either `DetermineBlockBoundaryOrder`
  (DESolver.m:705-728) or `CalcTaylor` (DESolver.m:752-790), so the
  C++ sorts re-routed which master ended up holding the free
  (unsolved) column after Gaussian elimination.  The two sorts
  compensated for each other on small / symmetric blocks (all 545
  pre-D8 gtests passed), so D8 only surfaced on
  `banana_4L_mixed`'s 21-master block in region 3 with monotone
  offsets `[0,0,1,2,2,3,3,3,3,4,4,4,4,4,5,5,5,5,6,6,7]`.  Both
  functions now return identity (concatenate `analyze_block`
  output without re-sorting).  `banana_4L_mixed` now matches MMA
  at rel < 1e-30 on all 20 sampled values (was 2/20); all 545
  gtests still pass.  Commit `c668f79`.
- **Audit divergence D9** — three coupled Kira-pipeline alignment
  fixes uncovered while debugging pentabox: (a) `src/ibp/kira_yaml.cpp`
  now substitutes `cfg.numeric_values` into `scalarproduct_rules` and
  propagator masses before writing `kinematics.yaml`, mirroring MMA's
  `SPToSTU /. IBPRule` (Kira/interface.m:53); without it, Masters-mode
  Kira sees symbolic SPs and over-enumerates masters (175 vs 172 for
  pentabox top).  (b) `src/ibp/kira_run.cpp` stops forwarding
  `-s<var>=<val>` for those numeric values now baked into the YAML,
  matching `FilterRules[IBPRule, Prepend[MassScale, ep]]`
  (Kira/interface.m:274); without it, Kira locks up waiting on a
  symbol that no longer exists.  (c) `src/ibp/reduce.cpp` reuses the
  *input* `preferred` / `jpreferred` file at Reduce-mode write time
  instead of the sorted preheat output, mirroring upstream
  `AnalyticReduction`.  Commit `c37f1a0`.
- **Audit divergence D10** — precision-mismatch in
  acb→fmpq rationalization: `acb_real_to_fmpq_local` (and three
  parallel call-sites in `src/ode/path.cpp`,
  `src/pipeline/amfsystem.cpp`, `src/numeric/matrix.cpp`) used
  `rationalize_pre` directly without verifying it fits inside the
  working precision.  When `rationalize_pre` (decimal digits)
  exceeded `working_pre * log10(2) - 5`, the rationalization
  captured binary-representation noise from the acb as a "real"
  rational, leaking 10⁻²⁵ residuals into the diagonal of m_pure and
  propagating into 10⁷-10²¹ errors on masters with all-zero BC.
  Now caps `rationalize_digits` defensively.  Commit `650e369`.
- **Audit divergence D11** — Jordan eigenvector normalization
  mismatch: `fmpq_mat_nullspace_exact` (used internally by
  `jordan_decomposition_exact`'s eigenvector search) returns FLINT's
  integer-cleared null-space vectors — e.g. an eigenvector that
  Mathematica would emit as `(6993/998, 1)` was returned as
  `(6993, 998)`, a 998× scaling.  Without rescaling, downstream
  shearing / leading-Jordan T blocks accumulated huge integer scale
  factors that cascaded through the off-diagonal Sylvester step in
  `to_fuchsian_global`, producing T entries up to 10³⁰⁰⁺ in the
  pentabox 76-master sub-system and overwhelming PSMapRuleS at any
  practical working precision (the corner value came out as
  `7.84 × 10¹² - 3.04 × 10¹³ i` vs MMA reference
  `2.24 - 150.51 i` on the 6-prop subsection that triggered
  isolation).  `fmpq_mat_nullspace_exact` now rescales each
  null-space basis vector so its last non-zero entry is 1, matching
  Mathematica's Eigenvectors / JordanDecomposition convention.
  Regression test
  `JordanTest.JordanEigenvectorsNormalizedToLastEntryOne`.
  Commit `f4f2aee`.
- **Audit divergence D12** — interleaved 2-mass doublebox traced
  2026-05-15 to insufficient Taylor expansion order, not a code
  bug.  The cross-loop mass placement gives the top-sector diffeq
  matrix a complex-conjugate pole pair (`η² + η + 12 = 0`,
  `|η| ≈ 3.46`) which tightens the Frobenius series convergence
  radius for the 4 top-sector masters.  At the doublebox defaults
  `working_pre=160, x_order=200, extra_x_order=240` C++ produced a
  sign-flipped Re and O(0.1) spurious Im; at `200, 400, 480` C++
  matches MMA to all printed digits.  Commit `3da778a` (audit
  rewrite + memo); regression bench
  `tools/bench/doublebox2m_eps001_*` committed in `11ff61a`.
  Audit table is now
  **86 🟢 / 0 🟡 / 12 🔴 (11 fixed + 1 D5 out of scope) / 17 ⚪**.

### Changed
- **D5 ComplexMode reclassified from "deferred indefinitely" to "out
  of scope"** (maintainer decision 2026-05-13).  Complex-valued numeric
  kinematics will not be implemented; the entry-point rejection
  remains, but ROADMAP no longer carries a "Phase 1C" placeholder for
  future work on this.  See [`docs/AUDIT_MMA_PARITY.md`](docs/AUDIT_MMA_PARITY.md) §D5.

### Removed
- ROADMAP "Phase 1C" section (D5 future-work placeholder).
- Stale "deferred" / D3-aborts / Trivial-not-ported limitation notes
  in `docs/FAQ.md` and `docs/USER_GUIDE.md` — D3 was fixed in v1.1.0
  and Trivial was ported in Phase 1A; only the D5 out-of-scope note
  remains.
- Root-level `AUDIT.md` v1.0 snapshot (now a thin pointer to the live
  `docs/AUDIT_MMA_PARITY.md`).

## [1.1.0] — 2026-05-12

Post-v1.0 audit-driven correctness pass + Phase 3 oracle diversity
expansion.  Together they close all silent-wrong-result paths
surfaced by the post-release line-level MMA parity audit
([`docs/AUDIT_MMA_PARITY.md`](docs/AUDIT_MMA_PARITY.md)) and open
the project's five "we-didn't-think-of-this" oracle diversity axes.

After v1.1.0 the audit table reads
**86 🟢 / 0 🟡 / 7 🔴 (6 fixed + 1 D5 deferred indefinitely) /
17 ⚪**, and the oracle suite covers all 5 diversity axes
(loop number L ≤ 4, ≥ 3 kinematic invariants, multi-cut Cutkosky,
mixed-mass, and ε-extremes) with rel ≤ ~10⁻³⁰ across the new
oracles.

### Added — Phase 3 oracle diversity expansion
- 3.F **L=4 banana oracle** (`tools/bench/banana_4loop_eps001_*`):
  4-loop equal-mass banana sunrise at `psq=-3, msq=1, eps=1/1000`.
  Surfaced and locked audit divergence **D7** (see Fixed below).
  Matches MMA at rel 7.2 × 10⁻³¹ / 1.6 × 10⁻³⁰ (C++ 399 s vs MMA
  385 s).
- 3.G **multi-invariant electroweak box oracle**
  (`tools/bench/ewbox_1loop_eps001_*`): 1-loop electroweak box with
  alternating W/Z masses and 4 distinct kinematic invariants
  {`s, t, mWsq, mZsq`}.  Matches at rel ≤ 3.8 × 10⁻³⁰ (C++ 39 s
  vs MMA 62 s).  WW bubble above threshold gives physical imag
  part, locking cross-mass branch-cut handling.
- 3.H **multi-cut Cutkosky oracle** (`tools/bench/cutbanana_4L_eps001_*`):
  4-loop massless cutbanana with all 5 internal lines on-shell
  (5-particle Cutkosky cut).  Matches at rel ≤ 2.2 × 10⁻³⁰
  (C++ 141 s vs MMA 149 s).  5-level subsystem recursion verified
  end-to-end.
- 3.I **mixed-mass oracle** (`tools/bench/bn3mix_eps001_*`):
  3-loop banana with 1 W-massive + 3 massless internal propagators
  — first 3-loop mixed-mass case.  Matches at rel ≤ 2.9 × 10⁻³⁰
  (C++ 56 s vs MMA 96 s).  Exercises scaleless-sub-sector
  detection under mixed mass.
- 3.J **ε-extremes oracles** (`tools/bench/cutbubble_1L_eps{2,10000}_*`):
  cutbubble at eps = 1/2 (D = 3) matches at rel 1.9 × 10⁻⁶⁴ on
  the exact value 1/8; at eps = 10⁻⁴ matches at rel 1.0 × 10⁻³².
  An initial run at eps = 1 (D = 2) revealed an upstream MMA
  AMFlow limitation (DESolver returns partially-symbolic
  `(1/2π) Im[DESolver\`Private\`variables[1, 1]]`); eps = 1/2 was
  used as the next-most-extreme rational that yields a clean
  numeric — documented in the bench source comment.

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
- **D7 — dual-Kira-call master-count divergence.**  `ibp::reduce`
  and `ibp::diffeq` now mirror upstream's `BlackBoxReduce` /
  `BlackBoxDiffeq` two-call pattern: a Masters-mode preheat call
  (sector-wide enumeration via `select_mandatory_recursively`)
  followed by a Reduce-mode call (`select_mandatory_list` for the
  specific targets), both at the same `(rank, dot)`.  Each call
  runs in its own subdirectory (`<work_dir>/masters_preheat/` and
  `<work_dir>/target_reduce/`) because the Kira 2.x release
  refuses to share a `$ReductionDirectory` between the two calls
  (Masters-mode `run_initiate: masters` doesn't register the `-s`
  numeric substitutions, and the subsequent Reduce-mode call
  aborts with `Kira::update_auxiliary_file: Last Kira run set 0
  variables to numeric values, this time you request N`).
  `ibp::diffeq` also no longer calls nested `reduce()` (which
  would re-floor `(rank, dot)` over the derivative integrals — the
  L=4 banana failure mode); it inlines the Reduce-mode Kira
  invocation at `opts_eff`'s `(rank, dot)`, mirroring upstream's
  `AnalyticReduction` which inherits `IBPRank`/`IBPDot` globals
  from `IBPSystem`.  Both functions add a SubsetQ guard on the
  Reduce-mode master file against the Masters-mode sector
  enumeration, mirroring upstream's
  `If[!SubsetQ[masters, str], Abort["inconsistent masters from
  Kira"]]`.  Locked by the L=4 banana oracle.  (Audit D7.)

### Added
- [`docs/AUDIT_MMA_PARITY.md`](docs/AUDIT_MMA_PARITY.md): line-level
  upstream-parity audit report.  86 🟢 verified / 0 🟡 unverified /
  7 🔴 (6 fully fixed, 1 D5 deferred indefinitely) / 17 ⚪ not
  ported.  (Started at 65 🟢 / 21 🟡 / 6 🔴 in v1.0.0; the post-v1.0
  pass closed every 🟡, fixed 6 of the 7 🔴, and surfaced D7 as
  the seventh 🔴 — now fixed too.)
- [`docs/PERFORMANCE.md`](docs/PERFORMANCE.md): wall-clock baseline
  vs MMA on the 12 oracle benchmarks.  Median speedup ~2× on
  Kira-light benches; converges to ~1× on Kira-heavy benches.
- [`docs/ROADMAP.md`](docs/ROADMAP.md): post-v1.0 development plan
  (Phase 1 implementation completeness, Phase 2 oracle expansion,
  Phase 3 ongoing diversification).
- [`tools/bench/run_perf_audit.sh`](tools/bench/run_perf_audit.sh):
  shell driver to reproduce the perf run.
- New oracle benchmarks under `tools/bench/`:
  - `tradcut_phase_2L_eps001_*` — first Tradition-with-cut parity
    test (D3 acceptance gate).
  - `banana_4loop_eps001_*` — Phase 3.F first L=4 oracle (D7 gate).
  - `ewbox_1loop_eps001_*` — Phase 3.G multi-invariant oracle.
  - `cutbanana_4L_eps001_*` — Phase 3.H 5-particle Cutkosky oracle.
  - `bn3mix_eps001_*` — Phase 3.I 3-loop mixed-mass oracle.
  - `cutbubble_1L_eps{2,10000}_*` — Phase 3.J ε-extremes oracles.
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

[1.1.0]: https://github.com/chang18/amflow-cpp/releases/tag/v1.1
[1.0.0]: https://github.com/chang18/amflow-cpp/releases/tag/v1.0
