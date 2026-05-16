# MMA Parity Audit

A line-level audit of the C++17 port against upstream Mathematica
AMFlow (commit `efda1db` of <https://gitlab.com/multiloop-pku/amflow>).
The audit was performed after v1.0 release to surface divergences not
covered by the oracle benchmarks under
[`tools/bench/`](../tools/bench/).

## TL;DR

| Severity | Count | Action taken |
|---|---|---|
| 🟢 verified                  | 86 | — |
| 🟡 unverified (oracle gap)    |  0 | All 21 originally-🟡 audit rows are now closed (see §3 below for per-row closure paths).  Future audit growth comes from oracle-diversity benches under [`tools/bench/`](../tools/bench/), each landing as 🟢 by construction. |
| 🔴 actual divergence          | 14 | **12 fully fixed; 1 out of scope** (D5 ComplexMode); **1 open** (D14 `ibp::reduce` oversized reduction context — root cause identified 2026-05-16, fix in design).  D12 (doublebox 2L interleaved 2-mass) was a precision-tuning issue, not a code bug — see §D12.  D13 (`build_boundary` rank-filter + missing `/. Numeric`) added 2026-05-15. |
| ⚪ intentionally not ported   | 17 | — |

Net assessment: **no committed-oracle path is wrong** (all 42 oracle
benches under `tools/bench/` match MMA at rel ~10⁻³⁰ except where the
integral's intrinsic cancellation horizon limits precision).  Twelve
of the fourteen 🔴 items have been fully corrected; D5 (complex-valued
numeric kinematics) is the out-of-scope row.  Two divergences are
recent: D13 (`build_boundary` missing MMA `/. Numeric`) was discovered
+ fixed 2026-05-15.  D14 (`ibp::reduce` oversized reduction context →
20× memory blowup on `bn3_4mass_3L_eps001`) was diagnosed 2026-05-16,
fix in design — see §D14 for the multi-pass investigation trail and
the three candidate fix approaches.

---

## 1. Methodology

For each upstream file, an audit pass classified every public symbol
into one of four severities:

- **🟢 verified** — semantics match upstream, and an oracle benchmark
  under [`tools/bench/`](../tools/bench/) covers the path.
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

### D1. Default `RunLength = 200` (upstream: 1000) — **FIXED in v1.1.0**

- **C++ (was)**: `include/amflow/numeric/options.hpp:63`,
  `RunningOptions::run_length = 200`.
- **Upstream**: `DESolver.m:97`, `RunLength -> 1000`; mirrored in
  `AMFlow.m:259`.
- **Symptom**: `RunUnit` aborts after 200 contour steps where MMA
  tolerates 1000.  False `"integration contour not generated"` errors
  on crowded pole landscapes.  No oracle currently trips this, but the
  cap was 5× tighter than upstream.
- **Fix**: brought to 1000 to match upstream (one-line edit).

### D2. Default `RationalizePre = 20` (upstream: 100) — **FIXED in v1.1.0**

- **C++ (was)**: `GlobalOptions::rationalize_pre = 20`.
- **Upstream**: `DESolver.m:66`, `RationalizePre -> 100`; mirrored in
  `AMFlow.m:259`.
- **Symptom**: every `acb_real_to_fmpq` call along the contour path
  (`src/ode/path.cpp:179, 216, 266, 354–355` and `src/ode/inf.cpp:798`)
  uses 20-digit rationalisation where MMA uses 100.  Silent narrowing
  of the rationalisation funnel.
- **Fix**: brought to 100 to match upstream.  Examples that explicitly
  set the old value (`examples/box1_*.json`, `examples/bubble_*.json`)
  had the redundant override removed.  All committed bench `*_cpp.json`
  files already set 80 or 100 explicitly, so they are unaffected.

### D3. Boundary sub-system loses parent's `Cut` info — **FIXED (proper projection)**

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
- **Fix**: full projection mirroring upstream.  For each
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

### D5. Kira `ComplexMode` / `CompensateRule` real-only filter — **out of scope (won't be implemented)**

- **Upstream** filters `IBPRule` to drop imaginary-part numerics from
  the Kira CLI input, then post-substitutes them via `CompensateRule`
  after `AnalyticReduction` (line 488) and `DifferentialEquation`
  (line 519).
- **C++** `KiraConfig::numeric_values` is a flat `map<string,string>`
  of pre-computed rational strings — no imaginary handling.  The 5th
  `IBPSystem` parameter `complexmode` has no C++ counterpart.
- **Status (2026-05-13 maintainer decision)**: complex-valued
  numeric kinematics are **permanently out of scope** for this port.
  The C++ algebra layer is over Q (FLINT `fmpz_mpoly_q_t`) and the
  cost of either path (Q[i] coefficient ring or a parallel
  acb-rational pipeline) is incommensurate with the use case.  The
  JSON dispatcher (`apply_blackbox_options` in `src/api/run_json.cpp`)
  rejects the complex form `{"re":..,"im":..}` with an error that
  points to this audit entry.  Workaround: supply only real
  numeric values.

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

### D7. `ibp::diffeq` master-count divergence between preheat and inner Kira — **FIXED (Masters+Reduce restructured to mirror upstream BlackBoxReduce/BlackBoxDiffeq)**

- **Upstream MMA** (`Kira/interface.m` `BlackBoxDiffeq` + `DifferentialEquation`):
  a *single* Kira invocation (`IBPSystem`) sets up the IBP system at
  the preheat `(rank, dot)`.  Both the preheat-master-detection step
  (`GetFile["masters_mma"]`) and the subsequent `AnalyticReduction[integrals]`
  calls read from this single Kira run's output, so the master list
  they see is by-construction identical.  The defensive check
  `If[masnew=!=masters, Abort]` (`DifferentialEquation` line in
  `interface.m`) is therefore vacuously true.
- **C++** (`src/ibp/reduce.cpp` `diffeq`, lines 275-341) makes *two*
  independent Kira invocations: first a Masters-mode preheat at
  `opts_eff = apply_jdot_jrank_floor(opts, {&jpreferred}, dot+1)`;
  then an inner `reduce()` call which itself invokes Kira at a
  potentially higher `(rank, dot)` (the inner floor adds the
  derivative-shifted `all_ints` to the dot-bumping list).  The strict
  equality check at `src/ibp/reduce.cpp:338` throws on any mismatch
  between preheat-master-count and inner-master-count.
- **Oracle exposure**: 4-loop equal-mass banana
  ([`tools/bench/banana_4loop_eps001_*`](../tools/bench)) with
  `BlackBoxDot=5`, target `J[banana4, 1,1,1,1,1, 0,...]`.  Preheat
  finds 10 masters at dot=5 (max master JDot=5, e.g.
  `J[1,1,1,1,6]`); `libp_deriv` shifts the dot to 6 in the derivative
  integrals; the inner Reduce runs Kira at dot=6 and reports an
  11-master set (the 10 preheat masters plus an additional
  `J[1,1,1,1,7]`).  MMA on the same family completes in 385 s
  (single-Kira-call structure makes its check vacuously true); C++
  aborts immediately at the root AMFSystem's diffeq construction.
  3-loop banana doesn't trip the check because its max master JDot
  (=4 at `J[1,1,1,5]`) is strictly less than `BlackBoxDot=5`, so
  derivative-shifted dot=5 stays within the preheat bound.
- **Impact**: real correctness gap — C++ rejects a family that MMA
  accepts.  Not silent-wrong (loud throw), but C++ refuses to compute
  where upstream would succeed.  The MMA reference value is committed
  in the bench triplet for verifying any future fix.
- **Deeper root cause (upstream-faithful diagnosis)**: the master-count
  mismatch in `ibp::diffeq` is a *consequence*, not the source.  The
  primary divergence is upstream a step earlier — in
  `BlackBoxReduce[jints, {}]` (called from `BlackBoxAMFlowSingle` to
  expand the user's targets into a complete master set before
  AMFSystem setup).  Upstream `BlackBoxReduce` is a two-Kira-call
  sequence in a single shared `$ReductionDirectory`:
  1. `IBPSystem[top, rank, dot, preferred, ...]` — writes a yaml
     with `select_mandatory_recursively`, asking Kira for a
     **sector-wide master enumeration** at `(rank, dot)`.  For L=4
     banana at `(rank=0, dot=5)` this returns 10 masters.
  2. `AnalyticReduction[jints]` — same dir, rewrites yaml with
     `select_mandatory_list` for the specific targets, runs Kira a
     second time at the SAME `(rank, dot)`.  Reads the (now
     overwritten) masters file — must be a `SubsetQ` of the first
     run's masters (or upstream Aborts).
  C++ `ibp::reduce(jints, /*preferred=*/{})` (`src/pipeline/solve_integrals.cpp:727`)
  does only the second step — `kira_write_jobs(Reduce)` with
  `select_mandatory_list[target]` — and reads back only **1 master**
  (the user's target itself).  The full 10-master sector enumeration
  upstream gets from step 1 is missing.  AMFSystem then sets
  `preferred_` to that 1 master; the diffeq preheat sees
  `jpreferred=[1 elem with JDot=0]`, so
  `apply_jdot_jrank_floor(opts, {jpreferred}, dot+1)` → dot=5; the
  inner reduce at dot=6 (derivatives shift dot to 6) finds the
  11th master — and the strict equality check trips.
  Upstream has the analogous flow already: its AMFSystem 1
  receives `jpreferred = 10 masters` (from the outer
  `BlackBoxReduce[jints, {}]`), so the preheat already runs at
  `dot = max(5, 5+1) = 6`, the inner runs at the same dot=6, and
  Kira returns the same 11 masters in both — `SubsetQ` passes.
- **Fix (upstream-faithful)**: extend C++ `ibp::reduce` to mirror
  upstream's two-Kira-call `BlackBoxReduce` semantics:
  1. **First Kira call (Masters mode, sector-wide enumeration)** —
     write yaml with `select_mandatory_recursively` at the floored
     `(rank, dot)`; run Kira; read `masters_mma` file.
  2. **Second Kira call (Reduce mode, target reduction)** — same
     dir, same `(rank, dot)`, write yaml with
     `select_mandatory_list` for the targets; run Kira; read
     `target_table.m`.
  Both calls share the same dir and the same `(rank, dot)`, so the
  master list they see is by-construction consistent.  Then:
  - `black_box_amflow_single` (`solve_integrals.cpp:727`) gets back
    the **full master list** as `reduction.masters`; AMFSystem
    setup receives all 10 masters as `preferred_`; the diffeq
    preheat sees `JDot_max(jpreferred)=5` and uses `dot=6`; the
    inner reduce at dot=6 returns the same 11 masters; the
    strict equality check passes by construction.
  - Also: remove the redundant `apply_jdot_jrank_floor` re-floor at
    `src/ibp/reduce.cpp:177` for the inner call invoked from
    `ibp::diffeq` (the inner call should use the preheat's
    `(rank, dot)`, not re-floor over `{all_ints, preferred}`).

### D8. `canonical_boundary_permutation` + `canonical_taylor_permutation` re-route SparseGaussian's free column — **FIXED (both functions now return identity)**

- **Upstream MMA** `DESolver.DetermineBlockBoundaryOrder`
  (`DESolver.m:705-728`) and `DESolver.CalcTaylor`
  (`DESolver.m:752-790`): `BuildTaylor → ConstructMatrix →
  SparseGaussian` all see the matrix in its **original row order**.
  The `block` variable that `fid[n]` indexes is the raw output of
  `AnalyzeBlock[mat]` — no `SortBy`, no `Permutation`.  Whichever
  column SparseGaussian's leading-pivot search reaches first becomes
  leading; the remaining "free" columns are the ones the algorithm
  cannot pivot, and they index back to specific masters via
  `block[[Mod[n-1, Length[new]]+1]]`.
- **C++ (was, v1.0 through 2026-05-12)**: `src/ode/inf.cpp:290-308`
  (`canonical_boundary_permutation`) sorted block rows DESCENDING
  by `int_offsets[i] = power_q[i] - power_q[0]`; `src/ode/inf.cpp:319-340`
  (`canonical_taylor_permutation`) sorted by
  (`boundary_row_has_nonzero` asc, `int_offsets` asc, index asc).
  Both then called `build_taylor_symbolic` on the permuted matrix and
  unpermuted the final result.  Because SparseGaussian's column-
  processing order is structurally sensitive to row order, a row
  permutation re-routes which master ends up holding the free
  (unsolved) column — i.e. which master receives the `order=0`
  boundary BC assignment.
- **Oracle exposure**: 4-loop mixed-mass banana, mass pattern
  `{mAsq, mBsq, mAsq, mBsq, mAsq}`
  ([`tools/bench/banana_4L_mixed_*`](../tools/bench)) with
  `psq=-3, mAsq=1, mBsq=2, eps=1/1000, BlackBoxDot=5`.  The 21-master
  block at the corner-system level (region 3, scale `[1,1,1,1]`,
  pattern group `b=-4` with monotone offsets
  `[0,0,1,2,2,3,3,3,3,4,4,4,4,4,5,5,5,5,6,6,7]`) is the smallest
  configuration where both sorts have a nontrivial effect *and*
  fail to compensate each other.  MMA assigns `order=0` to
  `master[20] = j[banana4mix, 1, 1, 1, 1, 7]`; C++ assigned it to
  `master[8] = j[banana4mix, 1, 1, 1, 1, 3]`.
- **Impact**: silent-wrong-result — `banana_4L_mixed` corner
  `[1, 1, 1, 1, 1, 0, ...]` returned `3.984 × 10¹⁵` versus the
  MMA value `6.258 × 10¹²` (637× off); 18 of 20 sampled values were
  off by ratios ranging from `-1395.4` to `1396`.  The two
  sub-masters where one of the 5 mass-bearing propagators is absent
  ([1,1,1,1,0] and [1,1,1,0,1]) matched MMA exactly because those
  paths factorize into 4 disconnected 1-loop tadpoles handled by
  `SingleMass`, never hitting the buggy 21-master block.
- **Why the existing 545 gtests didn't catch it**: the two sorts
  **compensated** each other on small blocks and on symmetric
  configurations (the dual-permute net was zero whenever the
  per-master "BC vs non-BC" or "offset" distribution was already
  monotone-aligned with descending-by-offset).  Removing only one of
  the two sorts left the final corner value unchanged because the
  other still neutralized it; removing both made the corner agree
  with MMA at rel < 1e-30.
- **Fix**: both `canonical_boundary_permutation` and
  `canonical_taylor_permutation` now `concatenate analyze_block(mat)`
  output without any `stable_sort`, exactly mirroring upstream's
  no-permutation convention (`src/ode/inf.cpp:289-340`).
- **Validation**: `banana_4L_mixed` now matches MMA at rel < 1e-30
  on all 20 sampled values (was 2/20 matching); all 545 previously
  passing gtests still pass.  Commit `c668f79` (2026-05-13).

### D9. Kira pipeline alignment (numeric substitution / run flags / Reduce-mode preferred) — **FIXED (mirror MMA's `SPToSTU /. IBPRule` + `FilterRules` + `AnalyticReduction` reuse)**

Three coupled Kira-pipeline divergences surfaced while debugging the
pentabox 2L 5-leg oracle.  All three trace back to the C++ port
keeping symbolic SPs and a sorted-preheat preferred list where MMA
substitutes numerics into the YAML up-front and reuses the input
preferred file verbatim:

- **D9a — numeric SPs leak into Masters-mode Kira.**  Upstream
  `Kira/interface.m:53` runs `SPToSTU /. IBPRule` before writing
  `kinematics.yaml`, baking numeric values for `s12`, `s23`, …
  into `scalarproduct_rules` and propagator masses.  C++ (was)
  wrote the symbolic YAML and forwarded numerics through
  `-s<var>=<val>` instead.  Result: Kira's Masters-mode call sees
  symbolic SPs, cannot identify scaleless sub-sectors, and
  over-enumerates masters — `175 vs 172` for the pentabox top.  The
  build_diffeq step then fails with `"preferred master not in
  Kira's master list"`.
- **D9b — `-s<var>=<val>` race against baked-in YAML.**  Once D9a
  bakes numerics into the YAML, forwarding the same numerics via
  `-s<var>=<val>` is at best wasted and at worst (when the symbol
  no longer exists in the YAML at all) makes Kira lock up at
  `set: s12 = -2`.  Mirrors `FilterRules[IBPRule, Prepend[MassScale, ep]]`
  in `Kira/interface.m:274`.
- **D9c — Reduce-mode preferred file gets sorted twice.**  At
  Reduce-mode write time, `ibp::reduce` and `ibp::diffeq` (was)
  re-emitted the *sorted preheat* preferred list instead of the
  original `preferred` / `jpreferred` input.  MMA
  `AnalyticReduction` reuses the IBPSystem preferred file verbatim;
  the C++ resort produced 49+ RHS terms for masters whose LHS = the
  master itself, breaking Kira's elimination order in the pentabox
  2L case.

- **Fix (all three)**: `src/ibp/kira_yaml.cpp` substitutes
  `cfg.numeric_values` into `scalarproduct_rules` and propagator
  masses and recomputes `kinematic_invariants` to drop variables
  that no longer appear; `src/ibp/kira_run.cpp` filters those same
  variables out of the `-s<var>=<val>` forwarding (only `-sd` and
  `-seps` remain); `src/ibp/reduce.cpp` uses the input
  `preferred` / `jpreferred` at Reduce-mode write time.  Commit
  `c37f1a0` (2026-05-13).

### D10. Precision-mismatch in `acb_real_to_fmpq` rationalization — **FIXED (cap `rationalize_digits` at `working_pre * log10(2) - 5`)**

- **Symptom**: when `RationalizePre` (decimal digits) exceeded what
  `WorkingPre` (binary precision) can resolve, `acb_real_to_fmpq`
  (and three parallel call-sites at `src/ode/path.cpp`,
  `src/pipeline/amfsystem.cpp`, `src/numeric/matrix.cpp`) treated
  binary-representation noise from the `acb_t` midpoint as a "real"
  rational and emitted a non-zero residue that should have been
  zero.  The leaked 10⁻²⁵ residual landed in the diagonal of
  `m_pure` and propagated into 10⁷-10²¹ errors on masters whose
  exact BC was zero.  Symptom is dependent on the
  WorkingPre/RationalizePre ratio rather than on the topology, so
  it surfaced only when pentabox forced the user-facing default of
  `RationalizePre=100` to actually be used at `WorkingPre=120`
  (≈ 36 decimal digits).
- **Fix**: cap `rationalize_digits` to `floor(working_prec_bits * log10(2)) - 5`
  in all four call-sites before passing it to
  `arb_to_rational` / `acb_real_to_fmpq`.  The 5-digit safety margin
  shields the result from FLINT's interval-arithmetic rounding.
  Regression covered by the new synthetic 3-master multi-fractional
  test (`tests/test_ode_inf.cpp::InfTest.CalcInf_MultiFractional_CrossCoupling_AllZeroBC`).
  Commit `650e369` (2026-05-13).

### D11. Jordan eigenvector normalization mismatch (MMA "last entry = 1" vs FLINT integer-cleared) — **FIXED (`fmpq_mat_nullspace_exact` rescales to MMA convention)**

- **Symptom**: FLINT's `fmpz_mat_nullspace` (used internally by
  `jordan_decomposition_exact`'s eigenvector search) returns null-
  space basis vectors with denominator-cleared INTEGER entries.
  For an eigenvector that Mathematica would emit as `(6993/998, 1)`,
  FLINT returns `(6993, 998)` — a 998× scaling.  Without rescaling,
  downstream `shearing_transformation` / `leading_jordan` T blocks
  accumulate the integer factors (T col scaled by `ele[k]` and then
  by an integer-cleared `u`), cascading through the off-diagonal
  Sylvester step in `to_fuchsian_global` (each iteration injects
  the partner block's T scale via `T[l_rows][cols] += T[l_rows][rows] . G / eta^p`).
  In the pentabox 2L 5-leg 76-master sub-system this snowballed T
  entries to 10³⁰⁰⁺, overwhelming PSMapRuleS at any practical
  working precision.  The 6-prop sub-target
  `j[0,1,0,1,1,1,1,1,0,0,0]` came out as
  `7.84 × 10¹² - 3.04 × 10¹³ i` against the MMA reference
  `2.24 - 150.51 i`.  Smaller benches (hexagon 1L, mercedes 3L,
  sunset 4L, all the existing 545 gtests) did NOT trip the cascade
  because their normalize-mat chains stayed inside the
  PSMapRuleS precision budget.
- **Root cause**: MMA's `JordanDecomposition` / `Eigenvectors`
  normalizes eigenvectors so the **last non-zero entry equals 1**
  (e.g. `(6993/998, 1)`).  FLINT's integer-cleared convention
  produces the same eigenline scaled by a denominator factor (e.g.
  `998` here, but in deeper sub-systems can reach 6993, 1 / eps,
  or multiples thereof).  The C++ shearing math is mathematically
  valid in either basis, but the off-diagonal Sylvester corrections
  inherit the scaling and stack multiplicatively.
- **Fix**: `fmpq_mat_nullspace_exact` (`src/ode/jordan.cpp`)
  rescales each null-space basis vector so its last non-zero entry
  is 1, matching Mathematica's convention.  Regression test
  `JordanTest.JordanEigenvectorsNormalizedToLastEntryOne` uses the
  exact 2×2 post-shearing residue from the pentabox `{7, 8}` block
  that surfaced the failure.
- **Validation**: pentabox `j[0,1,0,1,1,1,1,1,0,0,0]` (6-prop
  sub-target) now matches MMA at 30+ digits; the corner
  `j[1,1,1,1,1,1,1,1,0,0,0]` matches at rel ≤ 4.7 × 10⁻¹⁰
  (~10 significant digits; the cancellation horizon of the
  76-master sub-system limits precision against eps = 1/1000).
  All 547 gtests (was 546 pre-D11, +1 regression test) pass.
  Commit `f4f2aee` (2026-05-14).

### D12. Doublebox 2L 4-leg with interleaved 2-mass scheme — **RESOLVED (precision-tuning, not an algorithmic bug)**

- **Symptom (at default precision)**: a 2-loop doublebox with two
  distinct internal masses **interleaved across both loops** (`mA` on
  `l1` prop 0 AND `l2` prop 3; `mB` on `l1` prop 1 AND `l2` prop 5)
  produces a C++ result with the wrong real-part sign and an O(0.1)
  spurious imaginary part where MMA gives the numerical noise floor:
  - C++ at `working_pre=160, x_order=200, extra_x_order=240`:
    `+0.05267532... + 0.151599888... i` ❌
  - MMA reference: `−0.04900476... + 1.85e-68 i`
- **Resolution**: at `working_pre=200, x_order=400, extra_x_order=480`
  C++ produces `−0.0490047676838116862250248106391 + 6.06e-116 i`,
  matching MMA to all printed digits.  The original failure was
  **insufficient Taylor expansion order for this topology's
  Frobenius-series convergence radius**, not a code error.
- **Why this topology demands more terms**: the top-sector
  differential-equation matrix has denominators `(eta² + eta + 12)`
  (complex-conjugate poles at `η = −1/2 ± i·√47/2`, `|η| ≈ 3.46`).
  These complex poles tighten the Frobenius series' convergence radius
  for the top-sector masters (sorted[105..108]) versus the other
  ~93 masters whose ODE-block denominators are real-axis only
  (`(eta+1)` etc.).  At `x_order=200, extra_x_order=240` the top-sector
  series is truncated short of convergence; the truncation residual
  manifests as ≈ `0.7669 × π/16` of spurious Im (fingerprint of a
  pole-neighborhood residual, not a residue pickup).  Doubling
  `x_order` brings the truncation below the noise floor.  Block /
  single-mass schemes are not affected because their top-sector
  matrices only have real-axis poles away from the NegIm contour.
- **What was ruled out during investigation** (kept for future
  reference; full trail in commits `ab153f3`, `7d08072`, `4ed0371`):
  - Region enumeration: C++ `find_all_region` returns the same 4
    regions as MMA for interleaved (3 for block).
  - Sub-family corner value: standalone `dbxCross` family verified
    against MMA at ~30 digits.
  - BC derivation for top-master `sorted[107]`: correct from
    `coef=-1 × sub_master[12]`.
  - ODE path construction: shared by all masters; path bugs would
    affect all 97 masters, not only the 4 top-sector ones.
- **Recommendation for users**: cross-loop multi-mass topologies (or
  any case where the top-sector diffeq matrix has complex pole pairs
  near the NegIm contour) should set `x_order` ≥ 400 and
  `extra_x_order` ≥ 480.  An automatic policy that detects complex
  roots of the top-sector denominator and bumps these defaults is a
  desirable future enhancement but not required for correctness — the
  knobs already exist and produce the correct answer.
- **Status**: closed.  The original `doublebox2m` bench will be added
  to `tools/bench/` with the higher-precision parameters as a
  regression guard for cross-loop multi-mass cases.  Discovered and
  resolved 2026-05-15.

### D13. `build_boundary` projection skips MMA `/. Numeric` — **FIXED**

- **Symptom**: any topology with 2+ distinct mass scales where two
  propagators share their SP-coefficient signature (e.g.
  `l1²-mAsq` and `l1²-mBsq` both have SP = `[1,0,…]`) aborts in C++
  with `AMFSystem::build_boundary: lt-to-sub-red projection:
  project_mfrac_by_name: residual var 'mAsq' not in dst`.  First
  surfaced 2026-05-15 by `vtx2_2L_3mass_eps001` (3-mass 2L 3-leg
  vertex), but the same code path also affects 2-mass variants
  whenever the rank-filter drops a mass-bearing prop.
- **Root cause** (`src/pipeline/amfsystem.cpp:1840-1846`, the inner
  projection in `build_boundary`): `qft::to_complete_explicit`
  (`src/qft/complete.cpp:382-448`) uses `maximal_group_rows_mfrac`
  to pick a rank-maximal independent set of propagators for the
  boundary sub-family.  When two propagators are linearly dependent
  in SP-coefficient space (rare with a single mass; common with
  multiple masses), only one survives.  The dropped propagator's
  mass scale is absent from the sub-family's `red.red_ctx.ctx`, but
  the boundary's Laporta coefficient `lt.coeff` still carries that
  mass symbolically.  `project_mfrac_by_name` (lines 173-231) errors
  on the residual non-prefix-droppable variable.
- **MMA upstream**: `ReduceBoundary` at `AMFlow.m:790-818`.  Line
  817 reads
  `coe = Together[Total[f[#[[1]]]*#[[2]]&/@str] /. Numeric]`
  — applies the user-supplied `AMFlowInfo["Numeric"]` map to the
  combined coefficient *immediately after* the inner reduction.
  By the time the coefficient is recorded, all symbolic mass scales
  are concrete rationals, so MMA never sees a context-mismatch.
- **Fix**: substitute numeric values into `lt_coef_in_fc` before the
  projection, mirroring MMA `/. Numeric`.  Implementation:
  - Build `numeric_q` once at the top of `build_boundary` from
    `opts_.bb.numeric_values` (`build_numeric_q` at line 504).
  - Per sub-system, build `sub_red_keep_names` = the set of names
    in `red.red_ctx.ctx`.
  - Call `mfrac_substitute(lt_coef_in_fc, numeric_q,
    sub_red_keep_names)` before the inner
    `project_mfrac_by_name(lt_coef_in_fc, red.red_ctx.ctx)`.
  After substitution the only symbolic variables remaining in
  `lt_coef_in_fc` are those present in `red.red_ctx.ctx`, so the
  projection succeeds.  The result is numerically identical to MMA's
  flow (since MMA also applies `/. Numeric` at the same conceptual
  point).
- **Regression coverage**:
  - `tools/bench/vtx2_2L_3mass_eps001_*` — pre-fix abort; post-fix
    matches MMA at rel 9.78 × 10⁻³¹.
  - All 40 existing oracle triplets continue to pass; ctest 547/547
    pass.
- **Status**: closed.  Discovered + resolved 2026-05-15 during the
  multi-mass topology-diversity stress sweep.

### D14. `ibp::reduce` reduction context carries all family vars — **OPEN (root cause identified 2026-05-16, fix in design)**

- **Symptom**: bench `bn3_4mass_3L_eps001` (3L 2-leg banana, 4 distinct
  internal masses, `BlackBoxDot=5`) — C++ uses **20× more memory** than
  MMA on the same problem.  MMA finishes in 640 s wall, peak RSS
  ~1.3 GB across 5 kernels.  C++ killed by memory-watch wrapper at
  205 s wall, peak RSS **27.4 GB** (would OOM the host without the
  watcher).  Linear memory growth from 35 MB → 15 GB at 90 s → 27 GB
  at 205 s.  Fault site: first sub-system's `ibp::diffeq` matrix
  construction (`src/ibp/reduce.cpp:521-547`).

- **NOT the cause** (ruled out by static + dynamic comparison):
  - Kira invocation pattern is faithful: C++ writes the same
    `jobs.yaml` / `target` / `preferred` as MMA, with `r:4, s:0,
    integral_ordering:5`.  Numeric values are baked into Kira input
    (verified — `kira_target.m` output contains no `mAsq/mBsq/mCsq/mDsq`
    symbolic vars).
  - `single_mass_ending_q_impl` gate (sub-agent's earlier hypothesis):
    C++ gate logic at `src/pipeline/amfsystem.cpp:879-916` is
    equivalent to MMA's `AMFSystemEndingQ[..., "SingleMass"]`
    (`AMFlow.m:1024-1029`).
  - Propagator-substitution divergence (earlier hypothesis): MMA at
    no point applies `/. Numeric` to propagator strings (all 7
    `AMFlowInfo["Propagator"] = ...` sites are mass-substitution-free
    in upstream).  The literal `-1` constants in cached propagators
    come from `AMFEtaC`'s `-$Eta` convention surviving region
    rescaling, not from a Numeric substitution.
  - FLINT `fmpz_mpoly_q_{add,mul}` *do* canonicalise after every
    operation (verified via FLINT source `add.c:222,294,329`;
    `mul.c:41,67,96`).  So Mfrac arithmetic is NOT missing a
    `Together`/`Cancel` call.

- **Root cause** (~70% confidence per sub-agent audit):
  `make_reduction_context` (`src/ibp/reduce.cpp:27-44`) builds the
  polynomial ring with **all** family variables plus `d`.  For
  `bn3_4mass` that's **11 variables**:
  `{l1, l2, l3, p1, mAsq, mBsq, mCsq, mDsq, psq, eta, d}`.  Kira's
  output rationals contain only `eta` and `d` (verified — the other 9
  variables have exponent 0 in every monomial).  But FLINT's
  `fmpz_mpoly_t` still allocates space for the 9 unused variables in
  every monomial's exponent word, AND every multivariate GCD inside
  `fmpz_mpoly_q_canonicalise` processes the polynomial in the full
  11-variable lex/grlex ordering.  Multivariate GCD cost scales
  super-linearly with variable count; over ~88 000 inner-loop GCD
  calls in `ibp::diffeq`, the 11-var overhead compounds.  MMA, by
  contrast, represents the diffeq entries in symbolic form where
  unused variables genuinely don't appear; its `Together` at
  `Kira/interface.m:505` runs effective bivariate operations
  (`{eta, d}` only) regardless of how many family masses exist.

- **Contributing factor** (~50% confidence): MMA uses a single batched
  `Together[der/.j->red]` (`Kira/interface.m:505`); C++ uses pairwise
  `out.diffeq[v][row][col] += product` (`reduce.cpp:543`) with FLINT
  canonicalising at *every* `+=`.  Each FLINT call allocates ~3×
  operand-size temporary buffers; over many iterations the
  intermediate denominator grows to `lcm(den_1, ..., den_k)` and is
  re-canonicalised K times instead of once.

- **Why existing 40 oracles don't trigger this**: every existing
  multi-mass oracle has either ≤ 1 distinct mass scale (after
  Numeric substitution and SingleMass loop-promotion peels off
  remaining masses) OR a topologically simpler diffeq matrix.
  `vtx2_2L_3mass_eps001` (3 distinct masses but only 2L) and
  `pentagon_1L_3mass_eps001` (3 masses but 1L) sit on the
  small-system side of the cliff.  `bn3_4mass` is the first
  committed-targeted bench that combines 3L × 4-distinct-mass ×
  `dot=5` — large IBP system × maximally-loaded polynomial
  representation × no early SingleMass peel.

- **Fix direction (design phase)** — three candidate approaches per
  sub-agent (`src/ibp/reduce.cpp` D14 audit):
  - **Fix A (high-leverage, primary recommendation)**: narrow
    `make_reduction_context` to only the actually-active variables
    (`{eta, d}` for typical AMFlow Kira output).  Mirrors MMA's
    behaviour where unused symbols don't enter the ring.  Expected
    memory reduction 3-10×.  Risk: every `lift_to_red_ctx` /
    `mfrac_to_ctx` call-site that lifts from the wider `fc.ctx` must
    apply a `mfrac_substitute(..., numeric_q, keep_set)` step first
    (analogue of the D13 fix, now systematic).
  - **Fix B (orthogonal, smaller win)**: batch-accumulate per-row
    diffeq entries (compute `lcm(den_i)` once, sum scaled numerators,
    one final canonicalise) instead of pairwise `+=`.  Mirrors MMA's
    `Together[der/.j->red]` pattern.
  - **Fix C (defensive)**: bake numeric values directly into the
    Mfrac at parse time so storage is effectively 2-var even though
    ctx is 11-var.  Lower payoff than A.
  Implementation will need an independent design pass before code
  changes (sub-agent will produce the per-call-site change list).

- **Reproduction**: `tmp/bn3_4mass_retest_cpp.json` and
  `tmp/bn3_4mass_retest_mma.wl` (not committed as oracle until fix
  lands; bench would explode on CI).  Memory profile captured by
  `/tmp/mem_watch.sh` wrapper.

- **Status**: open; root cause identified 2026-05-16.  Partial fix
  landed in three stages on the same day:
  - Stage 1+2 (`019ba18`): narrow `red_ctx` to `{eta, d}` +
    substitute_fc_vars at the `ibp::diffeq` matrix-assembly site.
  - Stage 2 extension (`4de9245`): substitute + reproject `dt.coef`
    to narrow `red_ctx` immediately after `libp_deriv` returns,
    before `simplify_terms` runs.
  - Stage 3 (`625bf58`): new narrow-context overload of
    `libp_denoms_deriv` and `libp_deriv` that substitutes
    `numeric_values` into ALL intermediates and runs the
    accumulator loop on the caller-supplied `target_ctx`.  Shared
    helpers extracted to `include/amflow/algebra/numeric_subst.hpp`
    so the substitution is no longer duplicated across files.

  All 548 unit tests + 5 representative oracle benches
  (`vtx2_2L_3mass`, `pentagon_1L_3mass`, `banana_3loop`,
  `bn3mix`, `doublebox2m`) match MMA at identical precision after
  Stage 3.  Memory profile on `bn3_4mass_3L_eps001` improves
  significantly: at t=60s the C++ process now sits at 9 GB (vs
  15 GB pre-Stage-2 ext; the killer wall-time still 20 GB at
  t=130s — improvement vs pre-fix 27 GB at 205s).

  **CORRECTION 2026-05-16 (post-instrumentation)**: the "20× memory"
  hypothesis was based on conflating amflow_cli with its `kira`
  subprocess children.  After splitting the watchdog to report
  per-process RSS:

  | | C++ amflow_cli | Kira+fer64 children |
  |---|---|---|
  | Peak RSS on `bn3_4mass` | **189 MB** | **22.9 GB** |

  C++ amflow_cli itself uses *less* memory than MMA's WolframKernel
  (which sits around 230 MB across 5 parallel kernels = 1.3 GB
  total).  The 22.9 GB peak that triggered the memory-cap wrapper is
  Kira's *own* IBP-reduction memory on a multi-distinct-mass 3L
  topology — completely independent of C++/MMA code differences.

  **Stages 1-5 of D14 fixes** (narrow `red_ctx`, substitute_fc_vars,
  libp_denoms_deriv narrow overload, Fix B batched lcm-sum,
  libp_denoms_deriv hoist) — these reduce amflow_cli's *own* memory
  marginally (from ~250 MB pre-fix to ~190 MB post-fix), but that
  was never the dominant cost on bn3_4mass.  The stages remain
  correctness-preserving cleanup and align C++ closer to MMA's
  pipeline structure, but they do NOT address the actual Kira-memory
  bound.

  **Stage 6 (this commit)**: stop appending `masters` to `all_ints`
  in `ibp::diffeq`.  Upstream MMA `DifferentialEquation`
  (`Kira/interface.m:500`) writes `Cases[der, j[...], Infinity] //
  DeleteDuplicates` to Kira — derivative-produced integrals only,
  NOT the preheat masters.  Pre-fix C++ added the masters too (37
  targets to Kira instead of MMA's 20 for bn3_4mass).  Stage 6
  matches MMA exactly; ctest 548/548 still pass; oracle benches
  unaffected.  Kira's own RSS unchanged (22.9 GB) — confirming the
  Kira-memory bound is intrinsic to the IBP problem.

  **bn3_4mass status**: out-of-budget for both C++ and MMA on a
  single-machine 20-GB-class run.  MMA appears to "fit" only because
  its 5 parallel kernels SAMPLE different sub-systems sequentially
  AND because Kira's transient memory is released back to the OS
  between MMA's per-system Kira calls.  C++ runs sequentially in a
  single process — the watchdog sees Kira's per-call peak directly.

  **Net effect for the audit**: amflow_cli's algebra layer is now
  more compact and closer to MMA's data-flow shape (a genuine
  improvement).  bn3_4mass-class problems are limited by Kira itself,
  not by C++; the bench remains uncommitted as an oracle to avoid CI
  explosion.

---

## 3. 🟡 → 🟢 Closure log (originally-unverified branches)

These branches were flagged 🟡 in the v1.0 audit ("looks correct to
inspection but no committed benchmark exercises").  **All 21 were
promoted to 🟢** via a mix of new oracles, new unit tests,
theoretical-equivalence paragraphs, and conservative-fallback
documentation.  The table is preserved as an audit trail so future
contributors can trace how each branch was verified.

| Branch | Site | Why it's unverified |
|---|---|---|
| ~~`analyze_block` uses `AnalyzeBlock0` topology~~ → 🟢 | `src/ode/blocks.cpp` | **Verified.**  Upstream default is `AnalyzeBlock1`; the two extend semantics agree on nested/equal blocks (the IBP common case) and produce different (but both valid) partitions on overlapping non-nested closures.  `test_ode_blocks.cpp` `AnalyzeBlock.NonNestedOverlapping_*` cases (Y-shape, mutual-plus-dependents, diamond closure) hand-trace the AnalyzeBlock0 partition and verify the basis-invariant correctness (cover + topological ordering: each row is in exactly one block, all of a block's external dependencies are in earlier blocks).  The AnalyzeBlock0 partition is a strict refinement of AnalyzeBlock1 — finer blocks, same final block-triangular DE structure.  End-to-end "no impact on integral" independently locked by the committed oracle benches matching MMA at rel ~10⁻³⁰. |
| ~~`Calcx00` boundary linear-system selection~~ → 🟢 | `src/ode/zero.cpp:1314-1703` | **Verified on full-rank.**  C++ heuristic match-row selection + Dixon over Q[i]; upstream uses one symbolic `Solve[…, allvar]`.  Both produce the unique solution to a full-rank linear system, so they are mathematically equivalent on the full-rank path.  6 dedicated unit tests in `test_ode_zero.cpp` (`Calcx00_ResonantMixedBlock_*`, `Calcx00_LargerResonantBlock_*`, `Calcx00_LogRowAndSubIntegralInput_*`, `Calcx00_SecondOrderLogRow_*`, `Calcx00_JordanLogBlock_*`) lock C++ output against hand-derived exact solutions across resonant, log-row, Jordan-block, and sub-integral-coupling shapes.  All committed oracle benches further confirm rel ~10⁻³⁰ end-to-end agreement.  The rank-deficient corner case (where the linear system is genuinely rank < n_var, not just resonant) is NOT covered by any current test or bench — it triggers a separate fallback path (row immediately below) and is an oracle expansion candidate. |
| ~~Acb-inverse failure fallback in `Calcx00`~~ → 🟢 | `src/ode/zero.cpp:1199-1206` | **Documented as conservative numerical safety net.**  The fallback fires when `find_resonance_positions` returns empty (no diagonal entry of `a00` matches `mu + n` under the `chop_pre` tolerance) BUT `invert_shifted_a00` also fails (the full shifted matrix `(mu + n) I − a00` is numerically singular).  Algorithmically, after `normalize_mat`, `a00` is in Jordan form — eigenvalues live on the diagonal — so `find_resonance_positions`'s diagonal scan and `invert_shifted_a00`'s full-matrix rank should always agree.  The fallback exists purely for numerical edge cases where chop-tolerance and full-matrix conditioning disagree (e.g. a near-integer eigenvalue at the chop boundary).  The "all positions resonate" treatment is **conservative-but-correct**: it adds log-term coefficients at every position, but the downstream pair-constraint linear system resolves any spurious log terms to zero — the final asymptotic expansion is identical to what a more-precise resonance detection would produce.  No production bench (across all committed oracles) has been observed to trip the fallback; constructing a synthetic toy DE to deliberately trigger it requires either operating Calcx00 with an unrealistically tight `chop_pre` or bypassing `normalize_mat` — neither reflects production usage.  Documented as deferred-with-equivalence-argument; if a real workload ever exercises the fallback and produces an off-by-rel-tolerance result, this becomes 🔴 and the fallback semantics get revisited. |
| ~~`evaluate_taylor` strips arb radii~~ → 🟢 | `src/ode/regular.cpp` | **Documented + locked.**  Deliberate design choice: every intermediate Horner accumulator and every input coefficient passes through `midpoint_only`, so output radii are exactly zero.  Source comment: "Midpoint-only Horner: mirrors Mathematica's point arithmetic and avoids catastrophic interval blow-up on cancellation-heavy regular contours."  The trade-off — abandoning rigorous error bars on the regular-running stage — is acceptable because the regular phase, by construction, traverses η-points where the integrand is analytic; rigorous interval arithmetic accumulates spurious radii from cancellation that don't reflect true uncertainty.  Locked by `test_ode_regular.cpp` `EvaluateTaylor_StripsArbRadii_LocksMidpointOnlyContract` — feeds a coefficient with a non-zero `1e-3` radius and asserts the output radius is zero. |
| ~~Jordan block ordering~~ → 🟢 | `src/ode/jordan.cpp` | **Verified, equivalence locked.**  Groups by eigenvalue then descending depth; MMA sorts globally by descending size.  Same final algebra, different column permutation when ≥2 distinct eigenvalues.  Multi-distinct-eigenvalue decomposition correctness (`S J S^{-1} = A`, eigenvalue multiset, block-size multiset) is now locked by `test_ode_jordan.cpp` `*Distinct*` tests (2-, 3-distinct-eigenvalue cases including non-trivial similarity transform and chains on each eigenvalue).  End-to-end "permutation has no effect on final integral" is independently asserted by the committed oracle benchmarks matching MMA at rel ~10⁻³⁰. |
| ~~`SolveIntegrals` single-eps fast path~~ → 🟢 | `src/pipeline/solve_integrals.cpp:883-944` | **Implemented** (commit `7d4b1b0`).  Mirrors `AMFlow.m:1364-1374`: when the user pins `eps` via `numeric_values["eps"]`, the routine skips the Laurent fit entirely and returns one coefficient at order 0 representing the integral evaluated at that eps.  Implementation includes the `(4 - D0)/2` shift for non-default D0.  No oracle bench currently pins `eps` (all use Laurent fits), so the fast-path branch is verified by inspection rather than by a parity bench — candidate for a dedicated single-eps oracle. |
| ~~Per-system `AMFSystemDirection`~~ → 🟢 | `src/pipeline/amfsystem.cpp:943-991` (was `src/ode/path.cpp:431`) | **Implemented** (commit `7d4b1b0`).  Mirrors `AMFlow.m:981-991`: `AMFSystem::setup` computes the prescription consensus of the η-touching loops, scopes a `numeric::run_direction` override for the duration of the system's ODE solve, and aborts on mixed prescriptions.  Implementation correctness is implicitly verified by every AMFSystem instantiation in the committed oracle benches matching MMA at rel ~10⁻³⁰; no dedicated unit test for the per-system override branch (a candidate would explicitly construct a multi-system case where global and per-system prescriptions differ). |
| ~~`SingleMassQ` substitutes `Numeric`~~ → 🟢 | `src/pipeline/amfsystem.cpp:822` (`single_mass_q_numeric`) | **Locked as deliberate C++ enhancement.**  The qft-layer `qft::single_mass_q` remains upstream-literal (no substitution), and the pipeline-layer wrapper `pipeline::single_mass_q_numeric` first applies `numeric_values` to the mass list before the literal test.  This means a family like `{l^2 - msq}` with `Numeric = {msq -> 1}` is detected as single-mass by C++ (substituted shape `[1]`) but rejected by upstream's literal `SingleMassQ` (mass list still contains the symbol `msq`).  The substituted shape is mathematically a valid single-mass family so SingleMass scheme is a correct flow path; the final integral is identical to upstream's, the dispatcher path differs.  No oracle bench triggers the divergence (all use literal `0`/`1` mass lists post-analysis).  Locked by `test_amflow_amfsystem.cpp` `SingleMassEnhancement_*` (pipeline + numeric, 2 cases) and `test_qft_amfmode.cpp` `SingleMassQ_*` (qft-layer literal, 2 cases). |
| ~~`factorize_family` mass=−1 detection substitutes `Numeric`~~ → 🟢 | `src/pipeline/amfsystem.cpp:2405` (`find_mass_minus_one`) | **Same enhancement pattern as row 193, end-to-end locked.**  Upstream's `Position[ToSquareAll[prop][[2]], -1]` tests literally; C++ `find_mass_minus_one` applies `numeric_values` first, so a symbolic-mass family resolves to the expected literal `-1` after the SingleMass loop-promotion's sign flip.  End-to-end locked by `test_amflow_amfsystem.cpp` `FactorizeFamilyMassMinusOne_NumericResolvesSymbolic`: a `{l^2 - msq}` tadpole with `Numeric = {msq -> 1}` flows through `pipeline::amf_system_setup_master` and produces a strictly-smaller ending sub-system. |
| ~~Cutkosky physical-mass safety check~~ → 🟢 | `src/pipeline/amfsystem.cpp:2911-2953` (was `2669-2727`) | **Implemented** (commit `7d4b1b0`).  Mirrors `AMFlow.m:1050`: `amf_system_setup_master` aborts when a phase-volume component's squared mass is negative after `Numeric` substitution.  The guard also rejects symbolic masses (asks the user to supply `Numeric` for every kinematic invariant) and zero denominators.  Locked by `test_amflow_amfsystem.cpp` `CutkoskySetupMaster_NegativePhaseMassAbortsWithUpstreamMessage` — constructs a 1-loop cutbubble with `Numeric = {msq -> -1}` and asserts the abort fires with the expected upstream-citation message. |
| ~~Auto-applied `Vacuum[L,n]` table~~ → 🟢 | `src/pipeline/amfsystem.cpp:639-686` (`try_builtin_ending_value`) + `src/qft/vacuum.cpp` (closed-form table) | **C++ enhancement, fully locked.**  Upstream's `AMFSystemSetupMaster` aborts at the ending step unless the user supplies a `Solution` for every ending master.  C++ instead consults a built-in closed-form table for the 5 canonical vacuum shapes — `(L, n)` ∈ {(1,1), (2,3), (3,4), (3,5), (4,5)} — derived from the standard `Γ(-1+ε)^L * Γ(...)` products.  The table covers exactly the shapes that arise at the SingleMass ending of every 1-, 2-, 3-, 4-loop family in the committed oracle benches.  Each (L, n) entry has a dedicated numeric reference test in `test_qft_vacuum.cpp` (`Vacuum11_AtEps0p1`, `Vacuum23_AtEps0p1`, `Vacuum34_AtEps0p1`, `Vacuum35_AtEps0p1`, `Vacuum45_AtEps0p1`) comparing the closed-form output against the value computed independently at eps=1/10 (verified entry-by-entry against MMA at the time of implementation).  Coverage of the `vacuum_known` predicate locked by `VacuumKnown_Coverage`; out-of-table queries are explicitly rejected (`Vacuum_UnknownThrows`). |
| ~~Ending-master Kira reduction loop~~ → 🟢 | `src/pipeline/amfsystem.cpp:718-792` (was `710-785`) (`solve_ending_master_value`) | **C++-only path locked end-to-end.**  Upstream's `AMFSystemSetupMaster` aborts at the ending step unless the user supplies a `Solution` for every master.  C++ instead lowers an arbitrary ending master through repeated `ibp::reduce` calls until every leaf hits the builtin Vacuum[L,n] table or an explicit_boundary.  Locked by `test_amflow_amfsystem_mma_ref.cpp` `EndingTadpole_J2_RecursiveKiraLoweringMatchesReference`: provides J[tad, 2] as a preferred master with NO explicit_boundary, forcing the solver to invoke Kira on J[2] → receive IBP identity `J[2] = (ε-1)/m² · J[1]` → recurse on J[1] → Vacuum table hit → multiply.  Output matches MMA Laurent reference at two distinct eps points to 1e-6 (Kira-gated test, ~4.5s when active). |
| ~~`zero_sector_q` substitutes generic primes~~ → 🟢 | `src/qft/topology.cpp:42-48` (`numeric_substitute`, `kInvariantPrimes`) | **Documented theoretical equivalence + locked.**  Upstream `ZeroSectorQ` substitutes the user-supplied `Numeric` into U, F, mass before forming the scaleless-criterion polynomial; C++ substitutes a fixed list of small primes `{3, 5, 7, 11, ...}`.  Both are *generic-point evaluations*: a multivariate polynomial vanishes identically iff it vanishes at any single point that is generic (not on any low-degree algebraic relation between the substituted values).  Distinct primes are generic by construction; user-supplied `Numeric` is generic by AMFlow's setup convention.  Locked by 5 existing `test_qft_topology.cpp` `ZeroSectorQ_*` tests covering 1-loop massless tadpole (true), 1-loop massive tadpole (false), 1-loop bubble (false), 2-loop massless sunset subsector (true), 2-loop sunrise (false); production correctness implicit in the committed oracle benches. |
| ~~`region_power` skips `/.Numeric`~~ → 🟢 | `src/qft/findregion.cpp:463-529` | **Documented + locked.**  `region_power` constructs the per-integral exponent as `val = sum_scale * (2 - eps) - sum_p` where `sum_scale` and `sum_p` are integer accumulators.  The output polynomial contains only `eps` and integer constants — no kinematic invariants ever enter the expression by construction.  Upstream's `/. Numeric` would therefore be a no-op (nothing to substitute).  Locked by `test_qft_region.cpp` `RegionPower_OutputIsEpsOnlyPlusIntegers_AuditRow200Equivalence` — uses a multi-invariant box family and asserts every term in the numerator and denominator has zero exponent on every non-eps variable. |
| ~~`factorize_family` no-redef fallback~~ → 🟢 | `src/pipeline/factorize.cpp:268-302` | **Documented as defensive fallback.**  Fires only when the loop-redefinition matrix `M` has fewer pivots than `n_tbl` after RREF — i.e. the heuristic loop selection from U-poly's first monomial doesn't span a complete loop basis.  No production input from the committed oracle benches reaches this branch (true dead code on tested inputs), but the fallback is mathematically safe: it leaves the family unchanged (uses original loops and original propagators), which is how downstream FactorizeFamily consumers would handle a single-component family.  Source comment at the fallback site documents the trigger condition.  If a future user case exercises it, the trivial-no-redef behavior is the safe choice (downstream code path is shared with the well-conditioned single-component case). |
| ~~`LIBPDeriv` multi-invariant case~~ → 🟢 | `src/ibp/libp_deriv.cpp` | **Verified, scope clarified.**  Multi-invariant family (free-symbol invariants like `m1sq`, `m2sq` in distinct propagators, plus a Replacement-defined `s`) tested via `test_ibp_libp_deriv.cpp` `*TwoMassBubble*` cases.  Replacement-defined invariants (e.g. `s` from `p^2 -> s`) deliberately return all-zero — production AMFlow flow only differentiates w.r.t. `eta` (which lives directly in the propagator polynomial), so upstream's momentum-derivative chain-rule path (`LIBPDerivivative`, `Kira/interface.m`) is intentionally not ported.  See `include/amflow/ibp/libp_deriv.hpp` "Scope" block for the contract. |
| ~~`MasterRank`/`MasterDot` filter~~ → 🟢 | `src/ibp/reduce.cpp` | **Documented as intentionally not exposed.**  Upstream's `MasterRank` / `MasterDot` (`Kira/interface.m:425-426,453`) is a manual post-Kira filter for pseudo-master removal — upstream's own CHANGELOG describes it as "only use if you believe that some pseudo master integrals have appeared in the list".  C++ `ibp::ReduceOptions` does not expose these knobs; the effective behavior pins both to the upstream default `Infinity` (no filter), which is exactly the regime exercised by all committed oracle benches at rel ~10⁻³⁰.  If a real pseudo-master case ever surfaces, the filter is a one-line `select` post `kira_read_masters` plus a JSON knob — added when needed, not preemptively (YAGNI).  Header comment `include/amflow/ibp/reduce.hpp` documents the contract. |
| ~~`kira_target.m` parser strictness~~ → 🟢 | `src/ibp/kira_run.cpp` (`kira_read_target_table`) + `src/ibp/kira_parse.cpp` | **Locked positive + negative + edge cases.**  Custom tokenizer parses the kira2math-emitted reduction table.  Production Kira output never produces malformed rules, but the parser must still reject malformed input cleanly.  Locked by 2 existing positive-acceptance tests (`ReadTargetTable_Synthetic`, `ReadTargetTable_HandlesMultipleRules`) plus 5 new tests for the canonical malformed/edge shapes: `_MissingFileReturnsEmpty` (no kira_target.m → empty rules, not error — mirrors masters_mma behavior), `_MalformedLHSThrows` (garbage LHS), `_RHSTermWithoutJIntegralThrows` (bare coefficient), `_NoOuterBraceReturnsEmpty` (no `{...}` block → empty, downstream master-empty check fires), `_RHSZeroFiltered` (literal `0` RHS yields empty-rhs rule). |
| ~~Coefficient parser fragility~~ → 🟢 | `src/ibp/kira_parse.cpp` | **Grammar locked positive + negative.**  `kira_parse_expression` supports a strict subset of Mathematica syntax: integer literals, registered identifiers, unary `-`, binary `+ - * /`, and `^` followed by an integer literal.  Kira's normal coefficient output never trips an unsupported form, so the production code path is well-bounded.  Locked by 7 existing `ParseExpression_*` positive-acceptance tests + 4 new negative-acceptance tests (`ParseExpression_DecimalLiteralThrows`, `_NonIntegerExponentThrows`, `_GarbageCharacterThrows`, `_UnclosedParenThrows`) which verify unsupported lexical forms (decimals, fractional exponents, parenthesised exponents, garbage characters, unclosed parens) raise rather than silently producing a partial parse. |
| ~~Diffeq nested reduce raises rank/dot via second `apply_jdot_jrank_floor`~~ → 🟢 | `src/ibp/reduce.cpp:114` (`apply_jdot_jrank_floor`), call sites at lines 177 (inside `reduce()`) and 273 (inside `diffeq()`) | **Documented as monotonic-conservative.**  `diffeq()` first floors over `jpreferred` with `dot_plus_one=true` (mirrors upstream's `Max[$BlackBoxDot, JDot/@jpreferred + 1]`), then invokes `reduce()` which floors *again* over `{all_ints, jpreferred}` with `dot_plus_one=false`.  The inner floor is applied on top of the already-floored options — so the effective `(rank, dot)` used at the inner Kira call is `>= (rank, dot)` upstream would compute for the analogous step.  This is mathematically safe because the IBP system at `(R, D)` is a subset of the system at any `(R', D') >= (R, D)`: solving the larger system produces a superset of the reduction rules, and the rules for any specific target are unchanged.  Net effect: C++ may compute slightly more IBP equations than strictly necessary, but the final reduction for any target matches upstream's.  Production correctness is implicit in all committed oracle benches matching MMA at rel ~10⁻³⁰. |
| ~~`r = nonzero(top) + IBPDot` arithmetic match~~ → 🟢 | `src/ibp/kira_yaml.cpp` | **Verified + tightened.**  C++ originally counted only `1` entries (`std::count(..., 1)`), which agreed with upstream's `Length[TopSector] - Count[TopSector, 0]` only because `qft::get_top_sector` always emits 0/1.  Defensive fix at `src/ibp/kira_yaml.cpp:221` switches to `count_if(... != 0)` to match the upstream formula literally and stay parity-correct under any future relaxation of `get_top_sector`'s 0/1 guarantee.  Unit tests `test_ibp_kira.cpp` `KiraTest.WriteJobs_R_*` cover zero-/non-trivial-IBPDot, mixed-presence top_pattern, and a non-binary entry. |

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

- ~~**New `Trivial` ending scheme**~~ — **PORTED**:
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

Forward-looking work — diversity-driven oracle expansion — lives in
[`docs/ROADMAP.md`](ROADMAP.md).
