# AUDIT — current parity status against Mathematica AMFlow

This file is the canonical home for live parity status against the
upstream Mathematica AMFlow files
(<https://gitlab.com/multiloop-pku/amflow>):
`AMFlow.m`, `diffeq_solver/DESolver.m`, and
`ibp_interface/Kira/interface.m`.

Rules: MMA AMFlow is the spec; a C++ mismatch on a ported branch is a
parity bug; current focus is sampled single-point parity at `eps = 1/1000`.

Last updated: v1.0 release (single-commit history; live status is the green CI badge on the README).

---

## 1. Snapshot

- `ctest --test-dir build -N` discovers the test set under target
  `amflow_tests` (508 tests across 77 suites).
- `ctest --test-dir build --output-on-failure -j 4` is green from a
  local build.
- `tools/bench/` holds 13 oracle triplets; all match Mathematica
  reference at `rel ~1e-30 ~ 1e-33` when re-run via `amflow_cli`.
- No known red harness in the committed test suite.

## 2. Domain scoreboard

| Domain | MMA scope | C++ status |
|---|---|---|
| `numeric` | precision / options / log / acb wrappers / rational helpers | Complete |
| `algebra` | mpoly / mpoly_matrix / context migration | Complete |
| `ode` | ODE solver layers (port of `DESolver.m`) | Complete; 12/12 oracle cases match |
| `qft` | family_config / family_uf / topology / amfmode / region / boundary / vacuum / jintegral / dense_q | Complete |
| `ibp` | libp_deriv / Kira yaml/run/parse / reduce | Complete (Kira-only by project constraint) |
| `pipeline` | `AMFSystem`, ending schemes, factorize, `BlackBoxAMFlow{,Single}`, `GenerateNumericalConfig`, `FitEps`, `SolveIntegrals` | Complete; bubble + box1 fitted green |
| `api` | JSON dispatcher (`run_json`) | Complete |
| `cli` | `amflow_cli` JSON entry point | Complete with three modes (`amflow`, `black_box_amflow`, `solve_integrals`) |

`SolveIntegralsGaugeLink`, HQET / SCET / Wilson lines are out of scope
for this port.

## 3. Validated surface

Workflows backed by green local tests or committed sampled benchmarks:

- 1×1 power-law core, bubble (direct + AMF + sampled + fitted),
  ending tadpole including auto-reduction of `J(2)`, sunrise, Cutkosky
  routing + `prefactor * Im[...]` mapping.
- `box1` sampled Layer-17 + CLI sampled + one live `solve_integrals()` target.
- `vtx2` sampled top + preferred masters at `eps = 1/1000` (all
  three masters_only targets match MMA to ~1e-30 after the dual
  fix in `boundary_laporta`; previously the
  `j[vtx2,0,0,1,1,0,1,0]` numerator master was off by a constant
  factor because the PFD numerator polynomial used the
  `mom^2 + shift` convention while Kira reduces against
  `mom^2 - shift`).
- `tt` double-box four-target sampled at `eps = 1/1000`.
- Smirnov-Veretin planar massless double-box sampled scalar at
  `eps = 1/1000`, `s = -3`, `t = -1`.
- `tt` higher-rank ISP sampled at `eps = 1/1000`: rank-5 `[…,-5,0]`
  / `[…,-4,-1]` and rank-6 `[…,-3,-2]` / `[…,-3,-3]`. Validates the
  Layer-14 strict-QFT contraction on negative ISP indices up to rank 6.
- `xbox` 2L non-planar crossed-box sampled scalar top at `eps = 1/1000`,
  `s = -3`, `t = -1`.
- `tt` `EndingScheme = {Tradition, Cutkosky, SingleMass}` sampled
  scalar top at `eps = 1/1000`. The Cutkosky entry is exercised at the
  protocol level; for this kinematics the Tradition scheme handles all
  boundary integrals, and the result agrees with the
  `EndingScheme = {Tradition}` `tt` baseline.
- 3L equal-mass banana sampled top at `eps = 1/1000`, `psq = -3`,
  `msq = 1`. Validates the Layer-7 boundary-matching path on a
  rank-deficient (over-complete master basis) system.
- L=1/2/3 Cutkosky path-firing series (`cutbubble_1L`,
  `cutsunrise_2L`, `cutbanana_3L`) at `eps = 1/100`. Each family has
  `cut = {1, ..., 1, 0, ..., 0}` shaped to satisfy `phase_volume_q`
  on the top sector, so `AMFSystemSetupMaster` actually picks
  `EndingScheme::Cutkosky` (verified via `AMFLOW_DEBUG_SCHEME=1`).
  Together they exercise the `(-1)^(L+1)` sign at L=1/2/3 (+/-/+)
  and confirm `phase_loop_num` propagates through subsystem recursion
  up to `system 1 -> 2 -> 3 -> trivial`.
- Layer 17 `solve_integrals` honours `GlobalOptions::d0` (mirrors
  AMFlow.m's `"D0"` option, default 4): samples the family at the
  internal eps grid `cfg.eps_samples + (4-D0)/2` (AMFlow.m:1342 / 1351)
  and fits the Laurent expansion against the original user-facing grid
  (AMFlow.m:1356). `BlackBoxAMFlow` itself does not shift, matching MMA.
  Unit-tested at `NumericOptionsResetFixture.D0Shift_*` and the
  options round-trip is covered by `ApiRunJsonTest.OptionsAreEchoedBack`.
  End-to-end parity
  bench at `D0 != 4` committed: `box1_d0_7_3_solve_integrals` (1-loop
  box at `s=100, t=-1, D0=7/3`) matches MMA on `j[box1, 1, 0, 1, 0]`
  to ~30 sig digits across orders 0..2.
- Layer 7: exact `fmpq` Jordan in `NormalizeMat`; FLINT `qqbar` for
  algebraic eigenvalues and `Floor[Re[lambda]]`; `fmpq` projector path
  in diagonal ToFuchsian for higher-rank reducible poles. `CalcZero`
  has an internal block-local algebraic Jordan fallback while public
  `normalize_mat()` keeps a rational-output API.

## 4. Benchmark inventory

| Benchmark | State | Files |
|---|---|---|
| `vtx2` sampled top + masters | Matches MMA at `eps = 1/1000` | [`tools/bench/vtx2_2loop_vertex_eps001_mma_reference.json`](tools/bench/vtx2_2loop_vertex_eps001_mma_reference.json), [`tools/bench/vtx2_2loop_vertex_masters_only_eps001_mma_reference.json`](tools/bench/vtx2_2loop_vertex_masters_only_eps001_mma_reference.json) |
| `tt` sampled double-box | Matches MMA on four targets at `eps = 1/1000`; max relative delta ~`4.4e-31` (real) / `4.9e-26` (imag); rerun `152.4s` | [`tools/bench/tt_2loop_box_eps001_mma_reference.json`](tools/bench/tt_2loop_box_eps001_mma_reference.json) |
| Smirnov-Veretin sampled double-box | Matches MMA on `J[1,1,1,1,1,1,1,0,0]` at `s=-3,t=-1`; rerun `113.75s` | [`tools/bench/doublebox_sv_eps001_mma_reference.json`](tools/bench/doublebox_sv_eps001_mma_reference.json), [`tools/bench/doublebox_sv_eps001_black_box_amflow_cpp.json`](tools/bench/doublebox_sv_eps001_black_box_amflow_cpp.json) |
| `tt` higher-rank ISP (rank 5/6) | Matches MMA on four ISP targets `[…,-5,0]`, `[…,-4,-1]`, `[…,-3,-2]`, `[…,-3,-3]` at `eps = 1/1000`; max relative delta `~6.2e-31` (real) / `~1.5e-30` (imag); MMA wallclock `876.6s` | [`tools/bench/tt_higher_rank_eps001_mma_reference.json`](tools/bench/tt_higher_rank_eps001_mma_reference.json), [`tools/bench/tt_higher_rank_eps001_black_box_amflow_cpp.json`](tools/bench/tt_higher_rank_eps001_black_box_amflow_cpp.json) |
| `xbox` 2L non-planar crossed-box | Matches MMA on `J[1,1,1,1,1,1,1,0,0]` at `s=-3,t=-1`; relative delta `~1e-30`; MMA wallclock `310.85s` | [`tools/bench/xbox_2loop_eps001_mma_reference.json`](tools/bench/xbox_2loop_eps001_mma_reference.json), [`tools/bench/xbox_2loop_eps001_black_box_amflow_cpp.json`](tools/bench/xbox_2loop_eps001_black_box_amflow_cpp.json) |
| `tt` Cutkosky-probe (`EndingScheme = {Tradition, Cutkosky, SingleMass}`) | Matches MMA on top scalar at `eps = 1/1000`; relative delta `~4.7e-31` (real) / `~4.0e-32` (imag); MMA wallclock `413.5s` | [`tools/bench/tt_cutkosky_probe_eps001_mma_reference.json`](tools/bench/tt_cutkosky_probe_eps001_mma_reference.json), [`tools/bench/tt_cutkosky_probe_eps001_black_box_amflow_cpp.json`](tools/bench/tt_cutkosky_probe_eps001_black_box_amflow_cpp.json) |
| 3L equal-mass banana | Matches MMA on `j[banana3,1,1,1,1,...]` at `psq=-3, msq=1, eps=1/1000`; max relative delta `~2.4e-30` (real) / `~2.2e-31` (imag); MMA wallclock `117.0s`, C++ wallclock `71.8s`. Requires `BlackBoxDot=5` for both Masters and Reduce calls and Kira Pak symmetry detection enabled. | [`tools/bench/banana_3loop_eps001_mma_reference.json`](tools/bench/banana_3loop_eps001_mma_reference.json), [`tools/bench/banana_3loop_eps001_black_box_amflow_cpp.json`](tools/bench/banana_3loop_eps001_black_box_amflow_cpp.json), [`tools/bench/banana_3loop_eps001_black_box_amflow_mma.wl`](tools/bench/banana_3loop_eps001_black_box_amflow_mma.wl) |
| `cutbubble_1L` Cutkosky probe (L=1) | Matches MMA on `j[cutbubble,1,1]` at `s=4, eps=1/100`; relative delta `<1e-29`; MMA wallclock `22s`, C++ wallclock `<1s`. `[scheme] Cutkosky fired: family=cutbubble phase_loop_num=1 ...` (verified via `AMFLOW_DEBUG_SCHEME=1`); leading eps^0 reduces to `1/(8*Pi)`. | [`tools/bench/cutbubble_1L_eps001_mma_reference.json`](tools/bench/cutbubble_1L_eps001_mma_reference.json), [`tools/bench/cutbubble_1L_eps001_black_box_amflow_cpp.json`](tools/bench/cutbubble_1L_eps001_black_box_amflow_cpp.json), [`tools/bench/cutbubble_1L_eps001_black_box_amflow_mma.wl`](tools/bench/cutbubble_1L_eps001_black_box_amflow_mma.wl) |
| `cutsunrise_2L` Cutkosky probe (L=2) | Matches MMA on `j[cutsunrise,1,1,1,0,0]` at `s=4, eps=1/100`; relative delta `~1.85e-31`; MMA wallclock `55s`, C++ wallclock `<1s`. `[scheme] Cutkosky fired: family=cutsunrise phase_loop_num=2 top_position={0,1,2}`; verifies the `(-1)^(L+1) = -1` sign flip at L=2. | [`tools/bench/cutsunrise_2L_eps001_mma_reference.json`](tools/bench/cutsunrise_2L_eps001_mma_reference.json), [`tools/bench/cutsunrise_2L_eps001_black_box_amflow_cpp.json`](tools/bench/cutsunrise_2L_eps001_black_box_amflow_cpp.json), [`tools/bench/cutsunrise_2L_eps001_black_box_amflow_mma.wl`](tools/bench/cutsunrise_2L_eps001_black_box_amflow_mma.wl) |
| `cutbanana_3L` Cutkosky probe (L=3) | Matches MMA on `j[cutbanana,1,1,1,1,0,0,0,0,0]` at `s=4, eps=1/100`; relative delta `<1e-29`; MMA wallclock `93s`, C++ wallclock `~90s`. `[scheme] Cutkosky fired: family=cutbanana phase_loop_num=3 top_position={0,1,2,3}`; sign returns to `+1` at L=3 and exercises subsystem recursion `1 -> 2 -> 3 -> trivial`. **Requires `BlackBoxDot=5`** (matching `banana_3loop`) — set explicitly in both the MMA `*_mma.wl` and C++ `*_cpp.json`; AMFlow.m's `Max[$BlackBoxDot, JDot/@...]` floor cannot raise it on its own when every positive index is 1 (`JDot = 0`). | [`tools/bench/cutbanana_3L_eps001_mma_reference.json`](tools/bench/cutbanana_3L_eps001_mma_reference.json), [`tools/bench/cutbanana_3L_eps001_black_box_amflow_cpp.json`](tools/bench/cutbanana_3L_eps001_black_box_amflow_cpp.json), [`tools/bench/cutbanana_3L_eps001_black_box_amflow_mma.wl`](tools/bench/cutbanana_3L_eps001_black_box_amflow_mma.wl) |
| `box1_d0_7_3` Laurent solve at `D0 = 7/3` | Matches MMA on `j[box1, 1, 0, 1, 0]` at `s=100, t=-1` to ~30 sig digits across orders 0..2 (`goal_digits=30, eps_order=4` → samples=12, internal `working_pre=182`, `x_order=364`). MMA wallclock `22.3s`, C++ wallclock `10.7s` (~2× speed ratio). Exercises the `(4-D0)/2` internal eps shift (AMFlow.m:1342, 1351) and `FitEps` Laurent fit in the user-facing eps (AMFlow.m:1356), the first end-to-end check that `GlobalOptions::d0` cancels cleanly in the user output. | [`tools/bench/box1_d0_7_3_solve_integrals_mma_reference.json`](tools/bench/box1_d0_7_3_solve_integrals_mma_reference.json), [`tools/bench/box1_d0_7_3_solve_integrals_cpp.json`](tools/bench/box1_d0_7_3_solve_integrals_cpp.json), [`tools/bench/box1_d0_7_3_solve_integrals_mma.wl`](tools/bench/box1_d0_7_3_solve_integrals_mma.wl) |
| `tradcut_phase_2L` Tradition-with-cut probe | Matches MMA on `j[phase, 1, 0, 1, 0, 1, 0, 0]` at `s=100, msq=1, eps=1/100` to relative delta `2.68e-30` (real) / `0` (imag); MMA wallclock `39.3s`. The first oracle to actually exercise the parent-cut → boundary-sub-family projection in `AMFSystem::build_boundary` (Phase 1B / D3 fix; mirror of AMFlow.m:790-803).  Family taken from upstream `examples/automatic_phasespace/run.wl`: 7 propagators, cut={1,0,1,0,1,0,0}, prescription={0,0}; cut topology is not phase_volume so Cutkosky does not apply and Tradition fires at an uncut prop. | [`tools/bench/tradcut_phase_2L_eps001_mma_reference.json`](tools/bench/tradcut_phase_2L_eps001_mma_reference.json), [`tools/bench/tradcut_phase_2L_eps001_black_box_amflow_cpp.json`](tools/bench/tradcut_phase_2L_eps001_black_box_amflow_cpp.json), [`tools/bench/tradcut_phase_2L_eps001_black_box_amflow_mma.wl`](tools/bench/tradcut_phase_2L_eps001_black_box_amflow_mma.wl) |
| Raw MMA caches | Regenerable; not tracked. See [`tools/bench/README.md`](tools/bench/README.md) §"Regenerating Mathematica caches" |

## 5. Workflow for sampled parity expansion

1. Generate or preserve MMA sampled values.
2. Run C++ on the same inputs.
3. Compare first direct subintegrals and child integrals when needed —
   not only final fitted outputs.
4. Lift matched results into committed tests or benchmark artifacts
   when runtime permits.

## 6. Practical definition of progress

A parity claim is credible only when:

- MMA and C++ agree on the target;
- no tolerance was weakened to manufacture agreement;
- compared inputs are the same family/replacements/numerics/`eps`/targets;
- MMA outputs are preserved in the repo for later reuse.
