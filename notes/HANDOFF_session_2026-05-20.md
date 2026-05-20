# Session Handoff — 2026-05-20

## Closed this session

**η-placement bug** — `src/pipeline/amfsystem.cpp::amf_system_setup_master`
previously built `fc_use = build_numeric_replacement_family(fc, numeric_q)`
for all amf-scheme decisions; the numeric substitution silently flipped
`info.vacQ` from false to true for benchmarks where a replacement-map
invariant resolved to zero, triggering Branch-mode fallback and
mis-placing η on a previously-massless propagator. MMA's `AMFPosition`
operates on symbolic `ReducedPropagator`. Removed `fc_use` indirection
completely; deleted orphan helpers.

Pushed to `origin/main`:
- `20883ed` — chore: untrack `mma_refs/` Kira cache + retire HANDOFF doc
- `bef11d9` — fix(amfsystem): align η-placement scheme decision with MMA's symbolic ReducedPropagator

## Validation
| Group | Result |
|---|---|
| 1L | 49/49 PASS |
| 2L | 45/46 PASS (1 pre-existing wzbox precision quirk, **unrelated to this fix**) |
| 3L | 61/61 PASS (incl. new `eta_c_3L_Mercedes_eps001` oracle) |
| 4L photon mass1 (pre-fix: hit rc=1) | ✓ 2h27m C++ vs 2h28m MMA (0.99×) |
| 4L photon mass2 (pre-fix: 13.5h hang) | ✓ 3h51m C++ vs 3h48m MMA (1.01×) |
| 4L batch | 59/59 PASS |
| ctest unit tests | 549/549 PASS |

## Open follow-ups

### 1. Mandatory-list 19,000× MMA-faithfulness gap (recorded, not actionable)
Even post-fix, C++ `amf/system_0_diffeq/masters_preheat` Kira call has
~336,760 mandatory list vs MMA ~18. Masters count and wallclock are now
MMA-equivalent (Fermat dominates), but Kira selection is doing 19,000×
more work than needed. Hypothesis: C++ passes a wider `(r, s, d)` to
Kira than MMA's `BlackBoxDiffeq` formula (AMFlow.m:1261-1273:
`dot = Max[$BlackBoxDot, JDot+1]`).
Investigation point: `src/pipeline/amfsystem.cpp::AMFSystem::build_diffeq`.
Recorded in `notes/future_optimization_proposals.md` (2026-05-20 entry).
**Not actionable until: (a) user confirms scope, (b) post-regression sweep,
(c) dedicated benchmark host.**

### 2. `cpp_verification_snapshot` injection — going forward only
- Helper at `/tmp/inject_cpp_snapshot.py` (one-shot, may need
  re-creation in new sessions).
- `bench_runner.sh` at `/tmp/bench_runner.sh` auto-calls it on PASS.
- Schema in `tools/bench/eta_c_form_factor/eta_c_3L_Mercedes_eps001_mma_reference.json`.
- **User's instruction**: save C++ snapshots only when a new run produces
  them; do NOT backfill historical reference.json that were committed
  pre-snapshot-convention.

### 3. Project artifact: `~/tmp/mass2_comparison_report.md`
Standalone report of mass2 C++ vs MMA numeric agreement + per-stage
timings + η-placement debug trace. User-requested 2026-05-20.

## Files of interest for future sessions

| Path | Note |
|---|---|
| `src/pipeline/amfsystem.cpp:2860-3120` | post-fix `amf_system_setup_master` (uses symbolic fc throughout) |
| `src/qft/amfmode.cpp:494-508` | `amf_candidate_component` cascade + Branch fallback |
| `reference/amflow-master/AMFlow.m:503-510` | MMA `AnalyzeTopology` — no /.Numeric (the canonical anchor for the fix) |
| `reference/amflow-master/AMFlow.m:1041-1090` | MMA `AMFSystemSetupMaster` / scheme dispatch |
| `notes/future_optimization_proposals.md` | open MMA-faithfulness gaps |
| `tools/bench/compare_sampled.py` | numeric comparator, uses `requested_targets` + `sampled_values_from_same_mma_run` |

## What was NOT changed (intentional)
- `tools/bench/*_mma_reference.json` for benches not re-run this session
  (per user: don't backfill).
- 4L benchmarks beyond the two photon SE staged-mass cases.
- 4L mass3+ (HANDOFF notes mass5 would need dot=9, infeasible).

## Memory updated
New entries added to user-level memory at
`/home/chang/.claude/projects/-home-chang-workspace-agent-trials-amflowRefactor/memory/`:
- `project_eta_placement_closed.md`
- `project_mandatory_list_19000x_pending.md`
- `feedback_direct_compare_no_rerun.md`
- `feedback_cpp_values_save_forward.md`

## Pickup recipe for new session

1. `git log --oneline -3` to confirm head is at `bef11d9` (or later if
   intervening commits).
2. `git status` should show working tree clean.
3. If asked about η-placement or 4L photon SE — read this handoff +
   the four memory entries above.
4. If launching new bench, use `/tmp/bench_runner.sh` (preserves
   cpp_verification_snapshot) and not the old in-bench runner.
5. If user wants to compare values, follow `feedback_direct_compare_no_rerun`:
   read saved data first, only launch a fresh amflow_cli if data
   genuinely missing.
