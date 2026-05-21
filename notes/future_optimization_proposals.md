# Future Optimization Proposals

This file collects performance / engineering improvement ideas that are
**not yet validated and not yet implemented**. Each entry is a candidate
for the user to evaluate (correctness, priority, scope) before any code
change is undertaken.

Entries are dated and include the discovery context. Once an entry is
either accepted-and-implemented or rejected, move it out of this file
(into `CHANGELOG.md`, a code comment, or just delete it).

---

## 2026-05-19 — Cache `sectormappings/` across `system_X_diffeq` runs of the same family

**Status**: unverified proposal, awaiting evaluation.

**Origin**: During a 4-loop benchmark (`photon_4L_SE_TwoBubbles` mass1/mass2
staged-mass batch), I initially misread Kira logs and claimed C++
amflow_cli was using a ~26,000× larger IBP than MMA on the
`amf/system_X_diffeq/masters_preheat/` step. An independent audit refuted
that — C++ and MMA are doing exactly the same Masters-preheat +
Reduce-step pair on each system_X_diffeq run, and the apparent
discrepancy was just MMA overwriting its own Step-1 log when Step-2
runs in the same directory. **No code divergence exists, no correctness
fix is needed.** See agent audit summary in chat history (2026-05-19).

But the same audit incidentally noted a **separate, real** optimization
opportunity unrelated to that false alarm:

> Kira's `sectormappings/` (the topology / sector-symmetry tables) are
> **η-independent**: they depend only on the propagator family
> structure, not on which propagator carries the auxiliary mass η.
> Therefore the sectormappings produced during the first
> `amf/system_0_diffeq/masters_preheat/` could be reused as input to
> `amf/system_1_diffeq/masters_preheat/`, `system_2_diffeq/`, etc.,
> instead of having each Kira invocation rebuild them from scratch.
>
> Kira spends a meaningful fraction of preheat wall-time on building
> these tables (symmetry/sector relations enumeration). For a 4-loop
> 14-slot family that's roughly ~30 min × N_systems of repeated work
> that, in principle, only needs to be done once per family.

**Where the change would land** (per the audit):
- C++ jobs.yaml writer: `src/ibp/kira_yaml.cpp:170-299` —
  `kira_write_config()` is where the sectormappings dir gets
  initialized per Kira call.
- C++ diffeq driver: `src/pipeline/amfsystem.cpp:1035-1046`
  (`AMFSystem::build_diffeq`) — this is the per-system loop where
  successive system_X_diffeq runs are launched.
- The fix would be: for system_0, materialize sectormappings as today;
  for system_1+, point Kira at the existing sectormappings directory
  (symlink, copy, or pass via Kira config), bypassing the rebuild.

**Open questions for the user to evaluate**:

1. **Correctness**: is the η-independence claim watertight? The
   propagator list passed to Kira at each `system_X_diffeq` differs
   only by which propagator has `-eta` added. Symmetries and sector
   relations are derived from the structural form of propagators (loop
   momenta + external momenta + masses), not the specific η placement.
   Plausible but worth double-checking with a Kira-internals reader
   before any code change.
2. **Verify MMA does NOT do this**: MMA might already cache
   sectormappings between calls (it reuses the same Kira working dir).
   If MMA already gets the speedup implicitly via dir-reuse, the C++
   sibling-dir pattern (introduced as a Kira-2.x bug workaround at
   `reduce.cpp:174-182`) might be the only thing preventing it on C++
   side — in which case the fix is much smaller (revisit the Kira-2.x
   workaround for whether it still applies on current Kira versions).
3. **Priority**: per the project's no-perf-claims-until-dedicated-host
   policy (see `docs/PERFORMANCE.md`), perf-only optimizations should
   wait until after correctness validation completes. This entry is
   recorded so it isn't forgotten, not as a request to act on now.

**Estimated savings if real**:
- 4L 14-slot family: ~30 min/system of preheat overlap × (N_systems − 1)
  systems saved. For mass1 with ~ N_systems = 5–10 typical: 2–5 hours.
- 3L 9-slot family: ~30 sec/system, marginal.

**Not actionable until**:
- (a) user evaluates and confirms; AND
- (b) correctness validation phase of the project completes (per
  `docs/PERFORMANCE.md` policy); AND
- (c) dedicated benchmark host is available to measure actual win.

---

## 2026-05-20 — C++ uses `select_mandatory_recursively` in `IBPSystem` stage where MMA uses `select_mandatory_list`

**Status (2026-05-21)**: **CLOSED — false alarm**.  The headline "19,000×
mandatory-list gap" was a misread of MMA's cached `kira.log`.  MMA cache
appends two Kira invocations (`IBPSystem` Masters-mode followed by
`AnalyticReduction` Reduce-mode in the same directory), and the
`length of mandatory list` numbers we tabulated (e.g. 18 for mass1,
136 for sunset_bubble) are the *Reduce-mode* output (literally equal
to the `target` file's integral count, which is the
`select_mandatory_list: [fam, target]` selector's mandatory list).
C++ instead reports the *Masters-mode* `select_mandatory_recursively`
mandatory list, which is the sector-wide enumeration size — a
fundamentally different quantity from the Reduce-mode target count.

The minimum-reproducer experiment that closed this:
take MMA's `config/` + `preferred` from any cached `<sysid>/diffeqsetup`,
write a fresh Masters-mode `jobs.yaml` with the same `(r, s, d)`
parameters MMA used, and run Kira from a clean dir → produces the
**same** Masters-mode mandatory list as C++.  Bench-wide validation
on sunset_bubble_4L_2leg_eqmass: 20/20 systems match C++ vs
fresh-MMA-reproduce numbers exactly (e.g. system_0 = 3360 vs
MMA-reproduce 3360; system_44 = 91 vs 91; all 20 line-by-line ✓).

What was actually fixed during the investigation (separately
warranted, see `58cb588`): three mirror-MMA-behavior alignment
changes in `src/ibp/kira_yaml.cpp` + `src/ibp/reduce.cpp` —
`[Propagator, 0]` yaml format, closed-form rank-1 propagator
reconstruction, and disabling the over-correcting `sort_integrals_like_amflow`.
These are clean alignment improvements but did not — and could not —
"fix" the 19,000× gap because the gap never existed.

Lesson recorded at user-level `feedback_no_wild_guess_use_minimum_reproducer.md`:
when third-party tool outputs disagree between C++ and MMA on the
same nominal inputs, isolate via minimum reproducer experiment
before patching; cosmetic input differences rarely matter to mature
parsers, and "same log keyword" can mean different stage outputs
across two parallel runs of the same binary.

