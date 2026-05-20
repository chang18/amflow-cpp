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

**Status**: confirmed MMA-faithfulness gap. Not a correctness bug
(numerical output matches oracle 30+ digits post η-placement fix), but
an algorithmic-implementation gap that inflates Kira intermediate work
by orders of magnitude on heavy 4L benchmarks.

**Origin**: While verifying the η-placement fix on
`photon_4L_SE_TwoBubbles_mass1` (2026-05-20), I diffed C++ vs MMA
`jobs.yaml` and Kira logs at the **initial** `reduce/masters_preheat`
stage (i.e. the first Kira call inside `BlackBoxReduce`, which MMA calls
`IBPSystem`). MMA's `jobs.yaml` at the same stage uses
`select_mandatory_list: [pseTB1, target]`, while C++ writes
`select_mandatory_recursively: {sectors:[2047], r:16, s:0, d:5}`.

The actual MMA reference behaviour (from `Kira/interface.m:177-189`) is:

```mathematica
target = Which[
  $ReductionMode === "Masters",
    "select_mandatory_recursively: {topologies:[`fam`], sectors:[`top`], r:`r`, s:`s`, d:`d`}",
  True,
    "select_mandatory_list: [`fam`, target]"
];
```

So MMA **also** uses `select_mandatory_recursively` for the `Masters`
(= IBPSystem) call. The `select_mandatory_list` in MMA's cached
`jobs.yaml` is from the **subsequent** `AnalyticReduction` call that
overwrites the same directory's `jobs.yaml`. **Not a C++ vs MMA
divergence at this stage.**

**But** mass1 Kira logs reveal a real, separate gap:

| Stage | C++ post-fix | MMA (HANDOFF) |
|---|---|---|
| `reduce/masters_preheat` mandatory list | — | — |
| `amf/system_0_diffeq/masters_preheat` mandatory list | **336,760** | **18** |
| `amf/system_0_diffeq/masters_preheat` masters count | **18** | **16** |
| `amf/system_0_diffeq/masters_preheat` wallclock | 1902 s | 1857 s |

The masters count and wallclock are now MMA-equivalent (η-placement
fix). But the mandatory-list size at `system_0_diffeq/masters_preheat`
is still ~19,000× MMA's. Wallclock parity is achieved because Kira's
solver is fast even with a bloated mandatory list, but the algorithm
asks Kira for vastly more equations than necessary.

**Hypothesis (root cause to investigate)**: at the `system_X_diffeq`
stage, the **rank/dot** that C++ passes to Kira via `(r, s, d)` may be
larger than what MMA passes. MMA's `BlackBoxDiffeq` (`AMFlow.m:1261-1273`):

```mathematica
rank = Max[$BlackBoxRank, JRank/@jpreferred];
dot  = Max[$BlackBoxDot, JDot/@jpreferred + 1];
IBPSystem[top, rank, dot, jpreferred, ...];
```

— note the **`+1` on dot** is intentional, and `jpreferred` for the
diffeq step is the η-injected master set, not the user-input preferred.
C++ may be computing these parameters differently, expanding the
`sectors/r/s/d` envelope Kira enumerates.

**Where the divergence likely lives**:
- `src/pipeline/amfsystem.cpp::AMFSystem::build_diffeq` — the per-system
  `ibp::ReduceOptions` construction where `ibp_rank`/`ibp_dot` are set
  before handing off to `ibp::reduce`.
- `src/ibp/kira_yaml.cpp:285-340` (the `reduce_sectors` / `r/s/d`
  envelope writer) — but this is just emitting the values it's handed.

**Open questions for the user**:

1. Diff C++ and MMA `system_X_diffeq/masters_preheat/jobs.yaml`'s
   `(r, s, d)` values for the same benchmark. If different, locate
   the C++ source line that computes them and align to MMA's
   `Max[..., JDot+1]` formula.
2. If C++ already mirrors MMA's formula, check whether `jpreferred`
   (η-injected master list) is enumerated differently — e.g. C++ may
   include masters that MMA's `IBPSystem` filters out via
   `$MasterRank` / `$MasterDot` (interface.m:453).
3. Confirm impact: even though mass1 wallclock now matches MMA
   (2h27 ≈ MMA 2h28), the 19,000× mandatory-list bloat suggests
   Kira's internal selection is doing redundant work that doesn't
   show up in the dominant Fermat-bound numeric reconstruction
   phase. For heavier benchmarks (mass2, future 5L+), this bloat
   may dominate.

**Estimated impact if rooted out**:
- mass1: minor (wallclock already at MMA parity, ~1:1).
- mass2 and similar heavy multi-mass: potential 10-20% wallclock
  reduction if the Kira selection phase becomes a meaningful fraction
  of total time (currently dominated by Fermat reconstruction).
- 5L+ projection: likely much larger (mandatory list grows faster
  than masters count with loop order, so this gap widens with scale).

**Not actionable until**:
- (a) user evaluates and confirms scope; AND
- (b) post mass2/3L/4L full regression sweep, when we know which
  benchmarks have headroom for this optimization; AND
- (c) dedicated benchmark host is available (per
  `docs/PERFORMANCE.md`).

**Diagnostic data point (mass1 2026-05-20, post η-fix)**:
```
/tmp/amflow_photon_4L_SE_TwoBubbles_mass1_p0_eps001_cpp/
  part_0/amf/system_0_diffeq/masters_preheat/kira.log:
    length of mandatory list: 336760
    Number of master integrals: 18
    Total time: 1902 s
```
vs MMA cached HANDOFF-era figures:
```
mandatory list ≈ 18
masters ≈ 16
Total time ≈ 1857 s
```

