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

## 2026-05-19 — C++ `system_0_diffeq/masters_preheat` is ~2.66× slower than MMA's analogous step

**Status**: root cause identified (η-injection strategy divergence
from MMA reference). NOT a numerical correctness bug. IS an algorithm
implementation faithfulness gap that should be corrected to align with
MMA, per project principle "MMA = correctness oracle, legacy = design
reference". Fix not yet applied, awaiting user decision.

**Origin**: During the same `photon_4L_SE_TwoBubbles mass1` batch that
seeded the previous entry, after MMA mass1 completed in `==TIME==
8872.571245s` (≈2h28m), C++ mass1 was still running. Drilling into
sub-stage times in C++'s `kira.log` files vs MMA's reported step times
revealed a real wallclock gap on the diffeq-preheat step (NOT on the
earlier full-family reduce steps):

| Phase | C++ measured (`kira.log` Total time) | MMA reported (`IBPSystem` line) | Ratio |
|---|---|---|---|
| 1. Initial `reduce/masters_preheat` | 1856.6 s | 1857 s | 1.00× |
| 2. `reduce/target_reduce` | 1864.8 s | 1853 s | 1.01× |
| 3. **`amf/system_0_diffeq/masters_preheat`** | **4932.8 s** | **1857 s** | **2.66×** |

Phases 1–2 are essentially identical (within seconds), so it is NOT a
general "C++ overhead per Kira call" issue. The slowdown is specific to
the system_X_diffeq stage.

**This contradicts the prior audit conclusion** that "C++ and MMA hand
Kira identical yaml at the diffeq stage". The audit (run earlier the
same day) did not diff the system_0_diffeq yaml specifically — it diffed
the initial-stage yaml. Either:

(a) the audit's identical-yaml claim was correct, and the slowdown is
    purely runtime/host-load (but the timing pattern doesn't fit:
    parallel competition was similar in both stages); OR
(b) the C++ amflow_cli writes a meaningfully different yaml/jobs.yaml
    at the system_X_diffeq stage than at the initial reduce stage —
    and that difference is what causes 2.66× slowdown.

**What's needed**: a subagent dive specifically into the system_0_diffeq
phase yaml on both sides, comparing:
- `tmp/amflow_photon_4L_SE_TwoBubbles_mass1_p0_eps001_cpp/part_0/amf/system_0_diffeq/masters_preheat/{integralfamilies.yaml, kinematics.yaml, jobs.yaml, preferred, config/*}`
- `tools/bench/photon_self_energy/mma_refs/photon_4L_SE_TwoBubbles_mass1_p0_eps001_mma_cache/1/diffeqsetup/{...same files...}`
- The C++ source code path that constructs these (`src/ibp/kira_yaml.cpp`
  with the call sites from `src/pipeline/amfsystem.cpp`).
- The MMA AMFlow code path that constructs these (`reference/amflow-master/AMFlow.m::BlackBoxDiffeq` →
  `Kira.m::IBPSystem`).

**Specific hypotheses to test**:
1. C++ asks Kira to compute over a larger sector set / larger rmax/dmax
   than MMA at this stage.
2. C++ doesn't pass `preferred_masters` (or passes a different list),
   causing Kira to re-discover masters from scratch.
3. C++ writes a `select_mandatory_recursively` job where MMA writes
   `select_mandatory_list` (or vice versa), changing Kira's algorithmic
   path.
4. C++ omits some option (e.g. `integral_ordering`, symmetry reuse from
   the previous stage's `sectormappings/`) that MMA uses to cut work.

**Don't conflate with previous entry**: the sectormappings-caching
opportunity (preceding entry) is a separate optimization across multiple
system_X runs. This entry is about a single system_0 run being 2.66×
slower in C++ than MMA on what is supposedly the same Kira job.

**Estimated impact if rooted out**: each `system_X_diffeq` stage in a 4L
14-slot family costs ~30 min in MMA vs ~80 min in C++. Across 5–10
systems per family, eliminating the 2.66× gap could save 4–8 hours per
4L mass-staged test. Most of the C++ vs MMA wallclock gap currently
observed in 4L benchmarks is concentrated in this stage.

**Root cause found (subagent audit, 2026-05-19, high confidence on yaml
diff, medium on exact source line)**:

**Important clarification — how to characterize this**:

This is **a strategy divergence between C++ and the MMA reference
implementation**, NOT a numerical correctness bug. Two distinct levels:

- **Numerical correctness** (the eps-expanded integral value at η=0):
  η-placement on different propagators all yield the same value via
  η-flow. C++ will still produce **numerically correct** output when it
  finishes — just takes 2.66× longer at this stage.
- **Algorithmic faithfulness to MMA**: MMA's `AMFPosition` deliberately
  places η on the propagator that minimizes Kira's intermediate work
  (the Branch-mode tie-breaker is part of MMA's algorithm design, not a
  performance afterthought). C++ deviates from this strategy and picks
  a worse propagator. Per the project's principle "MMA = correctness
  oracle, legacy = design reference", **this is an implementation
  faithfulness gap and should be brought into alignment with MMA**.

So this is NOT a "perf nice-to-have to defer until dedicated-host
benchmarking". It is "C++ implements a different algorithm than the
MMA reference; align to MMA". The fact that the numerical answer is
the same is reassuring (it means the gap can be closed without
worrying about regression on existing oracle tests) — but it does not
demote the priority.

C++ and MMA's `amf_position` choose **different propagators** for the η
insertion in this single-mass-with-symbolic-msq case:

| Side | η placement | mandatory list | masters | wallclock |
|---|---|---|---|---|
| MMA | prop 0 (l₁², combines with the existing m²) | 18 | 16 | 1,857 s |
| C++ | prop 2 (l₂², a previously massless propagator) | 477,748 | 154 | 4,932.8 s |

Putting η on a previously-massless propagator breaks zero-sector
identification and inflates the IBP system 10× in masters, 26,000× in
mandatory-list, leading to the 2.66× wallclock penalty.

**Where both differ in source**:

| Step | C++ | MMA |
|---|---|---|
| Compute η position | `qft::amf_position(*fc_use, top, opts.amf_modes)` at `src/pipeline/amfsystem.cpp:3067` | `AMFEtaC@AMFPosition[GetTopPosition[preferred], $AMFMode]` at `reference/amflow-master/AMFlow.m:1041` |
| Branch-mode fallback (where both modes converge here, since `single_mass_q` returns false for symbolic msq) | `src/qft/amfmode.cpp:239-283` (Branch mode), with `same_branch_with` at `:178-199` | `reference/amflow-master/AMFlow.m:585-592` (`findbranch[v_]`) |
| Tie-break sort comparator | `src/qft/amfmode.cpp:262-268`: stable_sort with strict `<` on branch size; ties broken by insertion order in `info.var` | MMA's `Sort[branches, Length@#1 <= Length@#2 &]` with MMA's monomial-iteration order |
| Insertion order origin | `src/qft/topology.cpp:342-397` `make_component`, building `ci.var = feynman_vars_in(u0_factor, first_x, n_x)`; depends on fmpz_mpoly's grlex monomial iteration | MMA's default monomial ordering on the same polynomial |

The C++ `info.var` ordering (FLINT grlex) and MMA's ordering yield
different insertion orders into `branches[]`, so when several branches
tie on size, C++ picks a different one than MMA. The picked branch
determines which propagator gets η — and putting η on a massless prop
is much worse for Kira.

**There is ONE algorithm to align to (MMA's), not "two strategies"**:

MMA's `findbranch` (`reference/amflow-master/AMFlow.m:585-592`)
implements a single algorithm. The fix is to faithfully reproduce
that algorithm in C++, not to invent a heuristic that "happens to
match on this case". Two candidate fix routes:

**Route A — faithful reproduction (correct alignment, harder)**:
Read MMA `findbranch[v_]` line-by-line, identify exactly what it
does (sort key, tie-break, iteration order over branches), then
make C++'s `same_branch_with` and the Branch-mode sort at
`src/qft/amfmode.cpp:178-199, 239-283` produce the same `pos` for
every input. This may require matching MMA's monomial iteration
order on the U polynomial (`src/qft/topology.cpp:342-397`
`make_component` `ci.var = feynman_vars_in(u0_factor, first_x, n_x)`
vs MMA's equivalent), since that's what feeds branch construction.
This is "true alignment with reference".

**Route B — secondary heuristic (NOT alignment, only workaround)**:
Add a tie-breaker like "prefer branch containing a massive
propagator" in `src/qft/amfmode.cpp:262-268`. This would happen to
pick the same propagator as MMA on the observed case, but is **not
a faithful reproduction of MMA's algorithm** — there's no such
explicit "prefer-mass" rule in MMA's source. It would diverge from
MMA on cases where MMA's monomial-ordering tie-break happens to
NOT pick the mass-bearing branch. So Route B should be considered
only as a stopgap if Route A is too costly.

Route A is the correct work for "implementation faithfulness". Route
B is a band-aid that may cause its own future divergences.

**Required preparatory step**: read MMA's `findbranch[v_]` source
end-to-end. The audit so far has cited the file:line but has not
unpacked the exact algorithm logic. That's the precondition for any
real fix.

A `AMFLOW_DEBUG_SCHEME=1` runtime trace re-run (see
`src/pipeline/amfsystem.cpp:3074-3092` — the trace already emits the
`pos={...}` decision) can confirm `pos = {2}` vs `pos = {0}` before
any fix is attempted.

**Test that would regress on the fix**: `AMFCandidateComponent_*` family
in the existing test suite. Note: since η-placement choice is
position-agnostic on correctness (different η placements all converge
to the same integral value via η-flow), no existing correctness test
should regress — only speed/intermediate-size tests would change.

**Workaround for users until fix lands**: pre-substitute `msq` in the
input cpp.json so `single_mass_q` succeeds and the Single-mass mode
fires (which DOES match MMA correctly). The benchmarks currently set
`numeric_values: {msq: "1"}` only at the `blackbox` level, leaving the
propagator-list version of `msq` symbolic at the time `amf_position` is
called.
