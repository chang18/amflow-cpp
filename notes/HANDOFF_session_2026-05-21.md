# Session Handoff — 2026-05-21

## TL;DR

Three mirror-MMA alignment fixes in the Kira yaml writer + a debunking
of the headline "19,000× mandatory-list gap" recorded in the
2026-05-20 entry of `notes/future_optimization_proposals.md`.  The
gap was a misread of MMA's cached `kira.log` (Reduce-mode target
count, not Masters-mode enumeration size).  Verified end-to-end on
two benches with fresh parallel C++/MMA runs and per-system
line-by-line intermediate-output diff.

## Commits

| Hash | Summary |
|---|---|
| `58cb588` | fix(ibp): align Kira yaml emission with MMA interface.m |
| `7c31bcc` | chore(cli,ibp): AMFLOW_STAGE_TIMING env-gated wallclock markers |
| `9407b34` | docs(notes): close 2026-05-20 19000x mandatory-list entry as false alarm |
| _(this commit)_ | docs+test+chore: cleanup + unit test + handoff |

## What changed in code

`src/ibp/kira_yaml.cpp`:
- Propagator yaml now emits `[Propagator, 0]` (full denominator + mass=0),
  matching MMA `Kira/interface.m:126` active code.  The old
  `[Momentum, -Mass]` form was the commented-out MMA `:125` legacy
  writer.
- New `try_render_closed_form` helper reconstructs `±(v·x)² + const`
  closed form from FLINT-expanded `fmpz_mpoly` via rank-1
  quadratic-form factoring (extracts coefficient matrix A from
  per-term coeffs, picks a non-zero diagonal pivot, derives v, then
  verifies A == sign·v·vᵀ).  Falls back to expanded form on
  non-rank-1.  Mathematica's symbolic engine preserves closed form
  natively; FLINT mpoly is expanded by construction.
- `kira_write_preferred` writes caller-supplied order verbatim
  (matches MMA `:164-170` `Preferred[preferred, dir]`, also no sort).

`src/ibp/reduce.cpp`:
- `sort_integrals_like_amflow` is now a no-op (body removed, signature
  kept with documenting comment).  Previous `stable_sort + group_reverse`
  was a net REVERSAL of Kira's natural ascending-by-sector output —
  fed the next Kira call a wrong-direction preferred file.  Kira's
  raw output is already MMA-compatible (verified line-by-line on
  sunset_bubble root reduce vs MMA `cache/0/results/masters_mma`).

`src/cli/main.cpp` + `src/ibp/kira_run.cpp`:
- Opt-in `AMFLOW_STAGE_TIMING` env var gates two stderr markers:
  `[kira_time]` per `kira_run`, `[total_time]` once at amflow_cli main
  exit.  Lets external tools split wallclock into Kira-subprocess
  time vs amflow-internal time without further code changes.
- Default-off → no behavior regression to silent_mode flow.

`tests/test_ibp_kira.cpp`:
- `KiraTest.WriteConfig_OneLoopBubble_HasExpectedYaml` extended to
  pin the new propagator yaml format (`[FullDenom, 0]` form + closed
  form `(l - p)^2 - 1`).  Prevents silent regression to the legacy
  `[Momentum, -Mass]` writer.

## Validation

| | Result |
|---|---|
| `ctest -j$(nproc)` | 549/549 PASS |
| `bubble_1L_diffmass_eps01` fresh parallel | C++ 30s vs MMA 49s; 2/2 systems all-quantity match; rel ~1.7e-30 |
| `sunset_bubble_4L_2leg_eqmass_eps001` fresh parallel | C++ 668s vs MMA 713s; **20/20 systems Masters-mode mandatory list match line-by-line** (C++ 3360 = MMA-fresh-reproduce 3360 for system_0, …); rel ~1.7e-30 |
| `AMFLOW_STAGE_TIMING` on bubble_1L | C++ amflow-only 0.08 s vs MMA amflow-only 13 s (≈153× faster non-Kira) — verified C++ has the expected algorithmic advantage outside Kira/Fermat |

## What was the "19,000× gap" really?

Misread of MMA's cached `kira.log`.  MMA writes `IBPSystem` (Masters
mode) + `AnalyticReduction` (Reduce mode) into the **same** directory;
the latter overwrites `jobs.yaml` and appends to `kira.log`.  Both
calls print a line of the form `length of mandatory list: N`, but the
N reported in the cached log is the *Reduce-mode* output (= target
file integral count, the `select_mandatory_list: [fam, target]`
selector's mandatory list), not the Masters-mode `select_mandatory_recursively`
sector-wide enumeration that C++ reports.

Minimum-reproducer that closed it:
1. Take MMA's `config/` + `preferred` from any
   `cache/<i>/diffeqsetup`.
2. Write a fresh Masters-mode `jobs.yaml` with the same `(r, s, d)`
   parameters MMA used (parse from MMA's `jobs.yaml`).
3. Run `kira` from a clean dir.
4. → Produces the **same** Masters-mode mandatory list as C++.

Bench-wide on sunset_bubble_4L_2leg_eqmass: 20/20 systems match.

## Files of interest for future sessions

| Path | Note |
|---|---|
| `src/ibp/kira_yaml.cpp:43-265` | `mfrac_to_kira` + `try_render_closed_form` rank-1 helper + propagator emit |
| `src/ibp/kira_yaml.cpp:605-625` | `kira_write_preferred` no-sort writer |
| `src/ibp/reduce.cpp:112-130` | `sort_integrals_like_amflow` no-op (with documenting comment) |
| `src/ibp/kira_run.cpp:142-154` | `[kira_time]` env-gated marker |
| `src/cli/main.cpp:50-59` | `[total_time]` env-gated marker |
| `tests/test_ibp_kira.cpp:170-200` | propagator yaml format pin |
| `notes/future_optimization_proposals.md` | 2026-05-20 entry closed; only `2026-05-19 sectormappings reuse` proposal remains open |

## Open follow-ups (not actionable this session)

1. **`photon_4L_SE_TwoBubbles_mass1` masters count 18 vs MMA 16**:
   HANDOFF-era figures showed C++ produced 2 spurious masters on this
   bench.  This session's fixes very likely fixed it (sub-system
   structure is now MMA-aligned), but a 2.5-hour re-run of mass1
   would be needed to confirm.  Numerical output already passes
   oracle, so this is a "verify reduction in spurious work" task,
   not a correctness bug.
2. **`sectormappings/` reuse across `system_X_diffeq` runs of the
   same family** (the still-open 2026-05-19 entry in
   `notes/future_optimization_proposals.md`).  η-independent topology
   tables could be reused, saving ~30 min/system on 4L families.
   Pending: correctness review of η-independence claim + dedicated
   benchmark host.

## How to reproduce the verification

```bash
# Build
cmake --build build -j$(nproc)
ctest -j$(nproc)   # expect 549/549 PASS

# Per-system intermediate-result match on sunset_bubble_4L_2leg_eqmass
# (~12 min parallel)
export AMF_BENCH_MMA_CACHE=/tmp/verify_mma_cache
CPP_WD=/tmp/verify_cpp_work
rm -rf "$AMF_BENCH_MMA_CACHE" "$CPP_WD"
python3 -c "import json; d=json.load(open('tools/bench/sunset_bubble_4L_2leg_eqmass_eps001_black_box_amflow_cpp.json')); d['work_dir']='$CPP_WD'; json.dump(d, open('/tmp/verify_cpp.json','w'))"
( wolframscript -file tools/bench/sunset_bubble_4L_2leg_eqmass_eps001_black_box_amflow_mma.wl > /tmp/verify_mma.out 2>&1 ) &
( ./build/src/cli/amflow_cli /tmp/verify_cpp.json /tmp/verify_cpp.out.json > /tmp/verify_cpp.out 2>&1 ) &
wait
./tools/bench/compare_sampled.py tools/bench/sunset_bubble_4L_2leg_eqmass_eps001_mma_reference.json /tmp/verify_cpp.out.json --rel-tol 1e-28

# Stage timing on bubble_1L_diffmass (~1 min parallel)
rm -rf /tmp/stage_cpp /tmp/stage_mma
python3 -c "import json; d=json.load(open('tools/bench/bubble_1L_diffmass_eps01_black_box_amflow_cpp.json')); d['work_dir']='/tmp/stage_cpp'; json.dump(d, open('/tmp/stage.json','w'))"
( AMF_BENCH_MMA_CACHE=/tmp/stage_mma wolframscript -file tools/bench/bubble_1L_diffmass_eps01_black_box_amflow_mma.wl > /tmp/stage_mma.out 2>&1 ) &
( AMFLOW_STAGE_TIMING=1 ./build/src/cli/amflow_cli /tmp/stage.json /tmp/stage_cpp.out.json > /tmp/stage_cpp.out 2> /tmp/stage_cpp.err ) &
wait
grep "(kira_time|total_time)" /tmp/stage_cpp.err
grep "in [0-9]+s|finished in" /tmp/stage_mma.out
```
