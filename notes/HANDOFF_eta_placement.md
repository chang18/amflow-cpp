# HANDOFF — C++ η-injection strategy divergence from MMA

**Created**: 2026-05-20
**Topic**: C++ `amflow_cli` is not faithfully reproducing MMA AMFlow's
η-placement algorithm; on 4-loop tests this manifests as either
catastrophic slowdown or outright failure.

This document is the single-page summary for a new session that will
attempt to align the C++ implementation with MMA's reference. The
detailed analysis lives in `notes/future_optimization_proposals.md`.

---

## TL;DR

In Branch-mode fallback of `qft::amf_position`, C++ chooses a
**different propagator** to attach η to than MMA does. For the staged
mass1/2 photon TwoBubbles 4L test, MMA puts η on the originally-massive
propagator D1 = l1²−msq while C++ puts η on a massless propagator
(D3 = l2²). The two are both numerically correct (η-flow converges to
the same integral value at η=0), but the C++ choice makes Kira's IBP
work explode by ~4 orders of magnitude.

This is **not a numerical correctness bug** but **is an algorithm
implementation faithfulness gap** — per project principle "MMA =
correctness oracle, legacy = design reference" (see memory note
`feedback_legacy_primary_mma_oracle`), C++ should faithfully reproduce
MMA's algorithm, including its Branch-mode tie-break.

---

## Evidence (from photon_4L_SE_TwoBubbles_mass1, dot=5)

| Side | η placement | Kira mandatory list | Kira masters | Kira wallclock | Final result |
|---|---|---|---|---|---|
| MMA `AMFPosition` | prop 0 (l1² with mass) | **18** | 16 | 1857 s | DONE 2h28m |
| C++ `qft::amf_position` | prop 2 (l2² massless) | **477,748** | 154 | 4932.8 s | C++ side hit rc=1 in η-flow later |

The 477,748 mandatory-list explosion is what tanks performance.

For `photon_4L_SE_TwoBubbles_mass2` (2 masses on D1+D2): mass2 C++
ran 13h30m at dot=5 and never finished a single sub-system iteration in
the AMFlow η-flow pipeline; same root cause.

---

## Where to look (audited 2026-05-19)

### C++ side

| Step | File:line | What it does |
|---|---|---|
| Top-level η-position call | `src/pipeline/amfsystem.cpp:3067` | `qft::amf_position(*fc_use, top, opts.amf_modes)` |
| Mode cascade | `src/qft/amfmode.cpp:497-507` | `amf_candidate_component` — falls through Prescription → Mass → Propagator → Branch |
| `single_mass_q` (Mass-mode gate) | `src/qft/amfmode.cpp` (search) | Returns false on symbolic msq → forces Branch fallback |
| **Branch mode (THE diverging code)** | `src/qft/amfmode.cpp:178-199, 239-283` | `same_branch_with` + stable_sort by branch size; tie-break decides which prop wins |
| `info.var` construction (drives tie-break) | `src/qft/topology.cpp:342-397` `make_component` | `ci.var = feynman_vars_in(u0_factor, first_x, n_x)` — FLINT fmpz_mpoly grlex iteration order |
| Debug trace switch | `src/pipeline/amfsystem.cpp:3074-3092` | `AMFLOW_DEBUG_SCHEME=1` emits `pos={...}` decision |

### MMA reference side

| Step | File:line | What it does |
|---|---|---|
| Top-level | `reference/amflow-master/AMFlow.m:1041` | `AMFEtaC@AMFPosition[GetTopPosition[preferred], $AMFMode]` |
| Mode cascade | `reference/amflow-master/AMFlow.m:608-611` | Vacuum / single_mass guard, falls back to Branch on symbolic-mass case |
| **Branch mode** | `reference/amflow-master/AMFlow.m:585-592` | `findbranch[v_]` — this is the function whose algorithm C++ must reproduce. **Needs careful read end-to-end as the precondition for any C++ fix.** |

---

## Two fix routes (only Route A is "alignment")

### Route A — faithful reproduction of MMA `findbranch`  ✅ recommended

Read `findbranch[v_]` line by line, identify what determines the
choice (sort key, iteration order over branches, tie-break), then make
C++'s `same_branch_with` + sort comparator at
`src/qft/amfmode.cpp:178-199, 239-283` produce **the same `pos` for
every input** that MMA does.

May require matching MMA's monomial iteration order on the U
polynomial (because that order seeds `info.var` insertion into the
branch list). If FLINT's grlex differs from MMA's default lexorder,
you may need to either re-sort `info.var` deterministically with
MMA-equivalent ordering, OR change the tie-break to a content-based
criterion that doesn't depend on iteration order.

### Route B — `prefer-mass-bearing-branch` heuristic  ❌ NOT alignment

A subagent recommended adding "prefer branches containing a
non-zero-mass propagator" as secondary tie-breaker. This **happens to
match MMA on the observed cases** but is NOT a faithful reproduction
of MMA's algorithm — `findbranch` has no such explicit rule. Route B
would close the wallclock gap on this test but could re-diverge on
other inputs. Use only as stopgap if Route A is too costly.

---

## Reproduction recipe

The two photon 4L tests are the cheapest reproducers (the η_c 3L tests
are too small to trigger the divergence visibly):

```bash
# 1. Run mass1 (smallest)
./build/src/cli/amflow_cli \
  tools/bench/photon_self_energy/photon_4L_SE_TwoBubbles_mass1_p0_eps001_black_box_amflow_cpp.json \
  /tmp/mass1_out.json
# Compare against MMA reference in
#   tools/bench/photon_self_energy/photon_4L_SE_TwoBubbles_mass1_p0_eps001_mma_reference.json
# Target value: Re ≈ -8.273e+7, Im ≈ -6.58e-61

# 2. Run mass2
./build/src/cli/amflow_cli \
  tools/bench/photon_self_energy/photon_4L_SE_TwoBubbles_mass2_p0_eps001_black_box_amflow_cpp.json \
  /tmp/mass2_out.json
# Target value: Re ≈ 1.658e+8, Im ≈ 1.32e-60
```

Before changing code: turn on `AMFLOW_DEBUG_SCHEME=1` and re-run mass1
**just to capture the `pos` decision** (don't need to wait for full
completion). Expected output: `pos={2}` (the bug). Goal of the fix:
`pos={0}` (matching MMA).

After the fix:
- mass1 should finish in roughly MMA's time (~2-3 hours, allowing for
  C++ vs MMA overhead variance).
- mass1 numerical output should match the reference to 30+ digits.
- mass2 should also finish, in roughly MMA's time (~3-4 hours).

---

## Other test cases (smaller, also useful for non-regression)

Once the fix is in, re-run these to confirm no regression — they
**already passed** before the η-placement code was the limiting factor
(3L tests are too small for the divergence to matter wallclock-wise,
but it could still affect correctness in edge cases):

- `tools/bench/electron_self_energy/electron_3L_SE_{LA3,CL3,TX}_euclid_eps001_*` (committed in `4963725`)
- `tools/bench/electron_self_energy/electron_3L_SE_Mercedes_onshell_eps001_*` (uncommitted as of handoff)
- `tools/bench/eta_c_form_factor/eta_c_3L_Mercedes_*` (committed)
- `tools/bench/ee_mumu_scalar_qed/ee_mumu_2L_{PB,NPB}_*` (committed)

---

## Project state at handoff

### Committed (latest commit `4963725`)
- 3L QED electron self-energy: LA3 + CL3 + TX (off-shell Euclidean)
- 2L ee→μμ scalar QED (PB + NPB physical)
- η_c → γγ 3L Mercedes form factor (staged-mass mass1-3)

### Untracked / uncommitted (to be committed in this handoff)
- `notes/future_optimization_proposals.md` — full η-placement analysis
- `notes/HANDOFF_eta_placement.md` — this document
- `tools/bench/electron_self_energy/electron_3L_SE_Mercedes_onshell_eps001_*` — 3 files (cpp.json, mma.wl, mma_reference.json)
- `tools/bench/photon_self_energy/photon_4L_SE_TwoBubbles_{mass1,mass2}_p0_eps001_*` — 6 files (cpp.json + mma.wl + mma_reference.json per mass case)
- `tools/bench/photon_self_energy/mma_refs/...` — MMA cache (gitignored)

### Deleted as no-oracle, not actionable
- `tools/bench/electron_self_energy/electron_4L_SE_Knot_onshell_eps001_*` (Knot dot=6 test was terminated mid-run; no MMA reference was produced)
- `tools/bench/photon_self_energy/photon_4L_SE_TwoBubbles_p0_eps001_*` (full-mass dot=6 test was terminated; no MMA reference)

### Currently running
- Nothing. Cron `333b0827` cancelled. mass2 C++ user-terminated. Host is idle.

---

## Recommended first actions for next session

1. **Read `findbranch[v_]` source**: `reference/amflow-master/AMFlow.m:585-592`. Until you understand what MMA actually does, no productive C++ fix can be written.
2. **Run mass1 with `AMFLOW_DEBUG_SCHEME=1`** to confirm `pos={2}` (the bug, you don't need to wait for completion — `pos` is decided early).
3. **Trace MMA on the same input** to get its `pos`. AMFlow's verbose mode (search `Print[` in AMFlow.m) will emit it.
4. **Diff the two `pos` decisions**; identify the deciding code path.
5. **Implement Route A** in `src/qft/amfmode.cpp`.
6. **Re-run mass1 → expect 30+-digit match with reference**.
7. **Re-run mass2 → expect 30+-digit match with reference**.
8. **Re-run 3L tests → expect no regression** (committed bench oracle values).
9. **Commit the fix with a regression test** (mass1 numerical match).

Per project memory `feedback_milestone_deep_review`, audit the diff
before commit via an independent agent (don't self-review code that
touches `qft::amf_position` semantics).
