# Roadmap

This document records the forward-looking development plan.  Historical
phase logs (Phase 1A/1B, Phase 2A/B/C, Phase 3.A–J, Phase 3 batch-2,
Phase 4.A Docker) are retired — see [`CHANGELOG.md`](../CHANGELOG.md)
for the per-release narrative and
[`docs/AUDIT_MMA_PARITY.md`](AUDIT_MMA_PARITY.md) for the line-level
parity audit.

---

## Active practice — diversity-driven oracle expansion

The audit is closed: **86 🟢 / 0 🟡 / 11 🔴 fixed / 1 🔴 out of scope
(D5, see below) / 17 ⚪ intentional**.  Going forward the project's
quality loop is **oracle diversification**: build new benchmarks
along axes the existing suite under-covers, and fix anything that
diverges from MMA inline.

The 37 committed oracle triplets currently cover:

| Diversity axis | Current coverage |
|---|---|
| Loop number          | L = 1, 2, 3, 4 |
| Kinematic invariants | 1 (`s`), 2 (`s, t`), 4 (`s, t, mWsq, mZsq`) |
| External legs        | 2, 3, 4, 5, **6** (pentagon 1L, pentagon 1L massive, pentabox 2L, hexagon 1L) |
| Mass configurations  | all-massless, all-equal, mixed at L = 1 / 2 / 3 / 4 |
| Cutkosky cuts        | 2-, 3-, 4-, 5-particle cuts |
| ε extremes           | `1/2`, `1/100`, `1/1000`, `1/10000` |
| Sector size          | ≤ 8 propagators (pentabox 2L corner) |

Target additions: higher-loop multi-invariant cases, larger sectors
(≥ 10 propagators), and any topology a downstream user reports.

Each new oracle is roughly one day of work (construct family → MMA
reference → C++ run → commit triplet).  When a new oracle exposes a
bug, fix it inline; the bug becomes a parity-validated 🟢 going
forward.

---

## Out of scope — won't be implemented

These are deliberate non-goals.  The JSON entry-point rejects them
with a clear error rather than silently mishandling.

- **`SolveIntegralsGaugeLink`**, HQET / SCET / Wilson-line workflows.
- **IBP backends other than Kira** — no FIRE / LiteRed / FiniteFlow /
  Blade.
- **Complex-valued numeric kinematics** (audit divergence D5).  The
  C++ algebra layer is over Q (FLINT `fmpz_mpoly_q_t`) and does not
  carry a complex coefficient ring.  The JSON dispatcher rejects the
  object form `{"re":..,"im":..}` in
  `amf_options.blackbox.numeric_values`.  Workaround: provide only
  purely-real numeric values.  See
  [`docs/FAQ.md`](FAQ.md) "What's *not* implemented?" and
  [`docs/AUDIT_MMA_PARITY.md`](AUDIT_MMA_PARITY.md) §D5.

---

## Update policy

- Each new oracle triplet (cpp.json + mma.wl + mma_reference.json)
  is added to `run_perf_audit.sh` rotation on landing.
- When a new oracle exposes a real divergence, raise it in
  [`docs/AUDIT_MMA_PARITY.md`](AUDIT_MMA_PARITY.md) §2 as a new D-row,
  fix in the same session if practical, and reference the closing
  commit.
- Do not silently delete items that turn out to be infeasible —
  document the reason in the audit instead.
