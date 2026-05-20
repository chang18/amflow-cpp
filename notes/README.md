# notes/ — MMA source maps + active work items

## MMA source maps

Symbol → file:line → semantics maps for the upstream Mathematica
AMFlow source (<https://gitlab.com/multiloop-pku/amflow>):

- `mma_amflow_map.md` — symbol map for `AMFlow.m`.
- `mma_desolver_map.md` — symbol map for `diffeq_solver/DESolver.m`.
- `mma_kira_interface_map.md` — symbol map for `ibp_interface/Kira/interface.m`.

These supplement [`docs/REFERENCE_MAP.md`](../docs/REFERENCE_MAP.md) by
giving line-resolved navigation inside each MMA file.

## Active work items

- [`HANDOFF_eta_placement.md`](HANDOFF_eta_placement.md) — single-page
  handoff for fixing the C++ η-placement divergence from MMA
  `findbranch`. The 4-loop photon TwoBubbles benchmark is blocked on
  this. Read first if picking up that fix.
- [`future_optimization_proposals.md`](future_optimization_proposals.md)
  — full analysis of (a) the η-placement strategy divergence (root
  cause + code locations + proposed fix routes), and (b) a secondary
  `sectormappings/` caching opportunity unrelated to that.

## Conventions

- `<symbol>` (file:line) — function definition site.
- `→ <symbol>` — calls into.
- `← <symbol>` — called by.
- `[mma:AMFlow.m:1031]` style cross-reference back into MMA.

Sections within each file go top-down through the source, so file:line
ranges always increase.
