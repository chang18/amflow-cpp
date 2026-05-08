# notes/ — MMA source maps

This directory holds the symbol → file:line → semantics maps for the
upstream Mathematica AMFlow source
(<https://gitlab.com/multiloop-pku/amflow>):

- `mma_amflow_map.md` — symbol map for `AMFlow.m`.
- `mma_desolver_map.md` — symbol map for `diffeq_solver/DESolver.m`.
- `mma_kira_interface_map.md` — symbol map for `ibp_interface/Kira/interface.m`.

These supplement [`docs/REFERENCE_MAP.md`](../docs/REFERENCE_MAP.md) by
giving line-resolved navigation inside each MMA file.

## Conventions

- `<symbol>` (file:line) — function definition site.
- `→ <symbol>` — calls into.
- `← <symbol>` — called by.
- `[mma:AMFlow.m:1031]` style cross-reference back into MMA.

Sections within each file go top-down through the source, so file:line
ranges always increase.
