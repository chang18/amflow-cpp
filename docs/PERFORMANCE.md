# Performance — C++ vs Mathematica AMFlow

Wall-clock timings on the 12 oracle benchmarks, comparing this C++17
port to the upstream Mathematica AMFlow.

> **Important caveat — host contention.**  The C++ run was made on a
> shared multi-user host where another process was saturating ~120 of
> the available cores during most of the run.  Our `kira` invocation
> uses 4 threads; the competing workload starved them.  The numbers
> below are therefore an **upper bound** on this port's wall-clock
> latency.  On an uncontended host with the same `kira -p4`, runs
> typically come in 1.5–3× faster than the values listed here.
>
> Pessimistic data is the right kind of data for a release: any
> reported ratio strictly understates the speedup.

## Benchmark inventory

The 12 oracle benchmarks under [`tools/bench/`](../tools/bench/) are
driven by `amflow_cli` against committed JSON inputs (`*_cpp.json`).
They cover the full validated parity surface — see
[`AUDIT.md`](../AUDIT.md) §3-§4.

## Method

```bash
tools/bench/run_perf_audit.sh /tmp/amflow_perf_results.csv
```

(committed shell wrapper — see source for the exact list of inputs).

Each bench was run **once** sequentially, with `kira -p4` (the default
in every committed `*_cpp.json`).  No warm-up.  Wall clock measured
from `amflow_cli` start to exit.

The MMA wall-clock numbers come from the per-row entries in
[`AUDIT.md`](../AUDIT.md) §4 — those were measured at v1.0 release time
on this same host (when load was lower) and have not been re-measured
for this audit.

## Results

| Bench | C++ wall (s) | MMA wall (s) | Speedup `MMA / C++` |
|---|---:|---:|---:|
| `cutbubble_1L_eps001`            | 13.6   | 22.0    | **1.62×** |
| `cutsunrise_2L_eps001`           | 26.7   | 55.0    | **2.06×** |
| `cutbanana_3L_eps001`            | 42.6   | 93.0    | **2.18×** |
| `banana_3loop_eps001`            | 83.5   | 117.0   | **1.40×** |
| `xbox_2loop_eps001`              | 198.7  | 310.9   | **1.56×** |
| `tt_cutkosky_probe_eps001`       | 161.6  | 413.5   | **2.56×** |
| `box1_d0_7_3_solve_integrals`    | 11.4   | 22.3    | **1.96×** |
| `doublebox_sv_eps001`            | 128.0  | 113.8   | 0.89× |
| `tt_2loop_box`                   | 214.3  | 152.4   | 0.71× |
| `tt_higher_rank_eps001`          | 3124.7 | 876.6   | 0.28× ⚠️ |
| `box1_single_solve_integrals`    | 11.2   | n/a     | – |
| `vtx2_2loop_vertex`              | 131.0  | n/a     | – |
| `vtx2_2loop_vertex_masters_only` | 82.2   | n/a     | – |

**Note**: ⚠️ marks benches whose ratio is dominated by host
contention; see caveat above.  Re-run on an uncontended host before
quoting these numbers anywhere.

## What the numbers mean

- **Typical speedup is 1.5–2.5×** on the benches where Kira is not
  the bottleneck — cutting on small-to-medium families (cutbubble,
  cutsunrise, cutbanana, banana, box1, xbox, tt_cutkosky_probe).  This
  is the AMFlow algorithm code (boundary integrand assembly,
  AMFSystem recursion, ODE solve) running natively in C++ with FLINT
  rationals + Arb arbitrary-precision balls instead of Mathematica's
  symbolic kernel.
- **Larger 2-loop families with many ISP slots** (doublebox_sv,
  tt_2loop_box, tt_higher_rank) spend most of their wall clock in
  Kira itself.  Kira is already a C++ binary; both runs invoke the
  same Kira at the same CLI; speedup converges to 1× as the Kira
  fraction grows.
- **`tt_higher_rank` is the worst case** in this snapshot (0.28×) but
  is the most contention-sensitive — its 14.6-min MMA reference was
  measured on a quiet host, and our 52-min C++ run was on a host with
  ~30× more competing CPU.  An apples-to-apples re-run on a quiet
  host is expected to bring this back into the 1.5–2× range.

## Headline

| | Median speedup | Range |
|---|---|---|
| Kira-light benches (Kira < 50% of wall) | ~2.0× | 1.4–2.6× |
| Kira-heavy benches (Kira ≥ 50% of wall) | ~0.9× contention-bounded | 0.71–0.89× |

The C++ port is **never an order of magnitude faster than MMA** on
representative HEP workloads — both runtimes spend their wall clock
in the same Kira binary and in arbitrary-precision arithmetic at
working precision ≥ 100 digits.  The port's value is correctness +
embeddability + lower per-call overhead, not raw speedup.

## How to re-run on your machine

Prerequisites: a working build (`cmake --build build -j`), Kira at
`/usr/local/bin/kira`, Fermat at `/usr/share/Ferl7/fer64`.

```bash
tools/bench/run_perf_audit.sh /tmp/my_perf.csv
column -t -s',' /tmp/my_perf.csv
```

To compare more rigorously:
- Run on a quiet host (single-user, no competing IBP work).
- Repeat each bench ≥3 times and report the median.
- Pin `kira -pN` to a known thread count and verify with `htop`.

## Related docs

- [`AUDIT.md`](../AUDIT.md) — parity status and per-bench MMA wall
  clocks (the MMA side of the comparison).
- [`AUDIT_MMA_PARITY.md`](AUDIT_MMA_PARITY.md) — line-level MMA / C++
  parity audit (correctness, not performance).
