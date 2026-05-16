#!/usr/bin/env bash
# Run all committed C++ oracle benchmarks once each, recording wall-clock time.
# Used for the v1.0 performance audit.  Outputs CSV: name,seconds,exit_code.

set -u
cd "$(dirname "$0")/../.."

OUT_CSV="${1:-/tmp/amflow_perf_results.csv}"
echo "name,seconds,exit_code" > "$OUT_CSV"

CLI="./build/src/cli/amflow_cli"
[ -x "$CLI" ] || { echo "missing $CLI; run cmake --build build -j first" >&2; exit 2; }

# Each bench gets a fresh tmp work_dir so Kira-cache state doesn't bias
# adjacent runs.  We respect the work_dir baked into the json input
# unless that is missing (which it shouldn't be for committed benches).

BENCHES=(
  tools/bench/banana_3loop_eps001_black_box_amflow_cpp.json
  tools/bench/banana_4L_mixed_black_box_amflow_cpp.json
  tools/bench/banana_4loop_eps001_black_box_amflow_cpp.json
  tools/bench/bn3_4mass_3L_eps001_black_box_amflow_cpp.json
  tools/bench/hexagon_1L_2mass_eps001_black_box_amflow_cpp.json
  tools/bench/hexagon_1L_3mass_eps001_black_box_amflow_cpp.json
  tools/bench/hexagon_1L_4mass_eps001_black_box_amflow_cpp.json
  tools/bench/pentagon_1L_2mass_eps001_black_box_amflow_cpp.json
  tools/bench/pentagon_1L_4mass_eps001_black_box_amflow_cpp.json
  tools/bench/pentagon_1L_5mass_eps001_black_box_amflow_cpp.json
  tools/bench/bn3mix_eps001_black_box_amflow_cpp.json
  tools/bench/box_1L_2mass_eps001_black_box_amflow_cpp.json
  tools/bench/box_1L_3mass_eps001_black_box_amflow_cpp.json
  tools/bench/box_1L_4mass_eps001_black_box_amflow_cpp.json
  tools/bench/box1_d0_7_3_solve_integrals_cpp.json
  tools/bench/box1_multitarget_black_box_amflow_cpp.json
  tools/bench/box1_single_solve_integrals_cpp.json
  tools/bench/bubble_1L_diffmass_black_box_amflow_cpp.json
  tools/bench/cutbanana_3L_eps001_black_box_amflow_cpp.json
  tools/bench/cutbanana_4L_eps001_black_box_amflow_cpp.json
  tools/bench/cutbubble_1L_eps001_black_box_amflow_cpp.json
  tools/bench/cutbubble_1L_eps2_black_box_amflow_cpp.json
  tools/bench/cutbubble_1L_eps10000_black_box_amflow_cpp.json
  tools/bench/cutsunrise_2L_eps001_black_box_amflow_cpp.json
  tools/bench/dotted_3L_banana_black_box_amflow_cpp.json
  tools/bench/doublebox_sv_eps001_black_box_amflow_cpp.json
  tools/bench/ewbox_1loop_eps001_black_box_amflow_cpp.json
  tools/bench/sunset_2L_3mass_eps001_black_box_amflow_cpp.json
  tools/bench/sunset_2L_2mass_eps001_black_box_amflow_cpp.json
  tools/bench/sunrise_2L_2mass_eps001_black_box_amflow_cpp.json
  tools/bench/sunrise_2L_threshold_black_box_amflow_cpp.json
  tools/bench/sunset_2L_onshell_black_box_amflow_cpp.json
  tools/bench/tadbubble_2L_black_box_amflow_cpp.json
  tools/bench/triangle_1L_2mass_black_box_amflow_cpp.json
  tools/bench/triangle_1L_3mass_black_box_amflow_cpp.json
  tools/bench/tt_2loop_box_black_box_amflow_cpp.json
  tools/bench/tt_cutkosky_probe_eps001_black_box_amflow_cpp.json
  tools/bench/tradcut_phase_2L_eps001_black_box_amflow_cpp.json
  tools/bench/tt_higher_rank_eps001_black_box_amflow_cpp.json
  tools/bench/vacuum_2L_sunrise_black_box_amflow_cpp.json
  tools/bench/vtx2_2L_2mass_eps001_black_box_amflow_cpp.json
  tools/bench/vtx2_2loop_vertex_black_box_amflow_cpp.json
  tools/bench/vtx2_2loop_vertex_masters_only_cpp.json
  tools/bench/xbox_2L_2mass_eps001_black_box_amflow_cpp.json
  tools/bench/xbox_2L_eps10_black_box_amflow_cpp.json
  tools/bench/xbox_2loop_eps001_black_box_amflow_cpp.json
)

for input in "${BENCHES[@]}"; do
    name="$(basename "$input" .json)"
    [ -f "$input" ] || { echo "  skip (missing): $name" >&2; continue; }
    out="/tmp/amflow_perf_${name}_out.json"
    echo "[$(date +%H:%M:%S)] running: $name" >&2
    start=$(date +%s.%N)
    "$CLI" "$input" "$out" >/dev/null 2>"/tmp/amflow_perf_${name}_err.log"
    rc=$?
    end=$(date +%s.%N)
    secs=$(python3 -c "print(f'{$end - $start:.2f}')")
    echo "$name,$secs,$rc" >> "$OUT_CSV"
    echo "  -> $secs s (exit $rc)" >&2
done

echo "" >&2
echo "results written to $OUT_CSV" >&2
echo "" >&2
column -t -s',' "$OUT_CSV"
