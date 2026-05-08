#!/usr/bin/env python3
"""Side-by-side diff of Layer-6 traces (MMA vs C++) for the 1x1 PowerLaw case.

Inputs:
  build/trace/mma_trace_minimal.jsonl
  build/trace/cpp_trace_minimal.jsonl

Output: a tabular per-step compare, then a highlight of the first divergence.
"""
import json, sys, pathlib, textwrap

ROOT = pathlib.Path(__file__).resolve().parents[2]
MMA = ROOT / "build/trace/mma_trace_minimal.jsonl"
CPP = ROOT / "build/trace/cpp_trace_minimal.jsonl"

def load(p):
    return [json.loads(line) for line in p.read_text().splitlines() if line.strip()]

mma = load(MMA)
cpp = load(CPP)

# Align by semantic meaning.  The trace drivers capture similarly-named
# points; define a mapping from a canonical label to the (fn, phase-or-substep)
# patterns used on each side.
alignment = [
    # label,                 mma filter,                         cpp filter
    ("de",                   [("AMFlow","in","de")],             [("calc_inf","in","de")]),
    ("deinf",                [("CalcInf","in","bcs")],           [("calc_inf","make_deinf","de_inf")]),
    ("bcs reverse",          [("ReverseBCS","out","result")],    [("calc_inf","bc_transform","bcinf_rev")]),
    ("bcs union",            [("UnionBCS","out","result")],      [("calc_inf","bc_transform","bcinf")]),
    ("read_bcs result",      [("ReadBCS","out","result")],       [("calc_inf","region_iter","ini")]),
    ("BuildTaylor output",   [("BuildTaylor","out","result")],   [("calc_taylor_with_ini","m_pure","m_pure")]),
    ("PoincareRank",         [("PoincareRank","out","result")],  [("prepare_taylor_system","in","poincare_rank_m_pure")]),
    ("NHEquations out",      [("NHEquations","out","result")],   [("prepare_taylor_system","symbolic_per_block",None)]),
    ("ConstructMatrix out",  [("ConstructMatrix","out","result")], [("ctpm","construct_matrix_out","cm")]),
    ("SparseGaussian out",   [("SparseGaussian","out","result")],   [("ctpm","sparse_gaussian_out","reduce")]),
    ("CalcTaylor out",       [("CalcTaylor","out","result")],    [("ctpm","out","f0_first5")]),
    ("AMFlow final",         [("AMFlow","out","result_N")],      [("calc_inf","out","first5_coeffs")]),
]

def find_record(records, fn, phase_or_sub, key):
    """Return the first record matching `fn` and (phase==phase_or_sub OR substep==phase_or_sub).
       `key` is optional: if given, return the value at that key; otherwise the full record."""
    for r in records:
        if r.get("fn") != fn: continue
        if r.get("phase") == phase_or_sub or r.get("substep") == phase_or_sub:
            if key is None:
                return r
            return r.get(key)
    return None

def show(label, val, max_len=110):
    if val is None: return "(missing)"
    s = val if isinstance(val, str) else json.dumps(val, ensure_ascii=False)
    if len(s) > max_len: s = s[:max_len] + f"…[+{len(s)-max_len}]"
    return s

print("=" * 120)
print(f"{'LABEL':<28} | {'MMA':<58} | {'C++':<58}")
print("=" * 120)

divergence_at = None
for label, mma_keys, cpp_keys in alignment:
    mma_val = None
    for (fn, ph, key) in mma_keys:
        mma_val = find_record(mma, fn, ph, key)
        if mma_val is not None: break
    cpp_val = None
    for (fn, ph, key) in cpp_keys:
        cpp_val = find_record(cpp, fn, ph, key)
        if cpp_val is not None: break
    ms = show(mma_val, 58); cs = show(cpp_val, 58)
    flag = "  "
    # simple equivalence test: string equality after light normalization
    if mma_val is not None and cpp_val is not None:
        # strip whitespace, lowercase, compare
        ns = lambda s: "".join(s.split()).lower() if isinstance(s, str) else json.dumps(s, sort_keys=True)
        if ns(str(mma_val)) != ns(str(cpp_val)):
            flag = "!!"
            if divergence_at is None:
                divergence_at = label
    print(f"{flag} {label:<25} | {ms:<58} | {cs:<58}")

print("=" * 120)
if divergence_at:
    print(f"\nFirst (possible) divergence at: {divergence_at!r}")
else:
    print("\nNo textual divergence found at the aligned checkpoints.")

# Separately show the key numeric diagnosis for §3 bug.
print("\n--- Key §3 diagnostics ---")
pr_mma = find_record(mma, "PoincareRank", "out", "result")
pr_cpp = find_record(cpp, "prepare_taylor_system", "in", "poincare_rank_m_pure")
print(f"  MMA PoincareRank(BuildTaylor(mat, ini)) = {pr_mma!r}")
print(f"  C++ poincare_rank(m_pure)                = {pr_cpp!r}")

bt_mma = find_record(mma, "BuildTaylor", "out", "result")
mp_cpp = find_record(cpp, "calc_taylor_with_ini", "m_pure", "m_pure")
print(f"  MMA BuildTaylor output mat              = {bt_mma!r}")
print(f"  C++ m_pure (no diagonal correction)     = {mp_cpp!r}")

nh_mma = find_record(mma, "NHEquations", "out", "result")
dxexp_cpp = None
for r in cpp:
    if r.get("fn") == "prepare_taylor_system" and r.get("substep") == "numeric_pre_corr":
        dxexp_cpp = r.get("dxexp")
        break
print(f"  MMA NHEquations(Taylor)                 = {nh_mma!r}")
print(f"  C++ numeric dxexp (pre correction)      = {dxexp_cpp!r}")

# Sparse structure
cm_mma = find_record(mma, "ConstructMatrix", "out", "result")
cm_cpp = find_record(cpp, "ctpm", "construct_matrix_out", "cm")
print(f"\n  MMA ConstructMatrix output (first 200 chars):\n    {show(cm_mma, 200)}")
print(f"  C++ ConstructMatrix output (cm total_columns, first rows):")
if isinstance(cm_cpp, dict):
    print(f"    total_columns={cm_cpp.get('total_columns')}")
    rows = cm_cpp.get('rows', [])[:4]
    for row in rows:
        print(f"    {row}")
