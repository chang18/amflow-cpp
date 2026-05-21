// SPDX-License-Identifier: MIT
// ibp::reduce — implementation (flattened BlackBox; Kira-only).
//

#include "amflow/ibp/reduce.hpp"

#include "amflow/algebra/numeric_subst.hpp"

#include <algorithm>
#include <filesystem>
#include <initializer_list>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_mpoly_q.h>

namespace amflow::ibp {

using algebra::Mfrac;
using algebra::Mpoly;
using algebra::MpolyContext;

namespace fs = std::filesystem;

ReductionContext make_reduction_context(const qft::FamilyConfig& /*fc*/) {
    // Narrow the reduction MpolyContext to `{eta, d}` — the only
    // variables Kira's output rationals carry symbolically (numeric
    // kinematics are baked into the yaml by `kira_yaml.cpp`'s
    // `NumericSubs`).  Declaring the full `fc.ctx` here drags 9+
    // unused symbols through every FLINT multivariate-GCD inside
    // `fmpz_mpoly_q_canonicalise`, ballooning memory on multi-mass
    // families.
    //
    // This mirrors MMA: upstream `DifferentialEquation`
    // (`Kira/interface.m:505`) stores the diffeq tensor as symbolic
    // Mathematica expressions where unused symbols never enter the
    // algebra.  Variables in `fc.ctx` but not in `red_ctx` (loop
    // momenta, kinematic invariants, masses) must be substituted to
    // numeric values before any Mfrac is lifted — see
    // `substitute_fc_vars`.
    std::vector<std::string> names = {"eta", "d"};
    auto rctx = std::make_shared<MpolyContext>(std::move(names));

    ReductionContext out;
    out.ctx = rctx;
    out.d_var = 1;
    return out;
}

namespace {

using algebra::substitute_fc_vars;
using algebra::mfrac_to_ctx_lenient;

inline Mfrac mfrac_to_ctx(const Mfrac& src,
                           const std::shared_ptr<MpolyContext>& dst_ctx) {
    return algebra::mfrac_to_ctx_lenient(src, dst_ctx);
}

Mfrac lift_to_red_ctx(const Mfrac& src, const ReductionContext& red_ctx) {
    return mfrac_to_ctx(src, red_ctx.ctx);
}

KiraConfig make_kira_config(const qft::FamilyConfig& fc,
                              const std::vector<int>& top_pattern,
                              const ReduceOptions& opts) {
    KiraConfig cfg;
    cfg.fc = &fc;
    cfg.top_pattern = top_pattern;
    cfg.ibp_rank = opts.ibp_rank;
    cfg.ibp_dot  = opts.ibp_dot;
    cfg.n_thread = opts.n_thread;
    cfg.integral_order = opts.integral_order;
    cfg.numeric_values = opts.numeric_values;
    cfg.kira_executable = opts.kira_executable;
    cfg.fermat_executable = opts.fermat_executable;
    return cfg;
}

ReduceOptions
apply_jdot_jrank_floor(const ReduceOptions& opts,
                       std::initializer_list<const std::vector<qft::JIntegral>*> lists,
                       bool dot_plus_one) {
    ReduceOptions out = opts;
    for (const auto* lst : lists) {
        for (const auto& j : *lst) {
            long d = j.n_dots() + (dot_plus_one ? 1 : 0);
            long r = j.rank();
            if (d > out.ibp_dot)  out.ibp_dot  = d;
            if (r > out.ibp_rank) out.ibp_rank = r;
        }
    }
    return out;
}

// Kira's raw master list is already in MMA-compatible ascending-by-
// sector order; MMA's `Preferred[preferred, dir]`
// (`Kira/interface.m:164-170`) passes it on unchanged.  This function
// is intentionally a no-op — kept as a hook so call sites stay
// visible if a future Kira release produces a non-MMA-compatible raw
// order.
void sort_integrals_like_amflow(std::vector<qft::JIntegral>& /*masters*/) {
}

}  // namespace

ReduceResult
reduce(const qft::FamilyConfig& fc,
       const std::vector<qft::JIntegral>& targets,
       const std::vector<qft::JIntegral>& preferred,
       const std::vector<int>& top_pattern,
       const ReduceOptions& opts) {
    ReduceResult out;
    out.red_ctx = make_reduction_context(fc);

    if (targets.empty()) return out;

    if (top_pattern.size() != fc.propagators_after_conservation.size()) {
        throw std::invalid_argument(
            "ibp::reduce: top_pattern length != #propagators");
    }

    ReduceOptions opts_eff =
        apply_jdot_jrank_floor(opts, {&targets, &preferred}, false);

    KiraConfig cfg = make_kira_config(fc, top_pattern, opts_eff);

    // Mirror upstream `BlackBoxReduce` (Kira/interface.m): a two-Kira-call
    // sequence at the same `(rank, dot)`.  First call (`IBPSystem`) is
    // `select_mandatory_recursively` for a sector-wide master
    // enumeration; second call (`AnalyticReduction`) is
    // `select_mandatory_list` for the specific targets.  Upstream
    // shares one `$ReductionDirectory` between the two calls; our
    // Kira 2.x version refuses to do that because the Masters-mode
    // `run_initiate: masters` doesn't register the `-s` numeric
    // substitutions, so the subsequent Reduce-mode call sees a stale
    // auxiliary state and Aborts with `Kira::update_auxiliary_file:
    // Last Kira run set 0 variables to numeric values, this time you
    // request N.`  We therefore run the two calls in *separate*
    // subdirectories (`masters_preheat/` and `target_reduce/`) — both
    // at the same `(rank, dot)`, so Kira returns consistent master
    // lists — and bridge them with the SubsetQ check below mirroring
    // upstream's
    // `If[!SubsetQ[masters, str], Abort["inconsistent masters from Kira"]]`.
    const std::string masters_dir = opts.work_dir + "/masters_preheat";
    const std::string reduce_dir  = opts.work_dir + "/target_reduce";

    // === Step 1: Masters-mode preheat (sector-wide enumeration) ===
    fs::remove_all(masters_dir);
    fs::create_directories(masters_dir);
    kira_write_config(cfg, masters_dir);
    kira_write_preferred(preferred, fc, masters_dir);
    kira_write_jobs(cfg, masters_dir, KiraReductionMode::Masters);
    kira_run(cfg, masters_dir, opts.log_file);

    out.masters = kira_read_masters(fc, masters_dir);
    if (out.masters.empty()) {
        // Empty masters_mma means Kira identified the sector as scaleless
        // (or otherwise had no integrals to enumerate at the chosen
        // (rank, dot)).  Mirror upstream `AnalyticReduction`
        // (`Kira/interface.m:461`): `If[Length[masters]===0, Return[{{},{}}]]`.
        // Caller (AMFSystem::build_boundary) treats empty `red.masters` as
        // "boundary integrand reduces to zero, skip family".
        out.rules.reserve(targets.size());
        for (std::size_t i = 0; i < targets.size(); ++i) out.rules.emplace_back();
        return out;
    }
    sort_integrals_like_amflow(out.masters);

    // Snapshot the canonical (preheat) master set for the SubsetQ check.
    auto sector_master_keys = std::set<std::string>();
    for (const auto& m : out.masters) {
        std::ostringstream key;
        key << m.family();
        for (long v : m.indices()) key << "|" << v;
        sector_master_keys.insert(key.str());
    }

    // === Step 2: Reduce-mode (target reduction in a sibling subdir) ===
    //
    // **Preferred file**: mirror MMA `AnalyticReduction`, which reuses
    // the preferred file IBPSystem wrote (the *input* `preferred`).
    // Do NOT write the sorted preheat output as preferred here; it
    // changes Kira's elimination order and breaks the resulting
    // reduction rule structure.
    fs::remove_all(reduce_dir);
    fs::create_directories(reduce_dir);
    kira_write_config(cfg, reduce_dir);
    kira_write_preferred(preferred, fc, reduce_dir);  // mirror MMA: reuse input preferred
    kira_write_targets(targets, fc, reduce_dir);
    kira_write_jobs(cfg, reduce_dir, KiraReductionMode::Reduce);
    kira_run(cfg, reduce_dir, opts.log_file);

    // SubsetQ guard: Reduce-mode masters file must be ⊆ Masters-mode
    // sector enumeration.  Both Kira calls ran at the same
    // `(rank, dot)`, so by construction this holds; firing the throw
    // indicates Kira-side nondeterminism.
    auto reduce_masters = kira_read_masters(fc, reduce_dir);
    for (const auto& m : reduce_masters) {
        std::ostringstream key;
        key << m.family();
        for (long v : m.indices()) key << "|" << v;
        if (sector_master_keys.find(key.str()) == sector_master_keys.end()) {
            throw std::runtime_error(
                "ibp::reduce: inconsistent masters from Kira "
                "(Reduce-mode master '" + m.to_string() + "' is absent "
                "from the Masters-mode sector enumeration at the same "
                "(rank, dot); mirrors upstream "
                "`AnalyticReduction: inconsistent masters from Kira`).");
        }
    }

    auto raw_rules = kira_read_target_table(fc, reduce_dir);

    std::map<std::string, std::size_t> master_index;
    for (std::size_t i = 0; i < out.masters.size(); ++i) {
        std::ostringstream key;
        key << out.masters[i].family();
        for (long v : out.masters[i].indices()) {
            key << "|" << v;
        }
        master_index[key.str()] = i;
    }

    out.rules.reserve(targets.size());
    for (std::size_t i = 0; i < targets.size(); ++i) {
        out.rules.emplace_back();
    }

    auto target_key = [](const qft::JIntegral& j) {
        std::ostringstream k;
        k << j.family();
        for (long v : j.indices()) k << "|" << v;
        return k.str();
    };
    std::map<std::string, std::size_t> target_index;
    for (std::size_t i = 0; i < targets.size(); ++i) {
        target_index[target_key(targets[i])] = i;
    }

    for (const auto& kr : raw_rules) {
        auto it = target_index.find(target_key(kr.lhs));
        if (it == target_index.end()) {
            continue;
        }
        std::size_t i = it->second;
        for (const auto& [coef_str, rhs_j] : kr.rhs) {
            // Mirror Kira/interface.m:486 (AnalyticReduction): every
            // J integral on the RHS of a Kira reduction must be in
            // the master list.  Upstream Aborts with "wrong reduction
            // results from Kira."; we throw with a more diagnostic
            // message so the offending integral is identifiable.
            if (master_index.find(target_key(rhs_j)) == master_index.end()) {
                std::ostringstream msg;
                msg << "ibp::black_box_reduce: Kira returned an RHS J "
                       "integral that is not in the master list (target="
                    << target_key(kr.lhs) << ", offending rhs="
                    << target_key(rhs_j)
                    << ").  This indicates an inconsistent Kira reduction.";
                throw std::runtime_error(msg.str());
            }
            DerivTerm dt;
            dt.coef = kira_parse_expression(out.red_ctx.ctx, coef_str);
            dt.integ = rhs_j;
            out.rules[i].push_back(std::move(dt));
        }
    }

    for (std::size_t i = 0; i < targets.size(); ++i) {
        if (!out.rules[i].empty()) continue;
        auto mit = master_index.find(target_key(targets[i]));
        if (mit != master_index.end()) {
            DerivTerm dt;
            dt.coef = Mfrac::one(out.red_ctx.ctx);
            dt.integ = targets[i];
            out.rules[i].push_back(std::move(dt));
        }
    }

    return out;
}

DiffeqResult
diffeq(const qft::FamilyConfig& fc,
       const std::vector<qft::JIntegral>& jpreferred,
       const std::vector<std::string>& vars,
       const std::vector<int>& top_pattern,
       const ReduceOptions& opts) {
    DiffeqResult out;
    out.red_ctx = make_reduction_context(fc);
    out.vars = vars;

    if (jpreferred.empty()) return out;

    ReduceOptions opts_eff =
        apply_jdot_jrank_floor(opts, {&jpreferred}, true);

    KiraConfig cfg = make_kira_config(fc, top_pattern, opts_eff);

    // Mirror upstream `BlackBoxDiffeq` (Kira/interface.m): a two-Kira-call
    // sequence at the same `(rank, dot)`.  Same separate-subdir
    // pattern as in `ibp::reduce` above (see the matching block for
    // the rationale).  Preheat is the sector-wide master enumeration;
    // the second call reduces the libp_deriv-derived integral list
    // against those masters.  Both Kira runs use `opts_eff`'s
    // `(rank, dot)`; re-flooring over `all_ints` for the inner call
    // would break master-count agreement with upstream.
    const std::string masters_dir = opts.work_dir + "/masters_preheat";
    const std::string reduce_dir  = opts.work_dir + "/target_reduce";

    // === Step 1: Masters-mode preheat (sector-wide enumeration) ===
    fs::remove_all(masters_dir);
    fs::create_directories(masters_dir);
    kira_write_config(cfg, masters_dir);
    kira_write_preferred(jpreferred, fc, masters_dir);
    kira_write_jobs(cfg, masters_dir, KiraReductionMode::Masters);
    kira_run(cfg, masters_dir, opts.log_file);

    auto masters = kira_read_masters(fc, masters_dir);
    if (masters.empty()) {
        throw std::runtime_error("ibp::diffeq: no masters from Kira");
    }
    sort_integrals_like_amflow(masters);
    out.sortedmasters = masters;

    auto int_key = [](const qft::JIntegral& j) {
        std::ostringstream k; k << j.family();
        for (long v : j.indices()) k << "|" << v;
        return k.str();
    };

    // Snapshot the preheat master set for the SubsetQ check below.
    std::set<std::string> preheat_master_keys;
    for (const auto& m : masters) preheat_master_keys.insert(int_key(m));

    // Compute libp_deriv per master per requested variable.  Mirrors
    // upstream `der = ComputeDerivative[masters, #]&/@vars`.
    //
    // Two memory/perf optimisations applied here:
    //   (a) reproject each per-term Mfrac to the narrow red_ctx
    //       ({eta, d}) immediately after libp_deriv returns and BEFORE
    //       `simplify_terms` runs.  Numeric substitution alone is not
    //       enough — without reprojection the Mfrac stays on `fc.ctx`
    //       (10+ variables for multi-mass 3L families) and FLINT runs
    //       full-ctx multivariate-GCD on every subsequent operation.
    //   (b) hoist `libp_denoms_deriv` out of the per-master loop: `dd`
    //       depends only on (fc, var), so compute it once per `var`
    //       instead of `masters.size()` times.  Per-master `libp_deriv`
    //       is then inlined using the shared `dd`, mirroring
    //       `src/ibp/libp_deriv.cpp::libp_deriv`'s body, and
    //       `neg_ak_cache[ak]` amortises `Mfrac::from_si` allocations.
    long N_props = (long)fc.propagators_after_conservation.size();
    std::vector<std::vector<std::vector<DerivTerm>>> der(vars.size());
    for (std::size_t v = 0; v < vars.size(); ++v) {
        auto dd = libp_denoms_deriv(fc, vars[v],
                                     opts.numeric_values,
                                     out.red_ctx.ctx);
        long n_d_local = (long)dd.completede.size();
        std::map<long, Mfrac> neg_ak_cache;

        der[v].reserve(masters.size());
        for (const auto& m : masters) {
            if ((long)m.n_indices() != N_props) {
                throw std::invalid_argument(
                    "ibp::diffeq: master indices length != #propagators");
            }
            std::vector<DerivTerm> raw;
            for (long k = 0; k < N_props; ++k) {
                long ak = m.indices()[(std::size_t)k];
                if (ak == 0) continue;

                auto cache_it = neg_ak_cache.find(ak);
                if (cache_it == neg_ak_cache.end()) {
                    cache_it = neg_ak_cache.emplace(
                        ak, Mfrac::from_si(out.red_ctx.ctx, -ak)).first;
                }
                const Mfrac& neg_ak = cache_it->second;

                for (long jp = 0; jp < n_d_local; ++jp) {
                    const Mfrac& cf = dd.coef[(std::size_t)k][(std::size_t)jp];
                    if (cf.is_zero()) continue;
                    std::vector<long> new_idx = m.indices();
                    new_idx[(std::size_t)k] += 1;
                    if (jp < N_props) new_idx[(std::size_t)jp] -= 1;
                    qft::JIntegral shifted(m.family(), new_idx);
                    DerivTerm dt;
                    dt.integ = shifted;
                    dt.coef = neg_ak.clone();
                    dt.coef *= cf;
                    raw.push_back(std::move(dt));
                }

                const Mfrac& ck = dd.constant[(std::size_t)k];
                if (!ck.is_zero()) {
                    std::vector<long> new_idx = m.indices();
                    new_idx[(std::size_t)k] += 1;
                    qft::JIntegral shifted(m.family(), new_idx);
                    DerivTerm dt;
                    dt.integ = shifted;
                    dt.coef = neg_ak.clone();
                    dt.coef *= ck;
                    raw.push_back(std::move(dt));
                }
            }
            auto simp = simplify_terms(std::move(raw));
            der[v].push_back(std::move(simp));
        }
    }

    // Collect the unique integrals appearing in any derivative term,
    // plus the masters themselves (so the Reduce-mode call's target
    // list covers the basis).  Mirrors upstream
    // `integrals = Cases[der, _?(Head[#]===Symbol["j"]&), Infinity] //
    // DeleteDuplicates`.
    std::vector<qft::JIntegral> all_ints;
    {
        std::set<std::string> seen;
        for (const auto& vec : der) {
            for (const auto& terms : vec) {
                for (const auto& t : terms) {
                    std::string k = int_key(t.integ);
                    if (seen.insert(k).second) all_ints.push_back(t.integ);
                }
            }
        }
        // Do not add `masters` to `all_ints`.  MMA's
        // `DifferentialEquation` (`Kira/interface.m:500`) only sends
        // `Cases[der, j[...], Infinity] // DeleteDuplicates` to
        // `AnalyticReduction`.  Adding preheat masters here inflates
        // the Reduce target list and Kira's memory footprint without
        // gain — identity rules for masters-of-themselves are filled
        // in by the loop below (after raw_rules is parsed).
    }

    // === Step 2: Reduce-mode (target reduction in a sibling subdir) ===
    // Same `(rank, dot)` as preheat — NO re-floor over `all_ints`.
    // The re-floor used to be applied implicitly by the nested
    // `reduce()` call's own `apply_jdot_jrank_floor(opts, {targets,
    // preferred}, false)`; bypassing the nested call and writing the
    // Reduce yaml directly here preserves `opts_eff`'s `(rank, dot)`.
    //
    // **Preferred file**: mirrors MMA `AnalyticReduction`
    // (Kira/interface.m:458-489), which **does not call `Preferred[]`**
    // -- it reuses the file IBPSystem wrote (the *input* `preferred`).
    // Writing the sorted-preheat output (228 masters) here instead of
    // the original input (172 masters) changes Kira's elimination
    // order and produces a different reduction rule (the dotted
    // self-referencing form with 49+ RHS masters seen on pentabox).
    // That mismatched rule structure cascades through the DE matrix
    // construction and breaks pentabox top-sector master values by
    // ~10^12.  Pentagon/wzbox/hexagon/mercedes/sunset aren't sensitive
    // because their preheat masters happen to coincide with the input
    // preferred ordering closely enough.
    fs::remove_all(reduce_dir);
    fs::create_directories(reduce_dir);
    kira_write_config(cfg, reduce_dir);
    kira_write_preferred(jpreferred, fc, reduce_dir);  // mirror MMA: reuse input preferred
    kira_write_targets(all_ints, fc, reduce_dir);
    kira_write_jobs(cfg, reduce_dir, KiraReductionMode::Reduce);
    kira_run(cfg, reduce_dir, opts.log_file);

    // SubsetQ guard mirrors upstream's
    // `If[!SubsetQ[masters, str], Abort["inconsistent masters from Kira"]]`.
    auto reduce_masters = kira_read_masters(fc, reduce_dir);
    for (const auto& m : reduce_masters) {
        if (preheat_master_keys.find(int_key(m)) == preheat_master_keys.end()) {
            throw std::runtime_error(
                "ibp::diffeq: inconsistent masters from Kira "
                "(Reduce-mode master '" + m.to_string() + "' is absent "
                "from the Masters-mode sector enumeration at the same "
                "(rank, dot); mirrors upstream "
                "`AnalyticReduction: inconsistent masters from Kira`).");
        }
    }

    auto raw_rules = kira_read_target_table(fc, reduce_dir);

    std::map<std::string, std::size_t> master_index;
    for (std::size_t i = 0; i < masters.size(); ++i) {
        master_index[int_key(masters[i])] = i;
    }

    // Parse `raw_rules` (target_table.m) into a lhs-keyed map of
    // DerivTerm lists.  Identity rule for any `all_ints` element that
    // is itself a master and didn't appear as an explicit rule LHS.
    std::map<std::string, std::vector<DerivTerm>> rule_by_lhs;
    for (const auto& kr : raw_rules) {
        std::vector<DerivTerm> dts;
        dts.reserve(kr.rhs.size());
        for (const auto& [coef_str, rhs_j] : kr.rhs) {
            // Mirror Kira/interface.m:486 (AnalyticReduction RHS-in-master
            // guard); also the C++ D6 fix.  Now that the SubsetQ guard
            // above passes, RHS-not-in-master can only fire on a Kira
            // internal inconsistency.
            if (master_index.find(int_key(rhs_j)) == master_index.end()) {
                std::ostringstream msg;
                msg << "ibp::diffeq: Kira returned an RHS J integral "
                       "that is not in the master list (target="
                    << int_key(kr.lhs) << ", offending rhs="
                    << int_key(rhs_j)
                    << ").  This indicates an inconsistent Kira reduction.";
                throw std::runtime_error(msg.str());
            }
            DerivTerm dt;
            dt.coef = kira_parse_expression(out.red_ctx.ctx, coef_str);
            dt.integ = rhs_j;
            dts.push_back(std::move(dt));
        }
        rule_by_lhs[int_key(kr.lhs)] = std::move(dts);
    }
    for (const auto& integ : all_ints) {
        auto k = int_key(integ);
        if (rule_by_lhs.count(k)) continue;
        auto mit = master_index.find(k);
        if (mit != master_index.end()) {
            DerivTerm dt;
            dt.coef = Mfrac::one(out.red_ctx.ctx);
            dt.integ = integ;
            // `vector<DerivTerm>{std::move(dt)}` would go through
            // initializer_list which copies — DerivTerm holds a
            // move-only Mfrac, so the copy is deleted.  Build the
            // single-element vector via push_back instead.
            std::vector<DerivTerm> identity_rule;
            identity_rule.push_back(std::move(dt));
            rule_by_lhs[k] = std::move(identity_rule);
        }
    }

    // Batched lcm-sum per matrix entry, mirroring MMA's
    // `Together[der/.j->red]` at `Kira/interface.m:505`.  Pairwise
    // `out.diffeq[v][row][col] += product` would canonicalise each
    // addition (every step ~3× operand-size temporaries inside FLINT);
    // collecting all (num_i, den_i) pairs first, computing
    // `L = lcm(den_i)` once, scaling, summing, then canonicalising
    // once cuts that to a single canonicalisation per entry.
    auto batched_sum = [](std::vector<Mfrac>& terms,
                          const std::shared_ptr<MpolyContext>& ctx) -> Mfrac {
        if (terms.empty()) return Mfrac::zero(ctx);
        if (terms.size() == 1) return std::move(terms.front());
        // L = lcm of all denominators.
        // lcm(a, b) = a / gcd(a, b) * b  (uses exact_divide for the
        // division step, since gcd divides a by definition).
        Mpoly L = terms[0].denominator();
        for (std::size_t i = 1; i < terms.size(); ++i) {
            Mpoly g = Mpoly::gcd(L, terms[i].denominator());
            Mpoly L_div_g(ctx);
            if (!Mpoly::exact_divide(L_div_g, L, g)) {
                throw std::runtime_error(
                    "ibp::diffeq batched_sum: gcd did not divide lcm "
                    "(should be impossible)");
            }
            L = L_div_g * terms[i].denominator();
        }
        // N_total = sum_i (num_i * (L / den_i)).
        Mpoly N_total = Mpoly::zero(ctx);
        for (auto& t : terms) {
            Mpoly scale(ctx);
            if (!Mpoly::exact_divide(scale, L, t.denominator())) {
                throw std::runtime_error(
                    "ibp::diffeq batched_sum: term denominator does not "
                    "divide lcm (should be impossible)");
            }
            N_total += t.numerator() * scale;
        }
        return Mfrac(std::move(N_total), std::move(L));  // canonicalises once
    };

    out.diffeq.reserve(vars.size());
    for (std::size_t v = 0; v < vars.size(); ++v) {
        out.diffeq.emplace_back();
        out.diffeq[v].reserve(masters.size());
        for (std::size_t row = 0; row < masters.size(); ++row) {
            // Collect all products per column for this row, deferring
            // canonicalisation to one batched_sum call per column.
            std::map<std::size_t, std::vector<Mfrac>> per_col_products;
            for (const auto& dt : der[v][row]) {
                auto rit = rule_by_lhs.find(int_key(dt.integ));
                if (rit == rule_by_lhs.end()) continue;
                for (const auto& mt : rit->second) {
                    auto mit = master_index.find(int_key(mt.integ));
                    if (mit == master_index.end()) continue;
                    Mfrac product = dt.coef.clone();
                    product *= mt.coef;
                    per_col_products[mit->second].push_back(std::move(product));
                }
            }
            out.diffeq[v].emplace_back();
            out.diffeq[v][row].reserve(masters.size());
            for (std::size_t col = 0; col < masters.size(); ++col) {
                auto it = per_col_products.find(col);
                if (it == per_col_products.end()) {
                    out.diffeq[v][row].push_back(Mfrac::zero(out.red_ctx.ctx));
                } else {
                    out.diffeq[v][row].push_back(
                        batched_sum(it->second, out.red_ctx.ctx));
                }
            }
        }
    }

    return out;
}

}  // namespace amflow::ibp
