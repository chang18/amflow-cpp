// SPDX-License-Identifier: MIT
// ibp::reduce — implementation (flattened BlackBox; Kira-only).
//

#include "amflow/ibp/reduce.hpp"

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

ReductionContext make_reduction_context(const qft::FamilyConfig& fc) {
    std::vector<std::string> names;
    names.reserve((std::size_t)fc.ctx->n_vars() + 1);
    for (long i = 0; i < fc.ctx->n_vars(); ++i) {
        names.push_back(fc.ctx->var_name(i));
    }
    names.emplace_back("d");
    auto rctx = std::make_shared<MpolyContext>(std::move(names));

    ReductionContext out;
    out.ctx = rctx;
    out.d_var = fc.ctx->n_vars();
    out.lift_gens.reserve((std::size_t)fc.ctx->n_vars());
    for (long i = 0; i < fc.ctx->n_vars(); ++i) {
        out.lift_gens.push_back(Mpoly::variable(rctx, i));
    }
    return out;
}

namespace {

Mfrac mfrac_to_ctx(const Mfrac& src,
                    const std::shared_ptr<MpolyContext>& dst_ctx) {
    if (src.ctx().get() == dst_ctx.get()) return src.clone();

    std::map<std::string, long> name_to_dst;
    for (long i = 0; i < dst_ctx->n_vars(); ++i) {
        name_to_dst[dst_ctx->var_name(i)] = i;
    }
    std::vector<long> src_to_dst((std::size_t)src.ctx()->n_vars(), -1);
    for (long i = 0; i < src.ctx()->n_vars(); ++i) {
        auto it = name_to_dst.find(src.ctx()->var_name(i));
        if (it == name_to_dst.end()) {
            throw std::runtime_error(
                "mfrac_to_ctx: variable '" + src.ctx()->var_name(i)
                + "' not in destination ctx");
        }
        src_to_dst[(std::size_t)i] = it->second;
    }

    auto walk = [&](const Mpoly& p) -> Mpoly {
        Mpoly out(dst_ctx);
        long len = fmpz_mpoly_length(p.raw(), src.ctx()->raw());
        std::vector<unsigned long> exp((std::size_t)src.ctx()->n_vars());
        std::vector<unsigned long> dexp((std::size_t)dst_ctx->n_vars(), 0);
        fmpz_t coeff;
        fmpz_init(coeff);
        for (long t = 0; t < len; ++t) {
            fmpz_mpoly_get_term_exp_ui(exp.data(), p.raw(), t,
                                        src.ctx()->raw());
            fmpz_mpoly_get_term_coeff_fmpz(coeff, p.raw(), t,
                                            src.ctx()->raw());
            std::fill(dexp.begin(), dexp.end(), 0);
            for (long i = 0; i < src.ctx()->n_vars(); ++i) {
                dexp[(std::size_t)src_to_dst[(std::size_t)i]]
                    = exp[(std::size_t)i];
            }
            fmpz_mpoly_set_coeff_fmpz_ui(out.raw(), coeff,
                                          dexp.data(), dst_ctx->raw());
        }
        fmpz_clear(coeff);
        return out;
    };
    return Mfrac(walk(src.numerator()), walk(src.denominator()));
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

void sort_integrals_like_amflow(std::vector<qft::JIntegral>& masters) {
    std::stable_sort(masters.begin(), masters.end());

    struct SectorGroup {
        std::vector<int>                  pattern;
        std::vector<qft::JIntegral>       members;
    };
    std::vector<SectorGroup> groups;
    for (auto& m : masters) {
        auto pat = m.sector_pattern();
        auto git = std::find_if(groups.begin(), groups.end(),
                                [&](const SectorGroup& g) {
                                    return g.pattern == pat;
                                });
        if (git == groups.end()) {
            groups.push_back({std::move(pat), {std::move(m)}});
        } else {
            git->members.push_back(std::move(m));
        }
    }
    for (auto& g : groups) std::reverse(g.members.begin(), g.members.end());
    std::reverse(groups.begin(), groups.end());

    masters.clear();
    for (auto& g : groups) {
        for (auto& m : g.members) masters.push_back(std::move(m));
    }
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

    fs::create_directories(opts.work_dir);
    kira_write_config(cfg, opts.work_dir);
    kira_write_jobs(cfg, opts.work_dir, KiraReductionMode::Reduce);
    kira_write_preferred(preferred, fc, opts.work_dir);
    kira_write_targets(targets, fc, opts.work_dir);

    kira_run(cfg, opts.work_dir, opts.log_file);

    out.masters = kira_read_masters(fc, opts.work_dir);
    auto raw_rules = kira_read_target_table(fc, opts.work_dir);

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

    {
        KiraConfig cfg = make_kira_config(fc, top_pattern, opts_eff);
        fs::create_directories(opts.work_dir);
        kira_write_config(cfg, opts.work_dir);
        kira_write_jobs(cfg, opts.work_dir, KiraReductionMode::Masters);
        kira_write_preferred(jpreferred, fc, opts.work_dir);
        kira_run(cfg, opts.work_dir, opts.log_file);
    }
    auto masters = kira_read_masters(fc, opts.work_dir);
    if (masters.empty()) {
        throw std::runtime_error("ibp::diffeq: no masters from Kira");
    }

    sort_integrals_like_amflow(masters);

    out.sortedmasters = masters;

    std::vector<std::vector<std::vector<DerivTerm>>> der(vars.size());
    std::set<std::string> all_int_keys;
    auto int_key = [](const qft::JIntegral& j) {
        std::ostringstream k; k << j.family();
        for (long v : j.indices()) k << "|" << v;
        return k.str();
    };
    for (std::size_t v = 0; v < vars.size(); ++v) {
        der[v].reserve(masters.size());
        for (const auto& m : masters) {
            auto raw = libp_deriv(fc, m, vars[v]);
            auto simp = simplify_terms(std::move(raw));
            for (const auto& t : simp) {
                all_int_keys.insert(int_key(t.integ));
            }
            der[v].push_back(std::move(simp));
        }
    }

    std::vector<qft::JIntegral> all_ints;
    {
        std::set<std::string> seen;
        for (auto& vec : der) {
            for (auto& terms : vec) {
                for (auto& t : terms) {
                    std::string k = int_key(t.integ);
                    if (seen.insert(k).second) all_ints.push_back(t.integ);
                }
            }
        }
        for (const auto& m : masters) {
            std::string k = int_key(m);
            if (seen.insert(k).second) all_ints.push_back(m);
        }
    }

    ReduceOptions reduce_opts = opts_eff;
    reduce_opts.work_dir = opts.work_dir + "/reduce";
    auto reduce_res = reduce(fc, all_ints, jpreferred,
                              top_pattern, reduce_opts);

    std::map<std::string, std::size_t> rule_index;
    for (std::size_t i = 0; i < all_ints.size(); ++i) {
        rule_index[int_key(all_ints[i])] = i;
    }

    if (reduce_res.masters.size() != masters.size()) {
        throw std::runtime_error(
            "ibp::diffeq: master count differs between Masters and Reduce calls");
    }

    std::map<std::string, std::size_t> master_index;
    for (std::size_t i = 0; i < masters.size(); ++i) {
        master_index[int_key(masters[i])] = i;
    }

    out.diffeq.reserve(vars.size());
    for (std::size_t v = 0; v < vars.size(); ++v) {
        out.diffeq.emplace_back();
        out.diffeq[v].reserve(masters.size());
        for (std::size_t row = 0; row < masters.size(); ++row) {
            out.diffeq[v].emplace_back();
            out.diffeq[v][row].reserve(masters.size());
            for (std::size_t col = 0; col < masters.size(); ++col) {
                out.diffeq[v][row].push_back(Mfrac::zero(out.red_ctx.ctx));
            }
            for (const auto& dt : der[v][row]) {
                std::string ik = int_key(dt.integ);
                auto rit = rule_index.find(ik);
                if (rit == rule_index.end()) {
                    continue;
                }
                const auto& rule = reduce_res.rules[rit->second];
                Mfrac coef_lifted = lift_to_red_ctx(dt.coef, out.red_ctx);
                for (const auto& mt : rule) {
                    auto mit = master_index.find(int_key(mt.integ));
                    if (mit == master_index.end()) continue;
                    Mfrac product = coef_lifted.clone();
                    Mfrac mt_coef_in_red = mfrac_to_ctx(mt.coef, out.red_ctx.ctx);
                    product *= mt_coef_in_red;
                    out.diffeq[v][row][mit->second] += product;
                }
            }
        }
    }

    return out;
}

}  // namespace amflow::ibp
