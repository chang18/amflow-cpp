// SPDX-License-Identifier: MIT
// qft::findregion — FindAllRegion + ZeroRegionQ + RegionPower.
//

#include "amflow/qft/region.hpp"

#include <algorithm>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>

#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_mpoly_q.h>
#include <flint/flint.h>

#include "amflow/numeric/log.hpp"
#include "amflow/qft/dense_q.hpp"
#include "amflow/qft/jintegral.hpp"
#include "amflow/qft/topology.hpp"

namespace amflow::qft {

using algebra::Mfrac;
using algebra::Mpoly;
using algebra::MpolyContext;

PowersContext make_powers_context(const FamilyConfig& fc) {
    std::vector<std::string> names;
    names.reserve((std::size_t)fc.ctx->n_vars() + 1);
    for (long i = 0; i < fc.ctx->n_vars(); ++i) {
        names.push_back(fc.ctx->var_name(i));
    }
    names.emplace_back("__amf_eps");
    auto p_ctx = std::make_shared<MpolyContext>(std::move(names));

    PowersContext pctx;
    pctx.ctx     = p_ctx;
    pctx.eps_var = fc.ctx->n_vars();
    pctx.lift_gens.reserve((std::size_t)fc.ctx->n_vars());
    for (long i = 0; i < fc.ctx->n_vars(); ++i) {
        pctx.lift_gens.push_back(Mpoly::variable(p_ctx, i));
    }
    return pctx;
}

namespace {

long max_degree_in_var(const Mpoly& p, long var) {
    auto ctx = p.ctx();
    long len = fmpz_mpoly_length(p.raw(), ctx->raw());
    std::vector<unsigned long> exp((std::size_t)ctx->n_vars());
    long d = 0;
    for (long t = 0; t < len; ++t) {
        fmpz_mpoly_get_term_exp_ui(exp.data(), p.raw(), t, ctx->raw());
        if ((long)exp[(std::size_t)var] > d) d = (long)exp[(std::size_t)var];
    }
    return d;
}

bool depends_on_var(const Mpoly& p, long var) {
    return max_degree_in_var(p, var) > 0;
}

Mfrac merge_eta_into_half_eta(const FamilyConfig& fc,
                              const RegionContext& rctx,
                              const Mfrac& raw) {
    long fc_eta_var = fc.ctx->var_index("eta");
    if (fc_eta_var < 0) return raw.clone();

    auto sub_eta = [&](const Mpoly& p) -> Mpoly {
        Mpoly out(rctx.ctx);
        long len = fmpz_mpoly_length(p.raw(), rctx.ctx->raw());
        std::vector<unsigned long> exp((std::size_t)rctx.ctx->n_vars());
        fmpz_t coeff;
        fmpz_init(coeff);
        for (long t = 0; t < len; ++t) {
            fmpz_mpoly_get_term_exp_ui(exp.data(), p.raw(), t, rctx.ctx->raw());
            fmpz_mpoly_get_term_coeff_fmpz(coeff, p.raw(), t, rctx.ctx->raw());
            long ek = (long)exp[(std::size_t)fc_eta_var];
            std::vector<unsigned long> e2 = exp;
            e2[(std::size_t)fc_eta_var] = 0;
            e2[(std::size_t)rctx.half_eta_var] += (unsigned long)(2 * ek);
            fmpz_mpoly_set_coeff_fmpz_ui(out.raw(), coeff, e2.data(), rctx.ctx->raw());
        }
        fmpz_clear(coeff);
        return out;
    };

    return Mfrac(sub_eta(raw.numerator()), sub_eta(raw.denominator()));
}

}  // namespace

std::vector<Region>
find_all_region(const FamilyConfig& fc,
                const RegionContext& rctx,
                const std::vector<std::size_t>& topposi) {
    long L = (long)fc.n_loops();
    std::size_t N = fc.propagators_after_conservation.size();

    auto branch0 = branch_momenta(fc, fc.propagators_after_conservation);

    std::vector<Mpoly> dedup;
    std::vector<std::string> dedup_keys;
    auto add_unique = [&](const Mpoly& b) {
        std::string k = b.to_string();
        for (const auto& dk : dedup_keys) {
            if (dk == k) return;
        }
        dedup_keys.push_back(k);
        dedup.push_back(b.clone());
    };
    for (std::size_t k = 0; k < topposi.size(); ++k) {
        add_unique(branch0[topposi[k]]);
    }

    std::set<std::size_t> cutposi;
    for (std::size_t i = 0; i < N; ++i) {
        if (i < fc.cut.size() && fc.cut[i] == 1) {
            std::string k = branch0[i].to_string();
            for (std::size_t d = 0; d < dedup_keys.size(); ++d) {
                if (dedup_keys[d] == k) cutposi.insert(d);
            }
        }
    }

    if (L <= 0 || dedup.empty()) {
        return {};
    }
    std::vector<std::pair<LoopTransform, std::vector<std::size_t>>> trans_list;
    {
        std::vector<std::size_t> idx_tuple((std::size_t)L, 0);
        std::function<void(long)> rec = [&](long depth) {
            if (depth == L) {
                std::vector<Mpoly> cand;
                cand.reserve((std::size_t)L);
                for (long i = 0; i < L; ++i) {
                    cand.push_back(dedup[idx_tuple[(std::size_t)i]].clone());
                }
                LoopTransform tr = branch_to_loop(fc, rctx, cand);
                if (tr.ok) {
                    trans_list.emplace_back(std::move(tr), idx_tuple);
                }
                return;
            }
            for (std::size_t i = 0; i < dedup.size(); ++i) {
                idx_tuple[(std::size_t)depth] = i;
                rec(depth + 1);
            }
        };
        rec(0);
    }

    struct RegionEntry {
        LoopTransform     tr;
        std::vector<int>  scale;
        std::vector<int>  bscale;
    };

    std::vector<RegionEntry> all_entries;

    auto& dedup_ref = dedup;
    for (auto& [tr, _] : trans_list) {
        std::vector<int> scale((std::size_t)L, 0);
        std::function<void(long)> rec_scale = [&](long depth) {
            if (depth == L) {
            std::vector<int> identity_scale((std::size_t)L, 0);
            std::vector<Mfrac> identity_rule = region_rule(rctx, tr,
                                                            identity_scale);
            std::vector<Mfrac> transformed;
            transformed.reserve(dedup_ref.size());
            for (auto& d : dedup_ref) {
                transformed.push_back(
                    apply_region_rule(fc, rctx, d, identity_rule));
            }
            std::vector<int> bs = branch_scale(fc, transformed, scale);

            bool ok = true;
            for (std::size_t cp : cutposi) {
                if (cp >= bs.size()) continue;
                if (bs[cp] != 0) { ok = false; break; }
            }
            if (!ok) return;

            RegionEntry e;
            e.tr     = LoopTransform();
            e.tr.ok  = true;
            for (auto& mfm : tr.map) e.tr.map.push_back(mfm.clone());
            e.scale  = scale;
            e.bscale = std::move(bs);
            all_entries.push_back(std::move(e));
                return;
            }
            for (int bit : {0, 1}) {
                scale[(std::size_t)depth] = bit;
                rec_scale(depth + 1);
            }
        };
        rec_scale(0);
    }

    std::vector<std::pair<std::vector<int>, std::vector<RegionEntry>>> by_bscale;
    for (auto& e : all_entries) {
        auto it = std::find_if(
            by_bscale.begin(), by_bscale.end(),
            [&](const auto& item) { return item.first == e.bscale; });
        if (it == by_bscale.end()) {
            by_bscale.emplace_back();
            by_bscale.back().first = e.bscale;
            it = std::prev(by_bscale.end());
        }
        it->second.push_back(std::move(e));
    }

    std::vector<int> pres0;
    pres0.reserve(dedup.size());
    for (auto& d : dedup) {
        auto p = fc.prescription_of_prop(d);
        if (!p.has_value()) {
            throw std::runtime_error(
                "find_all_region: deduped branch has conflicting prescription");
        }
        pres0.push_back(p.value());
    }

    auto check_pres_match = [&](const RegionEntry& e) -> bool {
        std::vector<int> identity_scale((std::size_t)L, 0);
        std::vector<Mfrac> identity_rule = region_rule(rctx, e.tr,
                                                        identity_scale);
        for (std::size_t i = 0; i < dedup.size(); ++i) {
            Mfrac t = apply_region_rule(fc, rctx, dedup[i], identity_rule);
            std::vector<int> pres_list;
            long n_loop = (long)fc.n_loops();
            for (long k = 0; k < n_loop; ++k) {
                if (depends_on_var(t.numerator(), k)
                        || depends_on_var(t.denominator(), k)) {
                    pres_list.push_back(fc.prescription_of_loop(k));
                }
            }
            std::optional<int> p_new;
            if (pres_list.empty()) {
                p_new = 0;
            } else {
                bool all_zero = true;
                for (int pp : pres_list) if (pp != 0) { all_zero = false; break; }
                if (all_zero) p_new = 0;
                else {
                    bool any_pos = false, any_neg = false;
                    for (int pp : pres_list) {
                        if (pp > 0) any_pos = true;
                        else if (pp < 0) any_neg = true;
                    }
                    if (any_pos && !any_neg) p_new = 1;
                    else if (!any_pos && any_neg) p_new = -1;
                }
            }
            if (!p_new.has_value()) return false;
            if (p_new.value() != pres0[i]) return false;
        }
        return true;
    };

    std::vector<Region> result;
    for (auto& [bs, entries] : by_bscale) {
        std::vector<RegionEntry> kept;
        for (auto& e : entries) {
            if (check_pres_match(e)) {
                kept.push_back(std::move(e));
            }
        }
        if (kept.empty()) {
            throw std::runtime_error(
                "find_all_region: some regions are prohibited by prescriptions");
        }
        Region r;
        r.scale = kept.front().scale;
        r.transform.ok = true;
        for (auto& mfm : kept.front().tr.map) r.transform.map.push_back(mfm.clone());
        result.push_back(std::move(r));
    }
    return result;
}

namespace {

Mpoly coeff_half_eta_mpoly(const Mfrac& f, long half_eta_var, long net_he) {
    auto ctx = f.ctx();

    long den_he = 0;
    {
        const Mpoly& den = f.denominator();
        long dlen = fmpz_mpoly_length(den.raw(), ctx->raw());
        std::vector<unsigned long> dexp((std::size_t)ctx->n_vars());
        bool first = true;
        for (long t = 0; t < dlen; ++t) {
            fmpz_mpoly_get_term_exp_ui(dexp.data(), den.raw(), t, ctx->raw());
            long d = (long)dexp[(std::size_t)half_eta_var];
            if (first) {
                den_he = d;
                first = false;
            } else if (d != den_he) {
                throw std::runtime_error(
                    "coeff_half_eta_mpoly: denominator is not a monomial in half_eta");
            }
        }
    }

    long target_num_he = net_he + den_he;
    Mpoly out = Mpoly::zero(ctx);
    const Mpoly& num = f.numerator();
    long nlen = fmpz_mpoly_length(num.raw(), ctx->raw());
    std::vector<unsigned long> nexp((std::size_t)ctx->n_vars());
    fmpz_t coeff;
    fmpz_init(coeff);
    for (long t = 0; t < nlen; ++t) {
        fmpz_mpoly_get_term_exp_ui(nexp.data(), num.raw(), t, ctx->raw());
        if ((long)nexp[(std::size_t)half_eta_var] != target_num_he) continue;
        fmpz_mpoly_get_term_coeff_fmpz(coeff, num.raw(), t, ctx->raw());
        std::vector<unsigned long> e2 = nexp;
        e2[(std::size_t)half_eta_var] = 0;
        fmpz_mpoly_set_coeff_fmpz_ui(out.raw(), coeff, e2.data(), ctx->raw());
    }
    fmpz_clear(coeff);
    return out;
}

Mpoly project_to_fc(const FamilyConfig& fc,
                     const RegionContext& rctx,
                     const Mpoly& p) {
    if (depends_on_var(p, rctx.eta_var) ||
        depends_on_var(p, rctx.half_eta_var)) {
        throw std::runtime_error(
            "project_to_fc: polynomial unexpectedly contains eta/half_eta");
    }
    Mpoly p1 = p.clone();
    fmpz_t zero;
    fmpz_init(zero);
    Mpoly tmp(rctx.ctx);
    fmpz_mpoly_evaluate_one_fmpz(tmp.raw(), p1.raw(),
                                  rctx.half_eta_var, zero, rctx.ctx->raw());
    p1 = std::move(tmp);
    Mpoly tmp2(rctx.ctx);
    fmpz_mpoly_evaluate_one_fmpz(tmp2.raw(), p1.raw(),
                                  rctx.eta_var, zero, rctx.ctx->raw());
    p1 = std::move(tmp2);
    fmpz_clear(zero);

    Mpoly out(fc.ctx);
    long len = fmpz_mpoly_length(p1.raw(), rctx.ctx->raw());
    std::vector<unsigned long> rexp((std::size_t)rctx.ctx->n_vars());
    std::vector<unsigned long> fexp((std::size_t)fc.ctx->n_vars(), 0);
    fmpz_t coeff;
    fmpz_init(coeff);
    for (long t = 0; t < len; ++t) {
        fmpz_mpoly_get_term_exp_ui(rexp.data(), p1.raw(),
                                    t, rctx.ctx->raw());
        fmpz_mpoly_get_term_coeff_fmpz(coeff, p1.raw(),
                                        t, rctx.ctx->raw());
        for (long i = 0; i < fc.ctx->n_vars(); ++i) {
            fexp[(std::size_t)i] = rexp[(std::size_t)i];
        }
        fmpz_mpoly_set_coeff_fmpz_ui(out.raw(), coeff,
                                      fexp.data(), fc.ctx->raw());
    }
    fmpz_clear(coeff);
    return out;
}

bool vanishes_at_loop_and_leg_zero(const FamilyConfig& fc, const Mpoly& p) {
    Mpoly q = p.clone();
    fmpz_t zero;
    fmpz_init(zero);
    long n_zero_vars = (long)fc.n_loops() + (long)fc.n_red_legs();
    for (long v = 0; v < n_zero_vars; ++v) {
        Mpoly tmp(fc.ctx);
        fmpz_mpoly_evaluate_one_fmpz(tmp.raw(), q.raw(), v, zero, fc.ctx->raw());
        q = std::move(tmp);
    }
    fmpz_clear(zero);
    return q.is_zero();
}

}  // namespace

bool zero_region_q(const FamilyConfig& fc,
                   const RegionContext& rctx,
                   const Region& region,
                   const std::vector<std::size_t>& topposi) {
    auto rule = region_rule(rctx, region.transform, region.scale);

    std::vector<Mpoly> expde;
    expde.reserve(topposi.size());
    for (std::size_t k = 0; k < topposi.size(); ++k) {
        std::size_t i = topposi[k];
        const Mpoly& prop = fc.propagators_after_conservation[i];
        Mfrac t = apply_region_rule(fc, rctx, prop, rule);
        t = merge_eta_into_half_eta(fc, rctx, t);
        bool has_he = depends_on_var(t.numerator(), rctx.half_eta_var)
                   || depends_on_var(t.denominator(), rctx.half_eta_var);
        Mfrac stripped = t.clone();
        if (has_he) {
            Mpoly he = Mpoly::variable(rctx.ctx, rctx.half_eta_var);
            Mpoly he2 = he * he;
            stripped = Mfrac(stripped.numerator(), stripped.denominator() * he2);
        }
        Mpoly lead = coeff_half_eta_mpoly(stripped, rctx.half_eta_var, 0);
        Mpoly proj = project_to_fc(fc, rctx, lead);
        expde.push_back(std::move(proj));
    }

    std::stable_sort(expde.begin(), expde.end(),
                     [&](const Mpoly& a, const Mpoly& b) {
                         bool ak = vanishes_at_loop_and_leg_zero(fc, a);
                         bool bk = vanishes_at_loop_and_leg_zero(fc, b);
                         return ak < bk;
                     });

    long M = (long)fc.sp_list.size();
    long rows = (long)expde.size();
    if (rows == 0) {
        return false;
    }
    std::vector<fmpq*> mat_flat((std::size_t)(rows * M), nullptr);
    for (long i = 0; i < rows * M; ++i) {
        mat_flat[(std::size_t)i] = (fmpq*)flint_malloc(sizeof(fmpq));
        fmpq_init(mat_flat[(std::size_t)i]);
    }
    {
        std::vector<fmpq*> row((std::size_t)M, nullptr);
        for (long i = 0; i < M; ++i) {
            row[(std::size_t)i] = (fmpq*)flint_malloc(sizeof(fmpq));
            fmpq_init(row[(std::size_t)i]);
        }
        for (long i = 0; i < rows; ++i) {
            coeff_over_splist(expde[(std::size_t)i], fc.sp_list, row);
            for (long k = 0; k < M; ++k) {
                fmpq_set(mat_flat[(std::size_t)(i * M + k)],
                          row[(std::size_t)k]);
            }
        }
        for (auto* p : row) { fmpq_clear(p); flint_free(p); }
    }

    auto group = maximal_group_rows((std::size_t)rows, (std::size_t)M, mat_flat);
    for (auto* p : mat_flat) { fmpq_clear(p); flint_free(p); }
    std::vector<Mpoly> pruned;
    pruned.reserve(group.size());
    for (std::size_t k = 0; k < group.size(); ++k) {
        pruned.push_back(expde[group[k]].clone());
    }

    if (pruned.empty()) {
        return true;
    }

    return zero_sector_q(fc, pruned);
}

std::vector<Mfrac>
region_power(const FamilyConfig& fc,
             const RegionContext& rctx,
             const PowersContext& pctx,
             const Region& region,
             const std::vector<JIntegral>& integrals) {
    std::size_t N = fc.propagators_after_conservation.size();

    auto rule = region_rule(rctx, region.transform, region.scale);
    std::vector<bool> factor_eta(N, false);
    for (std::size_t i = 0; i < N; ++i) {
        Mfrac t = apply_region_rule(fc, rctx,
            fc.propagators_after_conservation[i], rule);
        t = merge_eta_into_half_eta(fc, rctx, t);
        Mpoly num = t.numerator();
        auto ctx = num.ctx();
        long len = fmpz_mpoly_length(num.raw(), ctx->raw());
        if (len == 0) { factor_eta[i] = false; continue; }
        std::vector<unsigned long> exp((std::size_t)ctx->n_vars());
        bool has_he = false;
        for (long t2 = 0; t2 < len; ++t2) {
            fmpz_mpoly_get_term_exp_ui(exp.data(), num.raw(),
                                        t2, ctx->raw());
            if (exp[(std::size_t)rctx.half_eta_var] > 0) {
                has_he = true;
                break;
            }
        }
        factor_eta[i] = has_he;
    }

    long sum_scale = 0;
    for (int s : region.scale) sum_scale += s;

    long eps = pctx.eps_var;
    Mpoly two_minus_eps = Mpoly::constant(pctx.ctx, 2);
    {
        Mpoly e = Mpoly::variable(pctx.ctx, eps);
        two_minus_eps -= e;
    }

    std::vector<Mfrac> out;
    out.reserve(integrals.size());
    for (const auto& integ : integrals) {
        if ((long)integ.indices().size() != (long)N) {
            throw std::invalid_argument(
                "region_power: integral indices length != n_propagators");
        }
        long sum_p = 0;
        for (std::size_t i = 0; i < N; ++i) {
            if (factor_eta[i]) sum_p += integ.indices()[i];
        }
        Mpoly val(pctx.ctx);
        {
            fmpz_t s;
            fmpz_init_set_si(s, sum_scale);
            fmpz_mpoly_scalar_mul_fmpz(val.raw(),
                                        two_minus_eps.raw(),
                                        s, pctx.ctx->raw());
            fmpz_clear(s);
        }
        Mpoly minus_sum = Mpoly::constant(pctx.ctx, -sum_p);
        val += minus_sum;
        out.push_back(Mfrac::from_mpoly(std::move(val)));
    }
    return out;
}

}  // namespace amflow::qft
