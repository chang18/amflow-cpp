// SPDX-License-Identifier: MIT
// ibp::libp_deriv — implementation.
//

#include "amflow/ibp/libp_deriv.hpp"

#include <map>
#include <set>
#include <stdexcept>
#include <utility>

#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_mpoly_q.h>

#include "amflow/algebra/numeric_subst.hpp"

namespace amflow::ibp {

using algebra::Mfrac;
using algebra::Mpoly;

namespace {

long find_var_index(const qft::FamilyConfig& fc, const std::string& name) {
    long idx = fc.ctx->var_index(name);
    if (idx < 0) {
        throw std::invalid_argument(
            "libp_deriv: invariant '" + name + "' not in fc.ctx");
    }
    return idx;
}

std::map<std::vector<unsigned long>, long>
build_sp_index(const qft::FamilyConfig& fc) {
    std::map<std::vector<unsigned long>, long> out;
    long n_loop = (long)fc.n_loops();
    long n_red  = (long)fc.n_red_legs();
    long n_loop_leg = n_loop + n_red;
    for (long k = 0; k < (long)fc.sp_list.size(); ++k) {
        std::vector<unsigned long> sig((std::size_t)n_loop_leg, 0);
        std::vector<unsigned long> exp((std::size_t)fc.ctx->n_vars(), 0);
        fmpz_mpoly_get_term_exp_ui(exp.data(), fc.sp_list[(std::size_t)k].raw(),
                                    0, fc.ctx->raw());
        for (long i = 0; i < n_loop_leg; ++i) {
            sig[(std::size_t)i] = exp[(std::size_t)i];
        }
        out[sig] = k;
    }
    return out;
}

struct SpDecomposition {
    std::vector<Mpoly> sp_coef;
    Mpoly              constant;
};

SpDecomposition decompose_into_splist(const qft::FamilyConfig& fc,
                                        const Mpoly& P) {
    long n_loop = (long)fc.n_loops();
    long n_red  = (long)fc.n_red_legs();
    long n_loop_leg = n_loop + n_red;
    long n_total = fc.ctx->n_vars();

    SpDecomposition out;
    out.constant = Mpoly::zero(fc.ctx);
    out.sp_coef.reserve(fc.sp_list.size());
    for (std::size_t k = 0; k < fc.sp_list.size(); ++k) {
        out.sp_coef.emplace_back(fc.ctx);
    }

    auto sp_index = build_sp_index(fc);

    long len = fmpz_mpoly_length(P.raw(), fc.ctx->raw());
    std::vector<unsigned long> exp((std::size_t)n_total);
    fmpz_t coeff;
    fmpz_init(coeff);
    for (long t = 0; t < len; ++t) {
        fmpz_mpoly_get_term_coeff_fmpz(coeff, P.raw(), t, fc.ctx->raw());
        fmpz_mpoly_get_term_exp_ui(exp.data(), P.raw(), t, fc.ctx->raw());

        std::vector<unsigned long> sig((std::size_t)n_loop_leg, 0);
        for (long i = 0; i < n_loop_leg; ++i) {
            sig[(std::size_t)i] = exp[(std::size_t)i];
        }
        std::vector<unsigned long> inv_exp = exp;
        for (long i = 0; i < n_loop_leg; ++i) inv_exp[(std::size_t)i] = 0;

        unsigned long deg_loopleg = 0;
        for (long i = 0; i < n_loop_leg; ++i) deg_loopleg += sig[(std::size_t)i];

        if (deg_loopleg == 0) {
            Mpoly mono = Mpoly::zero(fc.ctx);
            fmpz_mpoly_set_coeff_fmpz_ui(mono.raw(), coeff,
                                          inv_exp.data(), fc.ctx->raw());
            out.constant += mono;
        } else if (deg_loopleg == 2) {
            auto it = sp_index.find(sig);
            if (it == sp_index.end()) {
                fmpz_clear(coeff);
                throw std::runtime_error(
                    "decompose_into_splist: monomial not in sp_list");
            }
            Mpoly inv_mono = Mpoly::zero(fc.ctx);
            fmpz_mpoly_set_coeff_fmpz_ui(inv_mono.raw(), coeff,
                                          inv_exp.data(), fc.ctx->raw());
            out.sp_coef[(std::size_t)it->second] += inv_mono;
        } else {
            fmpz_clear(coeff);
            throw std::runtime_error(
                "decompose_into_splist: monomial has loop+leg degree " +
                std::to_string(deg_loopleg) + ", expected 0 or 2");
        }
    }
    fmpz_clear(coeff);
    return out;
}

struct SpToDMap {
    std::vector<Mpoly>                completede;
    std::vector<std::vector<Mfrac>>   alpha;
    std::vector<Mfrac>                beta;
};

SpToDMap build_sp_to_d_map(const qft::FamilyConfig& fc) {
    auto completede = qft::to_complete_explicit(fc, fc.propagators_after_conservation);
    auto dctx = qft::make_dlist_context(fc, completede.size());
    auto sp_to_d = qft::sp_list_to_dlist_symbol(fc, dctx, completede);

    SpToDMap out;
    out.completede = std::move(completede);
    long n_d = dctx.n_d;
    long first_d = dctx.first_d_var;
    long n_total_d = dctx.ctx->n_vars();

    out.alpha.reserve(sp_to_d.size());
    out.beta.reserve(sp_to_d.size());

    auto extract_in_fc = [&](const Mpoly& p) -> Mpoly {
        Mpoly out_p(fc.ctx);
        long len = fmpz_mpoly_length(p.raw(), dctx.ctx->raw());
        std::vector<unsigned long> dexp((std::size_t)n_total_d);
        std::vector<unsigned long> fexp((std::size_t)fc.ctx->n_vars(), 0);
        fmpz_t cf;
        fmpz_init(cf);
        for (long t = 0; t < len; ++t) {
            fmpz_mpoly_get_term_exp_ui(dexp.data(), p.raw(), t,
                                        dctx.ctx->raw());
            fmpz_mpoly_get_term_coeff_fmpz(cf, p.raw(), t,
                                            dctx.ctx->raw());
            for (long k = 0; k < n_d; ++k) {
                if (dexp[(std::size_t)(first_d + k)] > 0) {
                    fmpz_clear(cf);
                    throw std::runtime_error("extract_in_fc: D dependence found");
                }
            }
            for (long j = 0; j < fc.ctx->n_vars(); ++j) {
                fexp[(std::size_t)j] = dexp[(std::size_t)j];
            }
            fmpz_mpoly_set_coeff_fmpz_ui(out_p.raw(), cf,
                                          fexp.data(), fc.ctx->raw());
        }
        fmpz_clear(cf);
        return out_p;
    };

    for (const auto& m : sp_to_d) {
        Mpoly num = m.numerator();
        Mpoly den = m.denominator();

        std::vector<Mpoly> alpha_num;
        alpha_num.reserve((std::size_t)n_d);
        for (long k_d = 0; k_d < n_d; ++k_d) alpha_num.emplace_back(fc.ctx);
        Mpoly beta_num(fc.ctx);

        long len = fmpz_mpoly_length(num.raw(), dctx.ctx->raw());
        std::vector<unsigned long> dexp((std::size_t)n_total_d);
        std::vector<unsigned long> rest_exp((std::size_t)n_total_d);
        fmpz_t cf;
        fmpz_init(cf);
        for (long t = 0; t < len; ++t) {
            fmpz_mpoly_get_term_exp_ui(dexp.data(), num.raw(), t,
                                        dctx.ctx->raw());
            fmpz_mpoly_get_term_coeff_fmpz(cf, num.raw(), t,
                                            dctx.ctx->raw());

            long which_d = -1;
            unsigned long which_d_exp = 0;
            for (long k = 0; k < n_d; ++k) {
                if (dexp[(std::size_t)(first_d + k)] > 0) {
                    if (which_d != -1) {
                        fmpz_clear(cf);
                        throw std::runtime_error(
                            "build_sp_to_d_map: monomial involves multiple D vars");
                    }
                    which_d = k;
                    which_d_exp = dexp[(std::size_t)(first_d + k)];
                }
            }
            if (which_d != -1 && which_d_exp != 1) {
                fmpz_clear(cf);
                throw std::runtime_error(
                    "build_sp_to_d_map: D variable appears with non-linear exponent");
            }

            rest_exp = dexp;
            for (long k = 0; k < n_d; ++k) {
                rest_exp[(std::size_t)(first_d + k)] = 0;
            }

            std::vector<unsigned long> fexp((std::size_t)fc.ctx->n_vars(), 0);
            for (long j = 0; j < fc.ctx->n_vars(); ++j) {
                fexp[(std::size_t)j] = rest_exp[(std::size_t)j];
            }
            Mpoly rest_mono = Mpoly::zero(fc.ctx);
            fmpz_mpoly_set_coeff_fmpz_ui(rest_mono.raw(), cf,
                                          fexp.data(), fc.ctx->raw());

            if (which_d == -1) {
                beta_num += rest_mono;
            } else {
                alpha_num[(std::size_t)which_d] += rest_mono;
            }
        }
        fmpz_clear(cf);

        Mpoly den_fc = extract_in_fc(den);

        std::vector<Mfrac> alpha_row;
        alpha_row.reserve((std::size_t)n_d);
        for (long k = 0; k < n_d; ++k) {
            alpha_row.push_back(Mfrac(std::move(alpha_num[(std::size_t)k]),
                                        den_fc.clone()));
        }
        out.alpha.push_back(std::move(alpha_row));
        out.beta.push_back(Mfrac(std::move(beta_num), den_fc.clone()));
    }
    return out;
}

}  // namespace

LibpDenomsDerivResult
libp_denoms_deriv(const qft::FamilyConfig& fc, const std::string& s_name) {
    long s_var = find_var_index(fc, s_name);

    SpToDMap sp_map = build_sp_to_d_map(fc);

    LibpDenomsDerivResult out;
    out.completede = std::move(sp_map.completede);
    long n_d = (long)out.completede.size();

    std::size_t N = fc.propagators_after_conservation.size();
    out.coef.reserve(N);
    out.constant.reserve(N);
    for (std::size_t k = 0; k < N; ++k) {
        Mpoly d_k = fc.propagators_after_conservation[k]
                        .derivative(s_var);
        SpDecomposition dec = decompose_into_splist(fc, d_k);

        std::vector<Mfrac> coef_row;
        coef_row.reserve((std::size_t)n_d);
        for (long j = 0; j < n_d; ++j) {
            Mfrac acc = Mfrac::zero(fc.ctx);
            for (std::size_t l = 0; l < dec.sp_coef.size(); ++l) {
                if (dec.sp_coef[l].is_zero()) continue;
                Mfrac t = Mfrac::from_mpoly(dec.sp_coef[l].clone());
                if (t.ctx().get() != sp_map.alpha[l][(std::size_t)j].ctx().get()) {
                    throw std::runtime_error(
                        "libp_denoms_deriv: ctx mismatch t vs alpha at l="
                        + std::to_string(l) + " j=" + std::to_string(j));
                }
                t *= sp_map.alpha[l][(std::size_t)j];
                acc += t;
            }
            coef_row.push_back(std::move(acc));
        }
        Mfrac const_term = Mfrac::from_mpoly(dec.constant.clone());
        for (std::size_t l = 0; l < dec.sp_coef.size(); ++l) {
            if (dec.sp_coef[l].is_zero()) continue;
            Mfrac t = Mfrac::from_mpoly(dec.sp_coef[l].clone());
            t *= sp_map.beta[l];
            const_term += t;
        }
        out.coef.push_back(std::move(coef_row));
        out.constant.push_back(std::move(const_term));
    }
    return out;
}

std::vector<DerivTerm>
libp_deriv(const qft::FamilyConfig& fc,
           const qft::JIntegral& j,
           const std::string& s_name) {
    auto dd = libp_denoms_deriv(fc, s_name);
    long n_d = (long)dd.completede.size();
    long N   = (long)fc.propagators_after_conservation.size();

    if ((long)j.n_indices() != N) {
        throw std::invalid_argument(
            "libp_deriv: J integral indices length != #propagators");
    }

    std::vector<DerivTerm> out;
    for (long k = 0; k < N; ++k) {
        long ak = j.indices()[(std::size_t)k];
        if (ak == 0) continue;

        Mfrac neg_ak = Mfrac::from_si(fc.ctx, -ak);

        for (long jp = 0; jp < n_d; ++jp) {
            const Mfrac& cf = dd.coef[(std::size_t)k][(std::size_t)jp];
            if (cf.is_zero()) continue;
            std::vector<long> new_idx = j.indices();
            new_idx[(std::size_t)k] += 1;
            if (jp < N) {
                new_idx[(std::size_t)jp] -= 1;
            }
            qft::JIntegral shifted(j.family(), new_idx);
            DerivTerm dt;
            dt.integ = shifted;
            dt.coef = neg_ak.clone();
            dt.coef *= cf;
            out.push_back(std::move(dt));
        }

        const Mfrac& ck = dd.constant[(std::size_t)k];
        if (!ck.is_zero()) {
            std::vector<long> new_idx = j.indices();
            new_idx[(std::size_t)k] += 1;
            qft::JIntegral shifted(j.family(), new_idx);
            DerivTerm dt;
            dt.integ = shifted;
            dt.coef = neg_ak.clone();
            dt.coef *= ck;
            out.push_back(std::move(dt));
        }
    }
    return out;
}

std::vector<DerivTerm>
compute_derivative(const qft::FamilyConfig& fc,
                   const std::vector<DerivTerm>& expr,
                   const std::string& x_name) {
    long x_var = find_var_index(fc, x_name);

    std::vector<DerivTerm> out;
    for (const auto& term : expr) {
        Mfrac dc = term.coef.derivative(x_var);
        if (!dc.is_zero()) {
            DerivTerm dt;
            dt.coef = std::move(dc);
            dt.integ = term.integ;
            out.push_back(std::move(dt));
        }
        auto dj = libp_deriv(fc, term.integ, x_name);
        for (auto& d : dj) {
            DerivTerm dt;
            dt.coef = term.coef.clone();
            dt.coef *= d.coef;
            dt.integ = std::move(d.integ);
            out.push_back(std::move(dt));
        }
    }
    return out;
}

std::vector<DerivTerm>
simplify_terms(std::vector<DerivTerm> terms) {
    std::map<std::string, std::vector<std::size_t>> groups;
    for (std::size_t i = 0; i < terms.size(); ++i) {
        std::string key = terms[i].integ.family();
        for (long v : terms[i].integ.indices()) {
            key += "|"; key += std::to_string(v);
        }
        groups[key].push_back(i);
    }
    std::vector<DerivTerm> out;
    for (auto& [_, idxs] : groups) {
        DerivTerm acc;
        acc.integ = terms[idxs.front()].integ;
        acc.coef = Mfrac::zero(terms[idxs.front()].coef.ctx());
        for (std::size_t i : idxs) acc.coef += terms[i].coef;
        if (!acc.coef.is_zero()) out.push_back(std::move(acc));
    }
    return out;
}

// ===========================================================================
// D14 (2026-05-16) narrow-context overloads.
//
// These mirror `libp_denoms_deriv` / `libp_deriv` above but substitute
// the user-supplied `numeric_values` into intermediates and reproject
// everything to a caller-provided narrow polynomial-ring context
// (`target_ctx`).  This avoids the FLINT multivariate-GCD overhead
// that dominates the wide-fc.ctx path on multi-distinct-mass 3L+
// topologies (audit §D14 — 27 GB blowup on `bn3_4mass_3L_eps001`).
//
// The wide-ctx version is preserved for tests and any caller that
// genuinely needs symbolic dependence on family variables.
// ===========================================================================

LibpDenomsDerivResult
libp_denoms_deriv(const qft::FamilyConfig& fc, const std::string& s_name,
                  const std::map<std::string, std::string>& numeric_values,
                  const std::shared_ptr<algebra::MpolyContext>& target_ctx) {
    (void)find_var_index(fc, s_name);  // sanity check: s_name is in fc.ctx

    // Build sp_map on fc.ctx (this is essentially free — sp_map's
    // Mfracs are small rational coefficients of the inverse of a
    // SP-coefficient matrix).
    SpToDMap sp_map = build_sp_to_d_map(fc);

    LibpDenomsDerivResult out;
    out.completede = std::move(sp_map.completede);
    long n_d = (long)out.completede.size();

    // Pre-substitute + reproject sp_map.alpha[l][j] and sp_map.beta[l]
    // onto `target_ctx`.  Done ONCE up front so the per-k loop below
    // does arithmetic only on the narrow ctx.
    const std::size_t n_sp = sp_map.alpha.size();
    std::vector<std::vector<Mfrac>> alpha_narrow(n_sp);
    std::vector<Mfrac> beta_narrow;
    beta_narrow.reserve(n_sp);
    for (std::size_t l = 0; l < n_sp; ++l) {
        alpha_narrow[l].reserve((std::size_t)n_d);
        for (long j = 0; j < n_d; ++j) {
            alpha_narrow[l].push_back(
                algebra::substitute_and_narrow(
                    sp_map.alpha[l][(std::size_t)j],
                    numeric_values, target_ctx));
        }
        beta_narrow.push_back(
            algebra::substitute_and_narrow(
                sp_map.beta[l], numeric_values, target_ctx));
    }

    std::size_t N = fc.propagators_after_conservation.size();
    long s_var = find_var_index(fc, s_name);
    out.coef.reserve(N);
    out.constant.reserve(N);
    for (std::size_t k = 0; k < N; ++k) {
        Mpoly d_k = fc.propagators_after_conservation[k]
                        .derivative(s_var);
        SpDecomposition dec = decompose_into_splist(fc, d_k);

        // Pre-substitute + reproject dec.sp_coef[l] (Mpoly on fc.ctx)
        // and dec.constant (Mpoly on fc.ctx) to target_ctx.
        std::vector<Mfrac> sp_coef_narrow;
        sp_coef_narrow.reserve(dec.sp_coef.size());
        for (auto& mp : dec.sp_coef) {
            sp_coef_narrow.push_back(
                algebra::substitute_and_narrow(
                    Mfrac::from_mpoly(std::move(mp)),
                    numeric_values, target_ctx));
        }
        Mfrac const_narrow =
            algebra::substitute_and_narrow(
                Mfrac::from_mpoly(std::move(dec.constant)),
                numeric_values, target_ctx);

        std::vector<Mfrac> coef_row;
        coef_row.reserve((std::size_t)n_d);
        for (long j = 0; j < n_d; ++j) {
            Mfrac acc = Mfrac::zero(target_ctx);
            for (std::size_t l = 0; l < sp_coef_narrow.size(); ++l) {
                if (sp_coef_narrow[l].is_zero()) continue;
                Mfrac t = sp_coef_narrow[l].clone();
                t *= alpha_narrow[l][(std::size_t)j];
                acc += t;
            }
            coef_row.push_back(std::move(acc));
        }
        Mfrac const_term = const_narrow.clone();
        for (std::size_t l = 0; l < sp_coef_narrow.size(); ++l) {
            if (sp_coef_narrow[l].is_zero()) continue;
            Mfrac t = sp_coef_narrow[l].clone();
            t *= beta_narrow[l];
            const_term += t;
        }
        out.coef.push_back(std::move(coef_row));
        out.constant.push_back(std::move(const_term));
    }
    return out;
}

std::vector<DerivTerm>
libp_deriv(const qft::FamilyConfig& fc,
           const qft::JIntegral& j,
           const std::string& s_name,
           const std::map<std::string, std::string>& numeric_values,
           const std::shared_ptr<algebra::MpolyContext>& target_ctx) {
    auto dd = libp_denoms_deriv(fc, s_name, numeric_values, target_ctx);
    long n_d = (long)dd.completede.size();
    long N   = (long)fc.propagators_after_conservation.size();

    if ((long)j.n_indices() != N) {
        throw std::invalid_argument(
            "libp_deriv: J integral indices length != #propagators");
    }

    std::vector<DerivTerm> out;
    for (long k = 0; k < N; ++k) {
        long ak = j.indices()[(std::size_t)k];
        if (ak == 0) continue;

        // neg_ak is a constant rational on target_ctx (no symbolic
        // dependence — it's just `-ak`).
        Mfrac neg_ak = Mfrac::from_si(target_ctx, -ak);

        for (long jp = 0; jp < n_d; ++jp) {
            const Mfrac& cf = dd.coef[(std::size_t)k][(std::size_t)jp];
            if (cf.is_zero()) continue;
            std::vector<long> new_idx = j.indices();
            new_idx[(std::size_t)k] += 1;
            if (jp < N) {
                new_idx[(std::size_t)jp] -= 1;
            }
            qft::JIntegral shifted(j.family(), new_idx);
            DerivTerm dt;
            dt.integ = shifted;
            dt.coef = neg_ak.clone();
            dt.coef *= cf;
            out.push_back(std::move(dt));
        }

        const Mfrac& ck = dd.constant[(std::size_t)k];
        if (!ck.is_zero()) {
            std::vector<long> new_idx = j.indices();
            new_idx[(std::size_t)k] += 1;
            qft::JIntegral shifted(j.family(), new_idx);
            DerivTerm dt;
            dt.integ = shifted;
            dt.coef = neg_ak.clone();
            dt.coef *= ck;
            out.push_back(std::move(dt));
        }
    }
    return out;
}

}  // namespace amflow::ibp
