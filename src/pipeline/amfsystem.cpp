// SPDX-License-Identifier: MIT
// Layer 16: AMFSystem pipeline -- implementation.
//
// See amflow/pipeline/amfsystem.hpp for the conceptual model.

#include "amflow/pipeline/amfsystem.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <filesystem>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_mpoly_q.h>

#include "amflow/qft/amfmode.hpp"
#include "amflow/ibp/reduce.hpp"
#include "amflow/qft/boundary.hpp"
#include "amflow/pipeline/factorize.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/family_uf.hpp"
#include "amflow/ode/inf.hpp"
#include "amflow/numeric/log.hpp"
#include "amflow/numeric/options.hpp"
#include "amflow/numeric/rational.hpp"
#include "amflow/qft/region.hpp"
#include "amflow/ode/amflow.hpp"
#include "amflow/qft/topology.hpp"
#include "amflow/qft/vacuum.hpp"

namespace amflow::pipeline {

namespace fs = std::filesystem;

namespace {
// Kira keeps per-directory auxiliary state (numeric substitution count).
// Boundary reduction reuses similar path stems; append a unique tag so
// each BlackBoxReduce run gets a fresh working directory.
std::atomic<long long> g_boundary_reduce_unique{0};
}  // namespace

// ===========================================================================
//  EndingScheme name
// ===========================================================================

const char* ending_scheme_name(EndingScheme s) noexcept {
    switch (s) {
        case EndingScheme::Tradition:  return "Tradition";
        case EndingScheme::Cutkosky:   return "Cutkosky";
        case EndingScheme::SingleMass: return "SingleMass";
    }
    return "?";
}

// ===========================================================================
//  Helpers
// ===========================================================================

namespace {

std::string mpoly_str(const algebra::Mpoly& p) { return p.to_string(); }

std::string mfrac_str(const algebra::Mfrac& f) {
    algebra::Mpoly num = f.numerator();
    algebra::Mpoly den = f.denominator();
    if (den.is_one()) return num.to_string();
    return "(" + num.to_string() + ")/(" + den.to_string() + ")";
}

std::vector<std::pair<std::string, std::string>>
extract_conservation_strings(const qft::FamilyConfig& fc) {
    std::vector<std::pair<std::string, std::string>> out;
    for (const auto& [key, mpoly] : fc.conservation) {
        out.emplace_back(key, mpoly.to_string());
    }
    return out;
}

std::vector<std::pair<std::string, std::string>>
extract_replacement_strings(const qft::FamilyConfig& fc) {
    std::vector<std::pair<std::string, std::string>> out;
    for (const auto& [key, val] : fc.reduced_replacement) {
        out.emplace_back(key, mfrac_str(val));
    }
    return out;
}

numeric::AcbValue midpoint_only(const numeric::AcbValue& src) {
    numeric::AcbValue out;
    arb_set_arf(acb_realref(out.raw()), arb_midref(acb_realref(src.raw())));
    arb_set_arf(acb_imagref(out.raw()), arb_midref(acb_imagref(src.raw())));
    return out;
}

std::vector<std::string>
extract_propagator_strings(const qft::FamilyConfig& fc) {
    std::vector<std::string> out;
    for (const auto& p : fc.propagators_after_conservation) {
        out.push_back(p.to_string());
    }
    return out;
}

std::string jintegral_key(const qft::JIntegral& j) {
    std::ostringstream oss;
    oss << j.family();
    for (long v : j.indices()) oss << "|" << v;
    return oss.str();
}

// Build the eta-injected family.  Mathematica:
//   AMFlowInfo["Propagator"] = ReducedPropagator + $Eta * etac
// We add an `eta` invariant to fc.ctx (assuming it isn't already there)
// and rewrite each propagator string.  If etac is all-zero, returns a
// fresh copy of fc unchanged.
std::unique_ptr<qft::FamilyConfig>
build_injected_fc(const qft::FamilyConfig& parent, const std::vector<int>& etac) {
    if (etac.size() != parent.propagators_after_conservation.size()) {
        throw std::invalid_argument(
            "build_injected_fc: etac size mismatch");
    }
    bool any = false;
    for (int v : etac) if (v != 0) { any = true; break; }
    if (!any) {
        return std::make_unique<qft::FamilyConfig>(qft::FamilyConfig::build(
            parent.family, parent.loops, parent.legs,
            extract_conservation_strings(parent),
            extract_replacement_strings(parent),
            extract_propagator_strings(parent),
            parent.cut, parent.prescription));
    }
    std::vector<std::string> prop_strs;
    prop_strs.reserve(parent.propagators_after_conservation.size());
    for (std::size_t i = 0; i < parent.propagators_after_conservation.size();
            ++i) {
        std::ostringstream oss;
        oss << "(" << mpoly_str(parent.propagators_after_conservation[i]) << ")";
        if (etac[i] != 0) {
            if (etac[i] > 0) oss << " + ";
            else             oss << " - ";
            long ae = std::abs(etac[i]);
            if (ae != 1) oss << ae << "*";
            oss << "eta";
        }
        prop_strs.push_back(oss.str());
    }
    return std::make_unique<qft::FamilyConfig>(qft::FamilyConfig::build(
        parent.family,
        parent.loops, parent.legs,
        extract_conservation_strings(parent),
        extract_replacement_strings(parent),
        prop_strs,
        parent.cut, parent.prescription));
}

// Project an algebra::Mfrac on a "wide" ctx down to dst_ctx by walking
// monomials.  Variables in src that aren't in dst:
//   * "__amf_*" / "__feyn_*" / "__zsq_*" -> drop (set to 0)
//   * everything else -> error.
algebra::Mfrac project_mfrac_by_name(const algebra::Mfrac& src,
                              const std::shared_ptr<algebra::MpolyContext>& dst_ctx) {
    auto src_ctx = src.ctx();
    if (src_ctx.get() == dst_ctx.get()) return src.clone();

    std::map<std::string, long> name_to_dst;
    for (long i = 0; i < dst_ctx->n_vars(); ++i) {
        name_to_dst[dst_ctx->var_name(i)] = i;
    }
    std::vector<long> src_to_dst((std::size_t)src_ctx->n_vars(), -2);
    for (long i = 0; i < src_ctx->n_vars(); ++i) {
        const std::string& name = src_ctx->var_name(i);
        auto it = name_to_dst.find(name);
        if (it != name_to_dst.end()) {
            src_to_dst[(std::size_t)i] = it->second;
        } else if (name.rfind("__amf_", 0) == 0
                || name.rfind("__feyn_", 0) == 0
                || name.rfind("__zsq_", 0) == 0) {
            src_to_dst[(std::size_t)i] = -1;   // drop
        }
    }
    auto walk = [&](const algebra::Mpoly& p) -> algebra::Mpoly {
        algebra::Mpoly out(dst_ctx);
        long len = fmpz_mpoly_length(p.raw(), src_ctx->raw());
        std::vector<unsigned long> exp((std::size_t)src_ctx->n_vars());
        std::vector<unsigned long> dexp((std::size_t)dst_ctx->n_vars(), 0);
        fmpz_t coeff;
        fmpz_init(coeff);
        for (long t = 0; t < len; ++t) {
            fmpz_mpoly_get_term_exp_ui(exp.data(), p.raw(), t, src_ctx->raw());
            fmpz_mpoly_get_term_coeff_fmpz(coeff, p.raw(), t, src_ctx->raw());
            std::fill(dexp.begin(), dexp.end(), 0);
            bool drop = false;
            for (long i = 0; i < src_ctx->n_vars(); ++i) {
                if (exp[(std::size_t)i] == 0) continue;
                long pi = src_to_dst[(std::size_t)i];
                if (pi == -1) { drop = true; break; }
                if (pi == -2) {
                    fmpz_clear(coeff);
                    throw std::runtime_error(
                        "project_mfrac_by_name: residual var '"
                        + src_ctx->var_name(i) + "' not in dst");
                }
                dexp[(std::size_t)pi] = exp[(std::size_t)i];
            }
            if (drop) continue;
            fmpz_mpoly_set_coeff_fmpz_ui(out.raw(), coeff,
                                          dexp.data(), dst_ctx->raw());
        }
        fmpz_clear(coeff);
        return out;
    };
    algebra::Mpoly num = walk(src.numerator());
    algebra::Mpoly den = walk(src.denominator());
    if (den.is_zero()) {
        throw std::runtime_error("project_mfrac_by_name: denominator vanished");
    }
    return algebra::Mfrac(std::move(num), std::move(den));
}

// FmpqHolder: small RAII wrapper.
class FmpqHolder {
public:
    FmpqHolder() { fmpq_init(q_); }
    FmpqHolder(const FmpqHolder& o) { fmpq_init(q_); fmpq_set(q_, o.q_); }
    FmpqHolder(FmpqHolder&& o) noexcept { fmpq_init(q_); fmpq_swap(q_, o.q_); }
    FmpqHolder& operator=(const FmpqHolder& o) {
        if (this != &o) fmpq_set(q_, o.q_);
        return *this;
    }
    FmpqHolder& operator=(FmpqHolder&& o) noexcept {
        if (this != &o) fmpq_swap(q_, o.q_);
        return *this;
    }
    ~FmpqHolder() { fmpq_clear(q_); }
    fmpq* raw()             { return q_; }
    const fmpq* raw() const { return q_; }
private:
    fmpq_t q_;
};

// Test if a variable v has positive exponent in any term of `p`.
bool mpoly_has_var(const algebra::Mpoly& p, long v) {
    auto ctx = p.ctx();
    long len = fmpz_mpoly_length(p.raw(), ctx->raw());
    std::vector<unsigned long> exp((std::size_t)ctx->n_vars());
    for (long t = 0; t < len; ++t) {
        fmpz_mpoly_get_term_exp_ui(exp.data(), p.raw(), t, ctx->raw());
        if (exp[(std::size_t)v] > 0) return true;
    }
    return false;
}

// Substitute every variable in `src` *except* those flagged `keep_names`
// to the rational value supplied via numeric_q (by name).  For "d",
// derive d = 4 - 2*eps when not given directly.  Internal "__amf_*"
// vars not in numeric_q default to 0.  Variables that don't actually
// appear in `src` (positive exponent in num or den) are skipped.
algebra::Mfrac mfrac_substitute(const algebra::Mfrac& src,
                          const std::map<std::string, FmpqHolder>& numeric_q,
                          const std::set<std::string>& keep_names) {
    auto ctx = src.ctx();
    algebra::Mfrac cur = src.clone();
    for (long v = 0; v < ctx->n_vars(); ++v) {
        const std::string& name = ctx->var_name(v);
        if (keep_names.count(name)) continue;
        // Skip variables that don't appear at all (so we don't demand
        // a numeric value for them).  Note: `cur` has already been
        // partially substituted for earlier `v'`s; recompute.
        if (!mpoly_has_var(cur.numerator(), v)
                && !mpoly_has_var(cur.denominator(), v)) {
            continue;
        }
        auto it = numeric_q.find(name);
        if (it == numeric_q.end()) {
            if (name == "d") {
                auto eit = numeric_q.find("eps");
                if (eit == numeric_q.end()) eit = numeric_q.find("__amf_eps");
                if (eit == numeric_q.end()) {
                    throw std::runtime_error(
                        "mfrac_substitute: missing 'd' and 'eps'");
                }
                fmpq_t d_val, two_eps;
                fmpq_init(d_val); fmpq_init(two_eps);
                fmpq_set_si(d_val, 4, 1);
                fmpq_set(two_eps, eit->second.raw());
                fmpq_mul_si(two_eps, two_eps, 2);
                fmpq_sub(d_val, d_val, two_eps);
                fmpq_clear(two_eps);
                cur = cur.substitute(v, d_val);
                fmpq_clear(d_val);
                continue;
            }
            if (name.rfind("__amf_", 0) == 0
                    || name.rfind("__feyn_", 0) == 0
                    || name.rfind("__zsq_", 0) == 0) {
                fmpq_t zero; fmpq_init(zero);
                cur = cur.substitute(v, zero);
                fmpq_clear(zero);
                continue;
            }
            throw std::runtime_error(
                "mfrac_substitute: no value for '" + name + "'");
        }
        cur = cur.substitute(v, it->second.raw());
    }
    return cur;
}

// Convert an algebra::Mfrac in a single variable (after all other substitutions)
// to a numeric::RationalFunction over Q.  `eta_var` < 0 means the result is a
// constant.
numeric::RationalFunction mfrac_to_eta_rational(const algebra::Mfrac& src, long eta_var) {
    auto ctx = src.ctx();
    long n_vars = ctx->n_vars();
    auto walk = [&](const algebra::Mpoly& p) -> numeric::FmpqPoly {
        numeric::FmpqPoly out;
        long len = fmpz_mpoly_length(p.raw(), ctx->raw());
        std::vector<unsigned long> exp((std::size_t)n_vars);
        fmpz_t coeff;
        fmpz_init(coeff);
        for (long t = 0; t < len; ++t) {
            fmpz_mpoly_get_term_exp_ui(exp.data(), p.raw(), t, ctx->raw());
            fmpz_mpoly_get_term_coeff_fmpz(coeff, p.raw(), t, ctx->raw());
            for (long v = 0; v < n_vars; ++v) {
                if (v == eta_var) continue;
                if (exp[(std::size_t)v] != 0) {
                    fmpz_clear(coeff);
                    throw std::runtime_error(
                        "mfrac_to_eta_rational: residual symbolic var '"
                        + ctx->var_name(v) + "'");
                }
            }
            long k = (eta_var >= 0) ? (long)exp[(std::size_t)eta_var] : 0;
            fmpq_t qc, qexist;
            fmpq_init(qc); fmpq_init(qexist);
            fmpq_set_fmpz(qc, coeff);
            out.coeff(k, qexist);
            fmpq_add(qc, qc, qexist);
            out.set_coeff_fmpq(k, qc);
            fmpq_clear(qc); fmpq_clear(qexist);
        }
        fmpz_clear(coeff);
        return out;
    };
    numeric::FmpqPoly num_p = walk(src.numerator());
    numeric::FmpqPoly den_p = walk(src.denominator());
    return numeric::RationalFunction(std::move(num_p), std::move(den_p));
}

// Approximate an acb value as an fmpq (RationalizePre digits).
void acb_to_fmpq_rational(fmpq_t out, const acb_t v, long prec) {
    int digits = numeric::rationalize_pre();
    if (digits <= 0) digits = 20;

    const long bits = std::max<long>(prec, numeric::decimal_digits_to_bits(digits) + 32);
    arb_t scale, scaled, xb;
    arb_init(scale);
    arb_init(scaled);
    arb_init(xb);
    arb_set_si(scale, 10);
    arb_pow_ui(scale, scale, static_cast<unsigned long>(digits), bits);
    arb_set(xb, acb_realref(v));
    arb_mul(scaled, xb, scale, bits);

    fmpz_t numer, pow10;
    fmpz_init(numer);
    fmpz_init(pow10);
    arf_get_fmpz(numer, arb_midref(scaled), ARF_RND_NEAR);
    fmpz_set_si(pow10, 10);
    fmpz_pow_ui(pow10, pow10, static_cast<unsigned long>(digits));
    fmpq_set_fmpz_frac(out, numer, pow10);
    fmpq_canonicalise(out);

    fmpz_clear(numer);
    fmpz_clear(pow10);
    arb_clear(scale);
    arb_clear(scaled);
    arb_clear(xb);
}

}  // namespace

// ===========================================================================
//  Counter for system IDs
// ===========================================================================

namespace {
std::atomic<long> g_amf_system_id_counter{0};
}

// ===========================================================================
//  AMFSystem ctor / dtor / move
// ===========================================================================

AMFSystem::AMFSystem(qft::FamilyConfig fc,
                       std::vector<qft::JIntegral> preferred,
                       std::vector<int>       etac,
                       EndingScheme           ending_used,
                       const AMFSystemOptions& opts)
    : preferred_(std::move(preferred)),
      etac_(std::move(etac)),
      ending_(ending_used),
      opts_(opts),
      system_id_(g_amf_system_id_counter.fetch_add(1)) {
    bool all_zero = std::all_of(etac_.begin(), etac_.end(),
                                  [](int v) { return v == 0; });
    is_ending_ = all_zero;
    fc_ = std::make_unique<qft::FamilyConfig>(std::move(fc));
}

AMFSystem::AMFSystem(AMFSystem&&) noexcept            = default;
AMFSystem& AMFSystem::operator=(AMFSystem&&) noexcept = default;
AMFSystem::~AMFSystem() = default;

// ===========================================================================
//  EndingQ
// ===========================================================================

// Forward decl; implementation below uses topology analysis.
namespace {
bool single_mass_ending_q_impl(const qft::FamilyConfig& fc,
                                  const std::vector<qft::JIntegral>& preferred,
                                  const AMFSystemOptions& opts);
bool cutkosky_ending_q_impl(const qft::FamilyConfig& fc,
                              const std::vector<qft::JIntegral>& preferred);
}

bool ending_q(const qft::FamilyConfig& fc,
              const std::vector<qft::JIntegral>& preferred,
              EndingScheme scheme,
              const AMFSystemOptions& opts) {
    if (preferred.empty()) return true;
    auto top = qft::get_top_position(preferred);
    if (top.empty()) return true;

    switch (scheme) {
        case EndingScheme::Tradition: {
            auto pos = qft::amf_position(fc, top, opts.amf_modes);
            return pos.empty();
        }
        case EndingScheme::SingleMass:
            return single_mass_ending_q_impl(fc, preferred, opts);
        case EndingScheme::Cutkosky:
            return cutkosky_ending_q_impl(fc, preferred);
    }
    return false;
}

namespace {

// Apply numeric substitution (from opts.bb.numeric_values) to an algebra::Mfrac
// living on a context that prefixes fc.ctx (so its variables map by
// position to fc.ctx variables).  Variables not in numeric_values are
// left symbolic.  Returns a new algebra::Mfrac on the same ctx.
algebra::Mfrac apply_numeric(const algebra::Mfrac& src,
                     const std::map<std::string, FmpqHolder>& numeric_q) {
    auto ctx = src.ctx();
    algebra::Mfrac cur = src.clone();
    for (long v = 0; v < ctx->n_vars(); ++v) {
        const std::string& name = ctx->var_name(v);
        auto it = numeric_q.find(name);
        if (it == numeric_q.end()) continue;
        cur = cur.substitute(v, it->second.raw());
    }
    return cur;
}

// Build numeric_q from opts.bb.numeric_values.
std::map<std::string, FmpqHolder>
build_numeric_q(const ibp::ReduceOptions& bb) {
    std::map<std::string, FmpqHolder> out;
    auto parse = [&](const std::string& s, fmpq_t v) -> bool {
        std::size_t slash = s.find('/');
        try {
            if (slash != std::string::npos) {
                long p = std::stol(s.substr(0, slash));
                long q = std::stol(s.substr(slash + 1));
                fmpq_set_si(v, p, q);
                return true;
            }
            std::size_t pos;
            long p = std::stol(s, &pos);
            if (pos == s.size()) {
                fmpq_set_si(v, p, 1);
                return true;
            }
            // Decimal.
            std::size_t dot = s.find('.');
            if (dot == std::string::npos) return false;
            std::string sint = s.substr(0, dot) + s.substr(dot + 1);
            long denom = 1;
            for (std::size_t i = dot + 1; i < s.size(); ++i) denom *= 10;
            long num = std::stol(sint);
            fmpq_set_si(v, num, denom);
            return true;
        } catch (...) { return false; }
    };
    for (const auto& [name, val_str] : bb.numeric_values) {
        FmpqHolder h;
        if (parse(val_str, h.raw())) out.emplace(name, std::move(h));
    }
    return out;
}

std::map<std::string, FmpqHolder>
build_family_numeric_q(const ibp::ReduceOptions& bb) {
    auto out = build_numeric_q(bb);
    out.erase("eps");
    out.erase("__amf_eps");
    out.erase("d");
    out.erase("eta");
    out.erase("__amf_eta");
    return out;
}

std::vector<std::pair<std::string, std::string>>
extract_replacement_strings_with_numeric(
    const qft::FamilyConfig& fc,
    const std::map<std::string, FmpqHolder>& numeric_q) {
    std::vector<std::pair<std::string, std::string>> out;
    for (const auto& [key, val] : fc.reduced_replacement) {
        out.emplace_back(key, mfrac_str(apply_numeric(val, numeric_q)));
    }
    return out;
}

qft::FamilyConfig build_numeric_replacement_family(
    const qft::FamilyConfig& fc,
    const std::map<std::string, FmpqHolder>& numeric_q) {
    return qft::FamilyConfig::build(
        fc.family,
        fc.loops,
        fc.legs,
        extract_conservation_strings(fc),
        extract_replacement_strings_with_numeric(fc, numeric_q),
        extract_propagator_strings(fc),
        fc.cut,
        fc.prescription);
}

std::map<std::string, FmpqHolder>
build_numeric_q_for_eps(const ibp::ReduceOptions& bb,
                        const numeric::AcbValue& eps,
                        long prec) {
    std::map<std::string, FmpqHolder> out;
    {
        FmpqHolder h;
        acb_to_fmpq_rational(h.raw(), eps.raw(), prec);
        out.emplace("eps", h);
        out.emplace("__amf_eps", h);
    }

    auto parse = [&](const std::string& s, fmpq_t v) -> bool {
        std::size_t slash = s.find('/');
        try {
            if (slash != std::string::npos) {
                long p = std::stol(s.substr(0, slash));
                long q = std::stol(s.substr(slash + 1));
                fmpq_set_si(v, p, q);
                return true;
            }
            std::size_t pos = 0;
            long p = std::stol(s, &pos);
            if (pos == s.size()) {
                fmpq_set_si(v, p, 1);
                return true;
            }
            std::stod(s, &pos);
            if (pos != s.size()) return false;
            std::size_t dot = s.find('.');
            if (dot == std::string::npos) return false;
            std::string sint = s.substr(0, dot) + s.substr(dot + 1);
            long denom = 1;
            for (std::size_t i = dot + 1; i < s.size(); ++i) denom *= 10;
            long num = std::stol(sint);
            fmpq_set_si(v, num, denom);
            return true;
        } catch (...) {
            return false;
        }
    };

    for (const auto& [name, val_str] : bb.numeric_values) {
        if (out.count(name)) continue;
        FmpqHolder h;
        if (!parse(val_str, h.raw())) {
            throw std::runtime_error(
                "build_numeric_q_for_eps: cannot parse numeric '"
                + val_str + "' for '" + name + "'");
        }
        out.emplace(name, std::move(h));
    }
    return out;
}

void evaluate_const_mfrac(const algebra::Mfrac& m,
                          const std::map<std::string, FmpqHolder>& numeric_q,
                          long prec,
                          numeric::AcbValue& out) {
    std::set<std::string> empty_keep;
    algebra::Mfrac sub = mfrac_substitute(m, numeric_q, empty_keep);
    algebra::Mpoly num = sub.numerator();
    algebra::Mpoly den = sub.denominator();
    fmpz_t nint, dint;
    fmpz_init(nint);
    fmpz_init(dint);
    std::vector<unsigned long> z((std::size_t)sub.ctx()->n_vars(), 0);
    fmpz_mpoly_get_coeff_fmpz_ui(nint, num.raw(), z.data(), sub.ctx()->raw());
    fmpz_mpoly_get_coeff_fmpz_ui(dint, den.raw(), z.data(), sub.ctx()->raw());
    fmpq_t q;
    fmpq_init(q);
    fmpq_set_fmpz_frac(q, nint, dint);
    acb_set_fmpq(out.raw(), q, prec);
    fmpq_clear(q);
    fmpz_clear(nint);
    fmpz_clear(dint);
}

bool try_builtin_ending_value(const qft::FamilyConfig& fc,
                              const qft::JIntegral& j,
                              const AMFSystemOptions& opts,
                              std::size_t eps_index,
                              const numeric::AcbValue& eps,
                              numeric::AcbValue& out) {
    const std::string key = jintegral_key(j);
    auto eit = opts.explicit_boundary.find(key);
    if (eit != opts.explicit_boundary.end() && eit->second) {
        if (eps_index >= eit->second->size()) {
            throw std::runtime_error(
                "AMFSystem::solve: explicit_boundary too short for '" + key + "'");
        }
        out = (*eit->second)[eps_index].clone();
        return true;
    }

    bool any_pos = false;
    for (long v : j.indices()) {
        if (v > 0) {
            any_pos = true;
            break;
        }
    }
    if (!any_pos) {
        acb_one(out.raw());
        return true;
    }

    long n_props = 0;
    bool all_one_or_zero = true;
    for (long v : j.indices()) {
        if (v == 1) {
            ++n_props;
        } else if (v != 0) {
            all_one_or_zero = false;
            break;
        }
    }
    if (all_one_or_zero && n_props > 0
            && qft::vacuum_known((long)fc.n_loops(), n_props)) {
        qft::vacuum((long)fc.n_loops(), n_props, eps.raw(), out.raw(),
               numeric::working_prec_bits());
        return true;
    }

    return false;
}

ibp::ReduceResult&
get_ending_reduction(const qft::FamilyConfig& fc,
                     const qft::JIntegral& target,
                     long system_id,
                     const AMFSystemOptions& opts,
                     std::map<std::string, ibp::ReduceResult>& cache) {
    const std::string key = jintegral_key(target);
    auto it = cache.find(key);
    if (it != cache.end()) return it->second;

    ibp::ReduceOptions reduce_opts = opts.bb;
    if (!opts.cache_root.empty()) {
        reduce_opts.work_dir = opts.cache_root + "/system_"
                             + std::to_string(system_id)
                             + "_ending_reduce_" + key;
    } else if (reduce_opts.work_dir.empty()) {
        reduce_opts.work_dir = (fs::temp_directory_path() /
            ("amflow_amfsys_" + std::to_string(system_id)
             + "_ending_reduce_" + key)).string();
    } else {
        reduce_opts.work_dir += "/ending_reduce_" + key;
    }
    fs::remove_all(reduce_opts.work_dir);

    auto inserted = cache.emplace(
        key, ibp::reduce(fc, {target}, /*preferred=*/{},
                              qft::get_top_sector({target}), reduce_opts));
    return inserted.first->second;
}

numeric::AcbValue solve_ending_master_value(
        const qft::FamilyConfig& fc,
        const qft::JIntegral& target,
        long system_id,
        const AMFSystemOptions& opts,
        std::size_t eps_index,
        const numeric::AcbValue& eps,
        const std::map<std::string, FmpqHolder>& numeric_q,
        std::map<std::string, ibp::ReduceResult>& reduction_cache,
        std::map<std::string, numeric::AcbValue>& value_cache,
        std::set<std::string>& active_keys,
        long prec) {
    const std::string key = jintegral_key(target);
    auto vit = value_cache.find(key);
    if (vit != value_cache.end()) return vit->second.clone();

    if (!active_keys.insert(key).second) {
        throw std::runtime_error(
            "AMFSystem::solve: cyclic ending reduction for '" + key + "'");
    }

    numeric::AcbValue out;
    if (try_builtin_ending_value(fc, target, opts, eps_index, eps, out)) {
        value_cache.emplace(key, midpoint_only(out));
        active_keys.erase(key);
        return out;
    }

    ibp::ReduceResult& red =
        get_ending_reduction(fc, target, system_id, opts, reduction_cache);
    if (red.rules.size() != 1 || red.rules[0].empty()) {
        active_keys.erase(key);
        throw std::runtime_error(
            "AMFSystem::solve: ending system master '" + key
            + "' has no Vacuum entry, no explicit_boundary value, "
            "and no usable reduction rule");
    }

    numeric::AcbValue rhs_sum;
    acb_zero(rhs_sum.raw());
    numeric::AcbValue self_coef;
    acb_zero(self_coef.raw());
    bool saw_non_self = false;
    for (const auto& term : red.rules[0]) {
        numeric::AcbValue coef;
        evaluate_const_mfrac(term.coef, numeric_q, prec, coef);
        if (jintegral_key(term.integ) == key) {
            acb_add(self_coef.raw(), self_coef.raw(), coef.raw(), prec);
            continue;
        }
        saw_non_self = true;
        numeric::AcbValue rhs = solve_ending_master_value(
            fc, term.integ, system_id, opts, eps_index, eps, numeric_q,
            reduction_cache, value_cache, active_keys, prec);
        numeric::AcbValue prod;
        acb_mul(prod.raw(), coef.raw(), rhs.raw(), prec);
        acb_add(rhs_sum.raw(), rhs_sum.raw(), prod.raw(), prec);
    }

    numeric::AcbValue denom;
    acb_one(denom.raw());
    acb_sub(denom.raw(), denom.raw(), self_coef.raw(), prec);

    if (!saw_non_self && acb_is_one(self_coef.raw())) {
        active_keys.erase(key);
        throw std::runtime_error(
            "AMFSystem::solve: ending system master '" + key
            + "' has no Vacuum entry, no explicit_boundary value, "
            "and reduction did not lower it");
    }

    acb_div(out.raw(), rhs_sum.raw(), denom.raw(), prec);
    value_cache.emplace(key, midpoint_only(out));
    active_keys.erase(key);
    return out;
}

// Test whether a qft::TopSectorComponentInfo, after applying numeric
// substitution to mass entries, has exactly one mass=1 and all others
// mass=0.
bool vacuum_q_numeric(const qft::TopSectorComponentInfo& info,
                      const std::map<std::string, FmpqHolder>&) {
    return qft::vacuum_q(info);
}

bool single_mass_q_numeric(const qft::TopSectorComponentInfo& info,
                              const std::map<std::string, FmpqHolder>& numeric_q) {
    if (!vacuum_q_numeric(info, numeric_q)) return false;
    int n_one = 0, n_zero = 0;
    for (const auto& m : info.mass) {
        algebra::Mfrac sub = apply_numeric(m, numeric_q);
        if (sub.is_zero()) ++n_zero;
        else if (sub.is_one()) ++n_one;
        else return false;   // some mass is neither 0 nor 1 -- not single-mass
    }
    return n_one == 1 && n_zero == (int)info.mass.size() - 1;
}

bool ending_q_numeric(const qft::TopSectorComponentInfo& info,
                       const std::map<std::string, FmpqHolder>& numeric_q) {
    return single_mass_q_numeric(info, numeric_q) || qft::phase_volume_q(info);
}

// Mirrors AMFSystemEndingQ[preferred_, "Cutkosky"] (.m line 996-1000):
//   ending = ! (all components are EndingQ and exactly one component is
//               PhaseVolumeQ)
bool cutkosky_ending_q_impl(const qft::FamilyConfig& fc,
                              const std::vector<qft::JIntegral>& preferred) {
    auto top = qft::get_top_position(preferred);
    if (top.empty()) return true;
    auto info = qft::analyze_top_sector(fc, top);
    if (info.empty()) return true;

    long phase_count = 0;
    for (const auto& comp : info) {
        if (!qft::ending_q(comp)) return true;
        if (qft::phase_volume_q(comp)) ++phase_count;
    }
    return phase_count != 1;
}

// Mirrors AMFSystemEndingQ[preferred_, "SingleMass"] (.m line 1002-1007):
//   ending = ! ( all components ending && cutcom==0 && Loop>0 && all idx>=0 )
// i.e. SingleMass is *applicable* (ending=False) only when all four
// conditions hold. In particular, a heterogeneous-mass vacuum is still
// ending in the Mathematica reference because EndingQ = SingleMassQ ||
// PhaseVolumeQ, so such a component fails condition 1 and never reaches
// the later FactorizeFamily setup branch.
bool single_mass_ending_q_impl(const qft::FamilyConfig& fc,
                                  const std::vector<qft::JIntegral>& preferred,
                                  const AMFSystemOptions& opts) {
    auto top = qft::get_top_position(preferred);
    if (top.empty()) return true;          // empty top -> already vacuum-like, ending
    auto info = qft::analyze_top_sector(fc, top);
    if (info.empty()) return true;         // no analyzed components -> SingleMass not applicable

    auto numeric_q = build_numeric_q(opts.bb);

    // Early exit for the same Mathematica condition-1 failure: a pure
    // vacuum component with masses other than {1,0,...,0} is VacuumQ but
    // not SingleMassQ, hence not EndingQ.
    for (const auto& comp : info) {
        if (vacuum_q_numeric(comp, numeric_q)
                && !single_mass_q_numeric(comp, numeric_q)) {
            return true;
        }
    }

    // condition 1: every component is "EndingQ" under numeric
    // (= single-mass-vacuum or phase-volume)
    for (const auto& comp : info) {
        if (!ending_q_numeric(comp, numeric_q)) return true;  // not ending -> SingleMass not applicable
    }
    // condition 2: no phase-volume (cut) component
    for (const auto& comp : info) {
        if (qft::phase_volume_q(comp)) return true;
    }
    // condition 3: Length[Loop] > 0
    if (fc.n_loops() == 0) return true;
    // condition 4: all preferred indices non-negative
    for (const auto& j : preferred) {
        for (long v : j.indices()) if (v < 0) return true;
    }
    // All four conditions satisfied -> SingleMass is applicable -> NOT ending.
    return false;
}

}  // namespace

bool ending_q_all(const qft::FamilyConfig& fc,
                  const std::vector<qft::JIntegral>& preferred,
                  const AMFSystemOptions& opts) {
    for (auto s : opts.ending_schemes) {
        if (!ending_q(fc, preferred, s, opts)) return false;
    }
    return true;
}

// ===========================================================================
//  Setup -- top level
// ===========================================================================
//
//   * If is_ending_: nothing to do.  solve() will look up Vacuum or
//     explicit_boundary.
//   * Else:
//       a. build fc_with_eta_ = inject(fc_, etac_);
//       b. build_diffeq() -- compute dM/dη on fc_with_eta_;
//       c. build_boundary() -- find regions, compute boundary
//          integrals, recursively setup sub-systems.

void AMFSystem::setup() {
    if (!is_ending_) {
        // We'll show after injection.
    }
    if (is_ending_) return;
    if (fc_with_eta_) return;
    if (opts_.current_depth >= opts_.max_recursion_depth) {
        throw std::runtime_error(
            "AMFSystem::setup: max recursion depth ("
            + std::to_string(opts_.max_recursion_depth)
            + ") reached.  This typically means SingleMass / "
              "FactorizeFamily is needed to terminate the boundary "
              "integrand recursion -- see AUDIT.md.");
    }
    fc_with_eta_ = build_injected_fc(*fc_, etac_);
    build_diffeq();
    build_boundary();
}

// ---------------------------------------------------------------------------
//  build_diffeq
// ---------------------------------------------------------------------------

void AMFSystem::build_diffeq() {
    ibp::ReduceOptions bb = opts_.bb;
    if (!opts_.cache_root.empty()) {
        bb.work_dir = opts_.cache_root + "/system_"
                     + std::to_string(system_id_) + "_diffeq";
    } else if (bb.work_dir.empty()) {
        bb.work_dir = (fs::temp_directory_path() /
            ("amflow_amfsys_" + std::to_string(system_id_) + "_diffeq")).string();
    }
    auto top_pattern = qft::get_top_sector(preferred_);
    diffeq_result_ = ibp::diffeq(*fc_with_eta_, preferred_,
                                        {"eta"}, top_pattern, bb);

    // Map preferred[i] -> sortedmasters index.
    std::map<std::string, long> sorted_index;
    for (std::size_t i = 0; i < diffeq_result_.sortedmasters.size(); ++i) {
        sorted_index[jintegral_key(diffeq_result_.sortedmasters[i])] = (long)i;
    }
    pref_to_sorted_.clear();
    pref_to_sorted_.reserve(preferred_.size());
    for (const auto& j : preferred_) {
        auto it = sorted_index.find(jintegral_key(j));
        if (it == sorted_index.end()) {
            throw std::runtime_error(
                "AMFSystem::build_diffeq: preferred master '"
                + jintegral_key(j) + "' not in Kira's master list");
        }
        pref_to_sorted_.push_back(it->second);
    }
    AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
        std::cerr << "[build_diffeq] family=" << fc_->family
                  << " preferred_to_sorted=" << pref_to_sorted_.size()
                  << std::endl;
        for (std::size_t i = 0; i < preferred_.size(); ++i) {
            const auto& pref = preferred_[i];
            const long si = pref_to_sorted_[i];
            std::cerr << "  preferred[" << i << "] = J[" << pref.family();
            for (long idx : pref.indices()) std::cerr << "," << idx;
            std::cerr << "] -> sorted[" << si << "]";
            if (si >= 0 && (std::size_t)si < diffeq_result_.sortedmasters.size()) {
                const auto& sm = diffeq_result_.sortedmasters[(std::size_t)si];
                std::cerr << " = J[" << sm.family();
                for (long idx : sm.indices()) std::cerr << "," << idx;
                std::cerr << "]";
            }
            std::cerr << std::endl;
        }
    }
}

// ---------------------------------------------------------------------------
//  build_boundary
// ---------------------------------------------------------------------------
//
//   Mirrors the inner Table loop in AMFSystemBoundaryCondition +
//   ReduceBoundary + AMFSystemsSetup (the recursive part).
//
//   For each non-zero region:
//     1. powers = RegionPower(integrals, region)
//     2. (boundary integrand families, terms) = BoundaryIntegrals(BoundaryIntegrands(...))
//     3. For each family (i.e. each sub-qft::FamilyConfig):
//          a. Build sub_fc by replacing parent's Propagator with the
//             family's prop list, keeping loops/legs/conservation/
//             replacement/numeric/cut/prescription from parent.
//          b. Collect all J integrals appearing in the family's terms.
//          c. Run BlackBoxReduce(those Js, {}) on sub_fc -> (sub_masters,
//             sub_rules).
//          d. Recursively run amf_systems_setup(sub_fc, sub_masters):
//             this either ends (vacuum) or further injects eta.
//          e. Compute the per-(master, term-order) "table" linking
//             boundary integrand terms to sub_masters via sub_rules.

// Helper: evaluate an algebra::Mfrac at ε = 10^-4 (the AMFlow convention) and
// any other vars at 0.  Returns the value as fmpq.  Used by
// build_boundary's adaptive border computation.
namespace {
void mfrac_eval_at_small_eps(const algebra::Mfrac& m, fmpq_t out) {
    auto ctx = m.ctx();
    long n_vars = ctx->n_vars();
    long eps_var = ctx->var_index("eps");
    if (eps_var < 0) eps_var = ctx->var_index("__amf_eps");

    auto eval_q = [&](const algebra::Mpoly& p, fmpq_t res) {
        fmpq_t accum, term, eps_pow, base;
        fmpq_init(accum); fmpq_init(term);
        fmpq_init(eps_pow); fmpq_init(base);
        std::vector<unsigned long> exp((std::size_t)n_vars);
        fmpz_t coeff;
        fmpz_init(coeff);
        long len = fmpz_mpoly_length(p.raw(), ctx->raw());
        for (long t = 0; t < len; ++t) {
            fmpz_mpoly_get_term_exp_ui(exp.data(), p.raw(),
                                        t, ctx->raw());
            fmpz_mpoly_get_term_coeff_fmpz(coeff, p.raw(), t, ctx->raw());
            // Skip if any non-eps variable has positive exp.
            bool drop = false;
            for (long v = 0; v < n_vars; ++v) {
                if (v == eps_var) continue;
                if (exp[(std::size_t)v] != 0) { drop = true; break; }
            }
            if (drop) continue;
            fmpq_set_fmpz(term, coeff);
            if (eps_var >= 0 && exp[(std::size_t)eps_var] > 0) {
                fmpq_set_si(eps_pow, 1, 1);
                fmpq_set_si(base,    1, 100);   // matches num_q
                for (unsigned long k = 0;
                        k < exp[(std::size_t)eps_var]; ++k) {
                    fmpq_mul(eps_pow, eps_pow, base);
                }
                fmpq_mul(term, term, eps_pow);
            }
            fmpq_add(accum, accum, term);
        }
        fmpz_clear(coeff);
        fmpq_set(res, accum);
        fmpq_clear(accum); fmpq_clear(term);
        fmpq_clear(eps_pow); fmpq_clear(base);
    };

    fmpq_t qn, qd;
    fmpq_init(qn); fmpq_init(qd);
    eval_q(m.numerator(),   qn);
    eval_q(m.denominator(), qd);
    if (fmpq_is_zero(qd)) {
        fmpq_zero(out);
    } else {
        fmpq_div(out, qn, qd);
    }
    fmpq_clear(qn); fmpq_clear(qd);
}
}  // namespace

void AMFSystem::build_boundary() {
    auto rctx = qft::make_region_context(*fc_with_eta_);
    auto top_pattern = qft::get_top_sector(preferred_);
    std::vector<std::size_t> top_posi;
    for (std::size_t i = 0; i < top_pattern.size(); ++i) {
        if (top_pattern[i]) top_posi.push_back(i);
    }
    if (top_posi.empty()) return;   // no boundary needed

    auto regions = qft::find_all_region(*fc_with_eta_, rctx, top_posi);
    std::vector<qft::Region> non_zero;
    for (auto& r : regions) {
        if (!qft::zero_region_q(*fc_with_eta_, rctx, r, top_posi)) {
            non_zero.push_back(std::move(r));
        }
    }
    if (non_zero.empty()) return;

    auto pctx = qft::make_powers_context(*fc_with_eta_);

    // We need a BC for *every* master in the diffeq matrix, not just
    // the user-supplied preferred[].  Mathematica builds the BC over
    // sortedmasters.  We use the same set.
    const auto& diff_masters = diffeq_result_.sortedmasters;

    // Step 1 (AMFSystemBoundaryOrder, .m line 922-926):
    //   Compute powers per region, then compute pattern via
    //   BoundaryPattern.  pattern is shape [group_id][integral_id].
    std::vector<std::vector<algebra::Mfrac>> all_powers;
    all_powers.reserve(non_zero.size());
    for (const auto& region : non_zero) {
        all_powers.push_back(qft::region_power(*fc_with_eta_, rctx, pctx,
                                            region, diff_masters));
    }
    auto pattern = qft::boundary_pattern(all_powers);
    // Stash the full pattern list; solve_one_eps will use it to append
    // (mu_pattern_i -> 0) entries to bc_sorted.  Mirrors MMA's final
    //   bc = Join[#1, Thread[(pattern/.epsrule) -> 0]] step.
    bc_pattern_.clear();
    bc_pattern_.reserve(pattern.size());
    for (const auto& pg : pattern) {
        std::vector<algebra::Mfrac> row;
        row.reserve(pg.size());
        for (const auto& m : pg) row.push_back(m.clone());
        bc_pattern_.push_back(std::move(row));
    }

    // Step 2 (AMFSystemBoundaryOrder wlscript template):
    //   Compute deinf = -η^-2 (diffeq /. eta -> 1/eta), evaluated at
    //   ε = 10^-4 (the AMFlow convention).
    //
    //   Then npattern = -pattern /. epsrule (per-pattern, per-integral).
    //   Then orders[g] = DetermineBoundaryOrder[deinf, npattern[g]] for
    //   each pattern group.
    long n_masters = (long)diff_masters.size();
    if ((long)diffeq_result_.diffeq.size() != 1) {
        throw std::runtime_error(
            "AMFSystem::build_boundary: expected single-variable diffeq");
    }
    long eta_in_red = diffeq_result_.red_ctx.ctx->var_index("eta");
    if (eta_in_red < 0) {
        throw std::runtime_error(
            "AMFSystem::build_boundary: red_ctx is missing 'eta'");
    }

    // Build deinf as a numeric::RationalMatrix(η) by:
    //   (a) for each (i, j), evaluate diffeq[i][j] at ε = 10^-4 and
    //       any numeric_values
    //   (b) substitute eta -> 1/eta (as a polynomial transform)
    //   (c) multiply by -1/η²
    //
    // For step (a) we need numeric values for any non-eta vars.
    // Use opts_.bb.numeric_values + ε = 1e-4 fallback.
    // Mathematica uses ε = 10^-4 in the wlscript; we use a smaller-
    // denominator fraction to avoid exact-rational arithmetic blowup
    // (Mathematica falls back to machine-precision floats internally).
    std::map<std::string, FmpqHolder> num_q;
    {
        FmpqHolder h;
        fmpq_set_si(h.raw(), 1, 100);
        num_q.emplace("eps", h);
        FmpqHolder h2;
        fmpq_set(h2.raw(), h.raw());
        num_q.emplace("__amf_eps", h2);
    }
    auto parse_str = [](const std::string& s, fmpq_t v) -> bool {
        std::size_t slash = s.find('/');
        try {
            if (slash != std::string::npos) {
                long p = std::stol(s.substr(0, slash));
                long q = std::stol(s.substr(slash + 1));
                fmpq_set_si(v, p, q); return true;
            }
            std::size_t pos;
            long p = std::stol(s, &pos);
            if (pos == s.size()) {
                fmpq_set_si(v, p, 1); return true;
            }
            std::size_t dot = s.find('.');
            if (dot == std::string::npos) return false;
            std::string sint = s.substr(0, dot) + s.substr(dot + 1);
            long denom = 1;
            for (std::size_t i = dot + 1; i < s.size(); ++i) denom *= 10;
            long num = std::stol(sint);
            fmpq_set_si(v, num, denom); return true;
        } catch (...) { return false; }
    };
    for (const auto& [name, val_str] : opts_.bb.numeric_values) {
        if (num_q.count(name)) continue;
        FmpqHolder h;
        if (parse_str(val_str, h.raw())) num_q.emplace(name, std::move(h));
    }

    std::set<std::string> keep_eta = {"eta"};
    numeric::RationalMatrix deinf((std::size_t)n_masters, (std::size_t)n_masters);
    const auto& dmat = diffeq_result_.diffeq[0];
    for (long i = 0; i < n_masters; ++i) {
        for (long j = 0; j < n_masters; ++j) {
            // Evaluate dmat[i][j] at numeric values (keep eta).
            algebra::Mfrac sub = mfrac_substitute(dmat[(std::size_t)i][(std::size_t)j],
                                            num_q, keep_eta);
            // Convert to numeric::RationalFunction in eta.
            numeric::RationalFunction r = mfrac_to_eta_rational(sub, eta_in_red);
            deinf(i, j) = std::move(r);
        }
    }
    // Apply  eta -> 1/eta  and multiply by -1/η²:
    //   deinf'(η) = -1/η^2 * deinf(1/η)
    //
    // numeric::RationalFunction::substitute_inverse_eta implements only r(1/η);
    // we must explicitly multiply by -1/η² afterwards (mirrors AMFlow's
    // `deinf = -eta^-2 * (de /. eta -> 1/eta)`).  Missing the -1/η² factor
    // caused DetermineBoundaryOrder to see a wrong matrix (e.g. a non-
    // proper ODE with an η in the numerator instead of a 1/η pole), which
    // in turn produced wrong `border[i]` values and left `rb.integrand_terms`
    // empty for the non-zero-sector regions.  See Stage B notes.
    fmpq_t neg_one;
    fmpq_init(neg_one);
    fmpq_set_si(neg_one, -1, 1);
    numeric::RationalFunction prefactor = numeric::RationalFunction::monomial(-2, neg_one);
    fmpq_clear(neg_one);
    for (std::size_t i = 0; i < (std::size_t)n_masters; ++i) {
        for (std::size_t j = 0; j < (std::size_t)n_masters; ++j) {
            deinf(i, j) = deinf(i, j).substitute_inverse_eta() * prefactor;
        }
    }


    // Step 3: For each pattern group g, compute orders[g] (length n_masters).
    //
    //   npattern[g][i] = -pattern[g][i] /. epsrule
    //
    //   orders[g] = ode::determine_boundary_order(deinf, npattern[g])
    std::vector<std::vector<long>> orders;
    orders.reserve(pattern.size());
    for (const auto& pat : pattern) {
        if ((long)pat.size() != n_masters) {
            throw std::runtime_error(
                "AMFSystem::build_boundary: pattern size mismatch");
        }
        // Build npat_q: per-integral ini value for DetermineBoundaryOrder.
        // We DO NOT round to integer — MMA accepts rational ini and gives
        // different answers for -1 vs -99/100 (see Stage B notes).
        std::vector<numeric::RationalFunction> npat_q;
        npat_q.reserve(pat.size());
        for (const auto& m : pat) {
            fmpq_t q;
            fmpq_init(q);
            mfrac_eval_at_small_eps(m, q);
            fmpq_neg(q, q);    // npattern = -pattern /. epsrule
            npat_q.push_back(numeric::RationalFunction::from_fmpq(q));
            fmpq_clear(q);
        }

        std::vector<long> ord = ode::determine_boundary_order(deinf, npat_q);
        AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
            std::cerr << "[build_boundary] family=" << fc_->family
                      << " pattern_group " << (orders.size())
                      << " npat_q=[";
            for (std::size_t i = 0; i < npat_q.size(); ++i) {
                if (i) std::cerr << ",";
                std::cerr << npat_q[i].to_string();
            }
            std::cerr << "] orders=[";
            for (std::size_t i = 0; i < ord.size(); ++i) {
                if (i) std::cerr << ",";
                std::cerr << ord[i];
            }
            std::cerr << "]" << std::endl;
        }
        orders.push_back(std::move(ord));
    }

    // Step 4: For each region, compute per-master border:
    //   k = first pattern group containing region's powers (matched by
    //       integer-shift on integral 0)
    //   border[r][i] = orders[k][i] - (pattern[k][i] - powers[r][i])
    //                 (intent: how many extra orders we need to expand
    //                  to capture all the boundary information for
    //                  master i in region r)
    //   border = Map[If[#<0, -1, #]&, border, {2}]
    //
    // We still need to find which group each region belongs to.
    // Local mfrac_is_rational: returns true iff `m` is a rational
    // constant (no algebra::Mpoly variables in num or den).  Sets `out` to the
    // rational value.
    auto mfrac_is_rat = [](const algebra::Mfrac& m, fmpq_t out) -> bool {
        algebra::Mpoly num = m.numerator();
        algebra::Mpoly den = m.denominator();
        auto ctx = m.ctx();
        std::vector<unsigned long> z((std::size_t)ctx->n_vars(), 0);
        // Check that num is purely the const term (no other monomials).
        long len_n = fmpz_mpoly_length(num.raw(), ctx->raw());
        long len_d = fmpz_mpoly_length(den.raw(), ctx->raw());
        if (len_n > 1 || len_d > 1) return false;
        if (len_n == 1) {
            std::vector<unsigned long> e((std::size_t)ctx->n_vars());
            fmpz_mpoly_get_term_exp_ui(e.data(), num.raw(), 0, ctx->raw());
            for (long v = 0; v < ctx->n_vars(); ++v) {
                if (e[(std::size_t)v] != 0) return false;
            }
        }
        if (len_d == 1) {
            std::vector<unsigned long> e((std::size_t)ctx->n_vars());
            fmpz_mpoly_get_term_exp_ui(e.data(), den.raw(), 0, ctx->raw());
            for (long v = 0; v < ctx->n_vars(); ++v) {
                if (e[(std::size_t)v] != 0) return false;
            }
        }
        fmpz_t cn, cd;
        fmpz_init(cn); fmpz_init(cd);
        fmpz_mpoly_get_coeff_fmpz_ui(cn, num.raw(), z.data(), ctx->raw());
        fmpz_mpoly_get_coeff_fmpz_ui(cd, den.raw(), z.data(), ctx->raw());
        if (fmpz_is_zero(cd)) {
            fmpz_clear(cn); fmpz_clear(cd);
            return false;
        }
        fmpq_set_fmpz_frac(out, cn, cd);
        fmpz_clear(cn); fmpz_clear(cd);
        return true;
    };

    auto find_group = [&](std::size_t r) -> long {
        const algebra::Mfrac& this_pow0 = all_powers[r][0];
        for (std::size_t g = 0; g < pattern.size(); ++g) {
            // Test IntegerQ[pattern[g][0] - powers[r][0]]
            algebra::Mfrac diff = pattern[g][0].clone();
            diff -= this_pow0;
            fmpq_t q;
            fmpq_init(q);
            bool ok = mfrac_is_rat(diff, q);
            bool is_int = ok && fmpz_is_one(fmpq_denref(q));
            fmpq_clear(q);
            if (is_int) return (long)g;
        }
        return -1;
    };

    std::vector<std::vector<long>> region_borders;
    region_borders.reserve(non_zero.size());
    for (std::size_t r = 0; r < non_zero.size(); ++r) {
        long k = find_group(r);
        if (k < 0) {
            throw std::runtime_error(
                "AMFSystem::build_boundary: region " + std::to_string(r)
                + " has no matching pattern group");
        }
        std::vector<long> b((std::size_t)n_masters, -1);
        for (long i = 0; i < n_masters; ++i) {
            // pattern[k][i] - powers[r][i] should be a rational number
            // (small integer + fractional from ε).  We round.
            algebra::Mfrac diff = pattern[(std::size_t)k][(std::size_t)i].clone();
            diff -= all_powers[r][(std::size_t)i];
            fmpq_t q;
            fmpq_init(q);
            mfrac_eval_at_small_eps(diff, q);
            // border = orders[k][i] - (pattern[k][i] - powers[r][i])
            fmpq_t order_q;
            fmpq_init(order_q);
            fmpq_set_si(order_q, orders[(std::size_t)k][(std::size_t)i], 1);
            fmpq_sub(order_q, order_q, q);
            // Round to nearest long.
            fmpz_t fz, abs_num, abs_den, two_abs_num, two_abs_den;
            fmpz_init(fz); fmpz_init(abs_num); fmpz_init(abs_den);
            fmpz_init(two_abs_num); fmpz_init(two_abs_den);
            fmpz_abs(abs_num, fmpq_numref(order_q));
            fmpz_abs(abs_den, fmpq_denref(order_q));
            fmpz_mul_si(two_abs_num, abs_num, 2);
            fmpz_add(two_abs_num, two_abs_num, abs_den);
            fmpz_mul_si(two_abs_den, abs_den, 2);
            fmpz_fdiv_q(fz, two_abs_num, two_abs_den);
            long bv = fmpz_get_si(fz);
            int bv_sgn = fmpq_sgn(order_q);
            if (bv_sgn < 0) bv = -bv;
            fmpz_clear(fz); fmpz_clear(abs_num); fmpz_clear(abs_den);
            fmpz_clear(two_abs_num); fmpz_clear(two_abs_den);
            fmpq_clear(order_q); fmpq_clear(q);
            // Map[If[#<0, -1, #]&, border, {2}]  (clamp negative to -1)
            if (bv < 0) bv = -1;
            b[(std::size_t)i] = bv;
        }
        AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
            std::cerr << "[build_boundary] family=" << fc_->family
                      << " region " << r << " border=[";
            for (std::size_t i = 0; i < b.size(); ++i) {
                if (i) std::cerr << ",";
                std::cerr << b[i];
            }
            std::cerr << "]" << std::endl;
        }
        region_borders.push_back(std::move(b));
    }

    // For each region: build boundary integrand decomposition + sub-system.
    for (std::size_t reg_i = 0; reg_i < non_zero.size(); ++reg_i) {
        const auto& region = non_zero[reg_i];

        // Use the per-master border computed above.
        std::vector<long> border = region_borders[reg_i];

        qft::BoundaryIntegrandsResult bi;
        std::vector<qft::BoundaryFamily> fams;
        try {
            bi = qft::boundary_integrands(*fc_with_eta_, rctx, diff_masters,
                                       border, region);
            fams = qft::boundary_integrals(*fc_with_eta_, rctx, bi);
        } catch (const std::exception& e) {
            throw std::runtime_error(
                "AMFSystem::build_boundary: region " + std::to_string(reg_i)
                + ": " + e.what());
        }
        AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
            std::cerr << "[build_boundary] family=" << fc_->family
                      << " reg_i=" << reg_i
                      << " border=[";
            for (std::size_t i = 0; i < border.size(); ++i) {
                if (i) std::cerr << ",";
                std::cerr << border[i];
            }
            std::cerr << "] fams.size=" << fams.size() << std::endl;
            for (std::size_t fi = 0; fi < fams.size(); ++fi) {
                std::cerr << "  fam[" << fi << "] prop.size="
                          << fams[fi].prop.size()
                          << " terms.size=" << fams[fi].terms.size()
                          << std::endl;
                for (std::size_t ti = 0; ti < fams[fi].terms.size(); ++ti) {
                    std::size_t nonempty = 0, total = 0;
                    for (const auto& tk : fams[fi].terms[ti]) {
                        total += tk.size();
                        if (!tk.empty()) nonempty++;
                    }
                    std::cerr << "    terms[" << ti << "] #ks="
                              << fams[fi].terms[ti].size()
                              << " nonempty_ks=" << nonempty
                              << " total_lts=" << total << std::endl;
                }
            }
        }

        auto make_region_boundary = [&]() {
            RegionBoundary rb;
            rb.powers.reserve(diff_masters.size());
            for (auto& m : all_powers[reg_i]) rb.powers.push_back(m.clone());

            rb.integrand_terms.resize(diff_masters.size());
            for (std::size_t i = 0; i < diff_masters.size(); ++i) {
                for (long oi = 0; oi <= border[i]; ++oi) {
                    rb.integrand_terms[i].emplace_back();
                }
            }
            return rb;
        };

        // Mathematica's ReduceBoundary returns one entry per boundary family:
        //   {powers, config, masters, table}
        // AMFSystemBoundaryCondition then flattens all entries across regions.
        // Store each family as its own RegionBoundary to keep that outer
        // multiplicity distinct from SingleMass multi-root systems.
        for (auto& fam : fams) {
            // 1. Build sub-fc.
            std::vector<std::string> prop_strs;
            prop_strs.reserve(fam.prop.size());
            for (auto& p : fam.prop) prop_strs.push_back(p.to_string());

            // The sub-fc inherits everything from parent except prop.
            // We give it a unique family name; users don't observe it.
            std::string sub_fam_name = fc_->family + "_b"
                + std::to_string(system_id_) + "_r"
                + std::to_string(reg_i) + "_f"
                + std::to_string(subsystems_.size());

            qft::FamilyConfig sub_fc = qft::FamilyConfig::build(
                sub_fam_name,
                fc_->loops, fc_->legs,
                extract_conservation_strings(*fc_),
                extract_replacement_strings(*fc_),
                prop_strs,
                /*cut*/ {},
                fc_->prescription);

            // 2. Collect all J integrals appearing in this family's terms.
            //    Translate from the parent's family name to sub_fam_name.
            std::set<std::string> seen_js;
            std::vector<qft::JIntegral> all_js;
            for (std::size_t i = 0; i < fam.terms.size(); ++i) {
                for (std::size_t k = 0; k < fam.terms[i].size(); ++k) {
                    for (const auto& lt : fam.terms[i][k]) {
                        qft::JIntegral j(sub_fam_name, lt.indices);
                        std::string key = jintegral_key(j);
                        if (seen_js.insert(key).second) all_js.push_back(j);
                    }
                }
            }

            AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
                std::cerr << "[build_boundary] parent=" << fc_->family
                          << " reg_i=" << reg_i
                          << " all_js.size=" << all_js.size() << std::endl;
                for (const auto& j : all_js) {
                    std::cerr << "    j=J[" << j.family();
                    for (auto x : j.indices()) std::cerr << "," << x;
                    std::cerr << "]" << std::endl;
                }
            }
            if (all_js.empty()) continue;

            // 3. Run BlackBoxReduce(all_js, {}) on sub_fc to find
            //    sub-masters + reduction rules.
            //
            //    The Kira call lets it choose its own masters (passing
            //    {} for preferred).
            ibp::ReduceOptions reduce_opts = opts_.bb;
            reduce_opts.work_dir = (fs::temp_directory_path() /
                ("amflow_amfsys_" + std::to_string(system_id_)
                  + "_b" + std::to_string(reg_i)
                  + "_f" + std::to_string(subsystems_.size())
                  + "_p" + std::to_string(getpid())
                  + "_u"
                  + std::to_string(
                        g_boundary_reduce_unique.fetch_add(1,
                            std::memory_order_relaxed)))).string();
            fs::remove_all(reduce_opts.work_dir);

            std::vector<int> top_pat = qft::get_top_sector(all_js);
            ibp::ReduceResult red;
            try {
                red = ibp::reduce(sub_fc, all_js, /*preferred=*/{},
                                          top_pat, reduce_opts);
            } catch (const std::exception& e) {
                throw std::runtime_error(
                    "AMFSystem::build_boundary: BlackBoxReduce failed: "
                    + std::string(e.what()));
            }

            if (red.masters.empty()) {
                // Boundary integrand reduces to zero -- skip family.
                continue;
            }

            RegionBoundary rb = make_region_boundary();

            // 4. Recursively setup the sub-system.  amf_systems_setup
            //    decides ending vs. inject based on red.masters.
            AMFSystemOptions sub_opts = opts_;
            sub_opts.current_depth = opts_.current_depth + 1;
            AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
                std::cerr << "[build_boundary] parent=" << fc_->family
                          << " region=" << reg_i
                          << " sub_fam=" << sub_fam_name
                          << " red.masters.size=" << red.masters.size()
                          << std::endl;
                for (const auto& m : red.masters) {
                    std::cerr << "  master J[" << m.family();
                    for (auto x : m.indices()) std::cerr << "," << x;
                    std::cerr << "]" << std::endl;
                }
                std::cerr << "  sub_fc loops=" << sub_fc.loops.size()
                          << " props=" << sub_fc.propagators_after_conservation.size()
                          << std::endl;
            }
            auto sub_systems = amf_systems_setup(sub_fc, red.masters, sub_opts);

            // Register every root system belonging to this boundary family.
            rb.sub_system_indices.reserve(sub_systems.size());
            for (auto& sub_sys : sub_systems) {
                long sub_idx = (long)subsystems_.size();
                subsystems_.push_back(std::move(sub_sys));
                rb.sub_system_indices.push_back(sub_idx);
            }

            // 5. For each (master, integrand_term), look up its
            //    sub-master coefficient via red.rules.
            //
            // red.rules[k] = list of ibp::DerivTerm{coef, integ} for all_js[k].
            // We build a key->index map for all_js, then use that to
            // resolve LaportaTerm indices; multiply LaportaTerm.coef
            // (in fc_with_eta_.ctx) by red.rules's per-master coef
            // (in red.red_ctx); accumulate into rb.integrand_terms.

            std::map<std::string, std::size_t> all_js_index;
            for (std::size_t k = 0; k < all_js.size(); ++k) {
                all_js_index[jintegral_key(all_js[k])] = k;
            }
            std::map<std::string, std::size_t> sub_master_idx;
            for (std::size_t m = 0; m < red.masters.size(); ++m) {
                sub_master_idx[jintegral_key(red.masters[m])] = m;
            }

            // After the 2024-12 fix to boundary_integrals, fam.terms is
            // shaped [input_i][order_k][laporta_term], with input_i
            // ranging over diff_masters and order_k ranging over
            // border[input_i]+1.  This means fam.terms.size() ==
            // rb.integrand_terms.size() and fam.terms[i].size() ==
            // rb.integrand_terms[i].size().  We assert to catch any
            // regression.
            if (fam.terms.size() != rb.integrand_terms.size()) {
                throw std::runtime_error(
                    "AMFSystem::build_boundary: fam.terms.size "
                    + std::to_string(fam.terms.size()) + " != "
                    + "rb.integrand_terms.size "
                    + std::to_string(rb.integrand_terms.size())
                    + " — boundary_integrals layout drift; see Stage B notes.");
            }
            for (std::size_t i = 0; i < fam.terms.size(); ++i) {
                if (fam.terms[i].size() != rb.integrand_terms[i].size()) {
                    throw std::runtime_error(
                        "AMFSystem::build_boundary: fam.terms[" +
                        std::to_string(i) + "].size " +
                        std::to_string(fam.terms[i].size()) + " != " +
                        "rb.integrand_terms[" + std::to_string(i) +
                        "].size " + std::to_string(rb.integrand_terms[i].size()));
                }
                for (std::size_t k = 0; k < fam.terms[i].size(); ++k) {
                    for (const auto& lt : fam.terms[i][k]) {
                        qft::JIntegral j(sub_fam_name, lt.indices);
                        std::string jkey = jintegral_key(j);
                        auto idx_it = all_js_index.find(jkey);
                        if (idx_it == all_js_index.end()) continue;

                        const auto& rule = red.rules[idx_it->second];
                        // rule is per-master: each entry is (coef, sub_J).
                        // For each entry, accumulate
                        //   contribution = lt.coef (project to fc_->ctx)
                        //                * rule_term.coef (already on red.red_ctx)
                        //                * sub_solution[sub_idx]
                        // We multiply lt.coef and rule_term.coef
                        // symbolically (after projecting both into a
                        // common ctx = sub-system's red_ctx).

                        // lt.coef is on bi.dctx (= parent fc + DListSymbol
                        // + eta + half_eta).  Project it down to
                        // fc_with_eta_.ctx (drop __amf_*).
                        algebra::Mfrac lt_coef_in_fc;
                        try {
                            lt_coef_in_fc = project_mfrac_by_name(
                                lt.coeff, fc_with_eta_->ctx);
                        } catch (const std::exception& e) {
                            throw std::runtime_error(
                                "AMFSystem::build_boundary: lt.coef projection: "
                                + std::string(e.what()));
                        }

                        for (const auto& rt : rule) {
                            auto mit = sub_master_idx.find(
                                jintegral_key(rt.integ));
                            if (mit == sub_master_idx.end()) continue;
                            // rt.coef is on red.red_ctx (= sub_fc.ctx + 'd').
                            // Project lt_coef_in_fc to that ctx.
                            algebra::Mfrac lt_in_sub;
                            try {
                                lt_in_sub = project_mfrac_by_name(
                                    lt_coef_in_fc, red.red_ctx.ctx);
                            } catch (const std::exception& e) {
                                throw std::runtime_error(
                                    "AMFSystem::build_boundary: lt-to-sub-red projection: "
                                    + std::string(e.what()));
                            }
                            algebra::Mfrac product = lt_in_sub.clone();
                            product *= rt.coef;
                            BoundaryTerm bt;
                            bt.coef = std::move(product);
                            bt.sub_master_idx = (long)mit->second;
                            bt.sub_master_key = jintegral_key(rt.integ);
                            rb.integrand_terms[i][k].push_back(std::move(bt));
                        }
                    }
                }
            }

            regions_.push_back(std::move(rb));
        }
    }
}

// ===========================================================================
//  Solve
// ===========================================================================

// Compute the SingleMass prefactor numerically:
//
//   prefactor = (-1)^(-n) * Gamma(2 - eps - m) * Gamma(-2 + eps + m + n) /
//               (Gamma(2 - eps) * Gamma(n))
//
// where m = m0 + m_eps_coef * eps.
namespace {
void compute_single_mass_prefactor(const AMFSystem::GammaParams& gp,
                                      const acb_t eps,
                                      acb_t out,
                                      long prec) {
    acb_t m, two_minus_eps, gamma1, gamma2, gamma3, gamma4;
    acb_t arg, sign;
    acb_init(m); acb_init(two_minus_eps);
    acb_init(gamma1); acb_init(gamma2);
    acb_init(gamma3); acb_init(gamma4);
    acb_init(arg); acb_init(sign);

    // m = m0 + m_eps_coef * eps
    acb_set_si(m, gp.m0);
    if (gp.m_eps_coef != 0) {
        acb_t tmp; acb_init(tmp);
        acb_set_si(tmp, gp.m_eps_coef);
        acb_mul(tmp, tmp, eps, prec);
        acb_add(m, m, tmp, prec);
        acb_clear(tmp);
    }

    // Gamma(2 - eps - m)
    acb_set_si(arg, 2);
    acb_sub(arg, arg, eps, prec);
    acb_sub(arg, arg, m, prec);
    acb_gamma(gamma1, arg, prec);

    // Gamma(-2 + eps + m + n)
    acb_set_si(arg, -2 + (long)gp.n);
    acb_add(arg, arg, eps, prec);
    acb_add(arg, arg, m, prec);
    acb_gamma(gamma2, arg, prec);

    // Gamma(2 - eps)
    acb_set_si(arg, 2);
    acb_sub(arg, arg, eps, prec);
    acb_gamma(gamma3, arg, prec);

    // Gamma(n)
    acb_set_si(arg, gp.n);
    acb_gamma(gamma4, arg, prec);

    // sign = (-1)^(-n) = (-1)^n
    acb_one(sign);
    if (gp.n % 2 != 0) acb_neg(sign, sign);

    // out = sign * gamma1 * gamma2 / (gamma3 * gamma4)
    acb_mul(out, sign, gamma1, prec);
    acb_mul(out, out, gamma2, prec);
    acb_div(out, out, gamma3, prec);
    acb_div(out, out, gamma4, prec);

    acb_clear(m); acb_clear(two_minus_eps);
    acb_clear(gamma1); acb_clear(gamma2);
    acb_clear(gamma3); acb_clear(gamma4);
    acb_clear(arg); acb_clear(sign);
}

void compute_cutkosky_prefactor(long phase_loop_num,
                                  const acb_t eps,
                                  acb_t out,
                                  long prec) {
    if (phase_loop_num < 0) {
        throw std::runtime_error(
            "compute_cutkosky_prefactor: phase loop number is negative");
    }

    acb_t pi, two_pi, exp1, exp2, four, pi_pow, two_pi_pow, base;
    acb_init(pi); acb_init(two_pi);
    acb_init(exp1); acb_init(exp2); acb_init(four);
    acb_init(pi_pow); acb_init(two_pi_pow);
    acb_init(base);

    acb_const_pi(pi, prec);
    acb_mul_2exp_si(two_pi, pi, 1);

    // Pi^(2 - eps)
    acb_set_si(exp1, 2);
    acb_sub(exp1, exp1, eps, prec);
    acb_pow(pi_pow, pi, exp1, prec);

    // (2 Pi)^(2 eps - 4)
    acb_mul_2exp_si(exp2, eps, 1);
    acb_set_si(four, 4);
    acb_sub(exp2, exp2, four, prec);
    acb_pow(two_pi_pow, two_pi, exp2, prec);

    acb_mul(base, pi_pow, two_pi_pow, prec);
    acb_pow_ui(out, base, static_cast<unsigned long>(phase_loop_num), prec);
    acb_mul_2exp_si(out, out, 1);
    if ((phase_loop_num + 1) % 2 != 0) acb_neg(out, out);

    acb_clear(pi); acb_clear(two_pi);
    acb_clear(exp1); acb_clear(exp2); acb_clear(four);
    acb_clear(pi_pow); acb_clear(two_pi_pow);
    acb_clear(base);
}
}  // namespace

void AMFSystem::solve(const std::vector<numeric::AcbValue>& epslist) {
    solutions_.clear();
    solutions_.resize(epslist.size());
    std::map<std::string, ibp::ReduceResult> ending_reduction_cache;

    if (is_ending_) {
        for (std::size_t k = 0; k < epslist.size(); ++k) {
            std::vector<numeric::AcbValue> values;
            values.reserve(preferred_.size());
            auto numeric_q =
                build_numeric_q_for_eps(opts_.bb, epslist[k], numeric::working_prec_bits());
            std::map<std::string, numeric::AcbValue> ending_value_cache;
            std::set<std::string> ending_active_keys;
            for (const auto& j : preferred_) {
                values.push_back(solve_ending_master_value(
                    *fc_, j, system_id_, opts_, k, epslist[k], numeric_q,
                    ending_reduction_cache, ending_value_cache,
                    ending_active_keys, numeric::working_prec_bits()));
            }
            solutions_[k].master_values = std::move(values);
        }
        // NOTE: do NOT return here -- we still need to apply SingleMass
        // prefactor below if sm_cfg_ is set (an ending SingleMass system
        // computes preferred = jpre, multiplies by Γ-prefactor, stores
        // result in global_values for the parent's boundary integrals).
    } else {
        // Non-ending: solve sub-systems first, then per-ε solve.
        for (auto& sub : subsystems_) {
            sub->solve(epslist);
        }
        for (std::size_t k = 0; k < epslist.size(); ++k) {
            solutions_[k] = solve_one_eps(k, epslist[k]);
        }
    }

    // Apply SingleMass prefactor to produce global_values.
    if (sm_cfg_.has_value()) {
        long prec = numeric::working_prec_bits();
        const auto& cfg = *sm_cfg_;
        for (std::size_t k = 0; k < epslist.size(); ++k) {
            auto& sol = solutions_[k];
            sol.global_preferred = cfg.original_preferred;
            sol.global_values.clear();
            sol.global_values.reserve(cfg.original_preferred.size());
            for (std::size_t i = 0; i < cfg.original_preferred.size(); ++i) {
                long jpre_i = cfg.jpre_idx[i];
                if (jpre_i < 0 || (std::size_t)jpre_i >= sol.master_values.size()) {
                    throw std::runtime_error(
                        "AMFSystem::solve: jpre_idx out of range");
                }
                numeric::AcbValue prefactor;
                compute_single_mass_prefactor(cfg.gamma_params[i],
                                                 epslist[k].raw(),
                                                 prefactor.raw(), prec);
                numeric::AcbValue gv;
                acb_mul(gv.raw(), prefactor.raw(),
                          sol.master_values[(std::size_t)jpre_i].raw(), prec);
                sol.global_values.push_back(std::move(gv));
            }
            AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
                std::cerr << "[single_mass_solve] family=" << fc_->family
                          << " eps_index=" << k
                          << " global_values.size=" << sol.global_values.size()
                          << std::endl;
                for (std::size_t i = 0; i < sol.global_values.size(); ++i) {
                    const auto& jp = sol.global_preferred[i];
                    std::cerr << "  J[" << jp.family();
                    for (long idx : jp.indices()) std::cerr << "," << idx;
                    std::cerr << "] = "
                              << sol.global_values[i].to_string(20) << std::endl;
                }
            }
        }
    }

    if (ck_cfg_.has_value()) {
        long prec = numeric::working_prec_bits();
        const auto& cfg = *ck_cfg_;
        for (std::size_t k = 0; k < epslist.size(); ++k) {
            auto& sol = solutions_[k];
            if (sol.master_values.size() < cfg.original_preferred.size()) {
                throw std::runtime_error(
                    "AMFSystem::solve: Cutkosky master value count mismatch");
            }

            numeric::AcbValue prefactor;
            compute_cutkosky_prefactor(cfg.phase_loop_num,
                                         epslist[k].raw(),
                                         prefactor.raw(), prec);

            sol.global_preferred = cfg.original_preferred;
            sol.global_values.clear();
            sol.global_values.reserve(cfg.original_preferred.size());
            for (std::size_t i = 0; i < cfg.original_preferred.size(); ++i) {
                numeric::AcbValue im_value;
                arb_set(acb_realref(im_value.raw()),
                        acb_imagref(sol.master_values[i].raw()));
                arb_zero(acb_imagref(im_value.raw()));

                numeric::AcbValue gv;
                acb_mul(gv.raw(), prefactor.raw(), im_value.raw(), prec);
                sol.global_values.push_back(std::move(gv));
            }
        }
    }
}

// ---------------------------------------------------------------------------
//  solve_one_eps
// ---------------------------------------------------------------------------

AMFSystemSolution
AMFSystem::solve_one_eps(std::size_t eps_index, const numeric::AcbValue& eps) {
    long prec = numeric::working_prec_bits();
    AMFSystemSolution out;

    // ---- Build numeric_q (ε + family numeric values) ----
    auto numeric_q = build_numeric_q_for_eps(opts_.bb, eps, prec);

    // ---- Step 1.  Build numeric::RationalMatrix(η) by evaluating diffeq. ----
    if (diffeq_result_.diffeq.size() != 1) {
        throw std::runtime_error(
            "AMFSystem::solve_one_eps: expected exactly one diffeq variable");
    }
    const auto& dmat = diffeq_result_.diffeq[0];
    std::size_t Nm = diffeq_result_.sortedmasters.size();
    long eta_in_red = diffeq_result_.red_ctx.ctx->var_index("eta");
    if (eta_in_red < 0) {
        throw std::runtime_error(
            "AMFSystem::solve_one_eps: red_ctx is missing 'eta'");
    }
    std::set<std::string> keep_eta = {"eta"};

    numeric::RationalMatrix de((std::size_t)Nm, (std::size_t)Nm);
    for (std::size_t i = 0; i < Nm; ++i) {
        for (std::size_t j = 0; j < Nm; ++j) {
            try {
                algebra::Mfrac sub = mfrac_substitute(dmat[i][j], numeric_q, keep_eta);
                de(i, j) = mfrac_to_eta_rational(sub, eta_in_red);
            } catch (const std::exception& e) {
                throw std::runtime_error(
                    "AMFSystem::solve_one_eps: diffeq("
                    + std::to_string(i) + "," + std::to_string(j) + "): "
                    + e.what());
            }
        }
    }

    // ---- Step 2.  Build ode::BoundarySpec for each sortedmaster. ----
    //
    //   For each region:
    //     For each sortedmaster i:
    //       evaluate mu(ε), then for each integrand-order k:
    //         evaluate Σ_term (coef(ε) * sub_solution_at_eps),
    //         append (μ - k, value) to BC[i].
    std::vector<ode::BoundarySpec> bc_sorted(Nm);
    struct CombinedSubValues {
        std::vector<numeric::AcbValue> values;
        std::map<std::string, std::size_t> key_to_index;
    };
    auto combined_sub_values =
        [&](const std::vector<long>& sub_indices) -> CombinedSubValues {
            CombinedSubValues out_vals;
            if (sub_indices.empty()) return out_vals;

            const auto& first_sys = *subsystems_[(std::size_t)sub_indices.front()];
            if (eps_index >= first_sys.solutions().size()) {
                throw std::runtime_error(
                    "AMFSystem::solve_one_eps: sub-system not solved at this ε");
            }
            const auto& first_sol = first_sys.solutions()[eps_index];

            // SingleMass multi-root children expose the physically relevant
            // boundary-family masters via global_values/global_preferred, and
            // the final answer is the product over roots.  This mirrors
            // Mathematica's AMFSystemsSolution / AMFSystemCombineSolution.
            if (!first_sol.global_values.empty()) {
                std::vector<numeric::AcbValue> combined;
                combined.reserve(first_sol.global_values.size());
                for (std::size_t i = 0; i < first_sol.global_values.size(); ++i) {
                    numeric::AcbValue one;
                    acb_one(one.raw());
                    combined.push_back(std::move(one));
                }
                for (long sub_idx : sub_indices) {
                    const auto& sub_sys = *subsystems_[(std::size_t)sub_idx];
                    if (eps_index >= sub_sys.solutions().size()) {
                        throw std::runtime_error(
                            "AMFSystem::solve_one_eps: sub-system not solved at this ε");
                    }
                    const auto& sub_sol = sub_sys.solutions()[eps_index];
                    if (sub_sol.global_values.size() != combined.size()) {
                        throw std::runtime_error(
                            "AMFSystem::solve_one_eps: multi-root sub-system "
                            "global_values size mismatch");
                    }
                    for (std::size_t i = 0; i < combined.size(); ++i) {
                        numeric::AcbValue tmp;
                        acb_mul(tmp.raw(), combined[i].raw(),
                                sub_sol.global_values[i].raw(), prec);
                        combined[i] = std::move(tmp);
                    }
                }
                for (std::size_t i = 0; i < combined.size(); ++i) {
                    combined[i] = midpoint_only(combined[i]);
                    if (i < first_sol.global_preferred.size()) {
                        out_vals.key_to_index[jintegral_key(first_sol.global_preferred[i])] = i;
                    }
                }
                out_vals.values = std::move(combined);
                return out_vals;
            }

            if (sub_indices.size() != 1) {
                throw std::runtime_error(
                    "AMFSystem::solve_one_eps: multi-root sub-system without "
                    "global_values is unsupported");
            }

            out_vals.values.reserve(first_sol.master_values.size());
            for (const auto& v : first_sol.master_values) {
                out_vals.values.push_back(midpoint_only(v));
            }
            for (std::size_t i = 0; i < first_sys.preferred().size(); ++i) {
                out_vals.key_to_index[jintegral_key(first_sys.preferred()[i])] = i;
            }
            return out_vals;
        };

    for (const auto& rb : regions_) {
        if (rb.sub_system_indices.empty()) continue;
        CombinedSubValues sub_vals = combined_sub_values(rb.sub_system_indices);
        AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
            std::cerr << "[solve_one_eps] family=" << fc_->family
                      << " sub_roots=[";
            for (std::size_t ridx = 0; ridx < rb.sub_system_indices.size(); ++ridx) {
                if (ridx) std::cerr << ",";
                const auto& sub_sys = *subsystems_[(std::size_t)rb.sub_system_indices[ridx]];
                std::cerr << sub_sys.family().family;
            }
            std::cerr << "] combined_sub_vals.size=" << sub_vals.values.size()
                      << std::endl;
            for (std::size_t z = 0; z < sub_vals.values.size(); ++z) {
                std::cerr << "    sub_vals[" << z << "] = "
                          << sub_vals.values[z].to_string(60) << std::endl;
            }
            std::cerr << "  rb.powers.size=" << rb.powers.size()
                      << " rb.integrand_terms.size=" << rb.integrand_terms.size()
                      << std::endl;
            for (std::size_t i = 0; i < Nm; ++i) {
                std::cerr << "  i=" << i << " integrand_terms[i].size="
                          << rb.integrand_terms[i].size() << std::endl;
                for (std::size_t k = 0; k < rb.integrand_terms[i].size(); ++k) {
                    std::cerr << "    k=" << k << " terms.size="
                              << rb.integrand_terms[i][k].size() << std::endl;
                    for (const auto& bt : rb.integrand_terms[i][k]) {
                        std::cerr << "      sub_master_idx=" << bt.sub_master_idx
                                  << " coef=" << bt.coef.to_string() << std::endl;
                    }
                }
            }
        }

        for (std::size_t i = 0; i < Nm; ++i) {
            // mu(ε)
            numeric::AcbValue mu_acb;
            try {
                evaluate_const_mfrac(rb.powers[i], numeric_q, prec, mu_acb);
            } catch (const std::exception& e) {
                throw std::runtime_error(
                    "AMFSystem::solve_one_eps: mu eval: " + std::string(e.what()));
            }
            for (std::size_t k = 0; k < rb.integrand_terms[i].size(); ++k) {
                const auto& terms = rb.integrand_terms[i][k];
                if (terms.empty()) continue;
                numeric::AcbValue value;
                acb_zero(value.raw());
                for (const auto& bt : terms) {
                    numeric::AcbValue coef_acb;
                    try {
                        evaluate_const_mfrac(bt.coef, numeric_q, prec, coef_acb);
                    } catch (const std::exception& e) {
                        throw std::runtime_error(
                            "AMFSystem::solve_one_eps: coef eval: "
                            + std::string(e.what()));
                    }
                    std::size_t sub_idx = (std::size_t)bt.sub_master_idx;
                    auto key_it = sub_vals.key_to_index.find(bt.sub_master_key);
                    if (key_it != sub_vals.key_to_index.end()) {
                        sub_idx = key_it->second;
                    }
                    if (sub_idx >= sub_vals.values.size()) {
                        throw std::runtime_error(
                            "AMFSystem::solve_one_eps: sub-master lookup out of range");
                    }
                    numeric::AcbValue prod;
                    acb_mul(prod.raw(), coef_acb.raw(),
                             sub_vals.values[sub_idx].raw(), prec);
                    acb_add(value.raw(), value.raw(), prod.raw(), prec);
                }
                ode::BoundaryEntry be;
                if (k == 0) {
                    acb_set(be.mu.raw(), mu_acb.raw());
                } else {
                    numeric::AcbValue offset;
                    acb_set_si(offset.raw(), -(long)k);
                    acb_add(be.mu.raw(), mu_acb.raw(), offset.raw(), prec);
                }
                acb_set(be.value.raw(), value.raw());
                bc_sorted[i].push_back(std::move(be));
            }
        }
    }

    // ---- Step 2.5. Append pattern-zero BC entries. ----
    //
    // Mirrors MMA's AMFSystemSolution final line:
    //   bc = MapThread[Join[#1,#2]&,
    //                  {bc, Transpose[Thread /@ Thread[(pattern/.epsrule) -> 0]]}];
    //
    // For every pattern group and every sortedmaster i, add a trivial
    // (mu = pattern[g][i] | eps, value = 0) boundary entry.  Without this,
    // a sub-system whose regions all have border=-1 ends up with an empty
    // bc_sorted and amflow produces 0 -- see Stage B notes.
    for (const auto& pg : bc_pattern_) {
        for (std::size_t i = 0; i < Nm && i < pg.size(); ++i) {
            numeric::AcbValue mu_acb;
            try {
                evaluate_const_mfrac(pg[i], numeric_q, prec, mu_acb);
            } catch (const std::exception& e) {
                throw std::runtime_error(
                    "AMFSystem::solve_one_eps: pattern mu eval: "
                    + std::string(e.what()));
            }
            ode::BoundaryEntry be;
            acb_set(be.mu.raw(), mu_acb.raw());
            acb_zero(be.value.raw());
            bc_sorted[i].push_back(std::move(be));
        }
    }

    // ---- Step 3. amflow ----
    AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
        std::cerr << "[solve_one_eps] family=" << fc_->family
                  << " sortedmasters=" << diffeq_result_.sortedmasters.size()
                  << " boundary_entries=" << regions_.size() << std::endl;
        for (std::size_t i = 0; i < diffeq_result_.sortedmasters.size(); ++i) {
            const auto& m = diffeq_result_.sortedmasters[i];
            std::cerr << "  smaster[" << i << "] = J[" << m.family();
            for (auto x : m.indices()) std::cerr << "," << x;
            std::cerr << "]" << std::endl;
        }
        for (std::size_t i = 0; i < bc_sorted.size(); ++i) {
            std::cerr << "  bc_sorted[" << i << "]: " << bc_sorted[i].size()
                      << " entries" << std::endl;
            for (const auto& be : bc_sorted[i]) {
                std::cerr << "    mu=" << be.mu.to_string(60)
                          << "  value=" << be.value.to_string(60) << std::endl;
            }
        }
    }
    AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
        std::cerr << "[solve_one_eps] family=" << fc_->family
                  << " de matrix " << de.rows() << "x" << de.cols() << ":\n";
        for (std::size_t i = 0; i < de.rows(); ++i) {
            for (std::size_t j = 0; j < de.cols(); ++j) {
                std::cerr << "    de[" << i << "," << j << "] = "
                          << de(i, j).to_string() << std::endl;
            }
        }
    }
    std::vector<numeric::AcbValue> sol_acb;
    try {
        sol_acb = ode::amflow(de, bc_sorted, prec);
    } catch (const std::exception& e) {
        std::ostringstream msg;
        msg << "AMFSystem::solve_one_eps: amflow failed"
            << "; system_id=" << system_id_
            << "; family=" << fc_->family
            << "; eps_index=" << eps_index
            << "; masters=" << Nm
            << ": " << e.what();
        throw std::runtime_error(msg.str());
    }
    AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
        std::cerr << "[solve_one_eps] family=" << fc_->family
                  << " sol_acb (sorted) values:" << std::endl;
        for (std::size_t i = 0; i < sol_acb.size(); ++i) {
            std::cerr << "    sol[" << i << "] = "
                      << sol_acb[i].to_string(60) << std::endl;
        }
    }

    out.master_values.reserve(preferred_.size());
    for (std::size_t i = 0; i < preferred_.size(); ++i) {
        long si = pref_to_sorted_[i];
        out.master_values.push_back(sol_acb[(std::size_t)si].clone());
    }
    AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
        std::cerr << "[solve_one_eps] family=" << fc_->family
                  << " preferred_values=" << out.master_values.size()
                  << std::endl;
        for (std::size_t i = 0; i < preferred_.size(); ++i) {
            const auto& pref = preferred_[i];
            std::cerr << "  preferred[" << i << "] = J[" << pref.family();
            for (long idx : pref.indices()) std::cerr << "," << idx;
            std::cerr << "] value="
                      << out.master_values[i].to_string(60) << std::endl;
        }
    }
    return out;
}

// ===========================================================================
//  Top-level orchestrators
// ===========================================================================

// ===========================================================================
//  SingleMass setup
// ===========================================================================
//
//   Mirrors AMFSystemSetupMaster[..., "SingleMass"] (.m line 1039-1061).
//
//   The single_mass_setup_master function:
//     1. Apply legs->0 to ReducedPropagator[topposi].
//     2. Call FactorizeFamily(props, [FromJ[i][topposi] for i in preferred]).
//     3. For each component (loop, prop, mas):
//          a. Find `drop` = position where mass = -1 (numerically).
//          b. droploop = the loop momentum involved in prop[drop].
//          c. Build sub-qft::FamilyConfig with droploop as leg, replacement
//             {droploop² = -1}, propagators = ToCompleteExplicit[Drop[prop, drop]].
//          d. jpre[i] = ToJ[PadRight[Drop[mas[i], {drop}], n_orig_props]]
//             on the new sub-fc family.
//          e. prefactor[i] = analytic Gamma function expression.
//          f. Recursively setup the sub-system on jpre.
//     4. Return list of AMFSystems (one per component).

namespace {
std::vector<std::unique_ptr<AMFSystem>>
single_mass_setup_master(const qft::FamilyConfig& fc,
                           const std::vector<qft::JIntegral>& preferred,
                           const AMFSystemOptions& opts);

// Find Position[mass_list, -1] (under numeric sub).  Returns -1 if not found.
long find_mass_minus_one(const std::vector<algebra::Mfrac>& masses,
                            const std::map<std::string, FmpqHolder>& numeric_q) {
    for (std::size_t k = 0; k < masses.size(); ++k) {
        algebra::Mfrac sub = apply_numeric(masses[k], numeric_q);
        // mass=-1 means num = -1 * den, i.e. (num + den).is_zero().
        algebra::Mpoly sum = sub.numerator() + sub.denominator();
        if (sum.is_zero()) return (long)k;
    }
    return -1;
}

// Convert ToSquareResult.masses (algebra::Mfrac on fc.ctx) to a vector and
// return.  We need this for FactorizeFamily output's prop list.
std::vector<algebra::Mfrac>
to_square_masses(const qft::FamilyConfig& fc, const std::vector<algebra::Mpoly>& props) {
    auto sq = to_square_all(fc, props);
    return std::move(sq.masses);
}

// Find the loop momentum involved in `prop` (the chosen squared-loop
// variable index).  Returns -1 if no loop appears.
long find_loop_in_prop(const qft::FamilyConfig& fc, const algebra::Mpoly& prop) {
    long n_loop = (long)fc.n_loops();
    for (long j = 0; j < n_loop; ++j) {
        if (!prop.coeff_of(j, 2).is_zero()) return j;
    }
    return -1;
}
}  // namespace

// SingleMass scheme implementation.
//
// For each FactorizeFamily component, build a new sub-qft::FamilyConfig
// where one loop momentum is "promoted to a leg" (with on-shell
// constraint loop² = -1) and one propagator is dropped.  The new
// system is on (n_loops - 1, n_props - 1) form.  An analytic
// Gamma-function prefactor links the original master to the new
// (sub-)master.
namespace {

std::vector<std::unique_ptr<AMFSystem>>
single_mass_setup_master(const qft::FamilyConfig& fc,
                           const std::vector<qft::JIntegral>& preferred,
                           const AMFSystemOptions& opts) {
    auto top = qft::get_top_position(preferred);
    if (top.empty()) return {};

    // Step 1: extract top-sector propagators with legs->0 substitution.
    std::vector<algebra::Mpoly> top_props;
    top_props.reserve(top.size());
    long n_loop_orig = (long)fc.n_loops();
    long n_red_orig  = (long)fc.n_red_legs();
    fmpz_t zero;
    fmpz_init(zero);
    for (std::size_t pos : top) {
        algebra::Mpoly p = fc.propagators_after_conservation[pos].clone();
        for (long k = 0; k < n_red_orig; ++k) {
            algebra::Mpoly out(fc.ctx);
            fmpz_mpoly_evaluate_one_fmpz(out.raw(), p.raw(),
                                          n_loop_orig + k, zero,
                                          fc.ctx->raw());
            p = std::move(out);
        }
        top_props.push_back(std::move(p));
    }
    fmpz_clear(zero);

    // Step 2: build patts (per preferred, indices at top positions).
    std::vector<std::vector<long>> patts;
    patts.reserve(preferred.size());
    for (const auto& j : preferred) {
        std::vector<long> p;
        p.reserve(top.size());
        for (std::size_t pos : top) {
            if (pos < (std::size_t)j.indices().size()) {
                p.push_back(j.indices()[pos]);
            } else {
                p.push_back(0);
            }
        }
        patts.push_back(std::move(p));
    }

    // Step 3: factorize.
    auto faminfo = factorize_family(fc, top_props, patts);
    if (faminfo.empty()) {
        throw std::runtime_error(
            "single_mass_setup_master: FactorizeFamily returned no components");
    }

    auto numeric_q = build_numeric_q(opts.bb);

    std::vector<std::unique_ptr<AMFSystem>> out;

    for (auto& fcomp : faminfo) {
        // FactorizeFamily output: fcomp has loops (subset of fc.loops),
        // propagators (legs-zeroed), patterns (per master, projected
        // to component positions), original_positions (where in
        // top_props each prop comes from).
        //
        // For SingleMass, we operate on the COMPONENT's propagators and
        // loops.  Build a temporary qft::FamilyConfig for this component.
        // Note: legs are empty (legs were zeroed in FactorizeFamily input).
        std::vector<std::string> comp_prop_strs;
        comp_prop_strs.reserve(fcomp.propagators.size());
        for (const auto& pp : fcomp.propagators) {
            comp_prop_strs.push_back(pp.to_string());
        }
        qft::FamilyConfig comp_fc = qft::FamilyConfig::build(
            fc.family + "_smcomp", fcomp.loops, /*legs=*/{},
            /*conservation=*/{},
            /*replacement=*/{},
            comp_prop_strs);

        // Re-parse propagators on comp_fc.ctx (string round-trip).
        std::vector<algebra::Mpoly> comp_props;
        comp_props.reserve(fcomp.propagators.size());
        for (const auto& s : comp_prop_strs) {
            comp_props.push_back(algebra::Mpoly::from_string(comp_fc.ctx, s));
        }

        // Find drop = position where mass = -1 numerically.
        auto comp_masses = to_square_masses(comp_fc, comp_props);
        long drop = find_mass_minus_one(comp_masses, numeric_q);
        if (drop < 0) {
            throw std::runtime_error(
                "single_mass_setup_master: no mass=-1 propagator after "
                "numeric substitution; check that AMFlowInfo[\"Numeric\"] "
                "sets exactly one propagator's mass to 1");
        }

        // droploop = loop variable in dropped propagator (use comp_props on comp_fc.ctx).
        long droploop = find_loop_in_prop(comp_fc, comp_props[(std::size_t)drop]);
        if (droploop < 0) {
            throw std::runtime_error(
                "single_mass_setup_master: cannot find loop in dropped propagator");
        }
        std::string droploop_name = fcomp.loops[(std::size_t)droploop];

        // Build sub-fc: loops = component loops minus droploop;
        //              legs = {droploop};
        //              replacement = {droploop² -> -1};
        //              propagators = Drop[fcomp.propagators, drop],
        //              then ToCompleteExplicit (added implicitly when fc is
        //              built? No -- we do it manually).
        std::vector<std::string> sub_loops;
        for (const auto& l : fcomp.loops) {
            if (l != droploop_name) sub_loops.push_back(l);
        }
        std::vector<std::string> sub_legs = {droploop_name};
        std::vector<std::pair<std::string, std::string>> sub_replacement = {
            {droploop_name + "^2", "-1"}
        };

        // Build sub propagator list: drop the `drop`-th propagator,
        // keep the others.  Use comp_props (already on comp_fc.ctx).
        std::vector<algebra::Mpoly> dropped_props;
        for (std::size_t i = 0; i < comp_props.size(); ++i) {
            if ((long)i == drop) continue;
            dropped_props.push_back(comp_props[i].clone());
        }

        // Apply ToCompleteExplicit to dropped_props.  We need a
        // qft::FamilyConfig with the new loop/leg structure to do this.
        qft::FamilyConfig pre_sub_fc = qft::FamilyConfig::build(
            fc.family + "_sm" + std::to_string(out.size()),
            sub_loops, sub_legs,
            /*conservation=*/{}, sub_replacement,
            /*propagators=*/{},
            /*cut=*/{}, /*prescription=*/{});
        // ToCompleteExplicit expects propagators on this fc's ctx.
        // dropped_props live on comp_fc.ctx; translate by variable name
        // into pre_sub_fc.ctx after promoting droploop to an external leg.
        //
        // Build a name->index map for pre_sub_fc.ctx.
        std::map<std::string, long> name_to_dst;
        for (long i = 0; i < pre_sub_fc.ctx->n_vars(); ++i) {
            name_to_dst[pre_sub_fc.ctx->var_name(i)] = i;
        }
        std::vector<algebra::Mpoly> dropped_in_sub;
        dropped_in_sub.reserve(dropped_props.size());
        for (const auto& dp : dropped_props) {
            auto src_ctx = dp.ctx();
            // Walk monomials and re-encode on pre_sub_fc.ctx.
            algebra::Mpoly out_p(pre_sub_fc.ctx);
            long len = fmpz_mpoly_length(dp.raw(), src_ctx->raw());
            std::vector<unsigned long> sexp((std::size_t)src_ctx->n_vars());
            std::vector<unsigned long> dexp(
                (std::size_t)pre_sub_fc.ctx->n_vars(), 0);
            fmpz_t coeff;
            fmpz_init(coeff);
            for (long t = 0; t < len; ++t) {
                fmpz_mpoly_get_term_exp_ui(sexp.data(), dp.raw(),
                                            t, src_ctx->raw());
                fmpz_mpoly_get_term_coeff_fmpz(coeff, dp.raw(),
                                                t, src_ctx->raw());
                std::fill(dexp.begin(), dexp.end(), 0);
                for (long sv = 0; sv < src_ctx->n_vars(); ++sv) {
                    if (sexp[(std::size_t)sv] == 0) continue;
                    const std::string& name = src_ctx->var_name(sv);
                    auto it = name_to_dst.find(name);
                    if (it == name_to_dst.end()) {
                        fmpz_clear(coeff);
                        throw std::runtime_error(
                            "single_mass_setup_master: variable '" + name
                            + "' not in SingleMass sub-family context");
                    }
                    dexp[(std::size_t)it->second] = sexp[(std::size_t)sv];
                }
                fmpz_mpoly_set_coeff_fmpz_ui(out_p.raw(), coeff,
                                              dexp.data(), pre_sub_fc.ctx->raw());
            }
            fmpz_clear(coeff);
            dropped_in_sub.push_back(std::move(out_p));
        }

        // ToCompleteExplicit on dropped_in_sub.
        auto completede = qft::to_complete_explicit(pre_sub_fc, dropped_in_sub);
        std::vector<std::string> sub_prop_strs;
        sub_prop_strs.reserve(completede.size());
        for (const auto& pp : completede) {
            sub_prop_strs.push_back(pp.to_string());
        }
        // Rebuild final sub-fc with the completed propagator list.
        qft::FamilyConfig sub_fc = qft::FamilyConfig::build(
            pre_sub_fc.family,
            sub_loops, sub_legs,
            /*conservation=*/{}, sub_replacement,
            sub_prop_strs);

        // jpre[i] = ToJ[PadRight[Drop[mas[i], {drop}], n_orig_props]]
        // mas[i] = fcomp.patterns[i] (length = top.size())
        // We drop position `drop` from mas[i], then PadRight to
        // n_orig_props (i.e. pad with zeros on the right).
        //
        // Wait -- Mathematica's PadRight is to Length[ReducedPropagator]
        // = original family size.  But in our sub-system, the J integral
        // lives on sub_fc which has Length(sub_prop_strs) propagators.
        // We use sub_prop_strs.size() as the index length.
        std::vector<qft::JIntegral> jpre;
        jpre.reserve(fcomp.patterns.size());
        for (const auto& patt : fcomp.patterns) {
            std::vector<long> idx;
            idx.reserve(sub_prop_strs.size());
            for (std::size_t k = 0; k < patt.size(); ++k) {
                if ((long)k == drop) continue;
                idx.push_back(patt[k]);
            }
            // Pad with zeros up to sub_prop_strs.size().
            while (idx.size() < sub_prop_strs.size()) idx.push_back(0);
            jpre.emplace_back(sub_fc.family, idx);
        }

        // Compute prefactor[i] symbolically as algebra::Mfrac on a context that
        // includes 'd' / 'eps' (for runtime substitution).
        //
        // prefactor = (-1)^(-n) * Gamma[2-eps-m] * Gamma[-2+eps+m+n] /
        //              (Gamma[2-eps] * Gamma[n])
        //
        // where n = patt[drop], m = Total[patt] - patt[drop] - (loop_num-1)*(2-eps)
        //
        // We can't represent Gamma symbolically with algebra::Mfrac (which is
        // polynomial-rational only).  Instead we evaluate prefactor at
        // runtime in solve_one_eps using acb_gamma.  We store the
        // *integer* parameters (n, total, loop_num) per master and
        // evaluate symbolically when ε is given.
        //
        // For now, we represent the prefactor as a symbolic placeholder
        // using a small data structure attached to AMFSystem via
        // SingleMassConfig.  The fields n_param, m_int, total_int,
        // loop_num let us compute the Gamma expression at runtime.
        //
        // algebra::Mfrac is only used for *coefficient* multiplication; for
        // Gamma we extend the representation.

        AMFSystem::SingleMassConfig sm_cfg;
        sm_cfg.original_preferred = preferred;
        sm_cfg.jpre_idx.reserve(preferred.size());
        sm_cfg.gamma_params.reserve(preferred.size());

        for (std::size_t i = 0; i < preferred.size(); ++i) {
            const auto& patt = fcomp.patterns[i];
            long n_param = patt[(std::size_t)drop];
            long total = 0;
            for (long v : patt) total += v;
            long loop_num = (long)fcomp.loops.size();
            // m = total - n - (loop_num-1) * (2 - eps)
            //   = (total - n - 2*(loop_num-1)) + (loop_num-1) * eps
            //   = m0 + m_eps_coef * eps
            AMFSystem::GammaParams gp;
            gp.n            = n_param;
            gp.m0           = total - n_param - 2 * (loop_num - 1);
            gp.m_eps_coef   = loop_num - 1;
            gp.loop_num     = loop_num;
            sm_cfg.gamma_params.push_back(gp);

            // Find which jpre index this preferred maps to.  Since
            // jpre is built per-preferred from fcomp.patterns,
            // jpre[i] (before dedup) corresponds to preferred[i].
            sm_cfg.jpre_idx.push_back((long)i);     // placeholder; fixed below
        }

        // DeleteDuplicates[jpre] and remap jpre_idx accordingly.
        std::map<std::string, long> jpre_dedup;
        std::vector<qft::JIntegral> jpre_unique;
        std::vector<long> orig_to_unique(jpre.size(), -1);
        for (std::size_t i = 0; i < jpre.size(); ++i) {
            std::string key = jintegral_key(jpre[i]);
            auto it = jpre_dedup.find(key);
            if (it == jpre_dedup.end()) {
                jpre_dedup[key] = (long)jpre_unique.size();
                jpre_unique.push_back(jpre[i]);
                orig_to_unique[i] = (long)jpre_unique.size() - 1;
            } else {
                orig_to_unique[i] = it->second;
            }
        }
        // Remap sm_cfg.jpre_idx to point into jpre_unique.
        for (std::size_t i = 0; i < sm_cfg.jpre_idx.size(); ++i) {
            sm_cfg.jpre_idx[i] = orig_to_unique[(std::size_t)sm_cfg.jpre_idx[i]];
        }

        // Compute etac on jpre_unique and DIRECTLY construct the
        // sub-system (mirrors Mathematica's AMFSystemSetup call at .m
        // line 1057 -- NOT AMFSystemSetupMaster, so the sub-system
        // does not re-enter scheme selection).
        AMFSystemOptions sub_opts = opts;
        sub_opts.current_depth = opts.current_depth + 1;
        // Sub-fc has empty Numeric (Mathematica .m line 1052 clears it).
        // Drop all numeric values; ε is supplied at solve-time.
        sub_opts.bb.numeric_values.clear();

        // Compute etac for the sub-system directly:
        //   etac = AMFEtaC[AMFPosition[GetTopPosition[jpre_unique], $qft::AMFMode]]
        // If GetTopPosition is empty (vacuum), etac is all zeros and
        // the sub-system is ending.
        auto sub_top = qft::get_top_position(jpre_unique);
        std::vector<int> sub_etac(
            sub_fc.propagators_after_conservation.size(), 0);
        if (!sub_top.empty()) {
            auto sub_pos = qft::amf_position(sub_fc, sub_top, opts.amf_modes);
            if (!sub_pos.empty()) {
                sub_etac = qft::amf_eta_c(sub_fc, sub_pos);
            }
        }

        // Build sub_fc copy (qft::FamilyConfig is move-only).
        qft::FamilyConfig sub_fc_copy = qft::FamilyConfig::build(
            sub_fc.family, sub_fc.loops, sub_fc.legs,
            extract_conservation_strings(sub_fc),
            extract_replacement_strings(sub_fc),
            extract_propagator_strings(sub_fc),
            sub_fc.cut, sub_fc.prescription);

        AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
            std::cerr << "[single_mass_setup] parent fc family=" << fc.family
                      << " parent_n_props="
                      << fc.propagators_after_conservation.size()
                      << " fcomp.loops.size=" << fcomp.loops.size()
                      << " fcomp.props.size=" << fcomp.propagators.size()
                      << " fcomp.patterns.size=" << fcomp.patterns.size()
                      << std::endl;
            for (std::size_t pi = 0; pi < fcomp.patterns.size(); ++pi) {
                std::cerr << "  patt[" << pi << "] = [";
                for (auto v : fcomp.patterns[pi]) std::cerr << v << ",";
                std::cerr << "]" << std::endl;
            }
            std::cerr << "  drop=" << drop
                      << " droploop_name=" << droploop_name
                      << " comp_props.size=" << comp_props.size()
                      << " dropped_props.size=" << dropped_props.size()
                      << " completede.size=" << completede.size()
                      << std::endl;
            std::cerr << "[single_mass_setup] component fam="
                      << sub_fc.family
                      << " sub_loops=" << sub_loops.size()
                      << " sub_legs=" << sub_legs.size()
                      << " props=" << sub_fc.propagators_after_conservation.size()
                      << " sub_top.size=" << sub_top.size()
                      << " sub_etac=[";
            for (auto v : sub_etac) std::cerr << v << ",";
            std::cerr << "]" << std::endl;
            std::cerr << "[single_mass_setup] jpre.size=" << jpre.size()
                      << " jpre_unique.size=" << jpre_unique.size() << std::endl;
            for (const auto& jp : jpre_unique) {
                std::cerr << "  J[" << jp.family();
                for (auto x : jp.indices()) std::cerr << "," << x;
                std::cerr << "]" << std::endl;
            }
            std::cerr << "[single_mass_setup] gamma_params:" << std::endl;
            for (std::size_t i = 0; i < sm_cfg.gamma_params.size(); ++i) {
                const auto& gp = sm_cfg.gamma_params[i];
                std::cerr << "  gp[" << i << "] n=" << gp.n
                          << " m0=" << gp.m0
                          << " m_eps=" << gp.m_eps_coef
                          << " loop_num=" << gp.loop_num << std::endl;
            }
        }
        auto sub_sys = std::make_unique<AMFSystem>(
            std::move(sub_fc_copy), jpre_unique, std::move(sub_etac),
            EndingScheme::SingleMass, sub_opts);
        sub_sys->set_single_mass_config(std::move(sm_cfg));
        out.push_back(std::move(sub_sys));
    }

    return out;
}

}  // namespace

std::vector<std::unique_ptr<AMFSystem>>
amf_system_setup_master(const qft::FamilyConfig& fc,
                          const std::vector<qft::JIntegral>& preferred,
                          const AMFSystemOptions& opts) {
    if (preferred.empty()) return {};

    auto numeric_q = build_family_numeric_q(opts.bb);
    std::unique_ptr<qft::FamilyConfig> numeric_fc;
    const qft::FamilyConfig* fc_use = &fc;
    if (!numeric_q.empty()) {
        numeric_fc = std::make_unique<qft::FamilyConfig>(
            build_numeric_replacement_family(fc, numeric_q));
        fc_use = numeric_fc.get();
    }

    // Mirrors AMFSystemSetupMaster (line 1010):
    //   if (all schemes ending) -> use etac=0, scheme = first
    //   else -> recurse with the first non-ending scheme
    std::vector<std::unique_ptr<AMFSystem>> out;

    if (ending_q_all(*fc_use, preferred, opts)) {
        std::vector<int> zero_etac(fc_use->propagators_after_conservation.size(), 0);
        qft::FamilyConfig fc_copy = qft::FamilyConfig::build(
            fc_use->family, fc_use->loops, fc_use->legs,
            extract_conservation_strings(*fc_use),
            extract_replacement_strings(*fc_use),
            extract_propagator_strings(*fc_use),
            fc_use->cut, fc_use->prescription);
        auto sys = std::make_unique<AMFSystem>(
            std::move(fc_copy), preferred, std::move(zero_etac),
            opts.ending_schemes.front(), opts);
        out.push_back(std::move(sys));
        return out;
    }

    // Find the first non-ending scheme.
    EndingScheme scheme = EndingScheme::Tradition;
    for (auto s : opts.ending_schemes) {
        if (!ending_q(*fc_use, preferred, s, opts)) { scheme = s; break; }
    }

    if (scheme == EndingScheme::SingleMass) {
        AMFLOW_TRACE("AMFLOW_DEBUG_SCHEME") {
            auto top = qft::get_top_position(preferred);
            std::cerr << "[scheme] SingleMass fired: family="
                      << fc_use->family
                      << " top_position={";
            for (std::size_t i = 0; i < top.size(); ++i) {
                if (i) std::cerr << ',';
                std::cerr << top[i];
            }
            std::cerr << "}" << std::endl;
        }
        return single_mass_setup_master(*fc_use, preferred, opts);
    }

    if (scheme == EndingScheme::Cutkosky) {
        auto top = qft::get_top_position(preferred);
        auto info = qft::analyze_top_sector(*fc_use, top);
        long phase_loop_num = -1;
        for (const auto& comp : info) {
            if (qft::phase_volume_q(comp)) {
                if (phase_loop_num >= 0) {
                    throw std::runtime_error(
                        "amf_system_setup_master: Cutkosky found multiple "
                        "phase-volume components");
                }
                phase_loop_num = comp.loopnum;
            }
        }
        if (phase_loop_num < 0) {
            throw std::runtime_error(
                "amf_system_setup_master: Cutkosky found no "
                "phase-volume component");
        }

        AMFLOW_TRACE("AMFLOW_DEBUG_SCHEME") {
            std::cerr << "[scheme] Cutkosky fired: family="
                      << fc_use->family
                      << " phase_loop_num=" << phase_loop_num
                      << " top_position={";
            for (std::size_t i = 0; i < top.size(); ++i) {
                if (i) std::cerr << ',';
                std::cerr << top[i];
            }
            std::cerr << "}" << std::endl;
        }

        qft::FamilyConfig fc_cutless = qft::FamilyConfig::build(
            fc_use->family, fc_use->loops, fc_use->legs,
            extract_conservation_strings(*fc_use),
            extract_replacement_strings(*fc_use),
            extract_propagator_strings(*fc_use),
            /*cut=*/{}, /*prescription=*/{});

        auto pos = qft::amf_position(fc_cutless, top, opts.amf_modes);
        if (pos.empty() && !top.empty()) {
            throw std::runtime_error(
                "amf_system_setup_master: Cutkosky could not find an "
                "eta-injection position after clearing cut/prescription");
        }
        auto etac = qft::amf_eta_c(fc_cutless, pos);

        AMFSystem::CutkoskyConfig ck_cfg;
        ck_cfg.original_preferred = preferred;
        ck_cfg.phase_loop_num = phase_loop_num;

        auto sys = std::make_unique<AMFSystem>(
            std::move(fc_cutless), preferred, std::move(etac),
            scheme, opts);
        sys->set_cutkosky_config(std::move(ck_cfg));
        out.push_back(std::move(sys));
        return out;
    }

    if (scheme != EndingScheme::Tradition) {
        throw std::runtime_error(
            "amf_system_setup_master: unexpected scheme " +
            std::string(ending_scheme_name(scheme)));
    }

    auto top = qft::get_top_position(preferred);
    auto pos = qft::amf_position(*fc_use, top, opts.amf_modes);
    if (pos.empty()) {
        // This shouldn't happen since ending_q(Tradition) returned false.
        throw std::runtime_error(
            "amf_system_setup_master: no eta-injection position found");
    }
    auto etac = qft::amf_eta_c(*fc_use, pos);
    qft::FamilyConfig fc_copy = qft::FamilyConfig::build(
        fc_use->family, fc_use->loops, fc_use->legs,
        extract_conservation_strings(*fc_use),
        extract_replacement_strings(*fc_use),
        extract_propagator_strings(*fc_use),
        fc_use->cut, fc_use->prescription);
    auto sys = std::make_unique<AMFSystem>(
        std::move(fc_copy), preferred, std::move(etac), scheme, opts);
    out.push_back(std::move(sys));
    return out;
}

std::vector<std::unique_ptr<AMFSystem>>
amf_systems_setup(const qft::FamilyConfig& fc,
                  const std::vector<qft::JIntegral>& preferred,
                  const AMFSystemOptions& opts) {
    auto roots = amf_system_setup_master(fc, preferred, opts);
    for (auto& r : roots) r->setup();
    return roots;
}

std::vector<AMFSystemSolution>
amf_systems_solution(std::vector<std::unique_ptr<AMFSystem>>& systems,
                     const std::vector<numeric::AcbValue>& epslist) {
    if (systems.empty()) return {};
    for (auto& sys : systems) {
        sys->solve(epslist);
    }

    // Mirror Mathematica's AMFSystemCombineSolution (.m line 1188-1197):
    //   For each sysid, read GlobalPreferred = (global, coe, local, op).
    //   Compute  values[sysid][i] = coe[i] * op(local[i] /. sol[sysid]).
    //   All sysids share the same `global` list (Keys); the final
    //   answer per global[i] is the PRODUCT over sysids:
    //       result[i] = Product_sysid values[sysid][i]
    //
    // Tradition systems have an empty SingleMassConfig, so their
    // global_values is empty and they directly contribute their
    // master_values for the global preferred.
    std::vector<AMFSystemSolution> out;
    const auto& src0 = systems.front()->solutions();
    out.reserve(src0.size());
    long prec = numeric::working_prec_bits();

    for (std::size_t k = 0; k < src0.size(); ++k) {
        AMFSystemSolution cs;
        // master_values: root's solution (kept for inspection).
        cs.master_values.reserve(src0[k].master_values.size());
        for (const auto& v : src0[k].master_values) {
            cs.master_values.push_back(v.clone());
        }
        // global_preferred: take from the FIRST system's
        // global_preferred (all sysids share the same global list).
        // If the first system has no global_preferred (i.e. Tradition
        // single root), use master_values directly as the answer.
        if (!systems.front()->solutions()[k].global_preferred.empty()) {
            cs.global_preferred = systems.front()->solutions()[k].global_preferred;
            // Initialize global_values to all 1s (multiplicative identity).
            cs.global_values.reserve(cs.global_preferred.size());
            for (std::size_t i = 0; i < cs.global_preferred.size(); ++i) {
                numeric::AcbValue one;
                acb_one(one.raw());
                cs.global_values.push_back(std::move(one));
            }
            // Multiply per-system contributions.
            for (auto& sys : systems) {
                const auto& sol_k = sys->solutions()[k];
                if (sol_k.global_values.size() != cs.global_values.size()) {
                    throw std::runtime_error(
                        "amf_systems_solution: global_values size mismatch "
                        "across roots (FactorizeFamily components must "
                        "agree on global count)");
                }
                for (std::size_t i = 0; i < cs.global_values.size(); ++i) {
                    numeric::AcbValue tmp;
                    acb_mul(tmp.raw(), cs.global_values[i].raw(),
                              sol_k.global_values[i].raw(), prec);
                    cs.global_values[i] = std::move(tmp);
                }
            }
        }
        out.push_back(std::move(cs));
    }
    return out;
}

}  // namespace amflow::pipeline
