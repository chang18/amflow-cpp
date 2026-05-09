// SPDX-License-Identifier: MIT
#include "amflow/pipeline/solve_integrals.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/arf.h>
#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly_factor.h>

#include "amflow/ibp/reduce.hpp"
#include "amflow/numeric/options.hpp"
#include "amflow/numeric/rational.hpp"
#include "amflow/ode/amflow.hpp"
#include "amflow/ode/zero.hpp"

namespace amflow::pipeline {

namespace {

namespace fs = std::filesystem;

class FmpqHolder {
public:
    FmpqHolder() { fmpq_init(q_); }
    FmpqHolder(const FmpqHolder&) = delete;
    FmpqHolder& operator=(const FmpqHolder&) = delete;
    FmpqHolder(FmpqHolder&& other) noexcept {
        fmpq_init(q_);
        fmpq_swap(q_, other.q_);
    }
    FmpqHolder& operator=(FmpqHolder&& other) noexcept {
        if (this != &other) fmpq_swap(q_, other.q_);
        return *this;
    }
    ~FmpqHolder() { fmpq_clear(q_); }

    fmpq* raw() { return q_; }
    const fmpq* raw() const { return q_; }

private:
    fmpq_t q_;
};

void set_fmpq_from_ull(fmpq_t out, unsigned long long num, unsigned long long den) {
    fmpz_set_ui(fmpq_numref(out), num);
    fmpz_set_ui(fmpq_denref(out), den);
    fmpq_canonicalise(out);
}

FmpqHolder decimal_string_to_fmpq(long double x) {
    FmpqHolder out;
    if (x == 0.0L) {
        fmpq_zero(out.raw());
        return out;
    }

    std::ostringstream oss;
    oss << std::scientific
        << std::setprecision(std::numeric_limits<long double>::max_digits10)
        << x;
    std::string s = oss.str();

    bool neg = false;
    std::size_t pos = 0;
    if (s[pos] == '+' || s[pos] == '-') {
        neg = (s[pos] == '-');
        ++pos;
    }

    const auto epos = s.find_first_of("eE", pos);
    if (epos == std::string::npos) {
        throw std::runtime_error("decimal_string_to_fmpq: malformed scientific literal");
    }

    std::string mant = s.substr(pos, epos - pos);
    int exp10 = std::stoi(s.substr(epos + 1));

    std::string digits;
    long frac_digits = 0;
    bool seen_dot = false;
    for (char ch : mant) {
        if (ch == '.') {
            seen_dot = true;
            continue;
        }
        digits.push_back(ch);
        if (seen_dot) ++frac_digits;
    }
    while (!digits.empty() && digits.front() == '0') {
        digits.erase(digits.begin());
    }
    if (digits.empty()) {
        fmpq_zero(out.raw());
        return out;
    }

    fmpz_t num, den, pow10;
    fmpz_init(num);
    fmpz_init(den);
    fmpz_init(pow10);

    if (fmpz_set_str(num, digits.c_str(), 10) != 0) {
        fmpz_clear(num);
        fmpz_clear(den);
        fmpz_clear(pow10);
        throw std::runtime_error("decimal_string_to_fmpq: failed to parse mantissa");
    }
    fmpz_one(den);

    const long scale10 = exp10 - frac_digits;
    if (scale10 >= 0) {
        fmpz_set_ui(pow10, 10);
        fmpz_pow_ui(pow10, pow10, static_cast<unsigned long>(scale10));
        fmpz_mul(num, num, pow10);
    } else {
        fmpz_set_ui(pow10, 10);
        fmpz_pow_ui(pow10, pow10, static_cast<unsigned long>(-scale10));
        fmpz_mul(den, den, pow10);
    }
    if (neg) fmpz_neg(num, num);

    fmpq_set_fmpz_frac(out.raw(), num, den);
    fmpq_canonicalise(out.raw());

    fmpz_clear(num);
    fmpz_clear(den);
    fmpz_clear(pow10);
    return out;
}

FmpqHolder rationalize_machine_real(long double x, long double tol) {
    if (x == 0.0L) {
        FmpqHolder out;
        fmpq_zero(out.raw());
        return out;
    }
    const bool neg = x < 0.0L;
    if (neg) x = -x;

    if (x < 1.0L) {
        const long double inv = 1.0L / x;
        if (inv <= static_cast<long double>(
                std::numeric_limits<unsigned long long>::max())) {
            const auto den = static_cast<unsigned long long>(std::llround(inv));
            if (den != 0) {
                const long double approx = 1.0L / static_cast<long double>(den);
                if (std::fabs(approx - x) <= tol) {
                    FmpqHolder out;
                    set_fmpq_from_ull(out.raw(), 1ULL, den);
                    if (neg) fmpq_neg(out.raw(), out.raw());
                    return out;
                }
            }
        }
    }

    long double cur = x;
    unsigned long long h_nm2 = 0, h_nm1 = 1;
    unsigned long long k_nm2 = 1, k_nm1 = 0;
    for (int iter = 0; iter < 64; ++iter) {
        const auto a = static_cast<unsigned long long>(std::floor(cur));
        if (a != 0) {
            if (h_nm1 > (std::numeric_limits<unsigned long long>::max() - h_nm2) / a ||
                k_nm1 > (std::numeric_limits<unsigned long long>::max() - k_nm2) / a) {
                break;
            }
        }
        const unsigned long long h = a * h_nm1 + h_nm2;
        const unsigned long long k = a * k_nm1 + k_nm2;
        if (k == 0) break;
        const long double approx =
            static_cast<long double>(h) / static_cast<long double>(k);
        if (std::fabs(approx - x) <= tol) {
            FmpqHolder out;
            set_fmpq_from_ull(out.raw(), h, k);
            if (neg) fmpq_neg(out.raw(), out.raw());
            return out;
        }

        const long double frac = cur - std::floor(cur);
        if (std::fabs(frac) <= std::numeric_limits<long double>::epsilon()) break;
        cur = 1.0L / frac;
        h_nm2 = h_nm1; h_nm1 = h;
        k_nm2 = k_nm1; k_nm1 = k;
    }

    return decimal_string_to_fmpq(neg ? -x : x);
}

numeric::AcbValue acb_from_fmpq(const fmpq_t q, long prec = numeric::working_prec_bits()) {
    numeric::AcbValue out;
    out.set_fmpq(q, prec);
    return out;
}

numeric::AcbValue midpoint_only(const numeric::AcbValue& src) {
    numeric::AcbValue out;
    arb_set_arf(acb_realref(out.raw()), arb_midref(acb_realref(src.raw())));
    arb_set_arf(acb_imagref(out.raw()), arb_midref(acb_imagref(src.raw())));
    return out;
}

numeric::AcbValue midpoint_only_relative_chop(const numeric::AcbValue& src) {
    numeric::AcbValue out = midpoint_only(src);
    int digits = numeric::chop_pre();
    if (digits <= 0) return out;

    arf_t re_abs, im_abs;
    arf_init(re_abs);
    arf_init(im_abs);
    arf_abs(re_abs, arb_midref(acb_realref(out.raw())));
    arf_abs(im_abs, arb_midref(acb_imagref(out.raw())));

    arb_t scale, tol;
    arb_init(scale);
    arb_init(tol);
    arb_one(scale);
    if (arf_cmp(re_abs, arb_midref(scale)) > 0) {
        arb_set_arf(scale, re_abs);
    }
    if (arf_cmp(im_abs, arb_midref(scale)) > 0) {
        arb_set_arf(scale, im_abs);
    }

    const long bits = numeric::decimal_digits_to_bits(digits) + 32;
    arb_set_si(tol, 10);
    arb_inv(tol, tol, bits);
    arb_pow_ui(tol, tol, static_cast<unsigned long>(digits), bits);
    arb_mul(tol, tol, scale, bits);

    if (arf_cmp(re_abs, arb_midref(tol)) < 0) arb_zero(acb_realref(out.raw()));
    if (arf_cmp(im_abs, arb_midref(tol)) < 0) arb_zero(acb_imagref(out.raw()));

    arb_clear(scale);
    arb_clear(tol);
    arf_clear(re_abs);
    arf_clear(im_abs);
    return out;
}

void acb_pow_si_exact(acb_ptr out, acb_srcptr base, long exp, long prec) {
    if (exp == 0) {
        acb_one(out);
        return;
    }

    long e = exp;
    if (e < 0) e = -e;

    numeric::AcbValue result;
    result.set_one();
    numeric::AcbValue factor;
    acb_set(factor.raw(), base);
    while (e > 0) {
        if ((e & 1L) != 0) {
            acb_mul(result.raw(), result.raw(), factor.raw(), prec);
        }
        e >>= 1L;
        if (e > 0) {
            acb_mul(factor.raw(), factor.raw(), factor.raw(), prec);
        }
    }

    if (exp < 0) acb_inv(result.raw(), result.raw(), prec);
    acb_set(out, result.raw());
}

std::string jintegral_key(const qft::JIntegral& j) {
    std::ostringstream oss;
    oss << j.family();
    for (long idx : j.indices()) oss << "|" << idx;
    return oss.str();
}

std::string resolve_work_root(const AMFSystemOptions& opts,
                              const std::string& requested) {
    if (!requested.empty()) return requested;
    if (!opts.cache_root.empty()) return opts.cache_root;
    if (!opts.bb.work_dir.empty()) return opts.bb.work_dir;

    static std::atomic<unsigned long> counter{0};
    const auto id = counter.fetch_add(1, std::memory_order_relaxed) + 1;
    return (fs::temp_directory_path()
            / ("amflow_v2_" + std::to_string(id))).string();
}

bool parse_numeric_text_to_fmpq(const std::string& s, fmpq_t out) {
    const auto slash = s.find('/');
    try {
        if (slash != std::string::npos) {
            long p = std::stol(s.substr(0, slash));
            long q = std::stol(s.substr(slash + 1));
            fmpq_set_si(out, p, q);
            return true;
        }

        std::size_t pos = 0;
        long p = std::stol(s, &pos);
        if (pos == s.size()) {
            fmpq_set_si(out, p, 1);
            return true;
        }

        std::stod(s, &pos);
        if (pos != s.size()) return false;

        const auto dot = s.find('.');
        if (dot == std::string::npos) return false;
        std::string sint = s.substr(0, dot) + s.substr(dot + 1);
        long denom = 1;
        for (std::size_t i = dot + 1; i < s.size(); ++i) denom *= 10;
        long num = std::stol(sint);
        fmpq_set_si(out, num, denom);
        return true;
    } catch (...) {
        return false;
    }
}

numeric::AcbValue parse_numeric_text_to_acb(const std::string& s, long prec) {
    FmpqHolder q;
    if (!parse_numeric_text_to_fmpq(s, q.raw())) {
        throw std::runtime_error(
            "solve_integrals: cannot parse numeric value '" + s + "'");
    }
    return acb_from_fmpq(q.raw(), prec);
}

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

bool mfrac_has_var(const algebra::Mfrac& src, long v) {
    algebra::Mpoly num = src.numerator();
    if (mpoly_has_var(num, v)) return true;
    algebra::Mpoly den = src.denominator();
    return mpoly_has_var(den, v);
}

std::map<std::string, numeric::AcbValue>
build_numeric_acb_map(const ibp::ReduceOptions& bb,
                      const numeric::AcbValue& eps,
                      long prec) {
    std::map<std::string, numeric::AcbValue> out;
    out.emplace("eps", eps.clone());
    out.emplace("__amf_eps", eps.clone());

    numeric::AcbValue d_val;
    numeric::AcbValue four;
    four.set_si(4);
    numeric::AcbValue two_eps;
    acb_mul_2exp_si(two_eps.raw(), eps.raw(), 1);
    acb_sub(d_val.raw(), four.raw(), two_eps.raw(), prec);
    out.emplace("d", std::move(d_val));

    for (const auto& [name, value] : bb.numeric_values) {
        if (out.find(name) != out.end()) continue;
        out.emplace(name, parse_numeric_text_to_acb(value, prec));
    }
    return out;
}

numeric::AcbValue evaluate_numeric_mfrac(const algebra::Mfrac& src,
                                const std::map<std::string, numeric::AcbValue>& numeric_values,
                                long prec) {
    auto ctx = src.ctx();
    numeric::AcbVector values((std::size_t)ctx->n_vars());
    for (long v = 0; v < ctx->n_vars(); ++v) {
        acb_zero(values.at((std::size_t)v));
    }

    for (long v = 0; v < ctx->n_vars(); ++v) {
        if (!mfrac_has_var(src, v)) continue;

        const std::string& name = ctx->var_name(v);
        auto it = numeric_values.find(name);
        if (it != numeric_values.end()) {
            acb_set(values.at((std::size_t)v), it->second.raw());
            continue;
        }

        if (name.rfind("__amf_", 0) == 0
                || name.rfind("__feyn_", 0) == 0
                || name.rfind("__zsq_", 0) == 0) {
            continue;
        }

        throw std::runtime_error(
            "solve_integrals: no numeric value for '" + name + "'");
    }

    numeric::AcbValue out;
    src.evaluate_acb(out.raw(), values.raw(), prec);
    if (out.is_chop_zero()) out.set_zero();
    return out;
}

void validate_targets(const qft::FamilyConfig& fc,
                      const std::vector<qft::JIntegral>& jints) {
    for (const auto& j : jints) {
        if (j.family() != fc.family) {
            throw std::invalid_argument(
                "solve_integrals: target family '" + j.family()
                + "' does not match FamilyConfig '" + fc.family + "'");
        }
        if (j.n_indices() != fc.n_propagators()) {
            throw std::invalid_argument(
                "solve_integrals: target " + j.to_string()
                + " has wrong index count for family '" + fc.family + "'");
        }
    }
}

const numeric::AcbValue* lookup_solution_value(
        const std::vector<std::unique_ptr<AMFSystem>>& systems,
        const AMFSystemSolution& sol,
        const qft::JIntegral& target) {
    if (!sol.global_preferred.empty()) {
        for (std::size_t i = 0; i < sol.global_preferred.size(); ++i) {
            if (sol.global_preferred[i] == target) return &sol.global_values[i];
        }
    }

    const auto& preferred = systems.front()->preferred();
    for (std::size_t i = 0; i < preferred.size(); ++i) {
        if (preferred[i] == target) return &sol.master_values[i];
    }
    return nullptr;
}

std::vector<SampledIntegralSolution>
reorder_solutions(const std::vector<qft::JIntegral>& requested,
                  const std::map<std::string, std::vector<numeric::AcbValue>>& by_key) {
    std::vector<SampledIntegralSolution> out;
    out.reserve(requested.size());
    for (const auto& j : requested) {
        auto it = by_key.find(jintegral_key(j));
        if (it == by_key.end()) {
            throw std::runtime_error(
                "solve_integrals: internal missing solution for " + j.to_string());
        }
        SampledIntegralSolution row{j, {}};
        row.values.reserve(it->second.size());
        for (const auto& value : it->second) row.values.push_back(value.clone());
        out.push_back(std::move(row));
    }
    return out;
}

FmpqHolder fmpq_power_as_fmpq(const fmpq_t base, long exp) {
    if (fmpq_sgn(base) <= 0) {
        throw std::invalid_argument("fmpq_power_as_fmpq: base must be positive");
    }

    FmpqHolder out;
    fmpz_t num_pow, den_pow;
    fmpz_init(num_pow);
    fmpz_init(den_pow);
    const unsigned long abs_exp =
        static_cast<unsigned long>(exp < 0 ? -exp : exp);
    fmpz_pow_ui(num_pow, fmpq_numref(base), abs_exp);
    fmpz_pow_ui(den_pow, fmpq_denref(base), abs_exp);
    if (exp >= 0) {
        fmpz_set(fmpq_numref(out.raw()), num_pow);
        fmpz_set(fmpq_denref(out.raw()), den_pow);
    } else {
        fmpz_set(fmpq_numref(out.raw()), den_pow);
        fmpz_set(fmpq_denref(out.raw()), num_pow);
    }
    fmpq_canonicalise(out.raw());
    fmpz_clear(num_pow);
    fmpz_clear(den_pow);
    return out;
}

numeric::AcbValue integral_normalization_scale(const qft::FamilyConfig& fc,
                                      const qft::JIntegral& j,
                                      long prec) {
    if (fc.propagator_denominators.empty()) {
        numeric::AcbValue one;
        one.set_one();
        return one;
    }
    if (fc.propagator_denominators.size() != j.n_indices()) {
        throw std::runtime_error(
            "integral_normalization_scale: propagator scale count mismatch");
    }

    FmpqHolder acc;
    fmpq_one(acc.raw());
    for (std::size_t i = 0; i < j.n_indices(); ++i) {
        const long power = j.indices()[i];
        if (power == 0) continue;
        const std::string& den_text = fc.propagator_denominators[i];
        if (den_text == "1") continue;

        fmpq_t base;
        fmpq_init(base);
        if (fmpz_set_str(fmpq_numref(base), den_text.c_str(), 10) != 0) {
            fmpq_clear(base);
            throw std::runtime_error(
                "integral_normalization_scale: invalid propagator scale '"
                + den_text + "'");
        }
        fmpz_one(fmpq_denref(base));
        fmpq_canonicalise(base);
        FmpqHolder pw = fmpq_power_as_fmpq(base, power);
        fmpq_mul(acc.raw(), acc.raw(), pw.raw());
        fmpq_clear(base);
    }

    return acb_from_fmpq(acc.raw(), prec);
}

}  // namespace

NumericalConfig
generate_numerical_config(long loop_count,
                          long goal_digits,
                          long eps_order) {
    if (goal_digits <= 0) {
        throw std::invalid_argument(
            "generate_numerical_config: goal_digits must be positive");
    }
    if (eps_order < 0) {
        throw std::invalid_argument(
            "generate_numerical_config: eps_order must be non-negative");
    }

    const long loop = (loop_count > 0) ? loop_count : 1;
    const long number = static_cast<long>(
        std::ceil(2.5L * static_cast<long double>(eps_order)
                  + 2.0L * static_cast<long double>(loop)));
    if (number > 100) {
        throw std::runtime_error(
            "generate_numerical_config: the given order is too large to evaluate");
    }

    const long double exponent =
        -0.5L * static_cast<long double>(loop)
        - static_cast<long double>(goal_digits)
            / static_cast<long double>(eps_order + 1);
    const long double eps0 = std::pow(10.0L, exponent);
    auto eps0_q = rationalize_machine_real(eps0, eps0 / 100.0L);

    const long double single_pre_ld =
        (static_cast<long double>(number + 2 * loop))
        * (0.5L * static_cast<long double>(loop)
           + static_cast<long double>(goal_digits)
                / static_cast<long double>(eps_order + 1));
    const long single_pre = std::max<long>(
        static_cast<long>(std::ceil(single_pre_ld)), 30L);

    NumericalConfig out;
    out.working_pre = 2 * single_pre;
    out.x_order = 4 * single_pre;
    out.eps_samples.reserve(static_cast<std::size_t>(number));
    const long eps_prec =
        numeric::decimal_digits_to_bits(std::max<long>(out.x_order, 2 * out.working_pre));

    fmpq_t q, scale;
    fmpq_init(q);
    fmpq_init(scale);
    for (long i = 1; i <= number; ++i) {
        fmpq_set(q, eps0_q.raw());
        fmpq_set_si(scale, 100 + i, 100);
        fmpq_mul(q, q, scale);
        fmpq_canonicalise(q);
        out.eps_samples.push_back(acb_from_fmpq(q, eps_prec));
    }
    fmpq_clear(q);
    fmpq_clear(scale);

    return out;
}

std::vector<numeric::AcbValue>
fit_eps_numeric(const std::vector<numeric::AcbValue>& eps_samples,
                const std::vector<numeric::AcbValue>& values,
                long leading_order,
                long prec) {
    const long n = static_cast<long>(eps_samples.size());
    const long fit_prec =
        std::max<long>(2 * prec,
                       prec + numeric::decimal_digits_to_bits(std::max<long>(64, 16 * n)));

    numeric::AcbVector xs((std::size_t)n);
    numeric::AcbVector ys((std::size_t)n);
    for (long i = 0; i < n; ++i) {
        const auto eps_mid = midpoint_only_relative_chop(eps_samples[(std::size_t)i]);
        const auto val_mid = midpoint_only_relative_chop(values[(std::size_t)i]);
        if (eps_mid.is_zero()) {
            throw std::runtime_error("fit_eps: epsilon sample is zero");
        }
        numeric::AcbValue eps_leading;
        acb_pow_si_exact(eps_leading.raw(), eps_mid.raw(), leading_order, fit_prec);
        acb_div(ys.at((std::size_t)i), val_mid.raw(), eps_leading.raw(), fit_prec);
        acb_set(xs.at((std::size_t)i), eps_mid.raw());
    }

    numeric::AcbPoly poly;
    acb_poly_interpolate_newton(poly.raw(), xs.raw(), ys.raw(), n, fit_prec);

    std::vector<numeric::AcbValue> out;
    out.reserve((std::size_t)n);
    for (long i = 0; i < n; ++i) {
        numeric::AcbValue coeff;
        acb_poly_get_coeff_acb(coeff.raw(), poly.raw(), i);
        if (coeff.is_chop_zero()) coeff.set_zero();
        acb_set_round(coeff.raw(), coeff.raw(), prec);
        out.push_back(std::move(coeff));
    }
    return out;
}

std::vector<numeric::AcbValue>
fit_eps(const std::vector<numeric::AcbValue>& eps_samples,
        const std::vector<numeric::AcbValue>& values,
        long leading_order,
        long prec) {
    if (eps_samples.size() != values.size() || eps_samples.empty()) {
        throw std::invalid_argument(
            "fit_eps: eps_samples and values must have the same non-zero size");
    }

    return fit_eps_numeric(eps_samples, values, leading_order, prec);
}

std::vector<numeric::AcbValue>
trim_leading_zero_coeffs(std::vector<numeric::AcbValue> coeffs,
                         long* leading_order) {
    std::size_t first_nz = 0;
    while (first_nz < coeffs.size() && coeffs[first_nz].is_chop_zero()) {
        coeffs[first_nz].set_zero();
        ++first_nz;
    }
    if (first_nz == coeffs.size()) {
        if (leading_order) *leading_order = 0;
        std::vector<numeric::AcbValue> out(1);
        out[0].set_zero();
        return out;
    }
    if (leading_order) *leading_order += static_cast<long>(first_nz);

    std::vector<numeric::AcbValue> out;
    out.reserve(coeffs.size() - first_nz);
    for (std::size_t i = first_nz; i < coeffs.size(); ++i) {
        out.push_back(coeffs[i].clone());
    }
    return out;
}

std::vector<numeric::AcbValue>
truncate_coeffs_to_highest_order(std::vector<numeric::AcbValue> coeffs,
                                 long leading_order,
                                 long highest_order,
                                 long* adjusted_leading_order) {
    if (coeffs.empty()) {
        if (adjusted_leading_order) *adjusted_leading_order = 0;
        std::vector<numeric::AcbValue> out(1);
        out[0].set_zero();
        return out;
    }

    if (highest_order < leading_order) {
        if (adjusted_leading_order) *adjusted_leading_order = 0;
        std::vector<numeric::AcbValue> out(1);
        out[0].set_zero();
        return out;
    }

    const long keep =
        std::min<long>(static_cast<long>(coeffs.size()),
                       highest_order - leading_order + 1);
    coeffs.resize(static_cast<std::size_t>(keep));
    if (adjusted_leading_order) *adjusted_leading_order = leading_order;
    return coeffs;
}

std::vector<SampledIntegralSolution>
black_box_amflow_single(const qft::FamilyConfig& fc,
                        const std::vector<qft::JIntegral>& jints,
                        const std::vector<numeric::AcbValue>& eps_samples,
                        const AMFSystemOptions& opts,
                        const std::string& work_dir) {
    if (jints.empty()) return {};
    if (eps_samples.empty()) {
        throw std::invalid_argument(
            "black_box_amflow_single: eps_samples must be non-empty");
    }
    validate_targets(fc, jints);

    const auto root = resolve_work_root(opts, work_dir);
    const auto reduce_dir = root + "/reduce";
    const auto amf_dir = root + "/amf";
    fs::remove_all(reduce_dir);
    fs::remove_all(amf_dir);

    ibp::ReduceOptions bb = opts.bb;
    bb.work_dir = reduce_dir;
    auto reduction = ibp::reduce(fc, jints, /*preferred=*/{},
                                      qft::get_top_sector(jints), bb);

    std::vector<std::vector<numeric::AcbValue>> master_grid(reduction.masters.size());
    if (!reduction.masters.empty()) {
        AMFSystemOptions amf_opts = opts;
        amf_opts.cache_root = amf_dir;
        auto systems = amf_systems_setup(fc, reduction.masters, amf_opts);
        const auto sols = amf_systems_solution(systems, eps_samples);
        if (sols.size() != eps_samples.size()) {
            throw std::runtime_error(
                "black_box_amflow_single: AMF solution grid size mismatch");
        }

        for (std::size_t m = 0; m < reduction.masters.size(); ++m) {
            master_grid[m].reserve(eps_samples.size());
            for (std::size_t k = 0; k < eps_samples.size(); ++k) {
                const numeric::AcbValue* got =
                    lookup_solution_value(systems, sols[k], reduction.masters[m]);
                if (got == nullptr) {
                    throw std::runtime_error(
                        "black_box_amflow_single: master "
                        + reduction.masters[m].to_string()
                        + " was not solved by AMFSystem");
                }
                master_grid[m].push_back(got->clone());
            }
        }
    }

    std::map<std::string, std::size_t> master_index;
    for (std::size_t i = 0; i < reduction.masters.size(); ++i) {
        master_index.emplace(jintegral_key(reduction.masters[i]), i);
    }

    std::vector<SampledIntegralSolution> out;
    out.reserve(jints.size());
    const long prec = numeric::working_prec_bits();

    for (std::size_t i = 0; i < jints.size(); ++i) {
        SampledIntegralSolution row{jints[i], {}};
        row.values.reserve(eps_samples.size());
        for (std::size_t k = 0; k < eps_samples.size(); ++k) {
            numeric::AcbValue sum;
            sum.set_zero();

            const auto numeric_values =
                build_numeric_acb_map(opts.bb, eps_samples[k], prec);
            for (const auto& term : reduction.rules[i]) {
                auto mit = master_index.find(jintegral_key(term.integ));
                if (mit == master_index.end()) {
                    throw std::runtime_error(
                        "black_box_amflow_single: reduction term "
                        + term.integ.to_string() + " is not in masters");
                }
                const auto coef = evaluate_numeric_mfrac(term.coef, numeric_values, prec);
                numeric::AcbValue prod;
                acb_mul(prod.raw(), coef.raw(),
                        master_grid[mit->second][k].raw(), prec);
                acb_add(sum.raw(), sum.raw(), prod.raw(), prec);
            }

            numeric::AcbValue normalization =
                integral_normalization_scale(fc, jints[i], prec);
            if (!normalization.is_one()) {
                acb_mul(sum.raw(), sum.raw(), normalization.raw(), prec);
            }
            if (sum.is_chop_zero()) sum.set_zero();
            row.values.push_back(std::move(sum));
        }
        out.push_back(std::move(row));
    }

    return out;
}

std::vector<SampledIntegralSolution>
black_box_amflow(const qft::FamilyConfig& fc,
                 const std::vector<qft::JIntegral>& jints,
                 const std::vector<numeric::AcbValue>& eps_samples,
                 const AMFSystemOptions& opts,
                 const std::string& work_dir) {
    if (jints.empty()) return {};
    validate_targets(fc, jints);

    const auto root = resolve_work_root(opts, work_dir);
    const auto parts = qft::split_target(jints);

    std::map<std::string, std::vector<numeric::AcbValue>> by_key;
    for (std::size_t i = 0; i < parts.size(); ++i) {
        auto part_sol = black_box_amflow_single(
            fc, parts[i], eps_samples, opts, root + "/part_" + std::to_string(i));
        for (auto& row : part_sol) {
            auto key = jintegral_key(row.integral);
            if (by_key.find(key) != by_key.end()) {
                throw std::runtime_error(
                    "black_box_amflow: duplicated target solution for "
                    + row.integral.to_string());
            }
            by_key.emplace(std::move(key), std::move(row.values));
        }
    }

    return reorder_solutions(jints, by_key);
}

// Apply the D0-driven eps shift `(4 - D0) / 2` on top of `user_eps`.
// Mirrors AMFlow.m:1342 / 1351:  epslist = epslist0 + (4-$D0)*1/2.
//
// `user_eps` is the user-facing grid (where eps -> 0 corresponds to d = D0);
// the returned grid is what gets handed to BlackBoxAMFlow so the family is
// evaluated in the internal d = 4 - 2*eps_internal frame.
static std::vector<numeric::AcbValue>
apply_d0_eps_shift(const std::vector<numeric::AcbValue>& user_eps, long prec) {
    fmpq_t shift;
    fmpq_init(shift);
    numeric::d0_eps_shift_fmpq(shift);

    std::vector<numeric::AcbValue> out;
    out.reserve(user_eps.size());

    if (fmpq_is_zero(shift)) {
        // Fast path: D0 == 4, no work to do.
        for (const auto& e : user_eps) out.push_back(e.clone());
        fmpq_clear(shift);
        return out;
    }

    numeric::AcbValue shift_acb;
    shift_acb.set_fmpq(shift, prec);
    fmpq_clear(shift);

    for (const auto& e : user_eps) {
        numeric::AcbValue s;
        acb_add(s.raw(), e.raw(), shift_acb.raw(), prec);
        out.push_back(std::move(s));
    }
    return out;
}

std::vector<LaurentIntegralSolution>
solve_integrals(const qft::FamilyConfig& fc,
                const std::vector<qft::JIntegral>& jints,
                long goal_digits,
                long eps_order,
                const AMFSystemOptions& opts,
                const std::string& work_dir) {
    if (goal_digits <= 0) {
        throw std::invalid_argument(
            "solve_integrals: goal_digits must be positive");
    }
    if (eps_order < 0) {
        throw std::invalid_argument(
            "solve_integrals: eps_order must be non-negative");
    }

    // Single-eps fast path (mirrors AMFlow.m:1364-1374): if the user
    // pinned `eps` to a numeric value via `numeric_values["eps"]`, we
    // skip the Laurent fit entirely and return one coefficient at
    // order 0 representing the integral evaluated at that eps.
    auto eps_pin_it = opts.bb.numeric_values.find("eps");
    if (eps_pin_it != opts.bb.numeric_values.end()) {
        // Working precision and truncation tighten with the goal.
        numeric::GlobalScope fast_scope;
        fast_scope.global.working_pre = static_cast<int>(2 * goal_digits);
        fast_scope.expansion.x_order  = static_cast<int>(4 * goal_digits);
        if (fast_scope.expansion.extra_x_order < fast_scope.expansion.x_order) {
            fast_scope.expansion.extra_x_order = fast_scope.expansion.x_order;
        }
        fast_scope.commit();

        // Parse user's `eps` and apply the D0 shift `(4 - D0)/2`.
        fmpq_t user_eps; fmpq_init(user_eps);
        if (fmpq_set_str(user_eps, eps_pin_it->second.c_str(), 10) != 0) {
            fmpq_clear(user_eps);
            throw std::runtime_error(
                "solve_integrals: failed to parse user `eps` rational '"
                + eps_pin_it->second + "'");
        }
        fmpq_t shift; fmpq_init(shift);
        numeric::d0_eps_shift_fmpq(shift);
        fmpq_t internal_eps; fmpq_init(internal_eps);
        fmpq_add(internal_eps, user_eps, shift);

        const long prec = numeric::working_prec_bits();
        numeric::AcbValue sample;
        arb_set_fmpq(acb_realref(sample.raw()), internal_eps, prec);
        arb_zero(acb_imagref(sample.raw()));

        fmpq_clear(user_eps);
        fmpq_clear(shift);
        fmpq_clear(internal_eps);

        std::vector<numeric::AcbValue> internal_samples;
        internal_samples.push_back(std::move(sample));

        const auto sampled =
            black_box_amflow(fc, jints, internal_samples, opts, work_dir);

        std::vector<LaurentIntegralSolution> out;
        out.reserve(sampled.size());
        for (const auto& row : sampled) {
            if (row.values.size() != 1) {
                throw std::runtime_error(
                    "solve_integrals fast path: expected exactly one "
                    "sampled value per integral");
            }
            LaurentIntegralSolution sol;
            sol.integral = row.integral;
            // Express the result as a degenerate one-term Laurent at
            // order 0 (the value AT the user-supplied eps).  Mirrors
            // upstream `Thread[Keys[sol] -> Values[sol][[All, 1]]]`
            // (AMFlow.m:1372).
            sol.leading_order = 0;
            sol.coefficients.push_back(row.values[0].clone());
            out.push_back(std::move(sol));
        }
        return out;
    }

    auto cfg = generate_numerical_config(
        static_cast<long>(fc.n_loops()), goal_digits, eps_order);

    numeric::GlobalScope scope;
    scope.global.working_pre = static_cast<int>(cfg.working_pre);
    // AMFlow.m's default RationalizePre is 100; keeping the lower project
    // baseline here loses too much information on tiny-eps SolveIntegrals grids.
    if (scope.global.rationalize_pre < 100) {
        scope.global.rationalize_pre = 100;
    }
    // The Mathematica path keeps exact rationals much deeper into the
    // normalization pipeline. Our numeric bridge later asks Layer 7 to
    // recover exact rationals to RationalizePre digits, so low-goal fits
    // need a modest working-precision floor above that threshold.
    const int min_working_pre = scope.global.rationalize_pre + 32;
    if (scope.global.working_pre < min_working_pre) {
        scope.global.working_pre = min_working_pre;
    }
    scope.expansion.x_order = static_cast<int>(cfg.x_order);
    if (scope.expansion.extra_x_order < scope.expansion.x_order) {
        scope.expansion.extra_x_order = scope.expansion.x_order;
    }
    scope.commit();

    // D0 != 4 sampling: evaluate at the internal eps grid, but still fit
    // the Laurent expansion against the user-facing grid (matches MMA).
    const auto internal_samples =
        apply_d0_eps_shift(cfg.eps_samples, numeric::working_prec_bits());

    const auto sampled =
        black_box_amflow(fc, jints, internal_samples, opts, work_dir);

    const long fit_leading_order = -2 * static_cast<long>(fc.n_loops());
    const long highest_order = eps_order + fit_leading_order;

    std::vector<LaurentIntegralSolution> out;
    out.reserve(sampled.size());
    const long prec = numeric::working_prec_bits();
    for (const auto& row : sampled) {
        if (row.values.size() != cfg.eps_samples.size()) {
            throw std::runtime_error(
                "solve_integrals: sampled value count does not match epsilon grid");
        }

        auto coeffs = fit_eps(cfg.eps_samples, row.values, fit_leading_order, prec);
        long actual_leading = fit_leading_order;
        coeffs = trim_leading_zero_coeffs(std::move(coeffs), &actual_leading);
        coeffs = truncate_coeffs_to_highest_order(
            std::move(coeffs), actual_leading, highest_order, &actual_leading);

        LaurentIntegralSolution sol;
        sol.integral = row.integral;
        sol.leading_order = actual_leading;
        sol.coefficients.reserve(coeffs.size());
        for (std::size_t i = 0; i < coeffs.size(); ++i) {
            sol.coefficients.push_back(coeffs[i].clone());
        }
        out.push_back(std::move(sol));
    }

    return out;
}

}  // namespace amflow::pipeline
