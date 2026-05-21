// SPDX-License-Identifier: MIT
// ibp::kira yaml writers — implementation.
//

#include "amflow/ibp/kira.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_mpoly_q.h>

#include "amflow/qft/region.hpp"

namespace amflow::ibp {

using algebra::Mfrac;
using algebra::Mpoly;

namespace {

namespace fs = std::filesystem;

void ensure_dir(const std::string& dir) {
    fs::create_directories(dir);
}

std::string list_to_string(const std::vector<std::string>& items) {
    std::ostringstream oss;
    for (std::size_t i = 0; i < items.size(); ++i) {
        if (i) oss << ", ";
        oss << items[i];
    }
    return oss.str();
}

std::string mfrac_to_kira(const Mfrac& f) {
    Mpoly num = f.numerator();
    Mpoly den = f.denominator();
    if (den.is_one()) {
        return num.to_string();
    }
    if (num.is_zero()) {
        return "0";
    }
    std::ostringstream oss;
    oss << "(" << num.to_string() << ")/(" << den.to_string() << ")";
    return oss.str();
}

// Try to factor `prop` into closed form `±(v·x)² + const_part`,
// where `x` are momentum variables (loops + reduced legs).  Returns
// nullopt when the polynomial isn't a rank-1 quadratic in momenta —
// caller then falls back to plain expanded form.
//
// MMA's `interface.m:126` emits `[Propagator, 0]` using Mathematica's
// symbolic-expression printer, which preserves `(linear_combo)^2` form
// because Mathematica doesn't auto-expand `Power[Plus[...], 2]`.  C++
// stores propagators as FLINT fmpz_mpoly (expanded by construction);
// without this reconstruction step, Kira sees the expanded monomial sum
// and walks a generic-quadratic IBP path that produces ~25× larger
// mandatory lists than the closed-form path MMA triggers.  This is the
// follow-up to the `[Propagator, 0]` yaml-format fix recorded in
// `notes/future_optimization_proposals.md` 2026-05-20.
//
// Math:
//   A physical propagator `(Σ v_j x_j)² + const` has quadratic-part
//   coefficient matrix A = ±v v^T (rank-1).  Given expanded form,
//   extract A from term-by-term coeffs, pick a non-zero diagonal pivot
//   A_{ii}, set sign = sgn(A_{ii}) and v_i = √|A_{ii}|, derive v_j =
//   A_{j,i}/(sign·v_i), then verify A == sign·vv^T.  Fails (→ fallback)
//   on: non-polynomial denominator, mixed momentum-mass coupling like
//   `eta·l1·l2`, linear-in-momentum term, non-perfect-square pivot, or
//   non-rank-1 A.
std::optional<std::string>
try_render_closed_form(const Mfrac& prop,
                       const std::vector<long>& mom_var_idx,
                       const std::vector<std::string>& mom_var_name) {
    Mpoly num = prop.numerator();
    Mpoly den = prop.denominator();
    if (!den.is_one()) return std::nullopt;
    if (num.is_zero())  return std::nullopt;

    auto ctx = prop.ctx();
    const long n_total = ctx->n_vars();
    const long n_mom   = (long)mom_var_idx.size();

    // mom_pos[v] = k if ctx-var v is the k-th momentum, else -1
    std::vector<long> mom_pos((std::size_t)n_total, -1);
    for (long k = 0; k < n_mom; ++k) mom_pos[mom_var_idx[(std::size_t)k]] = k;

    std::vector<long> A_diag((std::size_t)n_mom, 0);
    std::map<std::pair<long, long>, long> A_off;  // {i<j} -> A_{ij}

    Mpoly const_part = Mpoly::zero(ctx);

    fmpz_mpoly_struct*       p_raw = num.raw();
    fmpz_mpoly_ctx_struct*   c_raw = ctx->raw();
    long n_terms = fmpz_mpoly_length(p_raw, c_raw);

    std::vector<unsigned long> exp((std::size_t)n_total);
    fmpz_t coeff;
    fmpz_init(coeff);

    auto fail = [&]() -> std::optional<std::string> {
        fmpz_clear(coeff);
        return std::nullopt;
    };

    for (long t = 0; t < n_terms; ++t) {
        fmpz_mpoly_get_term_coeff_fmpz(coeff, p_raw, t, c_raw);
        fmpz_mpoly_get_term_exp_ui(exp.data(), p_raw, t, c_raw);

        long mom_sum = 0;
        long mom_max = 0;
        long total_deg = 0;
        for (long v = 0; v < n_total; ++v) {
            total_deg += (long)exp[(std::size_t)v];
            if (mom_pos[(std::size_t)v] >= 0) {
                mom_sum += (long)exp[(std::size_t)v];
                if ((long)exp[(std::size_t)v] > mom_max)
                    mom_max = (long)exp[(std::size_t)v];
            }
        }

        if (mom_sum == 0) {
            Mpoly term = Mpoly::monomial(ctx, coeff, exp);
            const_part += term;
            continue;
        }
        if (mom_sum != 2)        return fail();   // linear or > quadratic
        if (total_deg != 2)      return fail();   // mass·momentum mixing
        if (!fmpz_fits_si(coeff)) return fail();
        long c = fmpz_get_si(coeff);

        if (mom_max == 2) {
            long i = -1;
            for (long k = 0; k < n_mom; ++k)
                if (exp[(std::size_t)mom_var_idx[(std::size_t)k]] == 2) {
                    i = k;
                    break;
                }
            A_diag[(std::size_t)i] += c;
        } else {
            long i = -1, j = -1;
            for (long k = 0; k < n_mom; ++k) {
                if (exp[(std::size_t)mom_var_idx[(std::size_t)k]] == 1) {
                    if      (i < 0) i = k;
                    else if (j < 0) { j = k; break; }
                }
            }
            if (i < 0 || j < 0) return fail();
            if (c % 2 != 0)     return fail();  // off-diag must be even
            A_off[{i, j}] += c / 2;
        }
    }
    fmpz_clear(coeff);

    // Find pivot.
    long pivot = -1;
    for (long i = 0; i < n_mom; ++i)
        if (A_diag[(std::size_t)i] != 0) { pivot = i; break; }
    if (pivot < 0) return std::nullopt;

    long pivot_val = A_diag[(std::size_t)pivot];
    long abs_pivot = pivot_val > 0 ? pivot_val : -pivot_val;
    int  sign      = pivot_val > 0 ? 1 : -1;

    long v_pivot = 0;
    for (long s = 1; s * s <= abs_pivot; ++s)
        if (s * s == abs_pivot) { v_pivot = s; break; }
    if (v_pivot == 0) return std::nullopt;

    std::vector<long> v((std::size_t)n_mom, 0);
    v[(std::size_t)pivot] = v_pivot;

    auto get_A = [&](long i, long j) -> long {
        if (i == j) return A_diag[(std::size_t)i];
        long a = std::min(i, j), b = std::max(i, j);
        auto it = A_off.find({a, b});
        return it != A_off.end() ? it->second : 0;
    };

    for (long j = 0; j < n_mom; ++j) {
        if (j == pivot) continue;
        long a_j_pivot = get_A(j, pivot);
        long denom = sign * v_pivot;
        if (a_j_pivot % denom != 0) return std::nullopt;
        v[(std::size_t)j] = a_j_pivot / denom;
    }

    // Verify A == sign·vv^T everywhere.
    for (long i = 0; i < n_mom; ++i) {
        if (A_diag[(std::size_t)i] != sign * v[(std::size_t)i] * v[(std::size_t)i])
            return std::nullopt;
    }
    for (long i = 0; i < n_mom; ++i) {
        for (long j = i + 1; j < n_mom; ++j) {
            long expected = sign * v[(std::size_t)i] * v[(std::size_t)j];
            if (get_A(i, j) != expected) return std::nullopt;
        }
    }

    // Render `(v·x)`.  Single-var with |v|=1 gets no parens (matches MMA's
    // `l1^2` rather than `(l1)^2`).
    std::vector<std::string> terms;
    for (long k = 0; k < n_mom; ++k) {
        long c = v[(std::size_t)k];
        if (c == 0) continue;
        std::ostringstream t;
        bool first = terms.empty();
        if (c == 1) {
            if (!first) t << " + ";
            t << mom_var_name[(std::size_t)k];
        } else if (c == -1) {
            t << (first ? "-" : " - ") << mom_var_name[(std::size_t)k];
        } else if (c > 0) {
            if (!first) t << " + ";
            t << c << "*" << mom_var_name[(std::size_t)k];
        } else {
            t << (first ? "-" : " - ") << (-c) << "*" << mom_var_name[(std::size_t)k];
        }
        terms.push_back(t.str());
    }
    if (terms.empty()) return std::nullopt;

    std::ostringstream out;
    if (sign < 0) out << "-";
    bool single_term = terms.size() == 1;
    if (single_term) {
        long k_only = -1;
        for (long k = 0; k < n_mom; ++k)
            if (v[(std::size_t)k] != 0) { k_only = k; break; }
        long c = v[(std::size_t)k_only];
        long abs_c = c > 0 ? c : -c;
        if (abs_c == 1) {
            out << mom_var_name[(std::size_t)k_only] << "^2";
        } else {
            out << "(" << abs_c << "*" << mom_var_name[(std::size_t)k_only] << ")^2";
        }
    } else {
        out << "(";
        for (const auto& s : terms) out << s;
        out << ")^2";
    }

    if (!const_part.is_zero()) {
        std::string c_str = const_part.to_string();
        // Cosmetic: keep "+ -" out by absorbing a leading minus into " - ".
        if (!c_str.empty() && c_str[0] == '-') {
            out << " - " << c_str.substr(1);
        } else {
            out << " + " << c_str;
        }
    }
    return out.str();
}

std::string top_pattern_to_int_str(const std::vector<int>& pat) {
    unsigned long long total = 0;
    for (std::size_t i = 0; i < pat.size(); ++i) {
        if (pat[i]) total |= (1ull << i);
    }
    return std::to_string(total);
}

std::string positions_where_one(const std::vector<int>& pat) {
    std::vector<std::string> out;
    for (std::size_t i = 0; i < pat.size(); ++i) {
        if (pat[i]) out.push_back(std::to_string(i + 1));
    }
    return list_to_string(out);
}

// Names of all symbolic scales appearing in fc (post-conservation
// variable context).  Used to compute kinematic_invariants BEFORE
// numeric substitution is applied; the actual yaml entry is then
// filtered to only those that REMAIN after substitution.
std::vector<std::string>
mass_scale_names(const KiraConfig& cfg) {
    std::set<std::string> names;
    auto fc = cfg.fc;
    for (long i = 0; i < fc->ctx->n_vars(); ++i) {
        const std::string& v = fc->ctx->var_name(i);
        long n_loop = (long)fc->n_loops();
        long n_red  = (long)fc->n_red_legs();
        if (i < n_loop) continue;
        if (i >= n_loop && i < n_loop + n_red) continue;
        names.insert(v);
    }
    for (const auto& [k, _] : cfg.numeric_values) {
        if (k != "eps") names.insert(k);
    }
    std::vector<std::string> sorted(names.begin(), names.end());
    return sorted;
}

// RAII fmpq with non-throwing move; used to hold pre-parsed numeric
// substitution values for the SubsExpr helper below.
struct FmpqOwned {
    fmpq_t q;
    FmpqOwned()  { fmpq_init(q); }
    ~FmpqOwned() { fmpq_clear(q); }
    FmpqOwned(const FmpqOwned&) = delete;
    FmpqOwned& operator=(const FmpqOwned&) = delete;
};

// Mirrors `SPToSTU /. IBPRule` (Kira/interface.m:53): apply the user's
// numeric replacements (e.g. {s12 -> -2, mWsq -> 1}) directly into the
// scalarproduct_rules and propagator masses we hand to Kira.  Without
// this, Kira's Masters-mode call enumerates symbolic SPs and cannot
// detect scaleless sub-sectors, over-including masters that MMA
// (which substitutes upstream) correctly drops.  This caused the
// pentabox_2L "preferred master not in Kira's master list" failure.
struct NumericSubs {
    std::vector<long> var_idx;
    std::vector<std::unique_ptr<FmpqOwned>> values;

    static NumericSubs build(const KiraConfig& cfg) {
        NumericSubs out;
        auto fc = cfg.fc;
        for (const auto& [name, val_str] : cfg.numeric_values) {
            if (name == "eps") continue;
            long idx = fc->ctx->var_index(name);
            if (idx < 0) continue;  // var unused by this fc
            auto holder = std::make_unique<FmpqOwned>();
            if (fmpq_set_str(holder->q, val_str.c_str(), 10) != 0) {
                throw std::runtime_error(
                    "kira_write_config: cannot parse numeric value '"
                    + val_str + "' for '" + name + "'");
            }
            out.var_idx.push_back(idx);
            out.values.push_back(std::move(holder));
        }
        return out;
    }

    algebra::Mfrac apply(const algebra::Mfrac& in) const {
        algebra::Mfrac r = in.clone();
        for (std::size_t i = 0; i < var_idx.size(); ++i) {
            r = r.substitute(var_idx[i], values[i]->q);
        }
        return r;
    }

    algebra::Mfrac apply_poly(const algebra::Mpoly& in) const {
        return apply(algebra::Mfrac::from_mpoly(in.clone()));
    }
};

// True iff `var_idx` appears with non-zero exponent in either the
// numerator or denominator of `f`.
bool mfrac_uses_var(const algebra::Mfrac& f, long var_idx) {
    return f.numerator().degree(var_idx) > 0 ||
           f.denominator().degree(var_idx) > 0;
}

std::string jintegral_to_kira(const qft::JIntegral& j,
                                const qft::FamilyConfig& fc) {
    std::ostringstream oss;
    oss << fc.family << "[";
    for (std::size_t i = 0; i < j.indices().size(); ++i) {
        if (i) oss << ", ";
        oss << j.indices()[i];
    }
    oss << "]";
    return oss.str();
}

}  // namespace

void kira_write_config(const KiraConfig& cfg, const std::string& dir) {
    auto fc = cfg.fc;
    if (!fc) {
        throw std::invalid_argument("kira_write_config: fc is null");
    }
    if (cfg.top_pattern.size() !=
            fc->propagators_after_conservation.size()) {
        throw std::invalid_argument(
            "kira_write_config: top_pattern length != #propagators");
    }

    std::string config_dir = dir + "/config";
    ensure_dir(config_dir);

    qft::ToSquareResult sq =
        qft::to_square_all(*fc, fc->propagators_after_conservation);

    std::vector<std::string> loop_names;
    for (const auto& l : fc->loops) loop_names.push_back(l);

    std::vector<std::string> leg_names;
    std::string aux_leg = fc->family + "AuxLeg";
    if (!fc->reduced_legs.empty()) {
        for (const auto& p : fc->reduced_legs) leg_names.push_back(p);
        leg_names.push_back(aux_leg);
    }

    std::string momcon;
    if (!fc->reduced_legs.empty()) {
        std::ostringstream oss;
        oss << aux_leg << ", -(";
        for (std::size_t i = 0; i < fc->reduced_legs.size(); ++i) {
            if (i) oss << " + ";
            oss << fc->reduced_legs[i];
        }
        oss << ")";
        momcon = oss.str();
    }

    std::string top = top_pattern_to_int_str(cfg.top_pattern);

    NumericSubs subs = NumericSubs::build(cfg);

    // Collect substituted SP rule RHSs and propagator masses; the
    // kinematic_invariants list is the set of mass_scale_names that
    // still appear after substitution (mirrors MMA's MassScale which
    // is computed from SPToSTU AFTER /. IBPRule).
    auto inv_names = mass_scale_names(cfg);
    std::vector<Mfrac> sub_sp_rhs;
    sub_sp_rhs.reserve(fc->reduced_replacement.size());
    for (const auto& [key, rhs] : fc->reduced_replacement) {
        sub_sp_rhs.push_back(subs.apply(rhs));
    }
    std::vector<Mfrac> sub_masses;
    sub_masses.reserve(sq.masses.size());
    for (const auto& m : sq.masses) {
        sub_masses.push_back(subs.apply(m));
    }
    // Propagator yaml: emit each propagator as a single complete
    // denominator expression with mass-field=0, mirroring MMA's current
    // `interface.m:126` (`[Propagator, 0]`).  The earlier
    // `[Momentum, -Mass]` form (still seen in the commented-out
    // `interface.m:125`) is mathematically equivalent but causes Kira
    // 2.x to walk a different IBP-enumeration path, inflating
    // mandatory-list size by ~10^4 and producing 2 spurious masters on
    // 4L photon_4L_SE_TwoBubbles mass1 (336760 vs 18; 18 vs 16 masters).
    //
    // Beyond the yaml-format fix, the polynomial expression is also
    // factored back to closed form `(v·x)² + const` via
    // `try_render_closed_form` so Kira sees `(l1 - l4)^2 - eta - 1`
    // instead of the FLINT-stored expanded `l1² - 2·l1·l4 + l4² - eta - 1`.
    // MMA preserves closed form natively (Mathematica's symbolic engine
    // doesn't auto-expand `Power[Plus[...], 2]`); FLINT mpoly is
    // expanded by construction, so we reconstruct.  Without this step
    // Kira still mismatches MMA by ~25× on sunset_bubble_4L because its
    // rank-1-quadratic-form symmetry detector relies on closed form to
    // identify single squared momenta.
    //
    // Momentum vars for the rank-1 factoring: loops + reduced (post-
    // conservation) legs.  Non-momentum vars (eta, mass scales,
    // kinematic invariants) flow into `const_part` and render as-is.
    std::vector<long> mom_var_idx;
    std::vector<std::string> mom_var_name;
    auto push_mom = [&](const std::string& n) {
        long vi = fc->ctx->var_index(n);
        if (vi >= 0) { mom_var_idx.push_back(vi); mom_var_name.push_back(n); }
    };
    for (const auto& l : fc->loops)        push_mom(l);
    for (const auto& l : fc->reduced_legs) push_mom(l);

    std::ostringstream props;
    for (const auto& p : fc->propagators_after_conservation) {
        Mfrac prop_full = Mfrac::from_mpoly(p.clone());
        Mfrac sub_prop = subs.apply(prop_full);
        auto closed = try_render_closed_form(sub_prop, mom_var_idx, mom_var_name);
        const std::string body =
            closed.has_value() ? *closed : mfrac_to_kira(sub_prop);
        props << "\n      - [ \"" << body << "\", 0 ]";
    }

    std::string cut = positions_where_one(fc->cut);

    std::ostringstream kin;
    for (const auto& v : inv_names) {
        long vidx = fc->ctx->var_index(v);
        if (vidx < 0) continue;
        bool used = false;
        for (const auto& f : sub_sp_rhs) {
            if (mfrac_uses_var(f, vidx)) { used = true; break; }
        }
        if (!used) {
            for (const auto& f : sub_masses) {
                if (mfrac_uses_var(f, vidx)) { used = true; break; }
            }
        }
        if (used) {
            kin << "\n    - [" << v << ", 2]";
        }
    }

    std::ostringstream rep;
    {
        std::size_t idx = 0;
        for (const auto& [key, rhs] : fc->reduced_replacement) {
            std::size_t sep = key.find('*');
            if (sep == std::string::npos) {
                throw std::runtime_error(
                    "kira_write_config: malformed leg-pair key '" + key + "'");
            }
            std::string p1 = key.substr(0, sep);
            std::string p2 = key.substr(sep + 1);
            rep << "\n    - [[" << p1 << "," << p2 << "], "
                << mfrac_to_kira(sub_sp_rhs[idx]) << "]";
            ++idx;
        }
    }

    std::ofstream f1(config_dir + "/integralfamilies.yaml");
    if (!f1) throw std::runtime_error(
        "kira_write_config: cannot open integralfamilies.yaml");
    f1 << "integralfamilies:\n"
       << "  - name: \"" << fc->family << "\"\n"
       << "    loop_momenta: [" << list_to_string(loop_names) << "]\n"
       << "    top_level_sectors: [" << top << "]\n"
       << "    propagators:" << props.str() << "\n"
       << "    cut_propagators: [" << cut << "]\n";
    f1.close();

    std::ofstream f2(config_dir + "/kinematics.yaml");
    if (!f2) throw std::runtime_error(
        "kira_write_config: cannot open kinematics.yaml");
    f2 << "kinematics:\n"
       << "  incoming_momenta: [" << list_to_string(leg_names) << "]\n"
       << "  outgoing_momenta: []\n"
       << "  momentum_conservation: [" << momcon << "]\n"
       << "  kinematic_invariants:" << kin.str() << "\n"
       << "  scalarproduct_rules:" << rep.str() << "\n";
    f2.close();
}

void kira_write_jobs(const KiraConfig& cfg,
                     const std::string& dir,
                     KiraReductionMode mode) {
    auto fc = cfg.fc;
    if (!fc) throw std::invalid_argument("kira_write_jobs: fc is null");

    std::string fam = fc->family;
    // Mirrors upstream `r = Length[TopSector] - Count[TopSector, 0] +
    // IBPDot` (Kira/interface.m yaml emit) — i.e. count of non-zero
    // entries in top_pattern, plus IBPDot.  In production
    // `qft::get_top_sector` returns a strict 0/1 vector
    // (`src/qft/jintegral.cpp:131`), so for the production-relevant
    // shape `count(==1)` and `count(!=0)` agree.  Use the
    // count-non-zero form here for internal consistency with
    // `top_pattern_to_int_str` (line 65) and `positions_where_one`
    // (line 73), both of which already treat any non-zero entry as
    // "present", and to match the upstream formula literally so the
    // yaml output stays parity-correct under any future relaxation
    // of `get_top_sector`'s 0/1 guarantee.
    long n_top = std::count_if(cfg.top_pattern.begin(),
                                cfg.top_pattern.end(),
                                [](int v) { return v != 0; });
    long r = n_top + cfg.ibp_dot;

    std::string top_int = top_pattern_to_int_str(cfg.top_pattern);

    std::ostringstream y;
    y << "jobs:\n"
      << " - reduce_sectors:\n"
      << "    reduce:\n"
      << "     - {topologies: [" << fam
      << "], sectors: [" << top_int
      << "], r: " << r << ", s: " << cfg.ibp_rank << "}\n";

    if (mode == KiraReductionMode::Masters) {
        y << "    select_integrals:\n"
          << "      select_mandatory_recursively:\n"
          << "        - {topologies: [" << fam
          << "], sectors: [" << top_int
          << "], r: " << r << ", s: " << cfg.ibp_rank
          << ", d: " << cfg.ibp_dot << "}\n";
        y << "    preferred_masters: preferred\n"
          << "    integral_ordering: " << cfg.integral_order << "\n";
        y << "    run_initiate: masters\n";
    } else {
        y << "    select_integrals:\n"
          << "      select_mandatory_list:\n"
          << "        - [" << fam << ", target]\n";
        y << "    preferred_masters: preferred\n"
          << "    integral_ordering: " << cfg.integral_order << "\n";
        y << "    run_initiate: true\n"
          << "    run_triangular: true\n"
          << "    run_back_substitution: true\n";
        y << " - kira2math:\n"
          << "    target:\n"
          << "     - [" << fam << ", target]\n";
    }

    std::ofstream out(dir + "/jobs.yaml");
    if (!out) throw std::runtime_error(
        "kira_write_jobs: cannot open jobs.yaml");
    out << y.str();
}

void kira_write_preferred(const std::vector<qft::JIntegral>& preferred,
                            const qft::FamilyConfig& fc,
                            const std::string& dir) {
    std::ofstream out(dir + "/preferred");
    if (!out) throw std::runtime_error(
        "kira_write_preferred: cannot open file");
    // Write in caller-provided order: MMA's `Preferred[preferred, dir]`
    // (Kira/interface.m:164-170) writes the list as-is, no sort.  Kira's
    // IBP enumeration is sensitive to the preferred order, so any
    // reordering here diverges from MMA's Kira-input pattern.  Upstream
    // master list is the Kira-raw output of the prior reduce/preheat
    // step, which is already in Kira's natural "ascending by
    // sector/complexity" order on both C++ and MMA when preferred is
    // empty (verified on sunset_bubble_4L).  See
    // `sort_integrals_like_amflow` for the previous over-correction.
    for (const auto& j : preferred) {
        out << jintegral_to_kira(j, fc) << "\n\n";
    }
}

void kira_write_targets(const std::vector<qft::JIntegral>& targets,
                         const qft::FamilyConfig& fc,
                         const std::string& dir) {
    std::ofstream out(dir + "/target");
    if (!out) throw std::runtime_error(
        "kira_write_targets: cannot open file");
    for (const auto& j : targets) {
        out << jintegral_to_kira(j, fc) << "\n\n";
    }
}

}  // namespace amflow::ibp
