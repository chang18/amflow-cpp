// SPDX-License-Identifier: MIT
// ibp::kira yaml writers — implementation.
//

#include "amflow/ibp/kira.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>

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
    std::vector<Mfrac> sub_momenta;
    sub_momenta.reserve(sq.momenta.size());
    for (const auto& m : sq.momenta) {
        sub_momenta.push_back(subs.apply(m));
    }

    std::ostringstream props;
    for (std::size_t i = 0; i < sub_momenta.size(); ++i) {
        Mfrac neg_mass = -sub_masses[i].clone();
        props << "\n      - [ \"" << mfrac_to_kira(sub_momenta[i])
              << "\", " << mfrac_to_kira(neg_mass) << " ]";
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
