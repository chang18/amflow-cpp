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

std::string mfrac_to_kira_neg(const Mfrac& f) {
    Mfrac neg = -f;
    return mfrac_to_kira(neg);
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

    std::ostringstream props;
    for (std::size_t i = 0; i < sq.momenta.size(); ++i) {
        props << "\n      - [ \"" << mfrac_to_kira(sq.momenta[i])
              << "\", " << mfrac_to_kira_neg(sq.masses[i]) << " ]";
    }

    std::string cut = positions_where_one(fc->cut);

    auto inv_names = mass_scale_names(cfg);
    std::ostringstream kin;
    for (const auto& v : inv_names) {
        kin << "\n    - [" << v << ", 2]";
    }

    std::ostringstream rep;
    for (const auto& [key, rhs] : fc->reduced_replacement) {
        std::size_t sep = key.find('*');
        if (sep == std::string::npos) {
            throw std::runtime_error(
                "kira_write_config: malformed leg-pair key '" + key + "'");
        }
        std::string p1 = key.substr(0, sep);
        std::string p2 = key.substr(sep + 1);
        rep << "\n    - [[" << p1 << "," << p2 << "], "
            << mfrac_to_kira(rhs) << "]";
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
    long n_top = std::count(cfg.top_pattern.begin(),
                              cfg.top_pattern.end(), 1);
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
