// SPDX-License-Identifier: MIT
// ibp::kira — Kira IBP-reduction interface.
//
//
// Mirrors AMFlow.m / ibp_interface/Kira/interface.m.
//
// Workflow
// --------
//   1. Caller constructs a `KiraConfig` from a FamilyConfig + JIntegral
//      preferred masters + IBP rank/dot/topsector + numeric values for
//      kinematic invariants.
//   2. yaml writers populate <dir>/config/{integralfamilies,kinematics}.yaml,
//      <dir>/jobs.yaml, <dir>/preferred, <dir>/target.
//   3. kira_run invokes the kira binary.
//   4. parsers read masters and target table.

#ifndef AMFLOW_IBP_KIRA_HPP
#define AMFLOW_IBP_KIRA_HPP

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "amflow/algebra/mpoly.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/jintegral.hpp"

namespace amflow::ibp {

struct KiraConfig {
    const qft::FamilyConfig* fc = nullptr;

    std::vector<int> top_pattern;

    long ibp_rank = 5;
    long ibp_dot  = 0;
    long n_thread = 1;
    long integral_order = 5;

    std::map<std::string, std::string> numeric_values;

    std::string kira_executable;
    std::string fermat_executable;
};

void kira_write_config(const KiraConfig& cfg, const std::string& dir);

enum class KiraReductionMode {
    Masters,
    Reduce,
};

void kira_write_jobs(const KiraConfig& cfg,
                     const std::string& dir,
                     KiraReductionMode mode);

void kira_write_preferred(const std::vector<qft::JIntegral>& preferred,
                            const qft::FamilyConfig& fc,
                            const std::string& dir);

void kira_write_targets(const std::vector<qft::JIntegral>& targets,
                         const qft::FamilyConfig& fc,
                         const std::string& dir);

double kira_run(const KiraConfig& cfg,
                const std::string& dir,
                const std::string& log_path = "");

std::vector<qft::JIntegral>
kira_read_masters(const qft::FamilyConfig& fc, const std::string& dir);

struct KiraReductionRule {
    qft::JIntegral lhs;
    std::vector<std::pair<std::string, qft::JIntegral>> rhs;
};

std::vector<KiraReductionRule>
kira_read_target_table(const qft::FamilyConfig& fc, const std::string& dir);

algebra::Mfrac
kira_parse_expression(const std::shared_ptr<algebra::MpolyContext>& ctx,
                       const std::string& expr);

}  // namespace amflow::ibp

#endif  // AMFLOW_IBP_KIRA_HPP
