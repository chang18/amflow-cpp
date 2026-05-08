// SPDX-License-Identifier: MIT
// ibp::reduce — top-level Kira-driven IBP reduction + diffeq.
//
//
// Mirrors AMFlow.m lines 1218-1247 plus ibp_interface/Kira/interface.m
// AnalyticReduction / DifferentialEquation (lines 454-505).
//
// The "BlackBox" abstraction (which originally existed to allow swapping
// reducers) has been flattened: this project is Kira-only per user
// constraint, so these helpers call Kira directly.

#ifndef AMFLOW_IBP_REDUCE_HPP
#define AMFLOW_IBP_REDUCE_HPP

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "amflow/algebra/mpoly.hpp"
#include "amflow/ibp/kira.hpp"
#include "amflow/ibp/libp_deriv.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/jintegral.hpp"

namespace amflow::ibp {

struct ReduceOptions {
    long ibp_rank = 5;
    long ibp_dot  = 0;
    long n_thread = 1;
    long integral_order = 5;
    std::map<std::string, std::string> numeric_values;
    std::string kira_executable;
    std::string fermat_executable;

    std::string work_dir;
    std::string log_file;
};

struct ReductionContext {
    std::shared_ptr<algebra::MpolyContext> ctx;
    long d_var;
    std::vector<algebra::Mpoly> lift_gens;
};

ReductionContext make_reduction_context(const qft::FamilyConfig& fc);

struct ReduceResult {
    std::vector<qft::JIntegral>             masters;
    std::vector<std::vector<DerivTerm>>     rules;
    ReductionContext                         red_ctx;
};

ReduceResult
reduce(const qft::FamilyConfig& fc,
       const std::vector<qft::JIntegral>& targets,
       const std::vector<qft::JIntegral>& preferred,
       const std::vector<int>& top_pattern,
       const ReduceOptions& opts);

struct DiffeqResult {
    std::vector<qft::JIntegral>                              sortedmasters;
    std::vector<std::string>                                 vars;
    std::vector<std::vector<std::vector<algebra::Mfrac>>>   diffeq;
    ReductionContext                                          red_ctx;
};

DiffeqResult
diffeq(const qft::FamilyConfig& fc,
       const std::vector<qft::JIntegral>& jpreferred,
       const std::vector<std::string>& vars,
       const std::vector<int>& top_pattern,
       const ReduceOptions& opts);

}  // namespace amflow::ibp

#endif  // AMFLOW_IBP_REDUCE_HPP
