// SPDX-License-Identifier: MIT
// ibp::libp_deriv — symbolic derivative of J integrals via chain rule.
//
//
// Mirrors AMFlow.m / Kira interface lines 282-405.
//
//   For J(a) = Integral[ Π_k 1 / D_k^{a_k} ]
//
//      ∂J(a)/∂s = Σ_k (∂D_k/∂s) · (-a_k · J(a + e_k))
//
//   where (∂D_k/∂s) is decomposed over the propagator basis as
//
//      ∂D_k/∂s = Σ_j coef_kj · D_j  +  const_k
//
//   so that  D_j · J(a + e_k) = J(a + e_k - e_j).

#ifndef AMFLOW_IBP_LIBP_DERIV_HPP
#define AMFLOW_IBP_LIBP_DERIV_HPP

#include <string>
#include <vector>

#include "amflow/algebra/mpoly.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/jintegral.hpp"
#include "amflow/qft/region.hpp"   // DListContext

namespace amflow::ibp {

struct DerivTerm {
    algebra::Mfrac    coef;
    qft::JIntegral    integ;
};

struct LibpDenomsDerivResult {
    std::vector<algebra::Mpoly>               completede;
    std::vector<std::vector<algebra::Mfrac>>  coef;
    std::vector<algebra::Mfrac>               constant;
};

LibpDenomsDerivResult
libp_denoms_deriv(const qft::FamilyConfig& fc, const std::string& s_name);

std::vector<DerivTerm>
libp_deriv(const qft::FamilyConfig& fc,
           const qft::JIntegral& j,
           const std::string& s_name);

std::vector<DerivTerm>
compute_derivative(const qft::FamilyConfig& fc,
                   const std::vector<DerivTerm>& expr,
                   const std::string& x_name);

std::vector<DerivTerm>
simplify_terms(std::vector<DerivTerm> terms);

}  // namespace amflow::ibp

#endif  // AMFLOW_IBP_LIBP_DERIV_HPP
