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
//
// Scope (audit divergence row 201)
// --------------------------------
//   Variables `s_name` for which differentiation is supported:
//
//     1. variables that appear directly as a polynomial-coefficient
//        in `fc.propagators_after_conservation[k]` — i.e. free
//        symbols (typical examples: `eta` after η-insertion, mass
//        parameters like `msq`, `m1sq`, `m2sq`).
//
//   For Replacement-defined kinematic invariants (e.g. `s` defined
//   via `p^2 -> s` in `fc.replacement`), this routine returns the
//   all-zero derivative because propagators are stored
//   pre-replacement (the symbol does not appear in the polynomial).
//   Upstream's `LIBPDeriv` provides an additional
//   momentum-derivative chain-rule path (`LIBPDerivivative` in
//   `Kira/interface.m`) that handles this case; the C++ port
//   intentionally does not implement it because the production
//   AMFlow flow (`pipeline::AMFSystem::setup`,
//   `src/pipeline/amfsystem.cpp:1012`) only ever calls
//   `libp_deriv(..., "eta")`, where `eta` is added directly to
//   propagators and lives in the polynomial.
//
//   Tests exercising the supported scope: `test_ibp_libp_deriv.cpp`
//   `*TwoMassBubble*` cases (multi-mass bubble with three
//   distinct invariants {m1sq, m2sq, s}).

#ifndef AMFLOW_IBP_LIBP_DERIV_HPP
#define AMFLOW_IBP_LIBP_DERIV_HPP

#include <map>
#include <memory>
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

// D14 (2026-05-16) overload: same as `libp_denoms_deriv(fc, s_name)`
// but applies `numeric_values` substitution to every intermediate
// Mfrac AND reprojects the result onto `target_ctx`.  Use this when
// you have a narrow target context (e.g. `{eta, d}` from
// `make_reduction_context`) and want to avoid FLINT's
// multivariate-GCD overhead on the wide `fc.ctx` polynomial ring.
//
// Preconditions: `target_ctx` must contain `s_name` as one of its
// variables, and `numeric_values` must supply rationals for every
// `fc.ctx` variable that is NOT in `target_ctx`.
//
// Returned `coef[k][j]` and `constant[k]` Mfracs live on `target_ctx`.
LibpDenomsDerivResult
libp_denoms_deriv(const qft::FamilyConfig& fc, const std::string& s_name,
                  const std::map<std::string, std::string>& numeric_values,
                  const std::shared_ptr<algebra::MpolyContext>& target_ctx);

// D14 (2026-05-16) overload of libp_deriv mirroring the above:
// returns DerivTerms whose `.coef` lives on `target_ctx`.
std::vector<DerivTerm>
libp_deriv(const qft::FamilyConfig& fc,
           const qft::JIntegral& j,
           const std::string& s_name,
           const std::map<std::string, std::string>& numeric_values,
           const std::shared_ptr<algebra::MpolyContext>& target_ctx);

std::vector<DerivTerm>
compute_derivative(const qft::FamilyConfig& fc,
                   const std::vector<DerivTerm>& expr,
                   const std::string& x_name);

std::vector<DerivTerm>
simplify_terms(std::vector<DerivTerm> terms);

}  // namespace amflow::ibp

#endif  // AMFLOW_IBP_LIBP_DERIV_HPP
