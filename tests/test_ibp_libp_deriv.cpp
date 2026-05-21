// SPDX-License-Identifier: MIT
// Tests for amflow::ibp::libp_deriv.

#include <gtest/gtest.h>

#include <algorithm>

#include "amflow/algebra/mpoly.hpp"
#include "amflow/ibp/libp_deriv.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/jintegral.hpp"

namespace alg = amflow::algebra;
namespace ibp = amflow::ibp;
namespace qft = amflow::qft;

TEST(LibpDerivTest, LibpDenomsDeriv_Bubble_Wrt_msq) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    auto res = ibp::libp_denoms_deriv(fc, "msq");
    ASSERT_EQ(res.coef.size(), 2u);
    for (auto& c : res.coef[0]) EXPECT_TRUE(c.is_zero()) << c.to_string();
    EXPECT_EQ(res.constant[0].to_string(), "-1");
    for (auto& c : res.coef[1]) EXPECT_TRUE(c.is_zero()) << c.to_string();
    EXPECT_EQ(res.constant[1].to_string(), "-1");
}

TEST(LibpDerivTest, LibpDenomsDeriv_Bubble_Wrt_s_AllZero) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    auto res = ibp::libp_denoms_deriv(fc, "s");
    for (std::size_t k = 0; k < res.coef.size(); ++k) {
        for (auto& c : res.coef[k]) EXPECT_TRUE(c.is_zero()) << c.to_string();
        EXPECT_TRUE(res.constant[k].is_zero()) << res.constant[k].to_string();
    }
}

TEST(LibpDerivTest, LibpDeriv_J11_Wrt_msq) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    qft::JIntegral j("bubble", {1, 1});
    auto terms = ibp::libp_deriv(fc, j, "msq");
    auto simplified = ibp::simplify_terms(std::move(terms));
    ASSERT_EQ(simplified.size(), 2u);

    std::sort(simplified.begin(), simplified.end(),
              [](const ibp::DerivTerm& a, const ibp::DerivTerm& b) {
                  return a.integ.indices() < b.integ.indices();
              });
    EXPECT_EQ(simplified[0].integ.indices(), (std::vector<long>{1, 2}));
    EXPECT_EQ(simplified[0].coef.to_string(), "1");
    EXPECT_EQ(simplified[1].integ.indices(), (std::vector<long>{2, 1}));
    EXPECT_EQ(simplified[1].coef.to_string(), "1");
}

TEST(LibpDerivTest, LibpDeriv_J_DegenerateZeroIndex) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    qft::JIntegral j("bubble", {0, 1});
    auto terms = ibp::libp_deriv(fc, j, "msq");
    auto s = ibp::simplify_terms(std::move(terms));
    ASSERT_EQ(s.size(), 1u);
    EXPECT_EQ(s[0].integ.indices(), (std::vector<long>{0, 2}));
    EXPECT_EQ(s[0].coef.to_string(), "1");
}

TEST(LibpDerivTest, LibpDeriv_J_HigherIndex) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    qft::JIntegral j("bubble", {2, 1});
    auto s = ibp::simplify_terms(ibp::libp_deriv(fc, j, "msq"));
    ASSERT_EQ(s.size(), 2u);
    std::sort(s.begin(), s.end(),
              [](const ibp::DerivTerm& a, const ibp::DerivTerm& b) {
                  return a.integ.indices() < b.integ.indices();
              });
    EXPECT_EQ(s[0].integ.indices(), (std::vector<long>{2, 2}));
    EXPECT_EQ(s[0].coef.to_string(), "1");
    EXPECT_EQ(s[1].integ.indices(), (std::vector<long>{3, 1}));
    EXPECT_EQ(s[1].coef.to_string(), "2");
}

TEST(LibpDerivTest, ComputeDerivative_LinearCombo) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    long msq_idx = fc.var_index("msq");
    alg::Mfrac two = alg::Mfrac::from_si(fc.ctx, 2);
    alg::Mfrac msq_mfrac = alg::Mfrac::from_mpoly(
        alg::Mpoly::variable(fc.ctx, msq_idx));

    std::vector<ibp::DerivTerm> expr;
    {
        ibp::DerivTerm t;
        t.coef = two.clone();
        t.integ = qft::JIntegral("bubble", {1, 1});
        expr.push_back(std::move(t));
    }
    {
        ibp::DerivTerm t;
        t.coef = msq_mfrac.clone();
        t.integ = qft::JIntegral("bubble", {0, 1});
        expr.push_back(std::move(t));
    }

    auto out = ibp::simplify_terms(ibp::compute_derivative(fc, expr, "msq"));

    std::sort(out.begin(), out.end(),
              [](const ibp::DerivTerm& a, const ibp::DerivTerm& b) {
                  return a.integ.indices() < b.integ.indices();
              });

    ASSERT_EQ(out.size(), 4u);
    EXPECT_EQ(out[0].integ.indices(), (std::vector<long>{0, 1}));
    EXPECT_EQ(out[0].coef.to_string(), "1");
    EXPECT_EQ(out[1].integ.indices(), (std::vector<long>{0, 2}));
    EXPECT_EQ(out[1].coef.to_string(), "msq");
    EXPECT_EQ(out[2].integ.indices(), (std::vector<long>{1, 2}));
    EXPECT_EQ(out[2].coef.to_string(), "2");
    EXPECT_EQ(out[3].integ.indices(), (std::vector<long>{2, 1}));
    EXPECT_EQ(out[3].coef.to_string(), "2");
}

// --- Multi-invariant cases -----
//
// The existing tests above all use the single-mass bubble family
// `{l^2 - msq, (l-p)^2 - msq}` — one mass parameter `msq` plus one
// Replacement-defined kinematic invariant `s` (`p^2 -> s`).  These
// tests use a two-mass bubble with three distinct invariants
// {m1sq, m2sq, s} so that LIBPDeriv is exercised on a family where:
//   - different propagators carry different free-symbol invariants
//     (m1sq only in p_0, m2sq only in p_1),
//   - dp_k/d(invariant) is nonzero on a strict subset of k's
//     (exercising the `if (!ck.is_zero())` short-circuit at
//     `src/ibp/libp_deriv.cpp:327`),
//   - simplify_terms / compute_derivative handle a coefficient
//     ring with multiple invariants flowing through.
//
// All expected outputs are hand-derived from the formula
//   d/dx (J[a_0, a_1]) = sum_k (-a_k) * (dp_k/dx as constant + sp
//                                          expansion) * J_shifted
// (mirror of upstream LIBPDeriv with the sign of upstream's `-ak`
// factor preserved in `libp_deriv.cpp:308` `neg_ak`).
//
// Scope note: differentiation w.r.t. a *Replacement-defined*
// kinematic invariant (e.g., `s` defined via `p^2 -> s`) is the
// `LibpDenomsDeriv_Bubble_Wrt_s_AllZero` case above — C++ returns
// all-zero, which is consistent with upstream's *naive* D path
// (`LIBPDenomsDeriv[topo, s_] := ... D[#, s] /. LIBPSpsToJ[topo]
// ...`, `Kira/interface.m`).  Upstream additionally provides a
// momentum-derivative chain-rule path via `LIBPDerivivative` that
// is precomputed for each `LIBPInvariants[topo]`; the C++ port
// does not implement that path because the production AMFlow flow
// only differentiates w.r.t. `eta` (added directly to propagators
// pre-replacement).  See
// `docs/AUDIT_MMA_PARITY.md`.

namespace {

qft::FamilyConfig make_two_mass_bubble() {
    return qft::FamilyConfig::build(
        "bub2m", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - m1sq", "(l - p)^2 - m2sq"});
}

}  // namespace

TEST(LibpDerivTest, LibpDenomsDeriv_TwoMassBubble_Wrt_m1sq) {
    // m1sq appears only in propagator 0 → only k=0 has a nonzero
    // constant term (-1); k=1 row is all zero.
    auto fc = make_two_mass_bubble();
    auto res = ibp::libp_denoms_deriv(fc, "m1sq");
    ASSERT_EQ(res.coef.size(), 2u);
    for (auto& c : res.coef[0]) EXPECT_TRUE(c.is_zero()) << c.to_string();
    EXPECT_EQ(res.constant[0].to_string(), "-1");
    for (auto& c : res.coef[1]) EXPECT_TRUE(c.is_zero()) << c.to_string();
    EXPECT_TRUE(res.constant[1].is_zero()) << res.constant[1].to_string();
}

TEST(LibpDerivTest, LibpDenomsDeriv_TwoMassBubble_Wrt_m2sq) {
    auto fc = make_two_mass_bubble();
    auto res = ibp::libp_denoms_deriv(fc, "m2sq");
    ASSERT_EQ(res.coef.size(), 2u);
    for (auto& c : res.coef[0]) EXPECT_TRUE(c.is_zero()) << c.to_string();
    EXPECT_TRUE(res.constant[0].is_zero()) << res.constant[0].to_string();
    for (auto& c : res.coef[1]) EXPECT_TRUE(c.is_zero()) << c.to_string();
    EXPECT_EQ(res.constant[1].to_string(), "-1");
}

TEST(LibpDerivTest, LibpDeriv_J11_TwoMassBubble_Wrt_m1sq) {
    // d/d(m1sq) J[1,1] = (-1)*(-1) * J[2,1] = +J[2,1].
    auto fc = make_two_mass_bubble();
    qft::JIntegral j("bub2m", {1, 1});
    auto s = ibp::simplify_terms(ibp::libp_deriv(fc, j, "m1sq"));
    ASSERT_EQ(s.size(), 1u);
    EXPECT_EQ(s[0].integ.indices(), (std::vector<long>{2, 1}));
    EXPECT_EQ(s[0].coef.to_string(), "1");
}

TEST(LibpDerivTest, LibpDeriv_J11_TwoMassBubble_Wrt_m2sq) {
    auto fc = make_two_mass_bubble();
    qft::JIntegral j("bub2m", {1, 1});
    auto s = ibp::simplify_terms(ibp::libp_deriv(fc, j, "m2sq"));
    ASSERT_EQ(s.size(), 1u);
    EXPECT_EQ(s[0].integ.indices(), (std::vector<long>{1, 2}));
    EXPECT_EQ(s[0].coef.to_string(), "1");
}

TEST(LibpDerivTest, LibpDeriv_J11_TwoMassBubble_Wrt_s_StaysZero) {
    // Replacement-defined invariant: C++ libp_denoms_deriv returns
    // all-zero (matches upstream's naive D path; differs from
    // upstream's momentum-derivative chain-rule path, which the
    // C++ port deliberately does not implement — see scope note
    // above).  Locks the
    // documented behavior: NO terms produced for s-derivatives.
    auto fc = make_two_mass_bubble();
    qft::JIntegral j("bub2m", {1, 1});
    auto s = ibp::simplify_terms(ibp::libp_deriv(fc, j, "s"));
    EXPECT_EQ(s.size(), 0u);
}

TEST(LibpDerivTest, LibpDeriv_J21_TwoMassBubble_Wrt_m1sq) {
    // For J[2,1]: a_0 = 2, neg_a0 = -2.  Contribution from k=0:
    // -2 * (-1) * J[3,1] = +2 * J[3,1].
    auto fc = make_two_mass_bubble();
    qft::JIntegral j("bub2m", {2, 1});
    auto s = ibp::simplify_terms(ibp::libp_deriv(fc, j, "m1sq"));
    ASSERT_EQ(s.size(), 1u);
    EXPECT_EQ(s[0].integ.indices(), (std::vector<long>{3, 1}));
    EXPECT_EQ(s[0].coef.to_string(), "2");
}

TEST(LibpDerivTest, LibpDeriv_TwoMassBubble_AllInvariantsCommute) {
    // Cross-check: d/d(m1sq) followed by d/d(m2sq) of J[1,1] equals
    // d/d(m2sq) followed by d/d(m1sq) (partial derivatives commute).
    // This exercises compute_derivative over a multi-invariant family.
    auto fc = make_two_mass_bubble();
    qft::JIntegral j("bub2m", {1, 1});

    auto first_m1_then_m2 = ibp::simplify_terms(
        ibp::compute_derivative(fc,
            ibp::libp_deriv(fc, j, "m1sq"),
            "m2sq"));
    auto first_m2_then_m1 = ibp::simplify_terms(
        ibp::compute_derivative(fc,
            ibp::libp_deriv(fc, j, "m2sq"),
            "m1sq"));

    auto sort_terms = [](std::vector<ibp::DerivTerm>& v) {
        std::sort(v.begin(), v.end(),
                  [](const ibp::DerivTerm& a, const ibp::DerivTerm& b) {
                      return a.integ.indices() < b.integ.indices();
                  });
    };
    sort_terms(first_m1_then_m2);
    sort_terms(first_m2_then_m1);

    ASSERT_EQ(first_m1_then_m2.size(), first_m2_then_m1.size());
    for (std::size_t i = 0; i < first_m1_then_m2.size(); ++i) {
        EXPECT_EQ(first_m1_then_m2[i].integ.indices(),
                  first_m2_then_m1[i].integ.indices());
        EXPECT_EQ(first_m1_then_m2[i].coef.to_string(),
                  first_m2_then_m1[i].coef.to_string());
    }
}

TEST(LibpDerivTest, SimplifyTerms_DropsZeroCoefAndCombinesLikeJ) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    std::vector<ibp::DerivTerm> terms;
    {
        ibp::DerivTerm t;
        t.coef = alg::Mfrac::from_si(fc.ctx, 1);
        t.integ = qft::JIntegral("bubble", {1, 1});
        terms.push_back(std::move(t));
    }
    {
        ibp::DerivTerm t;
        t.coef = alg::Mfrac::from_si(fc.ctx, -1);
        t.integ = qft::JIntegral("bubble", {1, 1});
        terms.push_back(std::move(t));
    }
    {
        ibp::DerivTerm t;
        t.coef = alg::Mfrac::from_si(fc.ctx, 5);
        t.integ = qft::JIntegral("bubble", {2, 1});
        terms.push_back(std::move(t));
    }

    auto s = ibp::simplify_terms(std::move(terms));
    ASSERT_EQ(s.size(), 1u);
    EXPECT_EQ(s[0].integ.indices(), (std::vector<long>{2, 1}));
    EXPECT_EQ(s[0].coef.to_string(), "5");
}
