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
