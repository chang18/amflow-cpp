// SPDX-License-Identifier: MIT
// Tests for amflow::pipeline::factorize.

#include <gtest/gtest.h>

#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>

#include "amflow/algebra/mpoly.hpp"
#include "amflow/pipeline/factorize.hpp"
#include "amflow/qft/family_config.hpp"

namespace alg     = amflow::algebra;
namespace pipeline = amflow::pipeline;
namespace qft     = amflow::qft;

namespace {

std::vector<alg::Mpoly>
zero_legs(const qft::FamilyConfig& fc) {
    std::vector<alg::Mpoly> props;
    long n_loop = (long)fc.n_loops();
    long n_red  = (long)fc.n_red_legs();
    fmpz_t zero;
    fmpz_init(zero);
    for (const auto& p : fc.propagators_after_conservation) {
        alg::Mpoly tmp = p.clone();
        for (long k = 0; k < n_red; ++k) {
            alg::Mpoly out(fc.ctx);
            fmpz_mpoly_evaluate_one_fmpz(out.raw(), tmp.raw(),
                                          n_loop + k, zero, fc.ctx->raw());
            tmp = std::move(out);
        }
        props.push_back(std::move(tmp));
    }
    fmpz_clear(zero);
    return props;
}

}  // namespace

TEST(FactorizeTest, OneLoopBubble_Identity) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    auto props = zero_legs(fc);
    std::vector<std::vector<long>> patts = {{1, 1}};
    auto fams = pipeline::factorize_family(fc, props, patts);

    ASSERT_EQ(fams.size(), 1u);
    EXPECT_EQ(fams[0].loops.size(), 1u);
    EXPECT_EQ(fams[0].loops[0], "l");
    EXPECT_EQ(fams[0].propagators.size(), 2u);
    EXPECT_EQ(fams[0].patterns.size(), 1u);
    EXPECT_EQ(fams[0].patterns[0].size(), 2u);
    EXPECT_EQ(fams[0].patterns[0][0], 1);
    EXPECT_EQ(fams[0].patterns[0][1], 1);
}

TEST(FactorizeTest, SortBy_PicksMonomialContainingMassOne) {
    auto fc = qft::FamilyConfig::build(
        "two_indep_loops", {"l1", "l2"}, {}, {}, {},
        {"l1^2 - 1", "l2^2"});
    std::vector<alg::Mpoly> props;
    for (const auto& p : fc.propagators_after_conservation) {
        props.push_back(p.clone());
    }
    std::vector<std::vector<long>> patts = {{1, 1}};
    auto fams = pipeline::factorize_family(fc, props, patts);
    ASSERT_EQ(fams.size(), 2u);
    for (auto& f : fams) {
        EXPECT_EQ(f.loops.size(), 1u);
        EXPECT_EQ(f.propagators.size(), 1u);
    }
}

TEST(FactorizeTest, TwoLoopSunrise_OneComponent) {
    auto fc = qft::FamilyConfig::build(
        "sunrise", {"l1", "l2"}, {"p"}, {}, {{"p^2", "s"}},
        {"l1^2 - m1sq", "l2^2 - m2sq", "(l1 + l2 - p)^2 - m3sq"});

    auto props = zero_legs(fc);
    std::vector<std::vector<long>> patts = {{1, 1, 1}};
    auto fams = pipeline::factorize_family(fc, props, patts);

    ASSERT_EQ(fams.size(), 1u);
    EXPECT_EQ(fams[0].loops.size(), 2u);
    EXPECT_EQ(fams[0].propagators.size(), 3u);
    EXPECT_EQ(fams[0].patterns.size(), 1u);
    EXPECT_EQ(fams[0].patterns[0].size(), 3u);
}
