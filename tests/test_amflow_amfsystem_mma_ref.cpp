// SPDX-License-Identifier: MIT
// Layer 16 MMA-reference differential tests, ported from
// tests/test_layer16_amfsystem.cpp.

#include <array>
#include <complex>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/fmpq.h>

#include "amflow/pipeline/amfsystem.hpp"
#include "amflow/pipeline/solve_integrals.hpp"
#include "amflow/numeric/acb_wrap.hpp"
#include "amflow/numeric/options.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/jintegral.hpp"

#include "util/math_ref_util.hpp"

namespace pipeline = amflow::pipeline;
namespace numeric = amflow::numeric;
namespace qft     = amflow::qft;
namespace fs      = std::filesystem;
using amflow::refactor_test::acb_from_complex;
using amflow::refactor_test::acb_imag_mid;
using amflow::refactor_test::acb_real_mid;
using amflow::refactor_test::j_to_math_ref_key;
using amflow::refactor_test::kira_available;
using amflow::refactor_test::load_math_ref_laurent;
using amflow::refactor_test::load_math_ref_values;
using amflow::refactor_test::lookup_solution_value;
using amflow::refactor_test::source_root;

TEST(AmfsystemMmaRefTest, EndingTadpole_J11_MatchesMathReferenceJson) {
    auto fc = qft::FamilyConfig::build(
        "tad", {"l"}, {}, {}, {}, {"l^2 - msq"});

    pipeline::AMFSystemOptions opts;
    opts.bb.numeric_values["msq"] = "1";

    std::vector<qft::JIntegral> preferred;
    preferred.push_back(qft::JIntegral("tad", {1}));
    std::vector<int> zero_etac(1, 0);

    pipeline::AMFSystem sys(std::move(fc), preferred, zero_etac,
                            pipeline::EndingScheme::Tradition, opts);
    EXPECT_TRUE(sys.is_ending());
    sys.setup();

    numeric::AcbValue eps;
    fmpq_t q;
    fmpq_init(q);
    fmpq_set_si(q, 1, 100);
    acb_set_fmpq(eps.raw(), q, 200);
    fmpq_clear(q);

    std::vector<numeric::AcbValue> epslist;
    epslist.push_back(std::move(eps));
    sys.solve(epslist);

    ASSERT_EQ(sys.solutions().size(), 1u);
    ASSERT_EQ(sys.solutions()[0].master_values.size(), preferred.size());

    const auto refs = load_math_ref_values(
        source_root() / "tests/data/math_ref/tadpole_1L.json", 1.0 / 100.0);
    for (const auto& target : preferred) {
        const auto it = refs.find(j_to_math_ref_key(target));
        ASSERT_NE(it, refs.end());

        const auto& got = lookup_solution_value(sys, sys.solutions()[0], target);
        EXPECT_NEAR(acb_real_mid(got), it->second.real(), 1e-6)
            << target.to_string() << " real part (MMA JSON ref)";
        EXPECT_NEAR(acb_imag_mid(got), it->second.imag(), 1e-6)
            << target.to_string() << " imag part (MMA JSON ref)";
    }
}

TEST(AmfsystemMmaRefTest, EndingTadpole_J2_UsesExplicitBoundaryJson) {
    auto fc = qft::FamilyConfig::build(
        "tad", {"l"}, {}, {}, {}, {"l^2 - msq"});

    pipeline::AMFSystemOptions opts;
    opts.bb.numeric_values["msq"] = "1";

    const auto refs_eps_100 = load_math_ref_values(
        source_root() / "tests/data/math_ref/tadpole_1L.json", 1.0 / 100.0);
    const auto refs_eps_50 = load_math_ref_values(
        source_root() / "tests/data/math_ref/tadpole_1L.json", 1.0 / 50.0);

    auto explicit_vals = std::make_shared<std::vector<numeric::AcbValue>>();
    explicit_vals->push_back(
        acb_from_complex(refs_eps_100.at("j[tad, 2]")));
    explicit_vals->push_back(
        acb_from_complex(refs_eps_50.at("j[tad, 2]")));
    opts.explicit_boundary["tad|2"] = explicit_vals;

    std::vector<qft::JIntegral> preferred;
    preferred.push_back(qft::JIntegral("tad", {1}));
    preferred.push_back(qft::JIntegral("tad", {2}));
    std::vector<int> zero_etac(1, 0);

    pipeline::AMFSystem sys(std::move(fc), preferred, zero_etac,
                            pipeline::EndingScheme::Tradition, opts);
    EXPECT_TRUE(sys.is_ending());
    sys.setup();

    std::vector<numeric::AcbValue> epslist;
    {
        numeric::AcbValue eps;
        fmpq_t q;
        fmpq_init(q);
        fmpq_set_si(q, 1, 100);
        acb_set_fmpq(eps.raw(), q, 200);
        fmpq_clear(q);
        epslist.push_back(std::move(eps));
    }
    {
        numeric::AcbValue eps;
        fmpq_t q;
        fmpq_init(q);
        fmpq_set_si(q, 1, 50);
        acb_set_fmpq(eps.raw(), q, 200);
        fmpq_clear(q);
        epslist.push_back(std::move(eps));
    }

    sys.solve(epslist);

    ASSERT_EQ(sys.solutions().size(), 2u);
    ASSERT_EQ(sys.solutions()[0].master_values.size(), preferred.size());
    ASSERT_EQ(sys.solutions()[1].master_values.size(), preferred.size());

    const std::array<std::map<std::string, std::complex<double>>, 2> refs = {
        refs_eps_100, refs_eps_50
    };
    for (std::size_t k = 0; k < refs.size(); ++k) {
        for (const auto& target : preferred) {
            const auto key = j_to_math_ref_key(target);
            const auto it = refs[k].find(key);
            ASSERT_NE(it, refs[k].end());

            const auto& got = lookup_solution_value(sys, sys.solutions()[k], target);
            EXPECT_NEAR(acb_real_mid(got), it->second.real(), 1e-6)
                << "eps slot " << k << " " << target.to_string()
                << " real part (vacuum + explicit_boundary mixed ending system)";
            EXPECT_NEAR(acb_imag_mid(got), it->second.imag(), 1e-6)
                << "eps slot " << k << " " << target.to_string()
                << " imag part (vacuum + explicit_boundary mixed ending system)";
        }
    }
}

// --- Ending-master Kira reduction loop -------
//
// `solve_ending_master_value` (`src/pipeline/amfsystem.cpp:718-792`)
// is the C++-only path that lowers an arbitrary ending master
// through repeated `ibp::reduce` calls until every leaf hits either
// (a) the builtin Vacuum[L,n] table or (b) an explicit_boundary
// value supplied by the user.  Upstream's `AMFSystemSetupMaster`
// aborts at the ending step unless the user supplies a `Solution`
// for every master; C++ relaxes that constraint.
//
// `EndingTadpole_J2_UsesExplicitBoundaryJson` (above) exercises the
// short-circuit path (try_builtin_ending_value finds an
// explicit_boundary for J[tad, 2]).  This test exercises the
// load-bearing recursive lowering path: it provides J[tad, 2] as
// a preferred master with NO explicit_boundary, forcing the
// solver to:
//   1. Vacuum table miss on J[2] (n_props = 0 in the all-ones
//      counter; vacuum_known(1, 0) is false).
//   2. Call get_ending_reduction → ibp::reduce on J[2].
//   3. Receive Kira's IBP identity J[2] = (ε-1)/m² * J[1].
//   4. Recurse on J[1] → Vacuum table hits (vacuum_known(1, 1)
//      is true), returning the closed-form Γ(ε-1) value.
//   5. Multiply by the IBP coefficient at ε, return J[2] numeric
//      value.
//
// Reference data is the same JSON used by
// `EndingTadpole_J2_UsesExplicitBoundaryJson`; the test asserts
// the recursive-lowering output matches the MMA Laurent
// evaluation at the same two eps points.

TEST(AmfsystemMmaRefTest, EndingTadpole_J2_RecursiveKiraLoweringMatchesReference) {
    if (!kira_available()) {
        GTEST_SKIP() << "Kira/Fermat not installed";
    }

    auto fc = qft::FamilyConfig::build(
        "tad", {"l"}, {}, {}, {}, {"l^2 - msq"});

    pipeline::AMFSystemOptions opts;
    opts.bb.kira_executable   = "/usr/local/bin/kira";
    opts.bb.fermat_executable = "/usr/share/Ferl7/fer64";
    opts.bb.numeric_values["msq"] = "1";
    opts.cache_root = (fs::temp_directory_path() /
                       "amflow_amfsys_mma_tad_recursive_lowering").string();
    fs::remove_all(opts.cache_root);

    // Deliberately do NOT set opts.explicit_boundary["tad|2"].
    // The solver must derive J[2] from J[1] via the recursive Kira
    // reduction loop.

    std::vector<qft::JIntegral> preferred;
    preferred.push_back(qft::JIntegral("tad", {1}));
    preferred.push_back(qft::JIntegral("tad", {2}));
    std::vector<int> zero_etac(1, 0);

    pipeline::AMFSystem sys(std::move(fc), preferred, zero_etac,
                            pipeline::EndingScheme::Tradition, opts);
    EXPECT_TRUE(sys.is_ending());
    sys.setup();

    const auto refs_eps_100 = load_math_ref_values(
        source_root() / "tests/data/math_ref/tadpole_1L.json", 1.0 / 100.0);
    const auto refs_eps_50 = load_math_ref_values(
        source_root() / "tests/data/math_ref/tadpole_1L.json", 1.0 / 50.0);

    std::vector<numeric::AcbValue> epslist;
    for (long denom : {100, 50}) {
        numeric::AcbValue eps;
        fmpq_t q;
        fmpq_init(q);
        fmpq_set_si(q, 1, denom);
        acb_set_fmpq(eps.raw(), q, 200);
        fmpq_clear(q);
        epslist.push_back(std::move(eps));
    }

    sys.solve(epslist);

    ASSERT_EQ(sys.solutions().size(), 2u);

    const std::array<std::map<std::string, std::complex<double>>, 2> refs = {
        refs_eps_100, refs_eps_50
    };
    for (std::size_t k = 0; k < refs.size(); ++k) {
        for (const auto& target : preferred) {
            const auto key = j_to_math_ref_key(target);
            const auto it = refs[k].find(key);
            ASSERT_NE(it, refs[k].end());

            const auto& got =
                lookup_solution_value(sys, sys.solutions()[k], target);
            EXPECT_NEAR(acb_real_mid(got), it->second.real(), 1e-6)
                << "eps slot " << k << " " << target.to_string()
                << " real part (recursive Kira-lowering path);"
                << " expected " << it->second.real()
                << ", got " << acb_real_mid(got);
            EXPECT_NEAR(acb_imag_mid(got), it->second.imag(), 1e-6)
                << "eps slot " << k << " " << target.to_string()
                << " imag part";
        }
    }
}

TEST(AmfsystemMmaRefTest, BubblePipeline_MatchesMathReferenceJson) {
    if (!kira_available()) {
        GTEST_SKIP() << "Kira/Fermat not installed";
    }

    auto fc = qft::FamilyConfig::build(
        "bubblefam", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    pipeline::AMFSystemOptions opts;
    opts.bb.kira_executable   = "/usr/local/bin/kira";
    opts.bb.fermat_executable = "/usr/share/Ferl7/fer64";
    opts.bb.numeric_values["s"]   = "5";
    opts.bb.numeric_values["msq"] = "1";
    opts.ending_schemes = {pipeline::EndingScheme::Tradition,
                            pipeline::EndingScheme::SingleMass};
    opts.cache_root = (fs::temp_directory_path() /
                         "amflow_v2_amfsys_mma_bubble").string();
    opts.max_recursion_depth = 16;
    fs::remove_all(opts.cache_root);

    std::vector<qft::JIntegral> preferred;
    preferred.push_back(qft::JIntegral("bubblefam", {1, 1}));

    auto systems = pipeline::amf_systems_setup(fc, preferred, opts);
    ASSERT_FALSE(systems.empty());

    numeric::AcbValue eps;
    fmpq_t q;
    fmpq_init(q);
    fmpq_set_si(q, 1, 100);
    acb_set_fmpq(eps.raw(), q, 200);
    fmpq_clear(q);

    std::vector<numeric::AcbValue> epslist;
    epslist.push_back(std::move(eps));

    auto sols = pipeline::amf_systems_solution(systems, epslist);
    ASSERT_EQ(sols.size(), 1u);
    ASSERT_FALSE(sols[0].master_values.empty());

    const auto refs = load_math_ref_values(
        source_root() / "tests/data/math_ref/bubble_1L.json", 1.0 / 100.0);
    const qft::JIntegral target("bubblefam", {1, 1});
    const auto it = refs.find(j_to_math_ref_key(target));
    ASSERT_NE(it, refs.end());

    const auto& got = lookup_solution_value(*systems[0], sols[0], target);
    EXPECT_NEAR(acb_real_mid(got), it->second.real(), 1e-6)
        << "bubble J[1,1] real part (MMA JSON ref)";
    EXPECT_NEAR(acb_imag_mid(got), it->second.imag(), 1e-6)
        << "bubble J[1,1] imag part (MMA JSON ref)";
}

TEST(AmfsystemMmaRefTest, SunrisePipeline_LeadingCoefficientsMatchMathReferenceJson) {
    if (!kira_available()) {
        GTEST_SKIP() << "Kira/Fermat not installed";
    }

    auto num_cfg = pipeline::generate_numerical_config(/*loop_count=*/2,
                                                       /*goal_digits=*/25,
                                                       /*eps_order=*/3);

    numeric::GlobalScope s;
    s.global.working_pre = static_cast<int>(num_cfg.working_pre);
    s.global.chop_pre = 20;
    s.global.rationalize_pre = 100;
    s.global.silent_mode = true;
    s.expansion.x_order = static_cast<int>(num_cfg.x_order);
    s.expansion.extra_x_order = static_cast<int>(num_cfg.x_order);
    s.expansion.learn_x_order = -1;
    s.expansion.test_x_order = 5;
    s.running.run_radius = 2;
    s.running.run_length = 200;
    s.running.run_candidate = 10;
    s.running.run_direction = numeric::RunningOptions::Direction::NegIm;
    s.commit();

    auto fc = qft::FamilyConfig::build(
        "sunrise", {"l1", "l2"}, {"p"}, {}, {{"p^2", "s"}},
        {"l1^2 - msq", "l2^2 - msq", "(l1 + l2 - p)^2 - msq",
         "(l1 + p)^2", "(l2 + p)^2"});

    pipeline::AMFSystemOptions opts;
    opts.bb.kira_executable   = "/usr/local/bin/kira";
    opts.bb.fermat_executable = "/usr/share/Ferl7/fer64";
    opts.bb.numeric_values["s"]   = "-3";
    opts.bb.numeric_values["msq"] = "1";
    opts.ending_schemes = {pipeline::EndingScheme::Tradition,
                            pipeline::EndingScheme::SingleMass};
    opts.cache_root = (fs::temp_directory_path() /
                         "amflow_v2_amfsys_mma_sunrise").string();
    opts.max_recursion_depth = 16;
    fs::remove_all(opts.cache_root);

    const qft::JIntegral target("sunrise", {1, 1, 1, 0, 0});
    std::vector<qft::JIntegral> preferred{target};

    auto systems = pipeline::amf_systems_setup(fc, preferred, opts);
    ASSERT_FALSE(systems.empty());

    auto& epslist = num_cfg.eps_samples;

    const auto sols = pipeline::amf_systems_solution(systems, epslist);
    ASSERT_EQ(sols.size(), epslist.size());

    std::vector<numeric::AcbValue> values;
    values.reserve(sols.size());
    for (std::size_t k = 0; k < sols.size(); ++k) {
        const auto& got = lookup_solution_value(*systems[0], sols[k], target);
        values.push_back(got.clone());
    }

    std::vector<numeric::AcbValue> fit_epslist;
    std::vector<numeric::AcbValue> fit_values;
    const std::size_t keep = 4;
    ASSERT_GE(epslist.size(), keep);
    fit_epslist.reserve(keep);
    fit_values.reserve(keep);
    for (std::size_t i = 0; i < keep; ++i) {
        fit_epslist.push_back(epslist[i].clone());
        fit_values.push_back(values[i].clone());
    }

    const auto coeffs = pipeline::fit_eps(fit_epslist, fit_values, -2,
                                          numeric::working_prec_bits());
    ASSERT_GE(coeffs.size(), 2u);
    const auto ref = load_math_ref_laurent(
        source_root() / "tests/data/math_ref/sunrise_2L.json",
        j_to_math_ref_key(target));

    const auto it_m2 = ref.find(-2);
    const auto it_m1 = ref.find(-1);
    ASSERT_NE(it_m2, ref.end());
    ASSERT_NE(it_m1, ref.end());

    EXPECT_NEAR(acb_real_mid(coeffs[0]), it_m2->second.real(), 1e-5)
        << "fitted coeff of eps^-2 for sunrise J[1,1,1,0,0]";
    EXPECT_NEAR(acb_imag_mid(coeffs[0]), it_m2->second.imag(), 1e-8)
        << "fitted imag coeff of eps^-2 for sunrise J[1,1,1,0,0]";
    EXPECT_NEAR(acb_real_mid(coeffs[1]), it_m1->second.real(), 1e-4)
        << "fitted coeff of eps^-1 for sunrise J[1,1,1,0,0]";
    EXPECT_NEAR(acb_imag_mid(coeffs[1]), it_m1->second.imag(), 1e-8)
        << "fitted imag coeff of eps^-1 for sunrise J[1,1,1,0,0]";
}
