// SPDX-License-Identifier: MIT
// Tests for amflow::pipeline::amfsystem.

#include <algorithm>
#include <filesystem>

#include <gtest/gtest.h>

#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/arf.h>
#include <flint/fmpq.h>

#include "amflow/pipeline/amfsystem.hpp"
#include "amflow/numeric/acb_wrap.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/jintegral.hpp"
#include "amflow/qft/topology.hpp"

namespace pipeline = amflow::pipeline;
namespace qft     = amflow::qft;
namespace numeric = amflow::numeric;
namespace fs      = std::filesystem;

namespace {

bool kira_available() {
    return fs::exists("/usr/local/bin/kira") &&
           fs::exists("/usr/share/Ferl7/fer64");
}

}  // namespace

TEST(AmfsystemTest, EndingQ_SingleMassScheme_GeneralVacuumStopsImmediately) {
    auto fc = qft::FamilyConfig::build(
        "vacuum2", {"l1", "l2"}, {}, {}, {},
        {"l1^2 - 1", "l2^2 - 2"});

    pipeline::AMFSystemOptions opts;
    std::vector<qft::JIntegral> preferred{qft::JIntegral("vacuum2", {1, 1})};

    EXPECT_TRUE(pipeline::ending_q(fc, preferred,
                                   pipeline::EndingScheme::SingleMass, opts));
}

TEST(AmfsystemTest, CutkoskyEndingQ_PhaseVolumeComponentSelectsScheme) {
    auto fc = qft::FamilyConfig::build(
        "cutbubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2", "(l+p)^2"},
        /*cut=*/{1, 1});

    pipeline::AMFSystemOptions opts;
    std::vector<qft::JIntegral> preferred{qft::JIntegral("cutbubble", {1, 1})};

    EXPECT_TRUE(pipeline::ending_q(fc, preferred,
                                   pipeline::EndingScheme::Tradition, opts));
    EXPECT_FALSE(pipeline::ending_q(fc, preferred,
                                    pipeline::EndingScheme::Cutkosky, opts));
}

TEST(AmfsystemTest, CutkoskySetupMaster_ClearsCutAndPrescriptionForSystem) {
    auto fc = qft::FamilyConfig::build(
        "cutbubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2", "(l+p)^2"},
        /*cut=*/{1, 1});

    pipeline::AMFSystemOptions opts;
    opts.ending_schemes = {
        pipeline::EndingScheme::Tradition,
        pipeline::EndingScheme::Cutkosky,
        pipeline::EndingScheme::SingleMass
    };
    std::vector<qft::JIntegral> preferred{qft::JIntegral("cutbubble", {1, 1})};

    auto systems = pipeline::amf_system_setup_master(fc, preferred, opts);
    ASSERT_EQ(systems.size(), 1u);

    const auto& root = *systems[0];
    EXPECT_FALSE(root.is_ending());
    EXPECT_EQ(root.ending_scheme(), pipeline::EndingScheme::Cutkosky);
    ASSERT_EQ(root.family().cut.size(), 2u);
    EXPECT_EQ(root.family().cut[0], 0);
    EXPECT_EQ(root.family().cut[1], 0);
    ASSERT_EQ(root.family().prescription.size(), 1u);
    EXPECT_EQ(root.family().prescription[0], 1);
    EXPECT_EQ(std::count(root.etac().begin(), root.etac().end(), -1), 1);
}

TEST(AmfsystemTest, SingleMassSetup_TadpoleStrictlyShrinksToEndingSystem) {
    auto fc = qft::FamilyConfig::build(
        "tad", {"l"}, {}, {}, {}, {"l^2 - 1"});

    pipeline::AMFSystemOptions opts;
    opts.ending_schemes = {pipeline::EndingScheme::SingleMass};
    std::vector<qft::JIntegral> preferred{qft::JIntegral("tad", {2})};

    auto systems = pipeline::amf_system_setup_master(fc, preferred, opts);
    ASSERT_EQ(systems.size(), 1u);

    const auto& child = *systems[0];
    EXPECT_TRUE(child.is_ending());
    EXPECT_LT(child.family().n_loops(), fc.n_loops());
    EXPECT_LT(qft::get_top_position(child.preferred()).size(),
              qft::get_top_position(preferred).size());
}

TEST(AmfsystemTest, SingleMassSetup_TranslatesComponentPropagatorsByOwnContext) {
    auto fc = qft::FamilyConfig::build(
        "smctx", {"l1", "l2", "l3"}, {"n"}, {}, {{"n^2", "-1"}},
        {"l1^2 - 1", "l2^2", "l3^2", "(l1 + l2 + l3)^2"});

    pipeline::AMFSystemOptions opts;
    opts.ending_schemes = {pipeline::EndingScheme::SingleMass};
    std::vector<qft::JIntegral> preferred{
        qft::JIntegral("smctx", {1, 1, 1, 1})};

    auto systems = pipeline::amf_system_setup_master(fc, preferred, opts);
    ASSERT_EQ(systems.size(), 1u);

    const auto& child_fc = systems[0]->family();
    EXPECT_EQ(child_fc.n_loops(), 2u);
    EXPECT_EQ(child_fc.n_red_legs(), 1u);
    EXPECT_EQ(child_fc.ctx->var_index("n"), -1);

    bool saw_promoted_leg = false;
    for (const auto& prop : child_fc.propagators_after_conservation) {
        const std::string s = prop.to_string();
        EXPECT_EQ(s.find("n"), std::string::npos);
        if (s.find("l1") != std::string::npos) saw_promoted_leg = true;
    }
    EXPECT_TRUE(saw_promoted_leg);
}

TEST(AmfsystemTest, EndingTadpole_AtEps0p1) {
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
    fmpq_set_si(q, 1, 10);
    acb_set_fmpq(eps.raw(), q, 200);
    fmpq_clear(q);

    std::vector<numeric::AcbValue> epslist;
    epslist.push_back(std::move(eps));

    sys.solve(epslist);

    ASSERT_EQ(sys.solutions().size(), 1u);
    ASSERT_EQ(sys.solutions()[0].master_values.size(), 1u);

    arf_t mid;
    arf_init(mid);
    arf_set(mid, arb_midref(acb_realref(
        sys.solutions()[0].master_values[0].raw())));
    double got = arf_get_d(mid, ARF_RND_NEAR);
    arf_clear(mid);

    EXPECT_NEAR(got, 10.5705641, 1e-6);
}

TEST(AmfsystemTest, BubblePipeline_RunsToCompletion) {
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
                         "amflow_v2_amfsys_bubble").string();
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
}
