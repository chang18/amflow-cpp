// SPDX-License-Identifier: MIT
// Tests for amflow::ibp::kira (yaml writers + parser + readers).
//
// kira_run is exercised separately in a (currently absent) integration test.

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <sstream>

#include "amflow/algebra/mpoly.hpp"
#include "amflow/ibp/kira.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/jintegral.hpp"

namespace alg = amflow::algebra;
namespace ibp = amflow::ibp;
namespace qft = amflow::qft;
namespace fs  = std::filesystem;

namespace {

std::string read_file(const std::string& path) {
    std::ifstream f(path);
    if (!f) return "";
    std::ostringstream oss;
    oss << f.rdbuf();
    return oss.str();
}

}  // namespace

TEST(KiraTest, ParseExpression_IntegerLiteral) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"d", "s"});
    EXPECT_EQ(ibp::kira_parse_expression(ctx, "5").to_string(), "5");
    EXPECT_EQ(ibp::kira_parse_expression(ctx, "-7").to_string(), "-7");
    EXPECT_EQ(ibp::kira_parse_expression(ctx, "(3)").to_string(), "3");
}

TEST(KiraTest, ParseExpression_SingleVariable) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"d", "s"});
    EXPECT_EQ(ibp::kira_parse_expression(ctx, "d").to_string(), "d");
    EXPECT_EQ(ibp::kira_parse_expression(ctx, "-d").to_string(), "-d");
}

TEST(KiraTest, ParseExpression_Arithmetic) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"d", "s"});
    EXPECT_EQ(ibp::kira_parse_expression(ctx, "d + 2").to_string(), "d+2");
    EXPECT_EQ(ibp::kira_parse_expression(ctx, "d - 4").to_string(), "d-4");
    EXPECT_EQ(ibp::kira_parse_expression(ctx, "2 * d").to_string(), "2*d");
}

TEST(KiraTest, ParseExpression_Power) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"d", "s"});
    EXPECT_EQ(ibp::kira_parse_expression(ctx, "d^2").to_string(), "d^2");
    EXPECT_EQ(ibp::kira_parse_expression(ctx, "(d+1)^3").to_string(),
              "d^3+3*d^2+3*d+1");
}

TEST(KiraTest, ParseExpression_NestedParens) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"d", "s"});
    auto r = ibp::kira_parse_expression(ctx, "(d-4)*(d-2)/(s+1)");
    alg::Mpoly num = r.numerator();
    alg::Mpoly den = r.denominator();
    EXPECT_FALSE(num.is_zero());
    EXPECT_FALSE(den.is_one());
}

TEST(KiraTest, ParseExpression_DivisionAndUnaryMinus) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"d"});
    auto r = ibp::kira_parse_expression(ctx, "-3/(d-4)");
    alg::Mfrac expected = alg::Mfrac::from_si(ctx, -3) /
                          (alg::Mfrac::from_mpoly(alg::Mpoly::variable(ctx, 0))
                             - alg::Mfrac::from_si(ctx, 4));
    EXPECT_EQ(r.to_string(), expected.to_string());
}

TEST(KiraTest, ParseExpression_UnknownIdentifierThrows) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"d"});
    EXPECT_THROW(ibp::kira_parse_expression(ctx, "x"), std::runtime_error);
}

TEST(KiraTest, ParseExpression_TrailingInputThrows) {
    auto ctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"d"});
    EXPECT_THROW(ibp::kira_parse_expression(ctx, "d ) extra"),
                 std::runtime_error);
}

TEST(KiraTest, WriteConfig_OneLoopBubble_HasExpectedYaml) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    ibp::KiraConfig cfg;
    cfg.fc = &fc;
    cfg.top_pattern = {1, 1};
    cfg.numeric_values["s"] = "100";
    cfg.numeric_values["msq"] = "1";
    cfg.kira_executable = "/usr/local/bin/kira";
    cfg.fermat_executable = "/usr/share/Ferl7/fer64";

    std::string dir =
        (fs::temp_directory_path() / "amflow_v2_kira_t1").string();
    fs::remove_all(dir);
    fs::create_directories(dir);

    ibp::kira_write_config(cfg, dir);

    auto fam_yaml = read_file(dir + "/config/integralfamilies.yaml");
    auto kin_yaml = read_file(dir + "/config/kinematics.yaml");

    EXPECT_NE(fam_yaml.find("name: \"bubble\""), std::string::npos)
        << fam_yaml;
    EXPECT_NE(fam_yaml.find("loop_momenta: [l]"), std::string::npos);
    EXPECT_NE(fam_yaml.find("top_level_sectors: [3]"), std::string::npos)
        << fam_yaml;
    EXPECT_NE(kin_yaml.find("incoming_momenta: [p, bubbleAuxLeg]"),
              std::string::npos) << kin_yaml;
    EXPECT_NE(kin_yaml.find("scalarproduct_rules:"), std::string::npos);
    EXPECT_NE(kin_yaml.find("[p,p]"), std::string::npos) << kin_yaml;

    fs::remove_all(dir);
}

TEST(KiraTest, WriteConfig_NoLegs_NoConservation) {
    auto fc = qft::FamilyConfig::build(
        "tad", {"l"}, {}, {}, {}, {"l^2 - msq"});

    ibp::KiraConfig cfg;
    cfg.fc = &fc;
    cfg.top_pattern = {1};
    cfg.kira_executable = "/usr/local/bin/kira";

    std::string dir =
        (fs::temp_directory_path() / "amflow_v2_kira_t2").string();
    fs::remove_all(dir);
    fs::create_directories(dir);

    ibp::kira_write_config(cfg, dir);

    auto kin_yaml = read_file(dir + "/config/kinematics.yaml");
    EXPECT_NE(kin_yaml.find("incoming_momenta: []"), std::string::npos);
    EXPECT_NE(kin_yaml.find("momentum_conservation: []"), std::string::npos);

    fs::remove_all(dir);
}

TEST(KiraTest, WriteJobs_MastersMode) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    ibp::KiraConfig cfg;
    cfg.fc = &fc;
    cfg.top_pattern = {1, 1};
    cfg.ibp_rank = 5;
    cfg.ibp_dot = 0;

    std::string dir =
        (fs::temp_directory_path() / "amflow_v2_kira_t3").string();
    fs::remove_all(dir);
    fs::create_directories(dir);

    ibp::kira_write_jobs(cfg, dir, ibp::KiraReductionMode::Masters);

    auto jobs = read_file(dir + "/jobs.yaml");
    EXPECT_NE(jobs.find("run_initiate: masters"), std::string::npos) << jobs;
    EXPECT_EQ(jobs.find("kira2math"), std::string::npos)
        << "Masters mode should not produce kira2math";

    fs::remove_all(dir);
}

TEST(KiraTest, WriteJobs_ReduceMode) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    ibp::KiraConfig cfg;
    cfg.fc = &fc;
    cfg.top_pattern = {1, 1};

    std::string dir =
        (fs::temp_directory_path() / "amflow_v2_kira_t4").string();
    fs::remove_all(dir);
    fs::create_directories(dir);

    ibp::kira_write_jobs(cfg, dir, ibp::KiraReductionMode::Reduce);

    auto jobs = read_file(dir + "/jobs.yaml");
    EXPECT_NE(jobs.find("run_initiate: true"), std::string::npos);
    EXPECT_NE(jobs.find("run_triangular: true"), std::string::npos);
    EXPECT_NE(jobs.find("run_back_substitution: true"), std::string::npos);
    EXPECT_NE(jobs.find("kira2math"), std::string::npos);

    fs::remove_all(dir);
}

TEST(KiraTest, WritePreferredAndTargets) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    std::string dir =
        (fs::temp_directory_path() / "amflow_v2_kira_t5").string();
    fs::remove_all(dir);
    fs::create_directories(dir);

    std::vector<qft::JIntegral> prefs;
    prefs.push_back(qft::JIntegral("bubble", {1, 1}));
    prefs.push_back(qft::JIntegral("bubble", {2, 1}));
    ibp::kira_write_preferred(prefs, fc, dir);

    auto pref_str = read_file(dir + "/preferred");
    EXPECT_NE(pref_str.find("bubble[1, 1]"), std::string::npos);
    EXPECT_NE(pref_str.find("bubble[2, 1]"), std::string::npos);

    std::vector<qft::JIntegral> targets;
    targets.push_back(qft::JIntegral("bubble", {1, 1}));
    ibp::kira_write_targets(targets, fc, dir);

    auto tgt_str = read_file(dir + "/target");
    EXPECT_NE(tgt_str.find("bubble[1, 1]"), std::string::npos);

    fs::remove_all(dir);
}

TEST(KiraTest, ReadMasters_Synthetic) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    std::string dir =
        (fs::temp_directory_path() / "amflow_v2_kira_t6").string();
    fs::remove_all(dir);
    fs::create_directories(dir + "/results/bubble");

    {
        std::ofstream f(dir + "/results/bubble/masters");
        f << "bubble(1,1)  some_extra_tag\n"
          << "\n"
          << "bubble(0,1)  another_tag\n";
    }

    auto masters = ibp::kira_read_masters(fc, dir);
    ASSERT_EQ(masters.size(), 2u);
    EXPECT_EQ(masters[0].family(), "bubble");
    EXPECT_EQ(masters[0].indices(), (std::vector<long>{1, 1}));
    EXPECT_EQ(masters[1].indices(), (std::vector<long>{0, 1}));

    fs::remove_all(dir);
}

TEST(KiraTest, ReadMasters_MissingFileReturnsEmpty) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});
    std::string dir =
        (fs::temp_directory_path() / "amflow_v2_kira_t7").string();
    fs::remove_all(dir);
    auto masters = ibp::kira_read_masters(fc, dir);
    EXPECT_TRUE(masters.empty());
}

TEST(KiraTest, ReadTargetTable_Synthetic) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    std::string dir =
        (fs::temp_directory_path() / "amflow_v2_kira_t8").string();
    fs::remove_all(dir);
    fs::create_directories(dir + "/results/bubble");

    {
        std::ofstream f(dir + "/results/bubble/kira_target.m");
        f << "{\n"
          << "  bubble(2,1) -> (d-3)/(2*s) * bubble(1,1) "
          << "- 1/(2*s) * bubble(0,1)\n"
          << "}\n";
    }

    auto rules = ibp::kira_read_target_table(fc, dir);
    ASSERT_EQ(rules.size(), 1u);
    EXPECT_EQ(rules[0].lhs.indices(), (std::vector<long>{2, 1}));
    ASSERT_EQ(rules[0].rhs.size(), 2u);
    EXPECT_EQ(rules[0].rhs[0].second.indices(),
              (std::vector<long>{1, 1}));
    EXPECT_EQ(rules[0].rhs[1].second.indices(),
              (std::vector<long>{0, 1}));

    auto cctx = std::make_shared<alg::MpolyContext>(
        std::vector<std::string>{"d", "s"});
    auto c0 = ibp::kira_parse_expression(cctx, rules[0].rhs[0].first);
    auto c1 = ibp::kira_parse_expression(cctx, rules[0].rhs[1].first);
    EXPECT_FALSE(c0.is_zero());
    EXPECT_FALSE(c1.is_zero());

    fs::remove_all(dir);
}

TEST(KiraTest, ReadTargetTable_HandlesMultipleRules) {
    auto fc = qft::FamilyConfig::build(
        "bubble", {"l"}, {"p"}, {}, {{"p^2", "s"}},
        {"l^2 - msq", "(l - p)^2 - msq"});

    std::string dir =
        (fs::temp_directory_path() / "amflow_v2_kira_t9").string();
    fs::remove_all(dir);
    fs::create_directories(dir + "/results/bubble");

    {
        std::ofstream f(dir + "/results/bubble/kira_target.m");
        f << "{\n"
          << "  bubble(2,1) -> 3 * bubble(1,1),\n"
          << "  bubble(1,2) -> -2 * bubble(0,1) + d * bubble(1,1)\n"
          << "}\n";
    }

    auto rules = ibp::kira_read_target_table(fc, dir);
    ASSERT_EQ(rules.size(), 2u);
    EXPECT_EQ(rules[1].rhs.size(), 2u);
    EXPECT_EQ(rules[1].rhs[1].first, "d");
}
