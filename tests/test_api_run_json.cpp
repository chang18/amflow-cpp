// SPDX-License-Identifier: MIT
// Tests for amflow::api::run_json.

#include <gtest/gtest.h>
#include <nlohmann/json.hpp>

#include "amflow/api/run_json.hpp"

namespace api = amflow::api;
using nlohmann::json;

TEST(ApiRunJsonTest, UnknownModeThrows) {
    json input = {{"mode", "no_such_mode"}};
    EXPECT_THROW(api::run_json(input), std::runtime_error);
}

TEST(ApiRunJsonTest, OptionsAreEchoedBack) {
    // Setting a recognised option in input should be visible in output.options.
    json input = {
        {"mode", "amflow"},
        {"options", {{"silent_mode", true}, {"working_pre", 60}}},
        // Provide a tiny 1x1 DE that solves the boundary mu=0,value=1 trivially.
        {"matrix", json::array({json::array({json{
            {"num", json::array({"0"})},
            {"den", json::array({"1"})},
        }})})},
        {"boundaries", json::array({json::array({json{
            {"mu",    json{{"re", "0"}, {"im", "0"}}},
            {"value", json{{"re", "1"}, {"im", "0"}}},
        }})})},
    };

    json out = api::run_json(input);
    EXPECT_EQ(out.at("mode"), "amflow");
    ASSERT_TRUE(out.contains("options"));
    EXPECT_EQ(out.at("options").at("silent_mode"), true);
    EXPECT_EQ(out.at("options").at("working_pre"), 60);
    ASSERT_TRUE(out.contains("result"));
    EXPECT_TRUE(out.at("result").is_array());
    EXPECT_EQ(out.at("result").size(), 1u);
}

TEST(ApiRunJsonTest, AmflowMissingMatrixThrows) {
    json input = {
        {"mode", "amflow"},
        {"boundaries", json::array()},
    };
    EXPECT_THROW(api::run_json(input), std::exception);
}

TEST(ApiRunJsonTest, AmflowBoundariesSizeMismatchThrows) {
    json input = {
        {"mode", "amflow"},
        {"matrix", json::array({json::array({json{
            {"num", json::array({"0"})},
            {"den", json::array({"1"})},
        }})})},
        // 0 boundaries for a 1x1 system.
        {"boundaries", json::array()},
    };
    EXPECT_THROW(api::run_json(input), std::runtime_error);
}

// Mirrors audit divergence D5: the upstream `IBPRule` /
// `CompensateRule` machinery for complex-valued numeric kinematics
// (`Kira/interface.m:50-57`) is not ported.  The dispatcher must
// reject the object form `{"re":..,"im":..}` in
// `amf_options.blackbox.numeric_values` with a clear error rather
// than silently truncating to the real part.  Locks down the
// loud-rejection contract so a future regression cannot silently
// re-enable a mishandled path.
TEST(ApiRunJsonTest, BlackboxComplexNumericValueIsRejected) {
    json input = {
        {"mode", "black_box_amflow"},
        {"family", {
            {"family", "fam"},
            {"loops", json::array({"l1"})},
            {"legs", json::array()},
            {"propagators", json::array({"l1*l1"})},
        }},
        {"targets", json::array({json::array({1})})},
        {"eps_samples", json::array({json{{"re", "1/1000"}, {"im", "0"}}})},
        {"amf_options", {{"blackbox", {{"numeric_values", {
            {"s", {{"re", "100"}, {"im", "1"}}},
        }}}}}},
    };
    try {
        api::run_json(input);
        FAIL() << "expected complex-form numeric_values to be rejected";
    } catch (const std::runtime_error& e) {
        const std::string what(e.what());
        EXPECT_NE(what.find("complex-numeric"), std::string::npos)
            << "error should name the unsupported feature, got: " << what;
        EXPECT_NE(what.find("D5"), std::string::npos)
            << "error should point at audit entry D5, got: " << what;
    }
}
