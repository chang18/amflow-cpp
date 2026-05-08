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
