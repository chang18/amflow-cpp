// SPDX-License-Identifier: MIT
// amflow_cli — argv parser + file I/O around api::run_json.
//
// All JSON dispatch logic lives in src/api/run_json.cpp; this binary is
// just the executable wrapper.

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include <nlohmann/json.hpp>

#include "amflow/api/run_json.hpp"
#include "amflow/numeric/options.hpp"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: amflow_cli input.json [output.json]\n";
        return 2;
    }

    std::string out_path = (argc >= 3) ? argv[2] : "";

    amflow::numeric::set_default_options();

    nlohmann::json input;
    try {
        std::ifstream in(argv[1]);
        if (!in) {
            std::cerr << "error: cannot open " << argv[1] << "\n";
            return 2;
        }
        in >> input;
    } catch (const std::exception& e) {
        std::cerr << "error: failed to parse JSON: " << e.what() << "\n";
        return 2;
    }

    try {
        auto out = amflow::api::run_json(input);

        if (out_path.empty()) {
            std::cout << out.dump(2) << "\n";
        } else {
            std::ofstream of(out_path);
            if (!of) {
                std::cerr << "error: cannot write to " << out_path << "\n";
                return 2;
            }
            of << out.dump(2) << "\n";
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
