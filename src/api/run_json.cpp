// SPDX-License-Identifier: MIT
// amflow::api::run_json — implementation.
//

#include "amflow/api/run_json.hpp"

#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>

#include "amflow/pipeline/amfsystem.hpp"
#include "amflow/pipeline/solve_integrals.hpp"
#include "amflow/ibp/reduce.hpp"
#include "amflow/numeric/acb_wrap.hpp"
#include "amflow/numeric/options.hpp"
#include "amflow/numeric/rational.hpp"
#include "amflow/ode/amflow.hpp"
#include "amflow/ode/inf.hpp"
#include "amflow/qft/amfmode.hpp"
#include "amflow/qft/family_config.hpp"
#include "amflow/qft/jintegral.hpp"

namespace amflow::api {

using nlohmann::json;
using algebra::Mfrac;
using algebra::Mpoly;
using numeric::AcbValue;
using numeric::FmpqPoly;
using numeric::RationalComplex;
using numeric::RationalFunction;
using numeric::RationalMatrix;

namespace {

std::string json_scalar_to_string(const json& value, const std::string& name) {
    if (value.is_string()) return value.get<std::string>();
    if (value.is_number_integer()) {
        return std::to_string(value.get<long long>());
    }
    if (value.is_number_unsigned()) {
        return std::to_string(value.get<unsigned long long>());
    }
    if (value.is_number_float()) {
        std::ostringstream os;
        os << value.get<double>();
        return os.str();
    }
    throw std::runtime_error(name + " must be a string or number");
}

std::vector<std::string> parse_string_array(const json& arr,
                                              const std::string& name) {
    if (!arr.is_array()) {
        throw std::runtime_error(name + " must be an array");
    }
    std::vector<std::string> out;
    out.reserve(arr.size());
    for (std::size_t i = 0; i < arr.size(); ++i) {
        if (!arr[i].is_string()) {
            throw std::runtime_error(name + "[" + std::to_string(i)
                                       + "] must be a string");
        }
        out.push_back(arr[i].get<std::string>());
    }
    return out;
}

std::vector<int> parse_int_array(const json& arr, const std::string& name) {
    if (!arr.is_array()) {
        throw std::runtime_error(name + " must be an array");
    }
    std::vector<int> out;
    out.reserve(arr.size());
    for (std::size_t i = 0; i < arr.size(); ++i) {
        if (!arr[i].is_number_integer()) {
            throw std::runtime_error(name + "[" + std::to_string(i)
                                       + "] must be an integer");
        }
        out.push_back(arr[i].get<int>());
    }
    return out;
}

std::vector<std::pair<std::string, std::string>>
parse_rule_list(const json& node, const std::string& name) {
    std::vector<std::pair<std::string, std::string>> out;
    if (node.is_null()) return out;

    if (node.is_object()) {
        for (auto it = node.begin(); it != node.end(); ++it) {
            out.emplace_back(it.key(),
                              json_scalar_to_string(it.value(),
                                                     name + "." + it.key()));
        }
        return out;
    }

    if (!node.is_array()) {
        throw std::runtime_error(name + " must be an object or array");
    }

    out.reserve(node.size());
    for (std::size_t i = 0; i < node.size(); ++i) {
        const auto& entry = node[i];
        if (entry.is_array()) {
            if (entry.size() != 2 || !entry[0].is_string()) {
                throw std::runtime_error(name + "[" + std::to_string(i)
                                           + "] must be [lhs, rhs]");
            }
            out.emplace_back(entry[0].get<std::string>(),
                              json_scalar_to_string(entry[1],
                                                     name + "[" + std::to_string(i) + "][1]"));
            continue;
        }
        if (!entry.is_object() || !entry.contains("lhs") || !entry.contains("rhs")
                || !entry.at("lhs").is_string()) {
            throw std::runtime_error(name + "[" + std::to_string(i)
                                       + "] must contain string lhs/rhs");
        }
        out.emplace_back(entry.at("lhs").get<std::string>(),
                          json_scalar_to_string(entry.at("rhs"),
                                                 name + "[" + std::to_string(i) + "].rhs"));
    }

    return out;
}

long parse_required_long(const json& obj,
                            const std::string& primary,
                            const std::string& alias = "") {
    if (obj.contains(primary)) {
        if (!obj.at(primary).is_number_integer()) {
            throw std::runtime_error(primary + " must be an integer");
        }
        return obj.at(primary).get<long>();
    }
    if (!alias.empty() && obj.contains(alias)) {
        if (!obj.at(alias).is_number_integer()) {
            throw std::runtime_error(alias + " must be an integer");
        }
        return obj.at(alias).get<long>();
    }
    throw std::runtime_error("missing required integer field '" + primary + "'");
}

qft::FamilyConfig parse_family_config(const json& obj) {
    if (!obj.is_object()) {
        throw std::runtime_error("family must be an object");
    }

    const std::string name = obj.contains("name")
        ? obj.at("name").get<std::string>()
        : obj.at("family").get<std::string>();
    const auto loops = parse_string_array(obj.at("loops"), "family.loops");
    const auto legs = parse_string_array(obj.at("legs"), "family.legs");
    const auto conservation = obj.contains("conservation")
        ? parse_rule_list(obj.at("conservation"), "family.conservation")
        : std::vector<std::pair<std::string, std::string>>{};
    const auto replacement = obj.contains("replacement")
        ? parse_rule_list(obj.at("replacement"), "family.replacement")
        : std::vector<std::pair<std::string, std::string>>{};
    const auto propagators =
        parse_string_array(obj.at("propagators"), "family.propagators");
    const auto cut = obj.contains("cut")
        ? parse_int_array(obj.at("cut"), "family.cut")
        : std::vector<int>{};
    const auto prescription = obj.contains("prescription")
        ? parse_int_array(obj.at("prescription"), "family.prescription")
        : std::vector<int>{};

    return qft::FamilyConfig::build(name, loops, legs, conservation, replacement,
                                       propagators, cut, prescription);
}

qft::JIntegral parse_jintegral(const json& entry,
                                  const qft::FamilyConfig& fc,
                                  const std::string& name) {
    std::string family = fc.family;
    const json* idx_node = &entry;

    if (entry.is_object()) {
        if (entry.contains("family")) {
            if (!entry.at("family").is_string()) {
                throw std::runtime_error(name + ".family must be a string");
            }
            family = entry.at("family").get<std::string>();
        }
        if (!entry.contains("indices")) {
            throw std::runtime_error(name + " must contain 'indices'");
        }
        idx_node = &entry.at("indices");
    }

    if (!idx_node->is_array()) {
        throw std::runtime_error(name + ".indices must be an array");
    }

    std::vector<long> indices;
    indices.reserve(idx_node->size());
    for (std::size_t i = 0; i < idx_node->size(); ++i) {
        if (!(*idx_node)[i].is_number_integer()) {
            throw std::runtime_error(name + ".indices[" + std::to_string(i)
                                       + "] must be an integer");
        }
        indices.push_back((*idx_node)[i].get<long>());
    }

    if (indices.size() != fc.n_propagators()) {
        throw std::runtime_error(name + " has " + std::to_string(indices.size())
                                   + " indices, expected "
                                   + std::to_string(fc.n_propagators()));
    }
    return qft::JIntegral(family, std::move(indices));
}

std::vector<qft::JIntegral> parse_jintegral_list(const json& obj,
                                                    const qft::FamilyConfig& fc) {
    const json* arr = nullptr;
    if (obj.contains("integrals")) arr = &obj.at("integrals");
    else if (obj.contains("jints")) arr = &obj.at("jints");
    else if (obj.contains("targets")) arr = &obj.at("targets");
    else throw std::runtime_error("missing required array 'integrals'");

    if (!arr->is_array()) {
        throw std::runtime_error("integrals must be an array");
    }

    std::vector<qft::JIntegral> out;
    out.reserve(arr->size());
    for (std::size_t i = 0; i < arr->size(); ++i) {
        out.push_back(parse_jintegral((*arr)[i], fc,
                                         "integrals[" + std::to_string(i) + "]"));
    }
    return out;
}

pipeline::EndingScheme parse_ending_scheme_string(const std::string& s) {
    if (s == "Tradition") return pipeline::EndingScheme::Tradition;
    if (s == "Cutkosky") return pipeline::EndingScheme::Cutkosky;
    if (s == "SingleMass") return pipeline::EndingScheme::SingleMass;
    if (s == "Trivial") return pipeline::EndingScheme::Trivial;
    throw std::runtime_error("unknown ending scheme '" + s + "'");
}

void apply_blackbox_options(const json& obj, ibp::ReduceOptions& bb) {
    if (obj.contains("ibp_rank")) bb.ibp_rank = obj.at("ibp_rank").get<long>();
    if (obj.contains("ibp_dot")) bb.ibp_dot = obj.at("ibp_dot").get<long>();
    if (obj.contains("n_thread")) bb.n_thread = obj.at("n_thread").get<long>();
    if (obj.contains("integral_order")) {
        bb.integral_order = obj.at("integral_order").get<long>();
    }
    if (obj.contains("kira_executable")) {
        bb.kira_executable = obj.at("kira_executable").get<std::string>();
    }
    if (obj.contains("fermat_executable")) {
        bb.fermat_executable = obj.at("fermat_executable").get<std::string>();
    }
    if (obj.contains("work_dir")) {
        bb.work_dir = obj.at("work_dir").get<std::string>();
    }
    if (obj.contains("log_file")) {
        bb.log_file = obj.at("log_file").get<std::string>();
    }
    if (obj.contains("numeric_values")) {
        const auto& nv = obj.at("numeric_values");
        if (!nv.is_object()) {
            throw std::runtime_error(
                "amf_options.blackbox.numeric_values must be an object");
        }
        for (auto it = nv.begin(); it != nv.end(); ++it) {
            bb.numeric_values[it.key()] =
                json_scalar_to_string(it.value(),
                                       "amf_options.blackbox.numeric_values." + it.key());
        }
    }
}

pipeline::AMFSystemOptions parse_amf_options(const json& obj) {
    if (!obj.is_object()) {
        throw std::runtime_error("amf_options must be an object");
    }

    pipeline::AMFSystemOptions out;
    out.bb.kira_executable = "/usr/local/bin/kira";
    out.bb.fermat_executable = "/usr/share/Ferl7/fer64";

    if (obj.contains("amf_modes")) {
        out.amf_modes.clear();
        for (const auto& mode : obj.at("amf_modes")) {
            if (!mode.is_string()) {
                throw std::runtime_error("amf_options.amf_modes entries must be strings");
            }
            out.amf_modes.push_back(qft::parse_amf_mode(mode.get<std::string>()));
        }
    }

    if (obj.contains("ending_schemes")) {
        out.ending_schemes.clear();
        for (const auto& scheme : obj.at("ending_schemes")) {
            if (!scheme.is_string()) {
                throw std::runtime_error(
                    "amf_options.ending_schemes entries must be strings");
            }
            out.ending_schemes.push_back(
                parse_ending_scheme_string(scheme.get<std::string>()));
        }
    }

    if (obj.contains("cache_root")) {
        out.cache_root = obj.at("cache_root").get<std::string>();
    }
    if (obj.contains("direction")) {
        out.direction = obj.at("direction").get<std::string>();
    }
    if (obj.contains("max_recursion_depth")) {
        out.max_recursion_depth = obj.at("max_recursion_depth").get<long>();
    }

    if (obj.contains("blackbox")) apply_blackbox_options(obj.at("blackbox"), out.bb);
    if (obj.contains("bb")) apply_blackbox_options(obj.at("bb"), out.bb);
    return out;
}

pipeline::AMFSystemOptions default_solve_integrals_options() {
    pipeline::AMFSystemOptions out;
    out.bb.kira_executable = "/usr/local/bin/kira";
    out.bb.fermat_executable = "/usr/share/Ferl7/fer64";
    return out;
}

void ensure_cli_blackbox_log_file(pipeline::AMFSystemOptions& opts) {
    if (opts.bb.log_file.empty()) {
        opts.bb.log_file = "/tmp/amflow_cli_kira.log";
    }
}

json jintegral_to_json(const qft::JIntegral& j) {
    return json{
        {"family", j.family()},
        {"indices", j.indices()},
    };
}

bool parse_fmpq(fmpq_t out, const std::string& s) {
    int rc = fmpq_set_str(out, s.c_str(), 10);
    return rc == 0;
}

void load_real_component_from_json(arb_t out,
                                     const json& value,
                                     const std::string& name,
                                     long prec) {
    const std::string text = json_scalar_to_string(value, name);

    fmpq_t q;
    fmpq_init(q);
    if (parse_fmpq(q, text)) {
        arb_set_fmpq(out, q, prec);
        fmpq_clear(q);
        return;
    }
    fmpq_clear(q);

    if (arb_set_str(out, text.c_str(), prec) == 0) {
        return;
    }
    throw std::runtime_error("failed to parse real value '" + text
                               + "' for " + name);
}

FmpqPoly parse_coeff_list(const json& arr) {
    FmpqPoly p;
    long n = 0;
    fmpq_t q;
    fmpq_init(q);
    for (const auto& c : arr) {
        std::string s;
        if (c.is_string()) s = c.get<std::string>();
        else if (c.is_number_integer()) s = std::to_string(c.get<long>());
        else if (c.is_number_float()) s = std::to_string(c.get<double>());
        else throw std::runtime_error("coefficient must be a string or integer");
        if (!parse_fmpq(q, s)) {
            fmpq_clear(q);
            throw std::runtime_error("failed to parse fmpq from '" + s + "'");
        }
        p.set_coeff_fmpq(n++, q);
    }
    fmpq_clear(q);
    return p;
}

RationalFunction parse_rational(const json& obj) {
    if (obj.contains("coeffs")) {
        return RationalFunction::from_polynomial(parse_coeff_list(obj.at("coeffs")));
    }
    if (!obj.contains("num") || !obj.contains("den")) {
        throw std::runtime_error(
            "rational entry must have 'num' and 'den' (or just 'coeffs')");
    }
    FmpqPoly num = parse_coeff_list(obj.at("num"));
    FmpqPoly den = parse_coeff_list(obj.at("den"));
    return RationalFunction(std::move(num), std::move(den));
}

RationalMatrix parse_matrix(const json& arr) {
    if (!arr.is_array() || arr.empty()) {
        throw std::runtime_error("matrix must be a non-empty array of rows");
    }
    std::size_t rows = arr.size();
    std::size_t cols = arr[0].size();
    if (cols != rows) {
        throw std::runtime_error("matrix must be square");
    }
    RationalMatrix out(rows, cols);
    for (std::size_t i = 0; i < rows; ++i) {
        const auto& row = arr[i];
        if (row.size() != cols) {
            throw std::runtime_error("matrix is not rectangular");
        }
        for (std::size_t j = 0; j < cols; ++j) {
            out(i, j) = parse_rational(row[j]);
        }
    }
    return out;
}

AcbValue parse_complex(const json& obj, long prec) {
    AcbValue out;
    if (obj.is_object()) {
        load_real_component_from_json(acb_realref(out.raw()), obj.at("re"),
                                         "complex.re", prec);
        if (obj.contains("im")) {
            load_real_component_from_json(acb_imagref(out.raw()), obj.at("im"),
                                             "complex.im", prec);
        } else {
            arb_zero(acb_imagref(out.raw()));
        }
        return out;
    }

    load_real_component_from_json(acb_realref(out.raw()), obj, "scalar", prec);
    arb_zero(acb_imagref(out.raw()));
    return out;
}

std::vector<AcbValue> parse_eps_samples(const json& obj, long prec) {
    const json* arr = nullptr;
    if (obj.contains("eps_samples")) arr = &obj.at("eps_samples");
    else if (obj.contains("epslist")) arr = &obj.at("epslist");
    else throw std::runtime_error("missing required array 'eps_samples'");

    if (!arr->is_array()) {
        throw std::runtime_error("eps_samples must be an array");
    }

    std::vector<AcbValue> out;
    out.reserve(arr->size());
    for (std::size_t i = 0; i < arr->size(); ++i) {
        out.push_back(parse_complex((*arr)[i], prec));
    }
    return out;
}

std::vector<ode::BoundarySpec> parse_boundaries(const json& arr, long prec) {
    if (!arr.is_array()) throw std::runtime_error("boundaries must be an array");
    std::vector<ode::BoundarySpec> out(arr.size());
    for (std::size_t i = 0; i < arr.size(); ++i) {
        const auto& list = arr[i];
        if (!list.is_array()) throw std::runtime_error("each boundary must be an array");
        for (const auto& entry : list) {
            ode::BoundaryEntry e;
            e.mu    = parse_complex(entry.at("mu"),    prec);
            e.value = parse_complex(entry.at("value"), prec);
            out[i].push_back(std::move(e));
        }
    }
    return out;
}

void apply_options(const json& obj) {
    if (obj.contains("working_pre"))    { auto g = numeric::global_options(); g.working_pre = obj.at("working_pre").get<int>(); numeric::set_global_options(g); }
    if (obj.contains("chop_pre"))       { auto g = numeric::global_options(); g.chop_pre    = obj.at("chop_pre").get<int>();    numeric::set_global_options(g); }
    if (obj.contains("rationalize_pre")){ auto g = numeric::global_options(); g.rationalize_pre = obj.at("rationalize_pre").get<int>(); numeric::set_global_options(g); }
    if (obj.contains("silent_mode"))    { auto g = numeric::global_options(); g.silent_mode = obj.at("silent_mode").get<bool>(); numeric::set_global_options(g); }
    if (obj.contains("d0")) {
        auto g = numeric::global_options();
        const auto& v = obj.at("d0");
        if (v.is_string())             g.d0 = v.get<std::string>();
        else if (v.is_number_integer()) g.d0 = std::to_string(v.get<long>());
        else throw std::runtime_error("d0 must be a rational string ('p/q') or integer");
        fmpq_t q; fmpq_init(q);
        const bool ok = parse_fmpq(q, g.d0);
        fmpq_clear(q);
        if (!ok) throw std::runtime_error("d0: failed to parse rational '" + g.d0 + "'");
        numeric::set_global_options(g);
    }

    if (obj.contains("x_order"))        { auto e = numeric::expansion_options(); e.x_order       = obj.at("x_order").get<int>();       numeric::set_expansion_options(e); }
    if (obj.contains("extra_x_order"))  { auto e = numeric::expansion_options(); e.extra_x_order = obj.at("extra_x_order").get<int>(); numeric::set_expansion_options(e); }
    if (obj.contains("learn_x_order"))  { auto e = numeric::expansion_options(); e.learn_x_order = obj.at("learn_x_order").get<int>(); numeric::set_expansion_options(e); }
    if (obj.contains("test_x_order"))   { auto e = numeric::expansion_options(); e.test_x_order  = obj.at("test_x_order").get<int>();  numeric::set_expansion_options(e); }

    if (obj.contains("run_radius"))     { auto r = numeric::running_options(); r.run_radius    = obj.at("run_radius").get<int>();    numeric::set_running_options(r); }
    if (obj.contains("run_length"))     { auto r = numeric::running_options(); r.run_length    = obj.at("run_length").get<int>();    numeric::set_running_options(r); }
    if (obj.contains("run_candidate"))  { auto r = numeric::running_options(); r.run_candidate = obj.at("run_candidate").get<int>(); numeric::set_running_options(r); }
    if (obj.contains("run_direction"))  {
        auto r = numeric::running_options();
        std::string s = obj.at("run_direction").get<std::string>();
        if      (s == "Re")    r.run_direction = numeric::RunningOptions::Direction::Re;
        else if (s == "Im")    r.run_direction = numeric::RunningOptions::Direction::Im;
        else if (s == "NegRe") r.run_direction = numeric::RunningOptions::Direction::NegRe;
        else if (s == "NegIm") r.run_direction = numeric::RunningOptions::Direction::NegIm;
        else throw std::runtime_error("unknown run_direction '" + s + "'");
        numeric::set_running_options(r);
    }
}

json acb_to_json(const AcbValue& v) {
    json o;
    char* re = arb_get_str(acb_realref(v.raw()), 30, 0);
    char* im = arb_get_str(acb_imagref(v.raw()), 30, 0);
    o["re"] = re ? re : "";
    o["im"] = im ? im : "";
    if (re) flint_free(re);
    if (im) flint_free(im);
    return o;
}

json options_to_json() {
    json o;
    o["working_pre"]     = numeric::global_options().working_pre;
    o["chop_pre"]        = numeric::global_options().chop_pre;
    o["rationalize_pre"] = numeric::global_options().rationalize_pre;
    o["silent_mode"]     = numeric::global_options().silent_mode;
    o["d0"]              = numeric::global_options().d0;
    o["x_order"]         = numeric::expansion_options().x_order;
    o["extra_x_order"]   = numeric::expansion_options().extra_x_order;
    o["learn_x_order"]   = numeric::expansion_options().learn_x_order;
    o["test_x_order"]    = numeric::expansion_options().test_x_order;
    o["run_radius"]      = numeric::running_options().run_radius;
    o["run_length"]      = numeric::running_options().run_length;
    o["run_candidate"]   = numeric::running_options().run_candidate;
    switch (numeric::running_options().run_direction) {
        case numeric::RunningOptions::Direction::Re:    o["run_direction"] = "Re"; break;
        case numeric::RunningOptions::Direction::Im:    o["run_direction"] = "Im"; break;
        case numeric::RunningOptions::Direction::NegRe: o["run_direction"] = "NegRe"; break;
        case numeric::RunningOptions::Direction::NegIm: o["run_direction"] = "NegIm"; break;
        case numeric::RunningOptions::Direction::Custom: o["run_direction"] = "Custom"; break;
    }
    return o;
}

json laurent_solutions_to_json(
        const std::vector<pipeline::LaurentIntegralSolution>& sol) {
    json out = json::array();
    for (const auto& row : sol) {
        json item;
        item["integral"] = jintegral_to_json(row.integral);
        item["leading_order"] = row.leading_order;
        item["coefficients"] = json::array();
        for (std::size_t i = 0; i < row.coefficients.size(); ++i) {
            item["coefficients"].push_back(json{
                {"order", row.leading_order + static_cast<long>(i)},
                {"value", acb_to_json(row.coefficients[i])},
            });
        }
        out.push_back(std::move(item));
    }
    return out;
}

json sampled_solutions_to_json(
        const std::vector<pipeline::SampledIntegralSolution>& sol,
        const std::vector<AcbValue>& eps_samples) {
    json out = json::array();
    for (const auto& row : sol) {
        if (row.values.size() != eps_samples.size()) {
            throw std::runtime_error(
                "sampled_solutions_to_json: epsilon/value grid size mismatch");
        }

        json item;
        item["integral"] = jintegral_to_json(row.integral);
        item["samples"] = json::array();
        for (std::size_t i = 0; i < row.values.size(); ++i) {
            item["samples"].push_back(json{
                {"eps", acb_to_json(eps_samples[i])},
                {"value", acb_to_json(row.values[i])},
            });
        }
        out.push_back(std::move(item));
    }
    return out;
}

}  // namespace

json run_json(const json& input) {
    if (input.contains("options")) apply_options(input.at("options"));
    std::string mode = input.value("mode", std::string("amflow"));

    json out;
    out["mode"] = mode;

    if (mode == "amflow") {
        long prec = numeric::working_prec_bits();
        RationalMatrix de = parse_matrix(input.at("matrix"));
        std::vector<ode::BoundarySpec> bcs =
            parse_boundaries(input.at("boundaries"), prec);
        if (bcs.size() != de.rows()) {
            throw std::runtime_error(
                "boundaries length must match matrix size; got "
                + std::to_string(bcs.size()) + " vs "
                + std::to_string(de.rows()));
        }

        auto sol = ode::amflow(de, bcs, prec);
        out["result"] = json::array();
        for (const auto& v : sol) out["result"].push_back(acb_to_json(v));
    } else if (mode == "solve_integrals") {
        auto fc = parse_family_config(input.at("family"));
        auto targets = parse_jintegral_list(input, fc);
        const long goal_digits =
            parse_required_long(input, "goal_digits", "goal");
        const long eps_order =
            parse_required_long(input, "eps_order", "order");
        pipeline::AMFSystemOptions opts = input.contains("amf_options")
            ? parse_amf_options(input.at("amf_options"))
            : default_solve_integrals_options();
        ensure_cli_blackbox_log_file(opts);
        const std::string work_dir =
            input.value("work_dir", std::string{});

        auto sol = pipeline::solve_integrals(fc, targets, goal_digits, eps_order,
                                              opts, work_dir);
        out["result"] = laurent_solutions_to_json(sol);
    } else if (mode == "black_box_amflow") {
        const long prec = numeric::working_prec_bits();
        auto fc = parse_family_config(input.at("family"));
        auto targets = parse_jintegral_list(input, fc);
        auto eps_samples = parse_eps_samples(input, prec);
        pipeline::AMFSystemOptions opts = input.contains("amf_options")
            ? parse_amf_options(input.at("amf_options"))
            : default_solve_integrals_options();
        ensure_cli_blackbox_log_file(opts);
        const std::string work_dir =
            input.value("work_dir", std::string{});

        auto sol = pipeline::black_box_amflow(fc, targets, eps_samples, opts, work_dir);
        out["result"] = sampled_solutions_to_json(sol, eps_samples);
    } else {
        throw std::runtime_error("unknown mode '" + mode + "'");
    }
    out["options"] = options_to_json();

    return out;
}

}  // namespace amflow::api
