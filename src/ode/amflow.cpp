// SPDX-License-Identifier: MIT
// ode::amflow — implementation.
//

#include "amflow/ode/amflow.hpp"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/arf.h>
#include <flint/fmpq.h>

#include "amflow/numeric/log.hpp"
#include "amflow/numeric/options.hpp"
#include "amflow/numeric/rational.hpp"
#include "amflow/ode/blocks.hpp"
#include "amflow/ode/inf.hpp"
#include "amflow/ode/path.hpp"
#include "amflow/ode/regular.hpp"
#include "amflow/ode/zero.hpp"

namespace amflow::ode {

using numeric::AcbValue;
using numeric::RationalComplex;
using numeric::RationalMatrix;

namespace {

bool debug_stage_enabled() {
    return numeric::log::trace_enabled("AMFLOW_DEBUG_STAGES");
}

void dump_acb_vector(const char* label,
                     const std::vector<AcbValue>& values) {
    if (!debug_stage_enabled()) return;
    std::cerr << "[system_stage] " << label
              << " size=" << values.size() << std::endl;
    for (std::size_t i = 0; i < values.size(); ++i) {
        std::cerr << "  [" << i << "] = "
                  << values[i].to_string(20) << std::endl;
    }
}

void dump_run(const std::vector<RationalComplex>& run) {
    if (!debug_stage_enabled()) return;
    std::cerr << "[system_stage] run size=" << run.size() << std::endl;
    for (std::size_t i = 0; i < run.size(); ++i) {
        std::cerr << "  run[" << i << "] = " << run[i].to_string()
                  << std::endl;
    }
}

int arb_mid_cmp(arb_srcptr a, arb_srcptr b) {
    return arf_cmp(arb_midref(a), arb_midref(b));
}

}  // namespace

System::System(RationalMatrix de,
               std::vector<BoundarySpec> bc,
               bool point_at_infinity)
    : de_(std::move(de)),
      bc_singular_(std::move(bc)),
      bc_kind_(BCKind::Singular),
      point_at_infinity_(point_at_infinity) {
    if (!numeric::silent_mode()) {
        numeric::log_line("LoadSystem: integrals=" + std::to_string(de_.rows())
                 + ", point=" + (point_at_infinity_ ? "infinity" : "regular"));
    }
}

System::System(RationalMatrix de,
               std::vector<AcbValue> bc_regular,
               RationalComplex point)
    : de_(std::move(de)),
      bc_regular_(std::move(bc_regular)),
      bc_kind_(BCKind::Regular),
      point_at_infinity_(false),
      point_(std::move(point)) {
    if (!numeric::silent_mode()) {
        numeric::log_line("LoadSystem: integrals=" + std::to_string(de_.rows())
                 + ", point=regular");
    }
}

void System::inf_to_regular(const RationalComplex& x0, long prec) {
    if (!point_at_infinity_) {
        throw std::runtime_error("inf_to_regular: current boundary is not at infinity");
    }
    if (bc_kind_ != BCKind::Singular) {
        throw std::runtime_error("inf_to_regular: bc is not in singular form");
    }

    auto asy = calc_inf(de_, bc_singular_, prec);
    asy_exp_ = std::move(asy);

    AcbValue inv_x0_acb;
    AcbValue x0_acb = x0.to_acb(prec);
    acb_inv(inv_x0_acb.raw(), x0_acb.raw(), prec);

    bc_regular_.clear();
    bc_regular_.reserve(asy_exp_.size());
    for (const auto& term : asy_exp_) {
        bc_regular_.push_back(evaluate_asy_expansion(term, inv_x0_acb.raw(), prec));
    }

    bc_kind_ = BCKind::Regular;
    point_at_infinity_ = false;
    point_ = x0;

    dump_acb_vector("after_inf_to_regular", bc_regular_);

    if (!numeric::silent_mode()) numeric::log_line("InfToRegular: switched to regular point");
}

void System::regular_run(const std::vector<RationalComplex>& run, long prec) {
    if (point_at_infinity_) {
        throw std::runtime_error("regular_run: current boundary is at infinity");
    }
    if (run.empty()) {
        throw std::invalid_argument("regular_run: empty run list");
    }
    if (!(run.front() == point_)) {
        throw std::runtime_error("regular_run: first run point does not match current boundary");
    }

    std::vector<AcbValue> run_acb;
    run_acb.reserve(run.size());
    for (const auto& p : run) run_acb.push_back(p.to_acb(prec));

    dump_run(run);
    dump_acb_vector("before_regular_run", bc_regular_);
    auto new_bc = calc_run(de_, bc_regular_, run_acb, prec);
    bc_regular_ = std::move(new_bc);
    point_ = run.back();
    dump_acb_vector("after_regular_run", bc_regular_);

    if (!numeric::silent_mode()) numeric::log_line("RegularRun: arrived at last point");
}

std::vector<std::vector<AcbValue>>
System::regular_interpolation(const std::vector<RationalComplex>& samples, long prec) {
    if (point_at_infinity_) {
        throw std::runtime_error("regular_interpolation: current boundary is at infinity");
    }

    const std::vector<AcbValue> saved_bc = [&]() {
        std::vector<AcbValue> tmp;
        tmp.reserve(bc_regular_.size());
        for (const auto& v : bc_regular_) tmp.push_back(v.clone());
        return tmp;
    }();
    const RationalComplex saved_point = point_;

    auto poles = numeric::get_poles(de_);

    auto radius = [&]() -> AcbValue {
        AcbValue minv;
        bool have = false;
        AcbValue p_acb = point_.to_acb(prec);
        for (const auto& pole : poles) {
            AcbValue diff;
            AcbValue pole_acb = pole.to_acb(prec);
            acb_sub(diff.raw(), p_acb.raw(), pole_acb.raw(), prec);
            AcbValue mag;
            acb_abs(acb_realref(mag.raw()), diff.raw(), prec);
            arb_zero(acb_imagref(mag.raw()));
            if (!have
                    || arb_mid_cmp(acb_realref(mag.raw()),
                                   acb_realref(minv.raw())) < 0) {
                minv = mag.clone();
                have = true;
            }
        }
        if (!have) {
            AcbValue out;
            out.set_zero();
            return out;
        }
        AcbValue out;
        arb_div_ui(acb_realref(out.raw()), acb_realref(minv.raw()),
                   static_cast<unsigned long>(numeric::run_radius()), prec);
        arb_zero(acb_imagref(out.raw()));
        return out;
    };

    auto findnum = [&](const std::vector<RationalComplex>& sam) -> std::size_t {
        std::size_t k = 0;
        AcbValue r = radius();
        AcbValue p_acb = point_.to_acb(prec);
        while (k < sam.size()) {
            AcbValue diff;
            AcbValue s_acb = sam[k].to_acb(prec);
            acb_sub(diff.raw(), p_acb.raw(), s_acb.raw(), prec);
            AcbValue d;
            acb_abs(acb_realref(d.raw()), diff.raw(), prec);
            arb_zero(acb_imagref(d.raw()));
            if (arb_mid_cmp(acb_realref(d.raw()),
                            acb_realref(r.raw())) > 0) break;
            ++k;
        }
        return k;
    };

    std::vector<std::vector<AcbValue>> vsets;
    std::vector<RationalComplex> residues = samples;

    if (!residues.empty() && residues.front() == point_) {
        std::vector<AcbValue> bc_clone;
        bc_clone.reserve(bc_regular_.size());
        for (const auto& v : bc_regular_) bc_clone.push_back(v.clone());
        vsets.push_back(std::move(bc_clone));
        residues.erase(residues.begin());
    }

    while (!residues.empty()) {
        std::size_t number = findnum(residues);
        if (number == 0) {
            auto seg = run_segment(poles, point_, residues.front(), prec);
            if (seg.empty()) {
                throw std::runtime_error(
                    "regular_interpolation: cannot lay a contour to next sample");
            }
            regular_run(seg, prec);
            std::vector<AcbValue> bc_clone;
            bc_clone.reserve(bc_regular_.size());
            for (const auto& v : bc_regular_) bc_clone.push_back(v.clone());
            vsets.push_back(std::move(bc_clone));
            residues.erase(residues.begin());
        } else {
            auto nheq  = nh_equations(de_, EquationMode::Regular);
            auto nheqn = nh_equations_num(nheq, prec);
            AcbValue point_acb = point_.to_acb(prec);
            auto rule = calcx1x2(nheqn, bc_regular_, point_acb.raw(), prec);
            for (std::size_t i = 0; i < number; ++i) {
                AcbValue dh;
                AcbValue s_acb = residues[i].to_acb(prec);
                acb_sub(dh.raw(), s_acb.raw(), point_acb.raw(), prec);
                vsets.push_back(evaluate_taylor(rule, dh.raw(), prec));
            }
            bc_regular_.clear();
            for (const auto& v : vsets.back()) bc_regular_.push_back(v.clone());
            point_ = residues[number - 1];
            residues.erase(residues.begin(), residues.begin() + number);
        }
    }

    bc_regular_.clear();
    for (const auto& v : saved_bc) bc_regular_.push_back(v.clone());
    point_ = saved_point;

    return vsets;
}

void System::solve_asy_exp(long prec) {
    if (point_at_infinity_) {
        throw std::runtime_error("solve_asy_exp: current boundary is at infinity");
    }
    AcbValue x0_acb = point_.to_acb(prec);
    asy_exp_ = calc_zero(de_, bc_regular_, x0_acb.raw(), prec);
    if (debug_stage_enabled()) {
        std::cerr << "[system_stage] after_solve_asy_exp size="
                  << asy_exp_.size() << std::endl;
        for (std::size_t i = 0; i < asy_exp_.size(); ++i) {
            std::cerr << "  asy[" << i << "] = "
                      << pick_zero_rule_s(asy_exp_[i], prec).to_string(20)
                      << std::endl;
        }
    }
}

std::vector<AcbValue> System::pick_zero_solution(long prec) const {
    std::vector<AcbValue> out;
    out.reserve(asy_exp_.size());
    for (const auto& a : asy_exp_) out.push_back(pick_zero_rule_s(a, prec));
    return out;
}

std::vector<AcbValue> amflow(const RationalMatrix& de,
                               const std::vector<BoundarySpec>& bcs,
                               long prec) {
    if (!numeric::silent_mode()) numeric::log_line("AMFlow: start.");

    std::vector<BoundarySpec> bcs_clone;
    bcs_clone.reserve(bcs.size());
    for (const auto& b : bcs) bcs_clone.push_back(clone_boundary(b));

    System sys(de, std::move(bcs_clone), /*point_at_infinity=*/true);

    if (!numeric::silent_mode()) numeric::log_line("AMFlow: building integration contour");
    auto poles = numeric::get_poles(de);
    auto run = run_eta(poles, prec);
    if (run.empty()) {
        throw std::runtime_error("amflow: integration contour not generated");
    }
    if (!numeric::silent_mode()) {
        numeric::log_line("AMFlow: contour has " + std::to_string(run.size()) + " regular points");
    }

    if (!numeric::silent_mode()) numeric::log_line("AMFlow: solving near eta = Infinity");
    sys.inf_to_regular(run.front(), prec);

    if (!numeric::silent_mode()) numeric::log_line("AMFlow: regular running");
    sys.regular_run(run, prec);

    if (!numeric::silent_mode()) numeric::log_line("AMFlow: solving near the last point");
    sys.solve_asy_exp(prec);

    if (!numeric::silent_mode()) numeric::log_line("AMFlow: picking zero values");
    auto sol = sys.pick_zero_solution(prec);

    if (!numeric::silent_mode()) numeric::log_line("AMFlow: finished.");
    return sol;
}

}  // namespace amflow::ode
