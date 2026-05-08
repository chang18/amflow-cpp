// SPDX-License-Identifier: MIT
// qft::vacuum — implementation.
//

#include "amflow/qft/vacuum.hpp"

#include <stdexcept>
#include <string>

#include <flint/acb.h>
#include <flint/arb.h>

namespace amflow::qft {

namespace {

void gamma_lin(acb_t out, long a, long b, const acb_t eps, long prec) {
    acb_t arg;
    acb_init(arg);
    acb_set(arg, eps);
    acb_mul_si(arg, arg, a, prec);
    acb_add_si(arg, arg, b, prec);
    acb_gamma(out, arg, prec);
    acb_clear(arg);
}

void acb_mul_pow_ui(acb_t out, const acb_t base, ulong p, long prec) {
    acb_t pw;
    acb_init(pw);
    acb_pow_ui(pw, base, p, prec);
    acb_mul(out, out, pw, prec);
    acb_clear(pw);
}

}  // namespace

bool vacuum_known(long loop_count, long prop_count) noexcept {
    if (loop_count == 1 && prop_count == 1) return true;
    if (loop_count == 2 && prop_count == 3) return true;
    if (loop_count == 3 && prop_count == 4) return true;
    if (loop_count == 3 && prop_count == 5) return true;
    if (loop_count == 4 && prop_count == 5) return true;
    return false;
}

void vacuum(long loop_count, long prop_count,
            const acb_t eps, acb_t out, long prec) {
    if (!vacuum_known(loop_count, prop_count)) {
        throw std::invalid_argument(
            "vacuum: unsupported (loop, prop) = ("
            + std::to_string(loop_count) + ", "
            + std::to_string(prop_count) + ")");
    }

    acb_t g, t;
    acb_init(g);
    acb_init(t);

    if (loop_count == 1 && prop_count == 1) {
        gamma_lin(out, 1, -1, eps, prec);
        acb_neg(out, out);
    }
    else if (loop_count == 2 && prop_count == 3) {
        acb_one(out);
        gamma_lin(g, -1, 1, eps, prec);
        acb_mul_pow_ui(out, g, 2, prec);
        gamma_lin(g, 1, 0, eps, prec);
        acb_mul(out, out, g, prec);
        gamma_lin(g, 2, -1, eps, prec);
        acb_mul(out, out, g, prec);
        gamma_lin(g, -1, 2, eps, prec);
        acb_div(out, out, g, prec);
        acb_neg(out, out);
    }
    else if (loop_count == 3 && prop_count == 4) {
        acb_one(out);
        gamma_lin(g, -1, 1, eps, prec);
        acb_mul_pow_ui(out, g, 3, prec);
        gamma_lin(g, 2, -1, eps, prec);
        acb_mul(out, out, g, prec);
        gamma_lin(g, 3, -2, eps, prec);
        acb_mul(out, out, g, prec);
        gamma_lin(g, -1, 2, eps, prec);
        acb_div(out, out, g, prec);
    }
    else if (loop_count == 3 && prop_count == 5) {
        acb_one(out);
        gamma_lin(g, -3, 2, eps, prec);
        acb_mul(out, out, g, prec);
        gamma_lin(g, -1, 1, eps, prec);
        acb_mul_pow_ui(out, g, 4, prec);
        gamma_lin(g, 1, 0, eps, prec);
        acb_mul_pow_ui(out, g, 2, prec);
        gamma_lin(g, 3, -1, eps, prec);
        acb_mul(out, out, g, prec);
        gamma_lin(g, -2, 2, eps, prec);
        acb_pow_ui(t, g, 2, prec);
        acb_div(out, out, t, prec);
        gamma_lin(g, -1, 2, eps, prec);
        acb_div(out, out, g, prec);
        acb_neg(out, out);
    }
    else if (loop_count == 4 && prop_count == 5) {
        acb_one(out);
        gamma_lin(g, -1, 1, eps, prec);
        acb_mul_pow_ui(out, g, 4, prec);
        gamma_lin(g, 3, -2, eps, prec);
        acb_mul(out, out, g, prec);
        gamma_lin(g, 4, -3, eps, prec);
        acb_mul(out, out, g, prec);
        gamma_lin(g, -1, 2, eps, prec);
        acb_div(out, out, g, prec);
        acb_neg(out, out);
    }

    acb_clear(g);
    acb_clear(t);
}

}  // namespace amflow::qft
