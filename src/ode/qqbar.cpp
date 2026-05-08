// SPDX-License-Identifier: MIT
// ode::qqbar — implementation.
//

#include "amflow/ode/qqbar.hpp"

#include <stdexcept>
#include <string>
#include <vector>

#include <flint/acb.h>
#include <flint/flint.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz.h>
#include <flint/qqbar.h>

namespace amflow::ode {

namespace {

long fmpz_to_long_exact(const fmpz_t x, const char* what) {
    if (!fmpz_fits_si(x)) {
        throw std::overflow_error(std::string(what) + " does not fit into long");
    }
    return fmpz_get_si(x);
}

}  // namespace

long QqbarValue::floor_real_si() const {
    qqbar_t re;
    qqbar_init(re);
    qqbar_re(re, handle_);

    fmpz_t floor_z;
    fmpz_init(floor_z);
    qqbar_floor(floor_z, re);
    long out = fmpz_to_long_exact(floor_z, "qqbar floor");

    fmpz_clear(floor_z);
    qqbar_clear(re);
    return out;
}

numeric::AcbValue QqbarValue::to_acb(long prec) const {
    numeric::AcbValue out;
    qqbar_get_acb(out.raw(), handle_, prec);
    return out;
}

std::string QqbarValue::to_string(long display_digits) const {
    char* raw = qqbar_get_str_nd(handle_, display_digits);
    std::string out = (raw != nullptr) ? raw : "<null>";
    if (raw != nullptr) flint_free(raw);
    return out;
}

std::vector<QqbarValue>
algebraic_roots_from_fmpq_poly(const fmpq_poly_t poly) {
    const slong deg = fmpq_poly_degree(poly);
    std::vector<QqbarValue> out;
    if (deg <= 0) return out;

    qqbar_ptr roots = _qqbar_vec_init(deg);
    qqbar_roots_fmpq_poly(roots, poly, 0);

    out.reserve(static_cast<std::size_t>(deg));
    for (slong i = 0; i < deg; ++i) {
        QqbarValue root;
        qqbar_set(root.raw(), roots + i);
        out.push_back(std::move(root));
    }

    _qqbar_vec_clear(roots, deg);
    return out;
}

std::vector<QqbarEigenvalue>
algebraic_eigenvalues_from_fmpq_mat(const fmpq_mat_t A) {
    const slong nrows = fmpq_mat_nrows(A);
    const slong ncols = fmpq_mat_ncols(A);
    if (nrows != ncols) {
        throw std::invalid_argument(
            "algebraic_eigenvalues_from_fmpq_mat: matrix must be square");
    }

    std::vector<QqbarEigenvalue> out;
    if (nrows == 0) return out;

    qqbar_ptr eigs = _qqbar_vec_init(nrows);
    qqbar_eigenvalues_fmpq_mat(eigs, A, 0);

    std::vector<bool> used(static_cast<std::size_t>(nrows), false);
    for (slong i = 0; i < nrows; ++i) {
        if (used[static_cast<std::size_t>(i)]) continue;

        QqbarEigenvalue ev;
        qqbar_set(ev.value.raw(), eigs + i);
        ev.multiplicity = 1;
        used[static_cast<std::size_t>(i)] = true;

        for (slong j = i + 1; j < nrows; ++j) {
            if (used[static_cast<std::size_t>(j)]) continue;
            if (qqbar_equal(eigs + i, eigs + j)) {
                used[static_cast<std::size_t>(j)] = true;
                ++ev.multiplicity;
            }
        }
        out.push_back(std::move(ev));
    }

    _qqbar_vec_clear(eigs, nrows);
    return out;
}

}  // namespace amflow::ode
