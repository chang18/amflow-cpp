// SPDX-License-Identifier: MIT
// Tests for amflow::qft::vacuum.

#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>

#include <flint/acb.h>

#include "amflow/numeric/acb_wrap.hpp"
#include "amflow/numeric/options.hpp"
#include "amflow/qft/vacuum.hpp"

namespace nm  = amflow::numeric;
namespace qft = amflow::qft;

namespace {

void set_eps_rational(acb_t out, long num, long den) {
    fmpq_t q;
    fmpq_init(q);
    fmpq_set_si(q, num, den);
    acb_set_fmpq(out, q, 200);
    fmpq_clear(q);
}

double acb_real_mid_double(const acb_t z) {
    arb_t r;
    arb_init(r);
    acb_get_real(r, z);
    double v = arf_get_d(arb_midref(r), ARF_RND_NEAR);
    arb_clear(r);
    return v;
}

::testing::AssertionResult close_to(double got, double expected, double tol) {
    double err = std::abs(got - expected);
    double scale = std::max(std::abs(expected), 1.0);
    if (err / scale <= tol) {
        return ::testing::AssertionSuccess();
    }
    return ::testing::AssertionFailure()
        << "got " << got << " vs expected " << expected
        << " (rel err " << (err / scale) << ", tol " << tol << ")";
}

}  // namespace

TEST(VacuumTest, VacuumKnown_Coverage) {
    EXPECT_TRUE(qft::vacuum_known(1, 1));
    EXPECT_TRUE(qft::vacuum_known(2, 3));
    EXPECT_TRUE(qft::vacuum_known(3, 4));
    EXPECT_TRUE(qft::vacuum_known(3, 5));
    EXPECT_TRUE(qft::vacuum_known(4, 5));
    EXPECT_FALSE(qft::vacuum_known(2, 1));
    EXPECT_FALSE(qft::vacuum_known(5, 5));
}

TEST(VacuumTest, Vacuum11_AtEps0p1) {
    nm::AcbValue eps;
    set_eps_rational(eps.raw(), 1, 10);
    nm::AcbValue v;
    qft::vacuum(1, 1, eps.raw(), v.raw(), 200);
    EXPECT_TRUE(close_to(acb_real_mid_double(v.raw()),
                         10.5705641096319264, 1e-10));
}

TEST(VacuumTest, Vacuum23_AtEps0p1) {
    nm::AcbValue eps;
    set_eps_rational(eps.raw(), 1, 10);
    nm::AcbValue v;
    qft::vacuum(2, 3, eps.raw(), v.raw(), 200);
    EXPECT_TRUE(close_to(acb_real_mid_double(v.raw()),
                         64.82276029909613, 1e-10));
}

TEST(VacuumTest, Vacuum34_AtEps0p1) {
    nm::AcbValue eps;
    set_eps_rational(eps.raw(), 1, 10);
    nm::AcbValue v;
    qft::vacuum(3, 4, eps.raw(), v.raw(), 200);
    EXPECT_TRUE(close_to(acb_real_mid_double(v.raw()),
                         -18.30483177772678, 1e-10));
}

TEST(VacuumTest, Vacuum35_AtEps0p1) {
    nm::AcbValue eps;
    set_eps_rational(eps.raw(), 1, 10);
    nm::AcbValue v;
    qft::vacuum(3, 5, eps.raw(), v.raw(), 200);
    EXPECT_TRUE(close_to(acb_real_mid_double(v.raw()),
                         549.3564627903148, 1e-10));
}

TEST(VacuumTest, Vacuum45_AtEps0p1) {
    nm::AcbValue eps;
    set_eps_rational(eps.raw(), 1, 10);
    nm::AcbValue v;
    qft::vacuum(4, 5, eps.raw(), v.raw(), 200);
    EXPECT_TRUE(close_to(acb_real_mid_double(v.raw()),
                         3.02927188829392, 1e-10));
}

TEST(VacuumTest, Vacuum_UnknownThrows) {
    nm::AcbValue eps;
    set_eps_rational(eps.raw(), 1, 10);
    nm::AcbValue v;
    EXPECT_THROW(qft::vacuum(2, 1, eps.raw(), v.raw(), 200),
                 std::invalid_argument);
}
