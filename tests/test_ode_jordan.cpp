// SPDX-License-Identifier: MIT
// Tests for amflow::ode::jordan (Layer 7 helper).

#include <gtest/gtest.h>

#include "amflow/ode/jordan.hpp"
#include "amflow/numeric/rational.hpp"

#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>

#include <stdexcept>
#include <vector>

namespace ode = amflow::ode;
namespace nm  = amflow::numeric;

namespace {

void set_si(fmpq_mat_t M, slong i, slong j, long v) {
    fmpq_set_si(fmpq_mat_entry(M, i, j), v, 1);
}

void expect_reconstructs(const fmpq_mat_t A,
                         const fmpq_mat_t S,
                         const fmpq_mat_t J,
                         const fmpq_mat_t Sinv) {
    slong n = fmpq_mat_nrows(A);
    fmpq_mat_t SJ, recon, I_check, I_ref;
    fmpq_mat_init(SJ, n, n);
    fmpq_mat_init(recon, n, n);
    fmpq_mat_init(I_check, n, n);
    fmpq_mat_init(I_ref, n, n);

    fmpq_mat_mul(SJ, S, J);
    fmpq_mat_mul(recon, SJ, Sinv);
    EXPECT_TRUE(fmpq_mat_equal(recon, A));

    fmpq_mat_mul(I_check, S, Sinv);
    fmpq_mat_one(I_ref);
    EXPECT_TRUE(fmpq_mat_equal(I_check, I_ref));

    fmpq_mat_clear(I_ref);
    fmpq_mat_clear(I_check);
    fmpq_mat_clear(recon);
    fmpq_mat_clear(SJ);
}

}  // namespace

TEST(JordanTest, SimilarityTransformedJordanBlock) {
    fmpq_mat_t Jref, P, Pinv, PJ, A, S, J, Sinv;
    fmpq_mat_init(Jref, 2, 2);
    fmpq_mat_init(P, 2, 2);
    fmpq_mat_init(Pinv, 2, 2);
    fmpq_mat_init(PJ, 2, 2);
    fmpq_mat_init(A, 2, 2);
    fmpq_mat_init(S, 2, 2);
    fmpq_mat_init(J, 2, 2);
    fmpq_mat_init(Sinv, 2, 2);

    set_si(Jref, 0, 0, 3);
    set_si(Jref, 0, 1, 1);
    set_si(Jref, 1, 1, 3);
    set_si(P, 0, 0, 1);
    set_si(P, 0, 1, 2);
    set_si(P, 1, 0, 3);
    set_si(P, 1, 1, 4);
    ASSERT_NE(fmpq_mat_inv(Pinv, P), 0);
    fmpq_mat_mul(PJ, P, Jref);
    fmpq_mat_mul(A, PJ, Pinv);

    std::vector<long> blocks;
    ode::jordan_decomposition_exact(S, J, Sinv, blocks, A);

    ASSERT_EQ(blocks.size(), 1u);
    EXPECT_EQ(blocks[0], 2);
    EXPECT_TRUE(fmpq_mat_equal(J, Jref));
    expect_reconstructs(A, S, J, Sinv);

    fmpq_mat_clear(Sinv);
    fmpq_mat_clear(J);
    fmpq_mat_clear(S);
    fmpq_mat_clear(A);
    fmpq_mat_clear(PJ);
    fmpq_mat_clear(Pinv);
    fmpq_mat_clear(P);
    fmpq_mat_clear(Jref);
}

TEST(JordanTest, MixedBlocksSameEigenvalue) {
    // A = diag(2, 2, 2) with one off-diagonal making a 2-block + 1-block.
    fmpq_mat_t A, S, J, Sinv;
    fmpq_mat_init(A, 3, 3);
    fmpq_mat_init(S, 3, 3);
    fmpq_mat_init(J, 3, 3);
    fmpq_mat_init(Sinv, 3, 3);

    set_si(A, 0, 0, 2);
    set_si(A, 0, 1, 1);
    set_si(A, 1, 1, 2);
    set_si(A, 2, 2, 2);

    std::vector<long> blocks;
    ode::jordan_decomposition_exact(S, J, Sinv, blocks, A);

    ASSERT_EQ(blocks.size(), 2u);
    EXPECT_EQ(blocks[0], 2);
    EXPECT_EQ(blocks[1], 1);
    expect_reconstructs(A, S, J, Sinv);

    fmpq_mat_clear(Sinv);
    fmpq_mat_clear(J);
    fmpq_mat_clear(S);
    fmpq_mat_clear(A);
}

TEST(JordanTest, IrrationalEigenvalueThrows) {
    // A = [[0,1],[2,0]] -> eigenvalues +-sqrt(2) are not in Q.
    fmpq_mat_t A, S, J, Sinv;
    fmpq_mat_init(A, 2, 2);
    fmpq_mat_init(S, 2, 2);
    fmpq_mat_init(J, 2, 2);
    fmpq_mat_init(Sinv, 2, 2);

    set_si(A, 0, 1, 1);
    set_si(A, 1, 0, 2);

    std::vector<long> blocks;
    EXPECT_THROW(ode::jordan_decomposition_exact(S, J, Sinv, blocks, A),
                 std::runtime_error);

    fmpq_mat_clear(Sinv);
    fmpq_mat_clear(J);
    fmpq_mat_clear(S);
    fmpq_mat_clear(A);
}

TEST(JordanTest, FmpqMatFromRationalResidue_PlainConstantMatrix) {
    // mat = [[2, 0], [0, 3]] (constant entries) -> shift = 0 -> just reads
    // num/den constant terms.
    nm::RationalMatrix m(2, 2);
    m(0, 0) = nm::RationalFunction::from_si(2);
    m(0, 1) = nm::RationalFunction();
    m(1, 0) = nm::RationalFunction();
    m(1, 1) = nm::RationalFunction::from_si(3);

    fmpq_mat_t out;
    fmpq_mat_init(out, 2, 2);
    ode::fmpq_mat_from_rational_residue(out, m, 0);

    fmpq_t two, three;
    fmpq_init(two);   fmpq_set_si(two, 2, 1);
    fmpq_init(three); fmpq_set_si(three, 3, 1);
    EXPECT_TRUE(fmpq_equal(fmpq_mat_entry(out, 0, 0), two));
    EXPECT_TRUE(fmpq_is_zero(fmpq_mat_entry(out, 0, 1)));
    EXPECT_TRUE(fmpq_is_zero(fmpq_mat_entry(out, 1, 0)));
    EXPECT_TRUE(fmpq_equal(fmpq_mat_entry(out, 1, 1), three));
    fmpq_clear(two);
    fmpq_clear(three);
    fmpq_mat_clear(out);
}

TEST(JordanTest, FmpqMatFromRationalResidue_ShiftRemovesPole) {
    // mat = [[5/eta]] -> shift 1 -> [[5]].
    nm::FmpqPoly num; num.set_coeff_si(0, 5);
    nm::FmpqPoly den; den.set_coeff_si(1, 1);
    nm::RationalMatrix m(1, 1);
    m(0, 0) = nm::RationalFunction(std::move(num), std::move(den));

    fmpq_mat_t out;
    fmpq_mat_init(out, 1, 1);
    ode::fmpq_mat_from_rational_residue(out, m, 1);

    fmpq_t five;
    fmpq_init(five); fmpq_set_si(five, 5, 1);
    EXPECT_TRUE(fmpq_equal(fmpq_mat_entry(out, 0, 0), five));
    fmpq_clear(five);
    fmpq_mat_clear(out);
}

TEST(JordanTest, FmpqMatFromRationalResidue_RemainingPoleThrows) {
    // mat = [[5/eta]] -> shift 0 leaves the pole.
    nm::FmpqPoly num; num.set_coeff_si(0, 5);
    nm::FmpqPoly den; den.set_coeff_si(1, 1);
    nm::RationalMatrix m(1, 1);
    m(0, 0) = nm::RationalFunction(std::move(num), std::move(den));

    fmpq_mat_t out;
    fmpq_mat_init(out, 1, 1);
    EXPECT_THROW(ode::fmpq_mat_from_rational_residue(out, m, 0),
                 std::runtime_error);
    fmpq_mat_clear(out);
}

TEST(JordanTest, RationalMatrixFromFmpqMat_RoundTrip) {
    fmpq_mat_t A;
    fmpq_mat_init(A, 2, 2);
    set_si(A, 0, 0, 7);
    fmpq_set_si(fmpq_mat_entry(A, 1, 1), 3, 2);   // 3/2
    auto rm = ode::rational_matrix_from_fmpq_mat(A);
    ASSERT_EQ(rm.rows(), 2u);
    ASSERT_EQ(rm.cols(), 2u);
    EXPECT_TRUE(rm(0, 1).is_zero());
    EXPECT_TRUE(rm(1, 0).is_zero());

    fmpq_mat_t out;
    fmpq_mat_init(out, 2, 2);
    ode::fmpq_mat_from_rational_residue(out, rm, 0);
    EXPECT_TRUE(fmpq_mat_equal(A, out));
    fmpq_mat_clear(A);
    fmpq_mat_clear(out);
}
