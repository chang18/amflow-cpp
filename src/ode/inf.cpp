// SPDX-License-Identifier: MIT
// ode::inf — implementation.
//
//
// An earlier draft carried a ~150-line JSONL trace block gated by
// AMFLOW_TRACE_LAYER6 for side-by-side comparison with Mathematica.  Those
// traces are dropped from this version; the short stderr AMFLOW_DEBUG_*
// gates are preserved and remain zero-cost when the env var is not set.

#include "amflow/ode/inf.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz.h>

#include "amflow/numeric/log.hpp"
#include "amflow/ode/asy.hpp"
#include "amflow/ode/regular.hpp"
#include "amflow/ode/sparse.hpp"

namespace amflow::ode {

using numeric::AcbValue;
using numeric::FmpqPoly;
using numeric::RationalFunction;
using numeric::RationalMatrix;
using numeric::acb_is_chop_zero;
using numeric::chop_pre;
using numeric::decimal_digits_to_bits;
using numeric::extra_x_order;
using numeric::log_line;
using numeric::rationalize_pre;
using numeric::silent_mode;
using numeric::working_prec_bits;
using numeric::x_order;

// ===========================================================================
//  BoundarySpec helpers
// ===========================================================================

BoundarySpec clone_boundary(const BoundarySpec& bc) {
    BoundarySpec out;
    out.reserve(bc.size());
    for (const auto& e : bc) {
        BoundaryEntry n;
        n.mu    = e.mu.clone();
        n.value = e.value.clone();
        out.push_back(std::move(n));
    }
    return out;
}

std::vector<BoundarySpec>
clone_boundaries(const std::vector<BoundarySpec>& bcs) {
    std::vector<BoundarySpec> out;
    out.reserve(bcs.size());
    for (const auto& b : bcs) out.push_back(clone_boundary(b));
    return out;
}

BoundarySpec reverse_bcs(const BoundarySpec& bc) {
    BoundarySpec out;
    out.reserve(bc.size());
    for (const auto& e : bc) {
        BoundaryEntry n;
        n.mu = e.mu.clone();
        acb_neg(n.mu.raw(), n.mu.raw());
        n.value = e.value.clone();
        out.push_back(std::move(n));
    }
    return out;
}

std::vector<BoundarySpec>
reverse_bcs(const std::vector<BoundarySpec>& bcs) {
    std::vector<BoundarySpec> out;
    out.reserve(bcs.size());
    for (const auto& b : bcs) out.push_back(reverse_bcs(b));
    return out;
}

namespace {

// Are two AcbValue close up to chop tolerance?
bool acb_close_enough(acb_srcptr a, acb_srcptr b, long prec) {
    AcbValue diff;
    acb_sub(diff.raw(), a, b, prec);
    return acb_is_chop_zero(diff.raw(), chop_pre());
}

}  // namespace

BoundarySpec union_bcs(const BoundarySpec& bc, long prec) {
    BoundarySpec out;
    std::vector<bool> done(bc.size(), false);
    for (std::size_t i = 0; i < bc.size(); ++i) {
        if (done[i]) continue;
        BoundaryEntry combined;
        combined.mu    = bc[i].mu.clone();
        combined.value = bc[i].value.clone();
        done[i] = true;
        for (std::size_t j = i + 1; j < bc.size(); ++j) {
            if (done[j]) continue;
            if (acb_close_enough(bc[i].mu.raw(), bc[j].mu.raw(), prec)) {
                acb_add(combined.value.raw(), combined.value.raw(),
                        bc[j].value.raw(), prec);
                done[j] = true;
            }
        }
        out.push_back(std::move(combined));
    }
    return out;
}

std::vector<BoundarySpec>
union_bcs(const std::vector<BoundarySpec>& bcs, long prec) {
    std::vector<BoundarySpec> out;
    out.reserve(bcs.size());
    for (const auto& b : bcs) out.push_back(union_bcs(b, prec));
    return out;
}

// ===========================================================================
//  ReadBCS
// ===========================================================================

namespace {

// Round Re(diff) to nearest integer; check residual within chop tolerance.
std::pair<bool, long> diff_is_integer(acb_srcptr a, acb_srcptr b, long prec) {
    AcbValue diff;
    acb_sub(diff.raw(), a, b, prec);

    int chop = chop_pre();
    if (chop <= 0) chop = 20;
    const double eps_double = std::pow(10.0, -chop);

    arf_t mid_im;
    arf_init(mid_im);
    arf_abs(mid_im, arb_midref(acb_imagref(diff.raw())));
    bool im_small = (arf_cmp_d(mid_im, eps_double) < 0);
    arf_clear(mid_im);
    if (!im_small) return {false, 0};

    fmpz_t z;
    fmpz_init(z);
    arf_get_fmpz(z, arb_midref(acb_realref(diff.raw())), ARF_RND_NEAR);

    arb_t resid;
    arb_init(resid);
    arb_set_fmpz(resid, z);
    arb_sub(resid, acb_realref(diff.raw()), resid, prec);

    arf_t mr;
    arf_init(mr);
    arf_abs(mr, arb_midref(resid));
    bool re_int = (arf_cmp_d(mr, eps_double) < 0);
    arf_clear(mr);
    arb_clear(resid);

    long val = re_int ? fmpz_get_si(z) : 0;
    fmpz_clear(z);
    return {re_int, val};
}

}  // namespace

ReadBcsResult read_bcs(const BoundarySpec& bc, acb_srcptr region, long prec) {
    ReadBcsResult out;

    std::vector<std::size_t> region_idx;
    std::vector<long>        region_off;
    region_idx.reserve(bc.size());
    region_off.reserve(bc.size());

    for (std::size_t i = 0; i < bc.size(); ++i) {
        auto [ok, off] = diff_is_integer(bc[i].mu.raw(), region, prec);
        if (ok) {
            region_idx.push_back(i);
            region_off.push_back(off);
        }
    }

    if (region_idx.empty()) {
        acb_set(out.ini.raw(), region);
        return out;
    }

    long min_off = region_off[0];
    for (long o : region_off) if (o < min_off) min_off = o;

    AcbValue add; add.set_si(min_off);
    acb_add(out.ini.raw(), region, add.raw(), prec);

    out.shifted.reserve(region_idx.size());
    for (std::size_t k = 0; k < region_idx.size(); ++k) {
        long shifted_off = region_off[k] - min_off;
        out.shifted.emplace_back(shifted_off, bc[region_idx[k]].value.clone());
    }
    return out;
}

// ===========================================================================
//  build_taylor_symbolic + permutation helpers
// ===========================================================================

namespace {

// Rationalize a *real* AcbValue to a fmpq_t.  Returns false if Im part is
// nonzero within rationalize_pre digits.
//
// PRECISION-MISMATCH FIX (pentabox-2L bug, 2026-05-13):
//   When rationalize_digits > floor(working_prec_bits * log10(2)) (i.e. the
//   user asks for more decimal digits than the working precision actually
//   carries), the rationalization captures the binary-representation noise
//   in the input AcbValue as a "real" rational like
//       18 * 10^29 + 41359
//       ────────────────── = 9/5 + 4.1e-26
//             10^30
//   for an input that is mathematically exact 9/5.  Downstream this
//   spurious 4.1e-26 residual flows into BuildTaylor's diagonal subtraction
//   `mat(i,i) - ini/eta`, leaks a 1e-25 entry into ConstructMatrix at the
//   "should be unsolved free coefficient" column, and SparseGaussian pivots
//   on it (scaling all other entries by 1/ε ≈ 1e25).  The amplified pivots
//   then propagate into back-substitution for masters with all-zero BC and
//   give the 1e7-1e21 ratios observed on pentabox 2L 5-leg.
//
//   Cap `rationalize_digits` at the largest value that working_prec can
//   resolve unambiguously.  This makes the rationalization round to the
//   nearest exact rational at the precision the user's working_prec
//   actually supports, instead of capturing the binary noise.
bool acb_real_to_fmpq_local(fmpq_t out, acb_srcptr x,
                             int rationalize_digits, long prec) {
    // Effective working precision in decimal digits, minus a small safety
    // margin so the *last* digit isn't sensitive to floor/ceiling of the
    // exact rational representation.
    long working_digits = static_cast<long>(prec * 0.301029995663981195L);
    if (working_digits >= 5) working_digits -= 5;
    if (working_digits < 1) working_digits = 1;
    if (rationalize_digits > working_digits) {
        rationalize_digits = static_cast<int>(working_digits);
    }

    arf_t mid_im;
    arf_init(mid_im);
    arf_abs(mid_im, arb_midref(acb_imagref(x)));
    bool im_zero = (arf_cmp_d(mid_im, std::pow(10.0, -rationalize_digits)) < 0);
    arf_clear(mid_im);
    if (!im_zero) return false;

    long bits = decimal_digits_to_bits(rationalize_digits) + 32;
    arb_t scale, scaled, xb;
    arb_init(scale); arb_init(scaled); arb_init(xb);
    arb_set_si(scale, 10);
    arb_pow_ui(scale, scale, static_cast<unsigned long>(rationalize_digits), bits);
    arb_set(xb, acb_realref(x));
    arb_mul(scaled, xb, scale, bits);
    fmpz_t numer, pow10;
    fmpz_init(numer); fmpz_init(pow10);
    arf_get_fmpz(numer, arb_midref(scaled), ARF_RND_NEAR);
    fmpz_set_si(pow10, 10);
    fmpz_pow_ui(pow10, pow10, static_cast<unsigned long>(rationalize_digits));
    fmpq_set_fmpz_frac(out, numer, pow10);
    fmpz_clear(numer); fmpz_clear(pow10);
    arb_clear(scale); arb_clear(scaled); arb_clear(xb);
    return true;
}

// Full BuildTaylor:
//   mat'[i, j] = eta^(-ini[i]) * mat[i, j] * eta^(ini[j])   if i != j
//              = mat[i, i] - ini[i] / eta                   if i == j
//
// In one fractional region all (ini[j] - ini[i]) are integers, so the
// off-diagonal powers are pure integer shifts taken from int_offsets
// (= ini[i] - ini[0]).  The diagonal subtraction needs the *full* rational
// value of ini[i], encoded in ini_rat as a degree-0 RationalFunction.
RationalMatrix build_taylor_symbolic(const RationalMatrix& mat,
                                     const std::vector<long>& int_offsets,
                                     const std::vector<RationalFunction>& ini_rat) {
    if (mat.rows() != int_offsets.size() || mat.rows() != ini_rat.size()) {
        throw std::invalid_argument("build_taylor_symbolic: size mismatch");
    }
    std::size_t N = mat.rows();
    RationalMatrix out(N, N);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            if (i != j) {
                long shift = int_offsets[j] - int_offsets[i];
                if (shift == 0) {
                    out(i, j) = mat(i, j);
                } else {
                    RationalFunction tmp = mat(i, j);
                    tmp.multiply_by_eta_power(shift);
                    out(i, j) = std::move(tmp);
                }
            } else {
                out(i, j) = mat(i, j);
                if (!ini_rat[i].is_zero()) {
                    RationalFunction corr = ini_rat[i];
                    corr.multiply_by_eta_power(-1);
                    out(i, j) -= corr;
                }
            }
        }
    }
    return out;
}

std::vector<std::size_t>
canonical_boundary_permutation(const std::vector<std::vector<std::size_t>>& blocks,
                               const std::vector<long>& int_offsets) {
    // MMA's DESolver `DetermineBlockBoundaryOrder` (DESolver.m:705-728) does
    // NOT permute the matrix before BuildTaylor/ConstructMatrix/SparseGaussian:
    // the column-processing order in SparseGaussian is the ORIGINAL master
    // order.  Sorting by offset (the previous behaviour) was a port-time
    // mistake: it changes which column SparseGaussian picks as the free
    // variable, so the boundary-order assignment can land on the wrong
    // master when a block contains masters with varied integer offsets.
    // Concretely, audit divergence D8 (banana_4L_mixed, region 3) attributed
    // the only order=0 in a 19-master block to master[8]=[1,1,1,1,3] under
    // the old descending sort, while MMA assigns it to master[20]=[1,1,1,1,7]
    // — propagating into a 637× error on the corner.
    std::vector<std::size_t> perm;
    perm.reserve(int_offsets.size());
    for (const auto& blk : blocks) {
        perm.insert(perm.end(), blk.begin(), blk.end());
    }
    return perm;
}

bool boundary_row_has_nonzero(const std::vector<std::pair<long, AcbValue>>& row) {
    for (const auto& [order, value] : row) {
        (void)order;
        if (!value.is_zero()) return true;
    }
    return false;
}

std::vector<std::size_t>
canonical_taylor_permutation(const std::vector<std::vector<std::size_t>>& blocks,
                             const std::vector<long>& int_offsets,
                             const TaylorRegionInput& bc) {
    // MMA's `DESolver.CalcTaylor` does NOT permute the matrix before
    // BuildTaylor/ConstructMatrix/SparseGaussian — see DESolver.m:752-790.
    // It iterates blocks straight from `AnalyzeBlock[mat]` and the `block`
    // variable retains the ORIGINAL master indices throughout SparseGaussian
    // and fid lookup.
    //
    // The previous version sorted within each block by
    // (boundary_row_has_nonzero asc, int_offsets asc, index asc) — a
    // port-time heuristic that re-routes Gauss-elimination's free column
    // to a different master.  See D8 root-cause analysis in ROADMAP.md.
    (void)int_offsets;
    (void)bc;
    std::vector<std::size_t> perm;
    perm.reserve(int_offsets.size());
    for (const auto& blk : blocks) {
        perm.insert(perm.end(), blk.begin(), blk.end());
    }
    return perm;
}

RationalMatrix permute_matrix(const RationalMatrix& mat,
                              const std::vector<std::size_t>& perm) {
    RationalMatrix out(mat.rows(), mat.cols());
    for (std::size_t i = 0; i < perm.size(); ++i) {
        for (std::size_t j = 0; j < perm.size(); ++j) {
            out(i, j) = mat(perm[i], perm[j]);
        }
    }
    return out;
}

std::vector<long> permute_long_vec(const std::vector<long>& src,
                                   const std::vector<std::size_t>& perm) {
    std::vector<long> out(perm.size(), 0);
    for (std::size_t i = 0; i < perm.size(); ++i) out[i] = src[perm[i]];
    return out;
}

std::vector<RationalFunction>
permute_rational_vec(const std::vector<RationalFunction>& src,
                     const std::vector<std::size_t>& perm) {
    std::vector<RationalFunction> out(perm.size());
    for (std::size_t i = 0; i < perm.size(); ++i) out[i] = src[perm[i]];
    return out;
}

std::vector<AcbValue> permute_acb_vec(const std::vector<AcbValue>& src,
                                      const std::vector<std::size_t>& perm) {
    std::vector<AcbValue> out(perm.size());
    for (std::size_t i = 0; i < perm.size(); ++i) out[i] = src[perm[i]].clone();
    return out;
}

TaylorRegionInput permute_region_input(const TaylorRegionInput& src,
                                       const std::vector<std::size_t>& perm) {
    TaylorRegionInput out(perm.size());
    for (std::size_t i = 0; i < perm.size(); ++i) {
        const auto& row = src[perm[i]];
        out[i].reserve(row.size());
        for (const auto& [order, value] : row) {
            out[i].emplace_back(order, value.clone());
        }
    }
    return out;
}

TaylorCoefficients unpermute_taylor_coeffs(TaylorCoefficients coeffs,
                                           const std::vector<std::size_t>& perm) {
    TaylorCoefficients out(coeffs.size());
    for (std::size_t i = 0; i < perm.size(); ++i) {
        out[perm[i]] = std::move(coeffs[i]);
    }
    return out;
}

std::vector<long> unpermute_long_vec(std::vector<long> values,
                                     const std::vector<std::size_t>& perm) {
    std::vector<long> out(values.size(), -1);
    for (std::size_t i = 0; i < perm.size(); ++i) out[perm[i]] = values[i];
    return out;
}

}  // namespace

// ===========================================================================
//  determine_block_boundary_order   /   determine_boundary_order
// ===========================================================================

namespace {

struct TaylorPreparedSystem {
    std::vector<BlockEquationNum> eqs;
    std::vector<long>             rank_per_block;
};

TaylorPreparedSystem
prepare_taylor_system(const RationalMatrix& m_pure,
                      const std::vector<AcbValue>& /*ini_per_integral*/,
                      long prec) {
    auto symbolic = nh_equations(m_pure, EquationMode::Taylor);
    auto numeric  = nh_equations_num(symbolic, prec);

    TaylorPreparedSystem out;
    out.eqs.reserve(numeric.size());
    out.rank_per_block.reserve(numeric.size());

    long global_rank = poincare_rank(m_pure);
    if (global_rank < 0) global_rank = -1;

    for (std::size_t k = 0; k < symbolic.size(); ++k) {
        out.eqs.push_back(std::move(numeric[k]));
        out.rank_per_block.push_back(global_rank);
    }
    return out;
}

std::vector<long>
determine_block_boundary_order_impl(const RationalMatrix& block_mat,
                                    const std::vector<long>& int_offsets,
                                    const std::vector<RationalFunction>& ini_rat) {
    std::size_t Nblock = block_mat.rows();
    if (Nblock == 0) return {};

    auto canonical_perm =
        canonical_boundary_permutation(analyze_block(block_mat), int_offsets);
    AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
        std::cerr << "[dbo] Nblock=" << Nblock << " perm=[";
        for (std::size_t i = 0; i < canonical_perm.size(); ++i) {
            if (i) std::cerr << ",";
            std::cerr << canonical_perm[i];
        }
        std::cerr << "]\n";
    }
    RationalMatrix block_mat_c = permute_matrix(block_mat, canonical_perm);
    auto int_offsets_c = permute_long_vec(int_offsets, canonical_perm);
    auto ini_rat_c     = permute_rational_vec(ini_rat, canonical_perm);

    auto m_pure = build_taylor_symbolic(block_mat_c, int_offsets_c, ini_rat_c);

    // Build ini_per_integral (AcbValue) from ini_rat; interface stability only.
    std::vector<AcbValue> ini_per_integral(Nblock);
    long prec_bits = working_prec_bits();
    for (std::size_t i = 0; i < Nblock; ++i) {
        const FmpqPoly& num = ini_rat_c[i].numerator();
        const FmpqPoly& den = ini_rat_c[i].denominator();
        if (num.length() <= 1 && den.length() <= 1) {
            fmpq_t n_q, d_q;
            fmpq_init(n_q); fmpq_init(d_q);
            if (num.length() == 1) num.coeff(0, n_q); else fmpq_zero(n_q);
            if (den.length() == 1) den.coeff(0, d_q); else fmpq_one(d_q);
            fmpq_t ratio; fmpq_init(ratio);
            if (!fmpq_is_zero(d_q)) fmpq_div(ratio, n_q, d_q);
            acb_set_fmpq(ini_per_integral[i].raw(), ratio, prec_bits);
            fmpq_clear(n_q); fmpq_clear(d_q); fmpq_clear(ratio);
        }
    }

    auto prepared = prepare_taylor_system(m_pure, ini_per_integral, prec_bits);
    long extra = static_cast<long>(extra_x_order());

    std::vector<long> result(Nblock, -1);

    for (std::size_t bk = 0; bk < prepared.eqs.size(); ++bk) {
        const auto& eq_num = prepared.eqs[bk];
        const auto& blk    = eq_num.block;
        const std::size_t Nblk = blk.size();

        auto cm = construct_matrix(eq_num.dxexp, eq_num.axexp, extra);
        std::vector<AcbValue> nh;
        auto sg = sparse_gaussian(cm, nh);

        long total_columns = cm.total_columns;

        std::vector<bool> is_solved(total_columns, false);
        for (const auto& r : sg.reduce) {
            if (!r.empty()) is_solved[r[0].col] = true;
        }
        AMFLOW_TRACE("AMFLOW_DEBUG_DBO") {
            std::cerr << "[DBO] block=" << bk
                      << " total_columns=" << total_columns
                      << " reduce.size=" << sg.reduce.size() << std::endl;
        }

        std::vector<std::pair<std::size_t, long>> fids;
        for (long n = 1; n <= total_columns; ++n) {
            if (is_solved[n - 1]) continue;
            std::size_t row = static_cast<std::size_t>((n - 1) % static_cast<long>(Nblk));
            long order = extra + 1 - (n - 1) / static_cast<long>(Nblk);
            fids.emplace_back(blk[row], order);
        }

        long order_cut = (extra + 1) / 2;
        bool found = false;
        for (; order_cut <= extra + 1; ++order_cut) {
            bool present = false;
            for (const auto& f : fids) {
                if (f.second == order_cut) { present = true; break; }
            }
            if (!present) { found = true; break; }
        }
        if (!found) {
            throw std::runtime_error(
                "determine_block_boundary_order: cannot find a cutoff order; "
                "increase ExtraXOrder.");
        }

        std::vector<std::pair<std::size_t, long>> kept;
        for (const auto& f : fids) if (f.second < order_cut) kept.push_back(f);

        for (std::size_t i = 0; i < Nblock; ++i) {
            long maxo = -1;
            for (const auto& f : kept) {
                if (f.first == i && f.second > maxo) maxo = f.second;
            }
            if (maxo > result[i]) result[i] = maxo;
        }
    }
    return unpermute_long_vec(std::move(result), canonical_perm);
}

std::vector<long>
int_offsets_from_rational_powers(const std::vector<RationalFunction>& power_q) {
    std::vector<long> out(power_q.size(), 0);
    if (power_q.empty()) return out;
    for (std::size_t i = 1; i < power_q.size(); ++i) {
        RationalFunction diff = power_q[i];
        diff -= power_q[0];
        const FmpqPoly& num = diff.numerator();
        const FmpqPoly& den = diff.denominator();
        if (num.length() > 1 || den.length() > 1 || !den.is_one()) {
            throw std::invalid_argument(
                "int_offsets_from_rational_powers: power[i] - power[0] is not an integer");
        }
        fmpq_t c;
        fmpq_init(c);
        if (num.length() == 1) num.coeff(0, c);
        if (!fmpz_is_one(fmpq_denref(c))) {
            fmpq_clear(c);
            throw std::invalid_argument(
                "int_offsets_from_rational_powers: power[i] - power[0] is not an integer");
        }
        out[i] = fmpz_get_si(fmpq_numref(c));
        fmpq_clear(c);
    }
    return out;
}

}  // namespace

std::vector<long>
determine_block_boundary_order(const RationalMatrix& block_mat,
                               const std::vector<long>& power) {
    std::vector<RationalFunction> ini_rat(power.size());
    for (std::size_t i = 0; i < power.size(); ++i) {
        fmpq_t q;
        fmpq_init(q);
        fmpq_set_si(q, power[i], 1);
        ini_rat[i] = RationalFunction::from_fmpq(q);
        fmpq_clear(q);
    }
    return determine_block_boundary_order_impl(block_mat, power, ini_rat);
}

std::vector<long>
determine_block_boundary_order(const RationalMatrix& block_mat,
                               const std::vector<RationalFunction>& power_q) {
    auto int_offsets = int_offsets_from_rational_powers(power_q);
    return determine_block_boundary_order_impl(block_mat, int_offsets, power_q);
}

std::vector<long>
determine_boundary_order(const RationalMatrix& mat,
                         const std::vector<long>& power_per_integral) {
    auto blocks = analyze_block(mat);
    std::vector<long> result(mat.rows(), -1);

    for (const auto& blk : blocks) {
        RationalMatrix sub = mat.submatrix(blk, blk);
        std::vector<long> sub_power;
        sub_power.reserve(blk.size());
        for (auto r : blk) sub_power.push_back(power_per_integral[r]);

        auto sub_result = determine_block_boundary_order(sub, sub_power);
        for (std::size_t i = 0; i < blk.size(); ++i) {
            result[blk[i]] = sub_result[i];
        }
    }
    return result;
}

std::vector<long>
determine_boundary_order(const RationalMatrix& mat,
                         const std::vector<RationalFunction>& power_q) {
    auto blocks = analyze_block(mat);
    std::vector<long> result(mat.rows(), -1);

    for (const auto& blk : blocks) {
        RationalMatrix sub = mat.submatrix(blk, blk);
        std::vector<RationalFunction> sub_power;
        sub_power.reserve(blk.size());
        for (auto r : blk) sub_power.push_back(power_q[r]);

        auto sub_result = determine_block_boundary_order(sub, sub_power);
        for (std::size_t i = 0; i < blk.size(); ++i) {
            result[blk[i]] = sub_result[i];
        }
    }
    return result;
}

// ===========================================================================
//  CalcTaylor
// ===========================================================================

namespace {

TaylorCoefficients
calc_taylor_with_pure_mat_and_ini(const RationalMatrix& m_pure,
                                  const std::vector<AcbValue>& ini_per_integral,
                                  const TaylorRegionInput& bc,
                                  long prec) {
    long xorder = static_cast<long>(x_order());
    long extra  = static_cast<long>(extra_x_order());
    long totalorder = xorder + extra;

    std::size_t n_int = bc.size();

    TaylorCoefficients f0(n_int);
    for (auto& row : f0) row.resize(totalorder + 1);

    auto prepared = prepare_taylor_system(m_pure, ini_per_integral, prec);

    auto get_boundary = [&](std::size_t intid, long order)
            -> std::pair<bool, AcbValue> {
        for (const auto& [o, v] : bc[intid]) {
            if (o == order) return {true, v.clone()};
        }
        return {false, AcbValue()};
    };

    AcbValue tmp;
    for (std::size_t k = 0; k < prepared.eqs.size(); ++k) {
        BlockEquationNum& eq = prepared.eqs[k];
        const auto& blk  = eq.block;
        const auto& sub  = eq.sub;
        const std::size_t Nblk = blk.size();

        // 1. nh = MapRationalExpansion(...)  per row of the block.
        std::vector<std::vector<AcbValue>> sub_expansions(sub.size());
        for (std::size_t i = 0; i < sub.size(); ++i) {
            sub_expansions[i].reserve(totalorder + 1);
            for (long n = 0; n <= totalorder; ++n) {
                sub_expansions[i].push_back(f0[sub[i]][n].clone());
            }
        }

        std::vector<std::vector<AcbValue>> nh_per_row =
            map_rational_expansion(eq.bxexpn, eq.bxexpd, sub_expansions, prec);

        for (auto& v : nh_per_row) {
            if (v.empty()) v.resize(totalorder + 1);
        }

        // Flatten: for order in totalorder..0:  for i in 0..Nblk-1: nh_flat <- nh_per_row[i][order]
        std::vector<AcbValue> nh_flat;
        nh_flat.reserve((totalorder + 1) * Nblk);
        for (long order = totalorder; order >= 0; --order) {
            for (std::size_t i = 0; i < Nblk; ++i) {
                AcbValue v;
                if (order < static_cast<long>(nh_per_row[i].size())) {
                    v = nh_per_row[i][order].clone();
                }
                nh_flat.push_back(std::move(v));
            }
        }

        // 2. Sparse linear solve.
        auto cm = construct_matrix(eq.dxexp, eq.axexp, totalorder);
        auto sg = sparse_gaussian(cm, nh_flat);

        // 3. fid mapping: 0-based, n in [0, total_cols).
        auto fid = [&](long n) -> std::pair<std::size_t, long> {
            std::size_t row = static_cast<std::size_t>(n % static_cast<long>(Nblk));
            long order = totalorder + 1 - n / static_cast<long>(Nblk);
            return {blk[row], order};
        };

        // 4. Identify unsolved variable columns and fetch boundary data.
        std::vector<bool> is_solved(cm.total_columns, false);
        for (const auto& r : sg.reduce) {
            if (!r.empty()) is_solved[r[0].col] = true;
        }
        std::vector<std::pair<std::size_t, long>> fids;
        for (long n = 0; n < cm.total_columns; ++n) {
            if (!is_solved[n]) {
                auto [intid, order] = fid(n);
                if (order <= xorder) fids.emplace_back(intid, order);
            }
        }
        AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
            std::cerr << "[calc_taylor] block " << k
                      << " Nblk=" << Nblk
                      << " total_columns=" << cm.total_columns
                      << " unsolved.fids.size=" << fids.size() << std::endl;
        }
        for (const auto& [intid, order] : fids) {
            auto [ok, val] = get_boundary(intid, order);
            if (!ok) {
                throw std::runtime_error(
                    "calc_taylor: unsolved coefficient encountered with no "
                    "boundary value supplied; supply more BCs or call "
                    "determine_boundary_order first.");
            }
            f0[intid][order] = std::move(val);
        }

        // 5. Back-substitute reduced rows in reverse to fill f0.
        for (auto it = sg.reduce.rbegin(); it != sg.reduce.rend(); ++it) {
            const SparseRow& rel = *it;
            if (rel.empty()) continue;
            auto [c_int, c_ord] = fid(rel[0].col);
            AcbValue acc;
            for (std::size_t kk = 1; kk < rel.size(); ++kk) {
                long col = rel[kk].col;
                if (col >= cm.total_columns) {
                    // Mathematica's back-sub treats the appended nh column as
                    // f0[..., -1] := -1.  We split RHS out, so subtract here.
                    acb_sub(acc.raw(), acc.raw(), rel[kk].val.raw(), prec);
                    continue;
                }
                auto [oi, oo] = fid(col);
                if (oo > totalorder || oo < 0) continue;
                acb_mul(tmp.raw(), rel[kk].val.raw(), f0[oi][oo].raw(), prec);
                acb_add(acc.raw(), acc.raw(), tmp.raw(), prec);
            }
            if (c_ord >= 0 && c_ord <= totalorder) {
                acb_neg(f0[c_int][c_ord].raw(), acc.raw());
            }
        }
    }

    // Return only the first (xorder + 1) coefficients per integral.
    TaylorCoefficients out(n_int);
    for (std::size_t i = 0; i < n_int; ++i) {
        out[i].reserve(xorder + 1);
        for (long n = 0; n <= xorder; ++n) {
            out[i].push_back(f0[i][n].clone());
        }
    }
    return out;
}

}  // namespace

TaylorCoefficients calc_taylor(const RationalMatrix& mat,
                               const TaylorRegionInput& bc,
                               long prec) {
    std::vector<AcbValue> zero_ini(mat.rows());
    return calc_taylor_with_pure_mat_and_ini(mat, zero_ini, bc, prec);
}

TaylorCoefficients
calc_taylor_with_ini(const RationalMatrix& de_inf,
                     const std::vector<AcbValue>& ini_per_integral,
                     const TaylorRegionInput& bc,
                     long prec) {
    std::vector<long> int_offsets(ini_per_integral.size(), 0);
    if (!ini_per_integral.empty()) {
        for (std::size_t i = 1; i < ini_per_integral.size(); ++i) {
            auto [ok, off] = diff_is_integer(ini_per_integral[i].raw(),
                                             ini_per_integral[0].raw(), prec);
            if (!ok) {
                throw std::invalid_argument(
                    "calc_taylor_with_ini: ini differences are not all integers; "
                    "the inputs do not belong to the same region");
            }
            int_offsets[i] = off;
        }
    }

    std::vector<RationalFunction> ini_rat(ini_per_integral.size());
    for (std::size_t i = 0; i < ini_per_integral.size(); ++i) {
        fmpq_t qi;
        fmpq_init(qi);
        bool ok = acb_real_to_fmpq_local(qi, ini_per_integral[i].raw(),
                                          rationalize_pre(), prec);
        if (!ok) {
            fmpq_clear(qi);
            throw std::runtime_error(
                "calc_taylor_with_ini: ini[" + std::to_string(i) +
                "] is not real-rational; symbolic BuildTaylor needs rational "
                "ini values.");
        }
        ini_rat[i] = RationalFunction::from_fmpq(qi);
        fmpq_clear(qi);
    }

    auto canonical_perm =
        canonical_taylor_permutation(analyze_block(de_inf), int_offsets, bc);
    RationalMatrix de_inf_c = permute_matrix(de_inf, canonical_perm);
    auto ini_per_integral_c = permute_acb_vec(ini_per_integral, canonical_perm);
    auto bc_c               = permute_region_input(bc, canonical_perm);
    auto int_offsets_c      = permute_long_vec(int_offsets, canonical_perm);
    auto ini_rat_c          = permute_rational_vec(ini_rat, canonical_perm);

    auto m_pure = build_taylor_symbolic(de_inf_c, int_offsets_c, ini_rat_c);
    auto coeffs_c = calc_taylor_with_pure_mat_and_ini(m_pure, ini_per_integral_c,
                                                       bc_c, prec);
    return unpermute_taylor_coeffs(std::move(coeffs_c), canonical_perm);
}

// ===========================================================================
//  CalcInf
// ===========================================================================

namespace {

// "deinf" = -1/eta^2 * de(1/eta).  Stays in fmpq.
RationalMatrix make_deinf(const RationalMatrix& de) {
    auto subbed = de.substituted_inverse_eta();
    fmpq_t neg_one;
    fmpq_init(neg_one);
    fmpq_set_si(neg_one, -1, 1);
    RationalFunction prefactor = RationalFunction::monomial(-2, neg_one);
    fmpq_clear(neg_one);

    RationalMatrix out(subbed.rows(), subbed.cols());
    for (std::size_t i = 0; i < subbed.rows(); ++i) {
        for (std::size_t j = 0; j < subbed.cols(); ++j) {
            out(i, j) = subbed(i, j) * prefactor;
        }
    }
    return out;
}

struct NormalEigen {
    AcbValue fractional;
    long     integer_part;
};

NormalEigen normal_eigen(acb_srcptr mu, long prec) {
    NormalEigen out;

    fmpz_t z;
    fmpz_init(z);
    arf_get_fmpz(z, arb_midref(acb_realref(mu)), ARF_RND_FLOOR);
    out.integer_part = fmpz_get_si(z);

    AcbValue floor_part;
    acb_set_fmpz(floor_part.raw(), z);
    fmpz_clear(z);

    acb_sub(out.fractional.raw(), mu, floor_part.raw(), prec);
    return out;
}

bool acb_already_in_list(acb_srcptr x, const std::vector<AcbValue>& list, long prec) {
    for (const auto& y : list) {
        AcbValue diff;
        acb_sub(diff.raw(), x, y.raw(), prec);
        if (acb_is_chop_zero(diff.raw(), chop_pre())) return true;
    }
    return false;
}

}  // namespace

std::vector<AsyExpansion>
calc_inf(const RationalMatrix& de,
         const std::vector<BoundarySpec>& bcs,
         long prec) {
    if (bcs.size() != de.rows()) {
        throw std::invalid_argument("calc_inf: bcs size does not match de.rows()");
    }
    std::size_t N = de.rows();

    // 1. Transform.
    RationalMatrix de_inf = make_deinf(de);
    auto bcinf_rev = reverse_bcs(bcs);
    auto bcinf     = union_bcs(bcinf_rev, prec);

    // 2. Distinct fractional regions.
    std::vector<AcbValue> regions;
    for (const auto& b : bcinf) {
        for (const auto& e : b) {
            auto ne = normal_eigen(e.mu.raw(), prec);
            if (!acb_already_in_list(ne.fractional.raw(), regions, prec)) {
                regions.push_back(std::move(ne.fractional));
            }
        }
    }

    // 3. Per region: assemble (ini, bc) per integral and call calc_taylor.
    std::vector<AsyExpansion> all_rule(N);

    for (auto& region : regions) {
        if (!silent_mode()) {
            log_line(std::string("CalcTaylor: current region -> ")
                     + region.to_string(20));
        }

        std::vector<AcbValue> ini_per_integral(N);
        TaylorRegionInput     bc_per_integral(N);

        for (std::size_t i = 0; i < N; ++i) {
            auto r = read_bcs(bcinf[i], region.raw(), prec);
            ini_per_integral[i] = std::move(r.ini);
            bc_per_integral[i]  = std::move(r.shifted);
        }

        auto coeffs = calc_taylor_with_ini(de_inf, ini_per_integral,
                                            bc_per_integral, prec);

        AMFLOW_TRACE("AMFLOW_DEBUG_BC") {
            for (std::size_t i = 0; i < N; ++i) {
                std::cerr << "[calc_inf] integral " << i
                          << " ini=" << ini_per_integral[i].to_string(15)
                          << " bc.size=" << bc_per_integral[i].size()
                          << " coeffs.size=" << coeffs[i].size() << std::endl;
            }
        }

        // 4. Append (ini -> coeffs) to all_rule[i] as a new AsyTerm.
        for (std::size_t i = 0; i < N; ++i) {
            AsyTerm term;
            term.mu = std::move(ini_per_integral[i]);
            term.exp.resize(1);
            term.exp[0] = std::move(coeffs[i]);
            all_rule[i].push_back(std::move(term));
        }
    }

    return all_rule;
}

}  // namespace amflow::ode
