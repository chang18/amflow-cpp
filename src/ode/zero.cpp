// SPDX-License-Identifier: MIT
// ode::zero — implementation.
//

#include "amflow/ode/zero.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/acb_poly.h>
#include <flint/arb.h>
#include <flint/ca.h>
#include <flint/ca_mat.h>
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>

#include "amflow/numeric/options.hpp"
#include "amflow/ode/acb_rational.hpp"
#include "amflow/ode/asy.hpp"
#include "amflow/ode/blocks.hpp"
#include "amflow/ode/jordan.hpp"
#include "amflow/ode/normalize.hpp"

namespace amflow::ode {

using numeric::AcbValue;
using numeric::FmpqPoly;
using numeric::RationalFunction;
using numeric::RationalMatrix;
using numeric::acb_is_chop_zero;
using numeric::chop_pre;
using numeric::decimal_digits_to_bits;
using numeric::learn_x_order;
using numeric::log_line;
using numeric::silent_mode;
using numeric::test_x_order;
using numeric::working_pre;
using numeric::working_prec_bits;
using numeric::x_order;

// Forward decl of the impl with optional rotations.
std::vector<AsyExpansion>
calcx00_impl(const std::vector<BlockEquation>&            nheq,
             const std::vector<BlockEquationNum>&         nheqn,
             const std::vector<AcbValue>&                 bc,
             acb_srcptr                                   x0,
             const AsymptoticBehaviorList&                behavior,
             const std::vector<BlockJordanRotation>*      rotations,
             long                                         prec);

// ===========================================================================
//  Local helpers
// ===========================================================================

namespace {

inline AcbValue zero_acb() { AcbValue v; return v; }

std::string truncate_for_error(std::string s) {
    constexpr std::size_t kMax = 1600;
    if (s.size() <= kMax) return s;
    s.resize(kMax);
    s += "...<truncated>";
    return s;
}

std::string acb_to_error_string(acb_srcptr x) {
    AcbValue v;
    acb_set(v.raw(), x);
    return v.to_string(40);
}

AcbValue evaluate_with_context(const RationalFunction& rf,
                               acb_srcptr x,
                               long prec,
                               const std::string& where) {
    try {
        return rf.evaluate(x, prec);
    } catch (const std::exception& e) {
        throw std::runtime_error(
            where + ": " + e.what()
            + "; x=" + acb_to_error_string(x)
            + "; rf=" + truncate_for_error(rf.to_string()));
    }
}

std::vector<AcbValue> clone_vec(const std::vector<AcbValue>& v) {
    std::vector<AcbValue> out;
    out.reserve(v.size());
    for (const auto& x : v) out.push_back(x.clone());
    return out;
}

bool acb_close(acb_srcptr a, acb_srcptr b, long prec) {
    AcbValue d;
    acb_sub(d.raw(), a, b, prec);
    return acb_is_chop_zero(d.raw(), chop_pre());
}

BlockBehavior union_behavior_local(const BlockBehavior& in, long prec);
BlockBehavior union_behavior2_local(const BlockBehavior& in, long prec);

double acb_mid_abs(acb_srcptr z) {
    arf_t mr, mi;
    arf_init(mr);
    arf_init(mi);
    arf_abs(mr, arb_midref(acb_realref(z)));
    arf_abs(mi, arb_midref(acb_imagref(z)));
    const double out = arf_get_d(mr, ARF_RND_UP)
                     + arf_get_d(mi, ARF_RND_UP);
    arf_clear(mr);
    arf_clear(mi);
    return out;
}

bool acb_close_to_zero(acb_srcptr z) {
    return acb_is_chop_zero(z, chop_pre());
}

long acb_poly_valuation_chop(const acb_poly_t poly) {
    const slong len = acb_poly_length(poly);
    for (slong i = 0; i < len; ++i) {
        AcbValue c;
        acb_poly_get_coeff_acb(c.raw(), poly, i);
        if (!acb_close_to_zero(c.raw())) return static_cast<long>(i);
    }
    return -1;
}

std::vector<AcbValue> acb_poly_coeffs(const acb_poly_t poly,
                                      long            len,
                                      long            prec) {
    std::vector<AcbValue> out;
    out.reserve(static_cast<std::size_t>(std::max<long>(0, len)));
    for (long i = 0; i < len; ++i) {
        AcbValue c;
        acb_poly_get_coeff_acb(c.raw(), poly, i);
        acb_set_round(c.raw(), c.raw(), prec);
        out.push_back(std::move(c));
    }
    return out;
}

std::vector<AcbValue> acb_rational_taylor_coeffs(
        AcbRationalFunction r,
        long                len,
        long                prec) {
    std::vector<AcbValue> out(static_cast<std::size_t>(std::max<long>(0, len)));
    if (len <= 0 || r.is_zero()) return out;
    r.strip_eta_prefix();

    const long vn = acb_poly_valuation_chop(r.num());
    const long vd = acb_poly_valuation_chop(r.den());
    if (vn < 0) return out;
    if (vd < 0) {
        throw std::runtime_error("acb_rational_taylor_coeffs: zero denominator");
    }
    if (vd != 0) {
        throw std::runtime_error(
            "acb_rational_taylor_coeffs: input has a pole at eta=0");
    }

    auto num = acb_poly_coeffs(r.num(), len, prec);
    auto den = acb_poly_coeffs(r.den(), len, prec);
    if (den.empty() || acb_close_to_zero(den[0].raw())) {
        std::ostringstream oss;
        oss << "acb_rational_taylor_coeffs: denominator still has zero "
            << "constant term after eta-prefix stripping"
            << " (vn=" << vn << ", vd=" << vd
            << ")";
        throw std::runtime_error(oss.str());
    }
    std::vector<AcbValue> unit(static_cast<std::size_t>(len));
    unit[0].set_one();
    return rational_expansion(num, den, unit, prec);
}

PowerSeries acb_rational_to_power_series(AcbRationalFunction r,
                                         long                prec) {
    PowerSeries ps;
    if (r.is_zero()) return ps;
    r.strip_eta_prefix();

    const long vn = acb_poly_valuation_chop(r.num());
    const long vd = acb_poly_valuation_chop(r.den());
    if (vn < 0) return ps;
    if (vd < 0) {
        throw std::runtime_error("acb_rational_to_power_series: zero denominator");
    }
    const long offset = vn - vd;
    r.multiply_by_eta_power(-offset);
    r.strip_eta_prefix();

    const long len = static_cast<long>(x_order()) + 1;
    auto num = acb_poly_coeffs(r.num(), len, prec);
    auto den = acb_poly_coeffs(r.den(), len, prec);
    if (den.empty() || acb_close_to_zero(den[0].raw())) {
        std::ostringstream oss;
        oss << "acb_rational_to_power_series: denominator still has zero "
            << "constant term after eta-prefix stripping"
            << " (vn=" << vn << ", vd=" << vd
            << ", offset=" << offset << ")";
        throw std::runtime_error(oss.str());
    }
    std::vector<AcbValue> unit(static_cast<std::size_t>(len));
    unit[0].set_one();
    ps.leading_offset = offset;
    ps.coeffs = rational_expansion(num, den, unit, prec);
    return ps;
}

std::vector<std::vector<PowerSeries>>
to_power_series_matrix_acb(const AcbRationalMatrix& m, long prec) {
    std::vector<std::vector<PowerSeries>> out(m.rows());
    for (std::size_t i = 0; i < m.rows(); ++i) {
        out[i].reserve(m.cols());
        for (std::size_t j = 0; j < m.cols(); ++j) {
            out[i].push_back(acb_rational_to_power_series(m(i, j), prec));
        }
    }
    return out;
}

AcbValue evaluate_acb_rational(const AcbRationalFunction& f,
                               acb_srcptr                 x,
                               long                       prec) {
    AcbValue num;
    AcbValue den;
    acb_poly_evaluate(num.raw(), f.num(), x, prec);
    acb_poly_evaluate(den.raw(), f.den(), x, prec);
    if (acb_close_to_zero(den.raw())) {
        throw std::runtime_error("evaluate_acb_rational: denominator vanished");
    }
    AcbValue out;
    acb_div(out.raw(), num.raw(), den.raw(), prec);
    return out;
}

std::vector<AcbValue> acb_flat_from_ca_mat_local(ca_mat_t mat,
                                                 ca_ctx_t ctx,
                                                 long     prec) {
    const slong rows = ca_mat_nrows(mat);
    const slong cols = ca_mat_ncols(mat);
    std::vector<AcbValue> out(static_cast<std::size_t>(rows * cols));
    for (slong i = 0; i < rows; ++i) {
        for (slong j = 0; j < cols; ++j) {
            ca_get_acb(out[static_cast<std::size_t>(i * cols + j)].raw(),
                       ca_mat_entry(mat, i, j), prec, ctx);
        }
    }
    return out;
}

std::vector<AcbValue> acb_flat_from_fmpq_mat_local(const fmpq_mat_t mat,
                                                   long             prec) {
    const slong rows = fmpq_mat_nrows(mat);
    const slong cols = fmpq_mat_ncols(mat);
    std::vector<AcbValue> out(static_cast<std::size_t>(rows * cols));
    for (slong i = 0; i < rows; ++i) {
        for (slong j = 0; j < cols; ++j) {
            out[static_cast<std::size_t>(i * cols + j)].set_fmpq(
                fmpq_mat_entry(mat, i, j), prec);
        }
    }
    return out;
}

std::vector<AcbValue> acb_flat_matmul_local(const std::vector<AcbValue>& A,
                                            const std::vector<AcbValue>& B,
                                            std::size_t                  n,
                                            long                         prec) {
    std::vector<AcbValue> out(n * n);
    AcbValue tmp;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            AcbValue acc;
            for (std::size_t k = 0; k < n; ++k) {
                acb_mul(tmp.raw(), A[i * n + k].raw(),
                        B[k * n + j].raw(), prec);
                acb_add(acc.raw(), acc.raw(), tmp.raw(), prec);
            }
            out[i * n + j] = std::move(acc);
        }
    }
    return out;
}

bool acb_flat_close_local(const std::vector<AcbValue>& A,
                          const std::vector<AcbValue>& B) {
    if (A.size() != B.size()) return false;
    for (std::size_t i = 0; i < A.size(); ++i) {
        AcbValue diff;
        acb_sub(diff.raw(), A[i].raw(), B[i].raw(), working_prec_bits());
        if (!acb_close_to_zero(diff.raw())) return false;
    }
    return true;
}

struct AcbJordanData {
    std::vector<AcbValue> u;
    std::vector<AcbValue> invu;
    std::vector<AcbValue> jor;
};

AcbJordanData acb_jordan_data_from_l0(const fmpq_mat_t L0, long prec) {
    const std::size_t n = static_cast<std::size_t>(fmpq_mat_nrows(L0));
    ca_ctx_t ctx;
    ca_ctx_init(ctx);

    ca_mat_t A, J_ca, P_ca, Pinv_ca;
    ca_mat_init(A, static_cast<slong>(n), static_cast<slong>(n), ctx);
    ca_mat_init(J_ca, static_cast<slong>(n), static_cast<slong>(n), ctx);
    ca_mat_init(P_ca, static_cast<slong>(n), static_cast<slong>(n), ctx);
    ca_mat_init(Pinv_ca, static_cast<slong>(n), static_cast<slong>(n), ctx);
    ca_mat_set_fmpq_mat(A, L0, ctx);
    const int jordan_status = ca_mat_jordan_form(J_ca, P_ca, A, ctx);
    const truth_t inv_status = ca_mat_inv(Pinv_ca, P_ca, ctx);
    if (jordan_status == 0 || inv_status != T_TRUE) {
        ca_mat_clear(Pinv_ca, ctx);
        ca_mat_clear(P_ca, ctx);
        ca_mat_clear(J_ca, ctx);
        ca_mat_clear(A, ctx);
        ca_ctx_clear(ctx);
        throw std::runtime_error(
            "acb_jordan_data_from_l0: calcium Jordan decomposition failed");
    }

    AcbJordanData out;
    out.u = acb_flat_from_ca_mat_local(P_ca, ctx, prec);
    out.invu = acb_flat_from_ca_mat_local(Pinv_ca, ctx, prec);
    out.jor = acb_flat_from_ca_mat_local(J_ca, ctx, prec);

    const auto A_flat = acb_flat_from_fmpq_mat_local(L0, prec);
    const auto lhs = acb_flat_matmul_local(
        acb_flat_matmul_local(out.invu, A_flat, n, prec), out.u, n, prec);
    if (!acb_flat_close_local(lhs, out.jor)) {
        const auto alt = acb_flat_matmul_local(
            acb_flat_matmul_local(out.u, A_flat, n, prec), out.invu, n, prec);
        if (acb_flat_close_local(alt, out.jor)) {
            std::swap(out.u, out.invu);
        } else {
            ca_mat_clear(Pinv_ca, ctx);
            ca_mat_clear(P_ca, ctx);
            ca_mat_clear(J_ca, ctx);
            ca_mat_clear(A, ctx);
            ca_ctx_clear(ctx);
            throw std::runtime_error(
                "acb_jordan_data_from_l0: calcium Jordan orientation check failed");
        }
    }

    ca_mat_clear(Pinv_ca, ctx);
    ca_mat_clear(P_ca, ctx);
    ca_mat_clear(J_ca, ctx);
    ca_mat_clear(A, ctx);
    ca_ctx_clear(ctx);
    return out;
}

long acb_floor_real_mid(acb_srcptr z) {
    return static_cast<long>(std::floor(
        arf_get_d(arb_midref(acb_realref(z)), ARF_RND_DOWN)));
}

std::vector<AsyExpansion>
calc_zero_acb_shearing_single_block(const NormalizationResult& norm,
                                    const std::vector<AcbValue>& bc,
                                    acb_srcptr x0,
                                    long prec) {
    const std::size_t n = norm.B.rows();
    if (n == 0 || norm.B.cols() != n || bc.size() != n) {
        throw std::runtime_error(
            "calc_zero_acb_shearing_single_block: inconsistent dimensions");
    }

    BlockPartition partition = build_partition(norm.B);
    AcbRationalMatrix T_alg(n, n);
    AcbRationalMatrix invT_alg(n, n);
    std::vector<long> global_floors(n, 0);
    std::vector<BlockBehavior> essential_behavior(partition.size());
    AsymptoticBehaviorList behavior(partition.size());
    bool any_shift = false;

    for (std::size_t bi = 0; bi < partition.size(); ++bi) {
        const auto& blk = partition.block(bi);
        const std::size_t nb = blk.size();
        RationalMatrix sub = norm.B.submatrix(blk, blk);

        fmpq_mat_t L0;
        fmpq_mat_init(L0, static_cast<slong>(nb), static_cast<slong>(nb));
        fmpq_mat_from_rational_residue(L0, sub, /*shift=*/1);
        AcbJordanData jd = acb_jordan_data_from_l0(L0, prec);
        fmpq_mat_clear(L0);

        std::vector<long> floors(nb, 0);
        for (std::size_t i = 0; i < nb; ++i) {
            floors[i] = acb_floor_real_mid(jd.jor[i * nb + i].raw());
            global_floors[blk[i]] = floors[i];
            if (floors[i] != 0) any_shift = true;
        }

        for (std::size_t a = 0; a < nb; ++a) {
            for (std::size_t b = 0; b < nb; ++b) {
                T_alg(blk[a], blk[b]) = AcbRationalFunction::monomial(
                    floors[b], jd.u[a * nb + b].raw());
                invT_alg(blk[a], blk[b]) = AcbRationalFunction::monomial(
                    -floors[a], jd.invu[a * nb + b].raw());
            }
        }

        BlockBehavior raw_beh;
        std::size_t start = 0;
        while (start < nb) {
            std::size_t block_size = 1;
            while (start + block_size < nb) {
                const auto& super = jd.jor[
                    (start + block_size - 1) * nb + start + block_size];
                if (acb_close_to_zero(super.raw())) break;
                ++block_size;
            }
            BehaviorEntry e;
            e.mu = jd.jor[start * nb + start].clone();
            AcbValue shift;
            shift.set_si(floors[start]);
            acb_sub(e.mu.raw(), e.mu.raw(), shift.raw(), prec);
            e.log_power = static_cast<long>(block_size) - 1;
            raw_beh.push_back(std::move(e));
            start += block_size;
        }
        essential_behavior[bi] = union_behavior_local(raw_beh, prec);
    }

    if (!any_shift) {
        throw std::runtime_error(
            "calc_zero_acb_shearing_single_block: no algebraic shearing needed");
    }

    for (std::size_t bi = 0; bi < partition.size(); ++bi) {
        BlockBehavior subbeh;
        for (auto dep : partition.lower_dependencies(bi)) {
            for (const auto& bh : behavior[dep]) {
                BehaviorEntry e;
                e.mu = bh.mu.clone();
                e.log_power = bh.log_power;
                subbeh.push_back(std::move(e));
            }
        }
        subbeh = union_behavior_local(subbeh, prec);

        BlockBehavior joined;
        for (const auto& e : essential_behavior[bi]) {
            joined.push_back(BehaviorEntry{e.mu.clone(), e.log_power});
        }
        for (const auto& e : subbeh) {
            joined.push_back(BehaviorEntry{e.mu.clone(), e.log_power});
        }
        behavior[bi] = union_behavior2_local(joined, prec);
    }

    AcbRationalMatrix B_alg = invT_alg
        .matmul(AcbRationalMatrix::from_rational(norm.B, prec), prec)
        .matmul(T_alg, prec);
    for (std::size_t i = 0; i < n; ++i) {
        if (global_floors[i] == 0) continue;
        AcbValue c;
        c.set_si(global_floors[i]);
        B_alg(i, i) -= AcbRationalFunction::monomial(-1, c.raw());
    }

    AcbRationalMatrix T_total =
        AcbRationalMatrix::from_rational(norm.T, prec).matmul(T_alg, prec);
    AcbRationalMatrix invT_total =
        invT_alg.matmul(AcbRationalMatrix::from_rational(norm.invT, prec), prec);

    std::vector<AcbValue> bcT(n);
    for (std::size_t i = 0; i < n; ++i) {
        AcbValue acc;
        for (std::size_t j = 0; j < n; ++j) {
            AcbValue tij = evaluate_acb_rational(invT_total(i, j), x0, prec);
            AcbValue tmp;
            acb_mul(tmp.raw(), tij.raw(), bc[j].raw(), prec);
            acb_add(acc.raw(), acc.raw(), tmp.raw(), prec);
        }
        bcT[i] = std::move(acc);
    }

    std::vector<BlockEquationNum> nheqn;
    nheqn.reserve(partition.size());
    const long series_len = static_cast<long>(x_order()) + 1;
    std::vector<AcbValue> den_one(1);
    den_one[0].set_one();
    for (std::size_t bi = 0; bi < partition.size(); ++bi) {
        const auto& blk = partition.block(bi);
        const auto sub = partition.sub_rows(bi);
        const std::size_t nb = blk.size();

        BlockEquationNum eq;
        eq.dxexp.resize(1);
        eq.dxexp[0].set_one();
        eq.axexp.resize(nb);
        eq.bxexpn.resize(nb);
        eq.bxexpd.resize(nb);
        for (std::size_t a = 0; a < nb; ++a) {
            eq.axexp[a].resize(nb);
            for (std::size_t b = 0; b < nb; ++b) {
                AcbRationalFunction r = B_alg(blk[a], blk[b]);
                r *= AcbRationalFunction::monomial(1);
                eq.axexp[a][b] = acb_rational_taylor_coeffs(
                    std::move(r), series_len, prec);
            }

            eq.bxexpn[a].resize(sub.size());
            eq.bxexpd[a].resize(sub.size());
            for (std::size_t s = 0; s < sub.size(); ++s) {
                AcbRationalFunction r = B_alg(blk[a], sub[s]);
                r *= AcbRationalFunction::monomial(1);
                eq.bxexpn[a][s] = acb_rational_taylor_coeffs(
                    std::move(r), series_len, prec);
                eq.bxexpd[a][s] = clone_vec(den_one);
            }
        }
        eq.block = blk;
        eq.sub = sub;
        nheqn.push_back(std::move(eq));
    }

    std::vector<BlockEquation> symbolic_unused;
    auto raw = calcx00_impl(symbolic_unused, nheqn, bcT, x0,
                            behavior, nullptr, prec);
    auto T_ps = to_power_series_matrix_acb(T_total, prec);
    return ps_map_rule_set(T_ps, raw, prec);
}

}  // namespace

// ===========================================================================
//  Calcx00 helpers
// ===========================================================================

namespace {

struct ChainStorage {
    std::vector<std::vector<std::vector<std::vector<AcbValue>>>> bb;
    std::vector<std::vector<std::vector<std::vector<AcbValue>>>> kk;

    long Nblk = 0;

    void resize(long n_behaviors, long max_p_plus_one, long n_eta) {
        bb.clear();
        kk.clear();
        bb.resize(n_behaviors);
        kk.resize(n_behaviors);
        for (long m = 0; m < n_behaviors; ++m) {
            bb[m].resize(max_p_plus_one);
            kk[m].resize(max_p_plus_one);
            for (long p = 0; p < max_p_plus_one; ++p) {
                bb[m][p].resize(n_eta);
                kk[m][p].resize(n_eta);
            }
        }
    }

    std::vector<AcbValue>& bb_at(long m, long p, long n) {
        ensure_size(m, p, n);
        auto& v = bb[m][p][n];
        if (v.empty()) v.resize(Nblk);
        return v;
    }
    std::vector<AcbValue>& kk_at(long m, long p, long n) {
        ensure_size(m, p, n);
        auto& v = kk[m][p][n];
        if (v.empty()) v.resize(Nblk * Nblk);
        return v;
    }

private:
    void ensure_size(long m, long p, long n) {
        if (m >= static_cast<long>(bb.size())) {
            bb.resize(m + 1);
            kk.resize(m + 1);
        }
        if (p >= static_cast<long>(bb[m].size())) {
            bb[m].resize(p + 1);
            kk[m].resize(p + 1);
        }
        if (n >= static_cast<long>(bb[m][p].size())) {
            bb[m][p].resize(n + 1);
            kk[m][p].resize(n + 1);
        }
    }
};

acb_srcptr mat_at(const std::vector<AcbValue>& v, long Nblk, long i, long j) {
    return v[i * Nblk + j].raw();
}

std::vector<AcbValue> matvec(const std::vector<AcbValue>& M,
                             const std::vector<AcbValue>& v,
                             long Nblk, long prec) {
    std::vector<AcbValue> out(Nblk);
    AcbValue tmp;
    for (long i = 0; i < Nblk; ++i) {
        AcbValue acc;
        for (long j = 0; j < Nblk; ++j) {
            acb_mul(tmp.raw(), mat_at(M, Nblk, i, j), v[j].raw(), prec);
            acb_add(acc.raw(), acc.raw(), tmp.raw(), prec);
        }
        out[i] = std::move(acc);
    }
    return out;
}

std::vector<AcbValue> matmul_flat(const std::vector<AcbValue>& M,
                                  const std::vector<AcbValue>& X,
                                  long Nblk, long prec) {
    std::vector<AcbValue> out(Nblk * Nblk);
    AcbValue tmp;
    for (long i = 0; i < Nblk; ++i) {
        for (long j = 0; j < Nblk; ++j) {
            AcbValue acc;
            for (long k = 0; k < Nblk; ++k) {
                acb_mul(tmp.raw(), mat_at(M, Nblk, i, k), mat_at(X, Nblk, k, j), prec);
                acb_add(acc.raw(), acc.raw(), tmp.raw(), prec);
            }
            out[i * Nblk + j] = std::move(acc);
        }
    }
    return out;
}

void vec_add_inplace(std::vector<AcbValue>& v, const std::vector<AcbValue>& w, long prec) {
    for (std::size_t i = 0; i < v.size(); ++i) {
        acb_add(v[i].raw(), v[i].raw(), w[i].raw(), prec);
    }
}

void vec_scale_inplace(std::vector<AcbValue>& v, acb_srcptr s, long prec) {
    AcbValue tmp;
    for (std::size_t i = 0; i < v.size(); ++i) {
        acb_mul(tmp.raw(), v[i].raw(), s, prec);
        acb_set(v[i].raw(), tmp.raw());
    }
}

std::vector<AcbValue> zero_vec(long n) { return std::vector<AcbValue>(n); }
std::vector<AcbValue> identity_mat_flat(long n) {
    std::vector<AcbValue> out(n * n);
    for (long i = 0; i < n; ++i) out[i * n + i].set_one();
    return out;
}

std::vector<AcbValue> compute_a00(
        const std::vector<AcbValue>& dxexp,
        const std::vector<std::vector<std::vector<AcbValue>>>& axexp,
        long Nblk, long prec) {
    if (dxexp.empty() || dxexp[0].is_zero()) {
        throw std::runtime_error(
            "calcx00: dxexp[0] is zero -- the matrix is not in the expected form.");
    }
    AcbValue inv_d0;
    acb_inv(inv_d0.raw(), dxexp[0].raw(), prec);

    std::vector<AcbValue> a00(Nblk * Nblk);
    for (long i = 0; i < Nblk; ++i) {
        for (long j = 0; j < Nblk; ++j) {
            AcbValue v;
            if (!axexp[i][j].empty()) {
                acb_set(v.raw(), axexp[i][j][0].raw());
            }
            acb_mul(a00[i * Nblk + j].raw(), v.raw(), inv_d0.raw(), prec);
        }
    }
    return a00;
}

std::vector<long> find_resonance_positions(acb_srcptr spbeh,
                                           long n,
                                           const std::vector<AcbValue>& a00,
                                           long Nblk,
                                           long prec) {
    std::vector<long> out;
    AcbValue lhs, rhs, diff;
    acb_set(lhs.raw(), spbeh);
    AcbValue n_acb; n_acb.set_si(n);
    acb_add(lhs.raw(), lhs.raw(), n_acb.raw(), prec);

    for (long p = 0; p < Nblk; ++p) {
        acb_set(rhs.raw(), mat_at(a00, Nblk, p, p));
        acb_sub(diff.raw(), lhs.raw(), rhs.raw(), prec);
        if (acb_is_chop_zero(diff.raw(), chop_pre())) out.push_back(p);
    }
    return out;
}

std::vector<AcbValue> compute_a00_from_mata(
        const std::vector<AcbValue>& dxexp,
        const std::vector<AcbValue>& mata0,
        long Nblk,
        long prec) {
    if (dxexp.empty() || dxexp[0].is_zero()) {
        throw std::runtime_error(
            "calcx00: dxexp[0] is zero -- the matrix is not in the expected form.");
    }
    AcbValue inv_d0;
    acb_inv(inv_d0.raw(), dxexp[0].raw(), prec);

    std::vector<AcbValue> a00(Nblk * Nblk);
    for (long i = 0; i < Nblk; ++i) {
        for (long j = 0; j < Nblk; ++j) {
            acb_mul(a00[i * Nblk + j].raw(),
                    mata0[i * Nblk + j].raw(), inv_d0.raw(), prec);
        }
    }
    return a00;
}

void rotate_series_rows_left(
        std::vector<std::vector<AcbValue>>& rows_by_n,
        const std::vector<AcbValue>&        left,
        long                                Nblk,
        long                                prec) {
    if (rows_by_n.empty()) return;

    const long n_eta = static_cast<long>(rows_by_n[0].size());
    std::vector<std::vector<AcbValue>> rotated(static_cast<std::size_t>(Nblk));
    for (long i = 0; i < Nblk; ++i) {
        rotated[static_cast<std::size_t>(i)].resize(static_cast<std::size_t>(n_eta));
    }

    for (long n = 0; n < n_eta; ++n) {
        std::vector<AcbValue> column(Nblk);
        for (long i = 0; i < Nblk; ++i) {
            if (n < static_cast<long>(rows_by_n[i].size())) {
                acb_set(column[i].raw(), rows_by_n[i][n].raw());
            }
        }
        std::vector<AcbValue> transformed = matvec(left, column, Nblk, prec);
        for (long i = 0; i < Nblk; ++i) {
            acb_set(rotated[i][n].raw(), transformed[i].raw());
        }
    }

    rows_by_n = std::move(rotated);
}

BlockBehavior union_behavior_local(const BlockBehavior& in, long prec) {
    BlockBehavior out;
    std::vector<bool> done(in.size(), false);
    for (std::size_t i = 0; i < in.size(); ++i) {
        if (done[i]) continue;
        BehaviorEntry combined;
        combined.mu = in[i].mu.clone();
        combined.log_power = in[i].log_power;
        done[i] = true;
        for (std::size_t j = i + 1; j < in.size(); ++j) {
            if (done[j]) continue;
            AcbValue diff;
            acb_sub(diff.raw(), in[i].mu.raw(), in[j].mu.raw(), prec);
            if (acb_is_chop_zero(diff.raw(), chop_pre())) {
                if (in[j].log_power > combined.log_power) {
                    combined.log_power = in[j].log_power;
                }
                done[j] = true;
            }
        }
        out.push_back(std::move(combined));
    }
    return out;
}

BlockBehavior union_behavior2_local(const BlockBehavior& in, long prec) {
    BlockBehavior out;
    std::vector<bool> done(in.size(), false);
    for (std::size_t i = 0; i < in.size(); ++i) {
        if (done[i]) continue;
        std::vector<long> log_powers;
        log_powers.push_back(in[i].log_power);
        AcbValue mu_keep = in[i].mu.clone();
        done[i] = true;
        for (std::size_t j = i + 1; j < in.size(); ++j) {
            if (done[j]) continue;
            AcbValue diff;
            acb_sub(diff.raw(), in[i].mu.raw(), in[j].mu.raw(), prec);
            if (acb_is_chop_zero(diff.raw(), chop_pre())) {
                log_powers.push_back(in[j].log_power);
                done[j] = true;
            }
        }
        BehaviorEntry e;
        e.mu = std::move(mu_keep);
        if (log_powers.size() == 1) {
            e.log_power = log_powers[0];
        } else {
            long s = 0;
            for (long lp : log_powers) s += lp;
            e.log_power = s + 1;
        }
        out.push_back(std::move(e));
    }
    return out;
}

BlockBehavior block_behavior_from_rotation(const BlockJordanRotation& rot,
                                           long                     prec) {
    const long Nblk = static_cast<long>(rot.block.size());
    BlockBehavior raw;
    for (long start = 0; start < Nblk; ) {
        long block_size = 1;
        while (start + block_size < Nblk) {
            AcbValue diff;
            acb_sub(diff.raw(),
                    rot.jor[(start + block_size) * Nblk + (start + block_size)].raw(),
                    rot.jor[start * Nblk + start].raw(), prec);
            if (!acb_is_chop_zero(diff.raw(), chop_pre())) break;

            AcbValue super = rot.jor[(start + block_size - 1) * Nblk + start + block_size].clone();
            if (acb_is_chop_zero(super.raw(), chop_pre())) break;
            ++block_size;
        }

        BehaviorEntry e;
        e.mu = rot.jor[start * Nblk + start].clone();
        e.log_power = block_size - 1;
        raw.push_back(std::move(e));
        start += block_size;
    }
    return union_behavior_local(raw, prec);
}

AsymptoticBehaviorList asymptotic_behavior_for_calc_zero(
        const RationalMatrix&                   mat,
        const std::vector<BlockJordanRotation>& rotations,
        long                                    prec) {
    auto partition = build_partition(mat);
    if (partition.size() != rotations.size()) {
        throw std::runtime_error(
            "asymptotic_behavior_for_calc_zero: block/rotation size mismatch");
    }

    AsymptoticBehaviorList behavior(partition.size());
    for (std::size_t i = 0; i < partition.size(); ++i) {
        BlockBehavior essbeh = block_behavior_from_rotation(rotations[i], prec);
        BlockBehavior subbeh;
        for (auto j : partition.lower_dependencies(i)) {
            for (const auto& bh : behavior[j]) {
                BehaviorEntry e;
                e.mu = bh.mu.clone();
                e.log_power = bh.log_power;
                subbeh.push_back(std::move(e));
            }
        }
        subbeh = union_behavior_local(subbeh, prec);

        BlockBehavior joined;
        for (auto& e : essbeh) joined.push_back(BehaviorEntry{e.mu.clone(), e.log_power});
        for (auto& e : subbeh) joined.push_back(BehaviorEntry{e.mu.clone(), e.log_power});
        behavior[i] = union_behavior2_local(joined, prec);
    }
    return behavior;
}

}  // namespace

// ===========================================================================
//  calcx00_impl
// ===========================================================================

std::vector<AsyExpansion>
calcx00_impl(const std::vector<BlockEquation>&            /*nheq*/,
             const std::vector<BlockEquationNum>&         nheqn,
             const std::vector<AcbValue>&                 bc,
             acb_srcptr                                   x0,
             const AsymptoticBehaviorList&                behavior,
             const std::vector<BlockJordanRotation>*      rotations,
             long                                         prec) {
    long xorder = static_cast<long>(x_order());

    std::vector<AsyExpansion> result(bc.size());

    for (std::size_t k = 0; k < nheqn.size(); ++k) {
        const auto& eq        = nheqn[k];
        const auto& dxexp     = eq.dxexp;
        const auto& axexp     = eq.axexp;
        const auto& bxexpn    = eq.bxexpn;
        const auto& bxexpd    = eq.bxexpd;
        const auto& block     = eq.block;
        const auto& sub       = eq.sub;
        const std::size_t Nblk = block.size();
        const auto& beh       = behavior[k];
        const BlockJordanRotation* rot = nullptr;
        if (rotations != nullptr) {
            if (k >= rotations->size()) {
                throw std::runtime_error(
                    "calcx00: rotation list shorter than block equation list");
            }
            rot = &(*rotations)[k];
            if (rot->block != block) {
                throw std::runtime_error(
                    "calcx00: block/rotation mismatch at block index "
                    + std::to_string(k));
            }
        }

        if (!silent_mode()) {
            std::string block_str = "[";
            for (auto v : block) { block_str += std::to_string(v) + " "; }
            block_str += "]";
            log_line("Calcx00: block " + block_str + " (behaviors=" +
                     std::to_string(beh.size()) + ")");
        }

        long maxd = static_cast<long>(dxexp.size()) - 1;
        long maxa = -1;
        for (const auto& row : axexp)
            for (const auto& vec : row)
                if (static_cast<long>(vec.size()) - 1 > maxa)
                    maxa = static_cast<long>(vec.size()) - 1;

        std::vector<std::vector<AcbValue>> mata(maxa < 0 ? 0 : (maxa + 1));
        for (long m = 0; m <= maxa; ++m) {
            mata[m].clear();
            mata[m].resize(Nblk * Nblk);
            for (long i = 0; i < (long)Nblk; ++i) {
                for (long j = 0; j < (long)Nblk; ++j) {
                    if (m < static_cast<long>(axexp[i][j].size())) {
                        acb_set(mata[m][i * Nblk + j].raw(), axexp[i][j][m].raw());
                    }
                }
            }
            if (rot != nullptr) {
                std::vector<AcbValue> tmp = matmul_flat(rot->invu, mata[m], Nblk, prec);
                mata[m] = matmul_flat(tmp, rot->u, Nblk, prec);
            }
        }

        std::vector<AcbValue> a00 = (rot != nullptr)
            ? clone_vec(rot->jor)
            : (maxa >= 0
                ? compute_a00_from_mata(dxexp, mata[0], Nblk, prec)
                : compute_a00(dxexp, axexp, Nblk, prec));

        auto lookup_f = [&](std::size_t intid, acb_srcptr spbeh, long p) -> std::vector<AcbValue> {
            for (const auto& term : result[intid]) {
                if (acb_close(term.mu.raw(), spbeh, prec)) {
                    if (p < static_cast<long>(term.exp.size())) {
                        return clone_vec(term.exp[p]);
                    } else {
                        return zero_vec(xorder + 1);
                    }
                }
            }
            return zero_vec(xorder + 1);
        };

        std::vector<std::vector<std::vector<std::vector<AcbValue>>>> nh_all(beh.size());
        for (std::size_t m = 0; m < beh.size(); ++m) {
            long logk = beh[m].log_power;
            nh_all[m].resize(logk + 1);
            for (long p = 0; p <= logk; ++p) {
                std::vector<std::vector<AcbValue>> sub_expansions(sub.size());
                for (std::size_t s = 0; s < sub.size(); ++s) {
                    sub_expansions[s] = lookup_f(sub[s], beh[m].mu.raw(), p);
                }
                std::vector<std::vector<AcbValue>> nh_per_row =
                    map_rational_expansion(bxexpn, bxexpd, sub_expansions, prec);
                for (auto& v : nh_per_row) {
                    if (v.empty()) v.resize(xorder + 1);
                }
                if (rot != nullptr) {
                    rotate_series_rows_left(nh_per_row, rot->invu, Nblk, prec);
                }
                nh_all[m][p] = std::move(nh_per_row);
            }
        }

        long max_logk = 0;
        for (const auto& b : beh) max_logk = std::max(max_logk, b.log_power);

        ChainStorage chain;
        chain.Nblk = Nblk;
        chain.resize(static_cast<long>(beh.size()), max_logk + 2, xorder + 1);

        std::vector<bool> is_essential(beh.size(), false);
        struct PairData {
            std::vector<AcbValue> kk_logp1;
            std::vector<AcbValue> bb_logp1;
            std::vector<long>     complement_positions;
            bool                  populated = false;
        };
        std::vector<PairData> pair_data(beh.size());
        std::vector<std::vector<std::pair<long, std::vector<long>>>> resonance_log(beh.size());

        AcbValue tmp, tmp2;

        auto rec1 = [&](bool is_bb, long m, long p, long n,
                        acb_srcptr spbeh) -> std::vector<AcbValue> {
            std::vector<AcbValue> acc(Nblk);

            long mm_max1 = std::min<long>(n, maxa);
            for (long mm = 0; mm <= mm_max1; ++mm) {
                std::vector<AcbValue>& src = is_bb
                    ? chain.bb_at(m, p, n - mm)
                    : chain.kk_at(m, p, n - mm);
                if (src.size() != static_cast<std::size_t>(is_bb ? Nblk : Nblk * Nblk)) {
                    src.resize(is_bb ? Nblk : Nblk * Nblk);
                }
                std::vector<AcbValue> contrib;
                if (is_bb) {
                    contrib = matvec(mata[mm], src, Nblk, prec);
                } else {
                    contrib = matmul_flat(mata[mm], src, Nblk, prec);
                }
                if (acc.size() != contrib.size()) { acc.clear(); acc.resize(contrib.size()); }
                vec_add_inplace(acc, contrib, prec);
            }

            if (is_bb) {
                if (p < static_cast<long>(nh_all[m].size())
                    && n < static_cast<long>(nh_all[m][p][0].size())) {
                    if (acc.size() != static_cast<std::size_t>(Nblk)) {
                        acc.clear();
                        acc.resize(Nblk);
                    }
                    for (long i = 0; i < (long)Nblk; ++i) {
                        if (n < static_cast<long>(nh_all[m][p][i].size())) {
                            acb_add(acc[i].raw(), acc[i].raw(),
                                    nh_all[m][p][i][n].raw(), prec);
                        }
                    }
                }
            }

            long mm_max3 = std::min<long>(n, maxd);
            for (long mm = 0; mm <= mm_max3; ++mm) {
                std::vector<AcbValue>& src = is_bb
                    ? chain.bb_at(m, p, n - mm)
                    : chain.kk_at(m, p, n - mm);
                AcbValue factor;
                acb_set(factor.raw(), spbeh);
                AcbValue add; add.set_si(n - mm);
                acb_add(factor.raw(), factor.raw(), add.raw(), prec);
                acb_mul(factor.raw(), factor.raw(), dxexp[mm].raw(), prec);
                acb_neg(factor.raw(), factor.raw());
                std::vector<AcbValue> scaled = clone_vec(src);
                vec_scale_inplace(scaled, factor.raw(), prec);
                if (acc.size() != scaled.size()) { acc.clear(); acc.resize(scaled.size()); }
                vec_add_inplace(acc, scaled, prec);
            }

            long mm_max4 = std::min<long>(n, maxd);
            for (long mm = 1; mm <= mm_max4; ++mm) {
                std::vector<AcbValue>& src = is_bb
                    ? chain.bb_at(m, p + 1, n - mm)
                    : chain.kk_at(m, p + 1, n - mm);
                AcbValue factor;
                AcbValue ppp; ppp.set_si(-(p + 1));
                acb_mul(factor.raw(), dxexp[mm].raw(), ppp.raw(), prec);
                std::vector<AcbValue> scaled = clone_vec(src);
                vec_scale_inplace(scaled, factor.raw(), prec);
                if (acc.size() != scaled.size()) { acc.clear(); acc.resize(scaled.size()); }
                vec_add_inplace(acc, scaled, prec);
            }

            AcbValue denom;
            acb_mul_si(denom.raw(), dxexp[0].raw(), p + 1, prec);
            AcbValue inv_denom;
            acb_inv(inv_denom.raw(), denom.raw(), prec);
            vec_scale_inplace(acc, inv_denom.raw(), prec);
            return acc;
        };

        auto rec2 = [&](bool is_bb, long m, long p, long n,
                        const std::vector<AcbValue>& inv_mat,
                        acb_srcptr spbeh) -> std::vector<AcbValue> {
            std::vector<AcbValue> acc(Nblk);

            long mm_max1 = std::min<long>(n, maxa);
            for (long mm = 1; mm <= mm_max1; ++mm) {
                std::vector<AcbValue>& src = is_bb
                    ? chain.bb_at(m, p, n - mm)
                    : chain.kk_at(m, p, n - mm);
                std::vector<AcbValue> contrib;
                if (is_bb) contrib = matvec(mata[mm], src, Nblk, prec);
                else       contrib = matmul_flat(mata[mm], src, Nblk, prec);
                if (acc.size() != contrib.size()) { acc.clear(); acc.resize(contrib.size()); }
                vec_add_inplace(acc, contrib, prec);
            }

            if (is_bb) {
                if (p < static_cast<long>(nh_all[m].size())
                    && n < static_cast<long>(nh_all[m][p][0].size())) {
                    if (acc.size() != static_cast<std::size_t>(Nblk)) {
                        acc.clear();
                        acc.resize(Nblk);
                    }
                    for (long i = 0; i < (long)Nblk; ++i) {
                        if (n < static_cast<long>(nh_all[m][p][i].size())) {
                            acb_add(acc[i].raw(), acc[i].raw(),
                                    nh_all[m][p][i][n].raw(), prec);
                        }
                    }
                }
            }

            long mm_max3 = std::min<long>(n, maxd);
            for (long mm = 1; mm <= mm_max3; ++mm) {
                std::vector<AcbValue>& src = is_bb
                    ? chain.bb_at(m, p, n - mm)
                    : chain.kk_at(m, p, n - mm);
                AcbValue factor;
                acb_set(factor.raw(), spbeh);
                AcbValue add; add.set_si(n - mm);
                acb_add(factor.raw(), factor.raw(), add.raw(), prec);
                acb_mul(factor.raw(), factor.raw(), dxexp[mm].raw(), prec);
                acb_neg(factor.raw(), factor.raw());
                std::vector<AcbValue> scaled = clone_vec(src);
                vec_scale_inplace(scaled, factor.raw(), prec);
                if (acc.size() != scaled.size()) { acc.clear(); acc.resize(scaled.size()); }
                vec_add_inplace(acc, scaled, prec);
            }

            long mm_max4 = std::min<long>(n, maxd);
            for (long mm = 0; mm <= mm_max4; ++mm) {
                std::vector<AcbValue>& src = is_bb
                    ? chain.bb_at(m, p + 1, n - mm)
                    : chain.kk_at(m, p + 1, n - mm);
                AcbValue factor;
                AcbValue ppp; ppp.set_si(-(p + 1));
                acb_mul(factor.raw(), dxexp[mm].raw(), ppp.raw(), prec);
                std::vector<AcbValue> scaled = clone_vec(src);
                vec_scale_inplace(scaled, factor.raw(), prec);
                if (acc.size() != scaled.size()) { acc.clear(); acc.resize(scaled.size()); }
                vec_add_inplace(acc, scaled, prec);
            }

            std::vector<AcbValue> result_vec;
            if (is_bb) {
                result_vec = matvec(inv_mat, acc, Nblk, prec);
            } else {
                result_vec = matmul_flat(inv_mat, acc, Nblk, prec);
            }
            return result_vec;
        };

        auto invert_shifted_a00 = [&](long n, acb_srcptr spbeh,
                                      std::vector<AcbValue>& inv_mat) -> bool {
            std::vector<AcbValue> M(static_cast<std::size_t>(Nblk * Nblk));
            AcbValue n_plus_s;
            acb_set(n_plus_s.raw(), spbeh);
            AcbValue add; add.set_si(n);
            acb_add(n_plus_s.raw(), n_plus_s.raw(), add.raw(), prec);
            for (long i = 0; i < (long)Nblk; ++i) {
                for (long j = 0; j < (long)Nblk; ++j) {
                    AcbValue v;
                    acb_neg(v.raw(), mat_at(a00, Nblk, i, j));
                    if (i == j) acb_add(v.raw(), v.raw(), n_plus_s.raw(), prec);
                    acb_set(M[i * Nblk + j].raw(), v.raw());
                }
            }
            acb_mat_t M_acb, inv_acb;
            acb_mat_init(M_acb, Nblk, Nblk);
            acb_mat_init(inv_acb, Nblk, Nblk);
            for (long i = 0; i < (long)Nblk; ++i)
                for (long j = 0; j < (long)Nblk; ++j)
                    acb_set(acb_mat_entry(M_acb, i, j), M[i * Nblk + j].raw());
            int ok = acb_mat_inv(inv_acb, M_acb, prec);
            if (ok) {
                AcbValue inv_d0;
                acb_inv(inv_d0.raw(), dxexp[0].raw(), prec);
                inv_mat.clear();
                inv_mat.resize(Nblk * Nblk);
                for (long i = 0; i < (long)Nblk; ++i)
                    for (long j = 0; j < (long)Nblk; ++j) {
                        acb_mul(inv_mat[i * Nblk + j].raw(),
                                acb_mat_entry(inv_acb, i, j),
                                inv_d0.raw(), prec);
                    }
            }
            acb_mat_clear(M_acb);
            acb_mat_clear(inv_acb);
            return ok != 0;
        };

        for (long m = 0; m < (long)beh.size(); ++m) {
            const long logk = beh[m].log_power;

            for (long n = 0; n <= xorder; ++n) {
                std::vector<long> posi = find_resonance_positions(
                    beh[m].mu.raw(), n, a00, Nblk, prec);
                std::vector<AcbValue> inv_mat;
                bool inv_ok = false;
                bool full_matrix_resonance = false;
                if (posi.empty()) {
                    inv_ok = invert_shifted_a00(n, beh[m].mu.raw(), inv_mat);
                    if (!inv_ok) {
                        full_matrix_resonance = true;
                        for (long p = 0; p < (long)Nblk; ++p) posi.push_back(p);
                    }
                }

                if (!posi.empty()) {
                    resonance_log[m].push_back({n, posi});
                    is_essential[m] = true;

                    chain.kk_at(m, 0, n) = identity_mat_flat(Nblk);
                    chain.bb_at(m, 0, n) = zero_vec(Nblk);

                    for (long p = 0; p <= logk; ++p) {
                        chain.bb_at(m, p + 1, n) = rec1(true,  m, p, n, beh[m].mu.raw());
                        chain.kk_at(m, p + 1, n) = rec1(false, m, p, n, beh[m].mu.raw());
                    }
                    pair_data[m].kk_logp1 = clone_vec(chain.kk_at(m, logk + 1, n));
                    pair_data[m].bb_logp1 = clone_vec(chain.bb_at(m, logk + 1, n));
                    std::vector<long> comp;
                    for (long p = 0; p < (long)Nblk; ++p) {
                        if (full_matrix_resonance
                            || std::find(posi.begin(), posi.end(), p) == posi.end()) {
                            comp.push_back(p);
                        }
                    }
                    pair_data[m].complement_positions = std::move(comp);
                    pair_data[m].populated = true;
                } else {
                    for (long p = logk; p >= 0; --p) {
                        chain.bb_at(m, p, n) = rec2(true,  m, p, n, inv_mat, beh[m].mu.raw());
                        chain.kk_at(m, p, n) = rec2(false, m, p, n, inv_mat, beh[m].mu.raw());
                    }
                }
            }
        }

        std::vector<long> essentialset;
        for (long m = 0; m < (long)beh.size(); ++m) {
            if (is_essential[m]) essentialset.push_back(m);
        }

        AcbValue x0_log;
        acb_log(x0_log.raw(), x0, prec);

        auto evaluate_chain_to_matrix = [&](long m, bool use_kk) -> std::vector<AcbValue> {
            const long logk = beh[m].log_power;
            std::vector<std::vector<AcbValue>> per_p_value;
            per_p_value.resize(logk + 1);

            for (long p = 0; p <= logk; ++p) {
                long item_size = use_kk ? Nblk * Nblk : Nblk;
                std::vector<AcbValue> val(item_size);
                for (long n = xorder; n >= 0; --n) {
                    for (auto& a : val) {
                        AcbValue tmp_h;
                        acb_mul(tmp_h.raw(), a.raw(), x0, prec);
                        acb_set(a.raw(), tmp_h.raw());
                    }
                    std::vector<AcbValue>& src = use_kk
                        ? chain.kk_at(m, p, n)
                        : chain.bb_at(m, p, n);
                    if (!src.empty()) {
                        if (val.size() != src.size()) { val.clear(); val.resize(src.size()); }
                        vec_add_inplace(val, src, prec);
                    }
                }
                per_p_value[p] = std::move(val);
            }

            std::vector<AcbValue> total(use_kk ? Nblk * Nblk : Nblk);
            AcbValue log_pow;
            log_pow.set_one();
            for (long p = 0; p <= logk; ++p) {
                if (p > 0) acb_mul(log_pow.raw(), log_pow.raw(), x0_log.raw(), prec);
                for (std::size_t i = 0; i < total.size(); ++i) {
                    AcbValue tmp_h;
                    acb_mul(tmp_h.raw(), per_p_value[p][i].raw(), log_pow.raw(), prec);
                    acb_add(total[i].raw(), total[i].raw(), tmp_h.raw(), prec);
                }
            }

            AcbValue x0_spbeh;
            acb_pow(x0_spbeh.raw(), x0, beh[m].mu.raw(), prec);
            for (auto& a : total) {
                AcbValue tmp_h;
                acb_mul(tmp_h.raw(), a.raw(), x0_spbeh.raw(), prec);
                acb_set(a.raw(), tmp_h.raw());
            }
            return total;
        };

        std::vector<std::vector<AcbValue>> finalmat;
        finalmat.reserve(essentialset.size());
        for (long m : essentialset) {
            finalmat.push_back(evaluate_chain_to_matrix(m, /*use_kk=*/true));
        }
        std::vector<AcbValue> bc_block(Nblk);
        for (long i = 0; i < (long)Nblk; ++i) {
            acb_set(bc_block[i].raw(), bc[block[i]].raw());
        }
        if (rot != nullptr) {
            bc_block = matvec(rot->invu, bc_block, Nblk, prec);
        }
        std::vector<AcbValue> finalvec = clone_vec(bc_block);
        for (long m = 0; m < (long)beh.size(); ++m) {
            std::vector<AcbValue> bb_eval = evaluate_chain_to_matrix(m, /*use_kk=*/false);
            for (long i = 0; i < (long)Nblk; ++i) {
                acb_sub(finalvec[i].raw(), finalvec[i].raw(), bb_eval[i].raw(), prec);
            }
        }

        long n_ess = static_cast<long>(essentialset.size());
        long n_var = Nblk * n_ess;
        long n_pair = 0;
        long compensate = 0;

        std::vector<AcbValue> sol_v(n_var);
        if (n_var > 0) {
            acb_mat_t Asys;
            acb_mat_t rhs;
            acb_mat_init(Asys, n_var, n_var);
            acb_mat_init(rhs,  n_var, 1);
            acb_mat_zero(Asys);
            acb_mat_zero(rhs);

            std::vector<std::vector<AcbValue>> match_rows(Nblk);
            std::vector<AcbValue> match_rhs(Nblk);
            for (long row = 0; row < (long)Nblk; ++row) {
                match_rows[row].resize(n_var);
                acb_set(match_rhs[row].raw(), finalvec[row].raw());
                long col_offset = 0;
                for (std::size_t mi = 0; mi < essentialset.size(); ++mi) {
                    for (long c = 0; c < (long)Nblk; ++c) {
                        acb_set(match_rows[row][col_offset + c].raw(),
                                finalmat[mi][row * Nblk + c].raw());
                    }
                    col_offset += Nblk;
                }
            }

            std::vector<std::vector<AcbValue>> pair_rows;
            std::vector<AcbValue> pair_rhs;
            pair_rows.reserve(n_var);
            pair_rhs.reserve(n_var);
            long col_offset = 0;
            for (std::size_t mi = 0; mi < essentialset.size(); ++mi) {
                long m = essentialset[mi];
                const auto& comp = pair_data[m].complement_positions;
                for (long ip : comp) {
                    pair_rows.emplace_back(n_var);
                    pair_rhs.emplace_back();
                    auto& prow = pair_rows.back();
                    for (long k_ = 0; k_ < (long)Nblk; ++k_) {
                        acb_set(prow[col_offset + k_].raw(),
                                pair_data[m].kk_logp1[ip * Nblk + k_].raw());
                    }
                    acb_neg(pair_rhs.back().raw(), pair_data[m].bb_logp1[ip].raw());
                }
                col_offset += Nblk;
            }

            auto rank_of_rows = [&](const std::vector<std::vector<AcbValue>>& rows,
                                    const std::vector<AcbValue>* extra) -> long {
                const long nr = static_cast<long>(rows.size() + (extra ? 1 : 0));
                if (nr == 0 || n_var == 0) return 0;
                acb_mat_t M;
                acb_mat_init(M, nr, n_var);
                for (long r = 0; r < static_cast<long>(rows.size()); ++r) {
                    for (long c = 0; c < n_var; ++c) {
                        acb_set(acb_mat_entry(M, r, c), rows[r][c].raw());
                    }
                }
                if (extra) {
                    const long r = static_cast<long>(rows.size());
                    for (long c = 0; c < n_var; ++c) {
                        acb_set(acb_mat_entry(M, r, c), (*extra)[c].raw());
                    }
                }
                const double tol = std::pow(10.0, -chop_pre());
                long rank = 0;
                AcbValue tmp_h, inv, factor;
                for (long c = 0; c < n_var && rank < nr; ++c) {
                    long pivot = -1;
                    double best = 0.0;
                    for (long r = rank; r < nr; ++r) {
                        double mag = acb_mid_abs(acb_mat_entry(M, r, c));
                        if (mag > best) {
                            best = mag;
                            pivot = r;
                        }
                    }
                    if (pivot < 0 || !(best > tol)) continue;
                    if (pivot != rank) {
                        for (long j = 0; j < n_var; ++j) {
                            acb_swap(acb_mat_entry(M, rank, j),
                                     acb_mat_entry(M, pivot, j));
                        }
                    }
                    acb_inv(inv.raw(), acb_mat_entry(M, rank, c), prec);
                    for (long j = c; j < n_var; ++j) {
                        acb_mul(acb_mat_entry(M, rank, j),
                                acb_mat_entry(M, rank, j), inv.raw(), prec);
                    }
                    for (long r = 0; r < nr; ++r) {
                        if (r == rank) continue;
                        acb_set(factor.raw(), acb_mat_entry(M, r, c));
                        if (acb_is_chop_zero(factor.raw(), chop_pre())) continue;
                        for (long j = c; j < n_var; ++j) {
                            acb_mul(tmp_h.raw(), factor.raw(),
                                    acb_mat_entry(M, rank, j), prec);
                            acb_sub(acb_mat_entry(M, r, j),
                                    acb_mat_entry(M, r, j), tmp_h.raw(), prec);
                        }
                    }
                    ++rank;
                }
                acb_mat_clear(M);
                return rank;
            };

            {
                auto raw_rows = std::move(pair_rows);
                auto raw_rhs = std::move(pair_rhs);
                pair_rows.clear();
                pair_rhs.clear();
                long current_pair_rank = 0;
                for (std::size_t rr = 0; rr < raw_rows.size(); ++rr) {
                    long trial_rank = rank_of_rows(pair_rows, &raw_rows[rr]);
                    if (trial_rank > current_pair_rank) {
                        pair_rows.push_back(std::move(raw_rows[rr]));
                        pair_rhs.push_back(std::move(raw_rhs[rr]));
                        current_pair_rank = trial_rank;
                    }
                }
                n_pair = static_cast<long>(pair_rows.size());
            }

            auto match_complexity = [&](long row) -> long {
                long count = 0;
                for (long c = 0; c < n_var; ++c) {
                    if (!acb_is_chop_zero(match_rows[row][c].raw(), chop_pre())) {
                        ++count;
                    }
                }
                if (!acb_is_chop_zero(match_rhs[row].raw(), chop_pre())) ++count;
                return count;
            };

            std::vector<long> match_order(Nblk);
            for (long i = 0; i < (long)Nblk; ++i) match_order[i] = i;
            std::sort(match_order.begin(), match_order.end(),
                      [&](long a, long b) {
                          long ca = match_complexity(a);
                          long cb = match_complexity(b);
                          if (ca != cb) return ca > cb;
                          return a < b;
                      });

            auto rank_for_selected = [&](const std::vector<long>& selected) -> long {
                const long rows = static_cast<long>(pair_rows.size() + selected.size());
                if (rows == 0 || n_var == 0) return 0;
                acb_mat_t M;
                acb_mat_init(M, rows, n_var);
                long out_row = 0;
                for (const auto& prow : pair_rows) {
                    for (long c = 0; c < n_var; ++c) {
                        acb_set(acb_mat_entry(M, out_row, c), prow[c].raw());
                    }
                    ++out_row;
                }
                for (long r : selected) {
                    for (long c = 0; c < n_var; ++c) {
                        acb_set(acb_mat_entry(M, out_row, c),
                                match_rows[r][c].raw());
                    }
                    ++out_row;
                }

                const double tol = std::pow(10.0, -chop_pre());
                long rank = 0;
                AcbValue tmp_h, inv, factor;
                for (long c = 0; c < n_var && rank < rows; ++c) {
                    long pivot = -1;
                    double best = 0.0;
                    for (long r = rank; r < rows; ++r) {
                        double mag = acb_mid_abs(acb_mat_entry(M, r, c));
                        if (mag > best) {
                            best = mag;
                            pivot = r;
                        }
                    }
                    if (pivot < 0 || !(best > tol)) continue;
                    if (pivot != rank) {
                        for (long j = 0; j < n_var; ++j) {
                            acb_swap(acb_mat_entry(M, rank, j),
                                     acb_mat_entry(M, pivot, j));
                        }
                    }
                    acb_inv(inv.raw(), acb_mat_entry(M, rank, c), prec);
                    for (long j = c; j < n_var; ++j) {
                        acb_mul(acb_mat_entry(M, rank, j),
                                acb_mat_entry(M, rank, j), inv.raw(), prec);
                    }
                    for (long r = 0; r < rows; ++r) {
                        if (r == rank) continue;
                        acb_set(factor.raw(), acb_mat_entry(M, r, c));
                        if (acb_is_chop_zero(factor.raw(), chop_pre())) continue;
                        for (long j = c; j < n_var; ++j) {
                            acb_mul(tmp_h.raw(), factor.raw(),
                                    acb_mat_entry(M, rank, j), prec);
                            acb_sub(acb_mat_entry(M, r, j),
                                    acb_mat_entry(M, r, j), tmp_h.raw(), prec);
                        }
                    }
                    ++rank;
                }
                acb_mat_clear(M);
                return rank;
            };

            compensate = n_var - rank_for_selected({});
            if (compensate < 0) compensate = 0;

            if (compensate > (long)Nblk) {
                std::ostringstream msg;
                msg << "calcx00: boundary matching needs more match equations "
                    << "than block rows"
                    << "; block=[";
                for (std::size_t bi = 0; bi < block.size(); ++bi) {
                    if (bi) msg << ",";
                    msg << block[bi];
                }
                msg << "]; behaviors=" << beh.size()
                    << "; n_var=" << n_var
                    << "; n_pair=" << n_pair
                    << "; compensate=" << compensate
                    << "; Nblk=" << Nblk;
                acb_mat_clear(Asys);
                acb_mat_clear(rhs);
                throw std::runtime_error(msg.str());
            }

            std::vector<long> selected_match;
            selected_match.reserve(compensate);
            long current_rank = rank_for_selected(selected_match);
            for (long row : match_order) {
                if ((long)selected_match.size() >= compensate) break;
                std::vector<long> trial = selected_match;
                trial.push_back(row);
                long trial_rank = rank_for_selected(trial);
                if (trial_rank > current_rank) {
                    selected_match.push_back(row);
                    current_rank = trial_rank;
                }
            }
            for (long row : match_order) {
                if ((long)selected_match.size() >= compensate) break;
                if (std::find(selected_match.begin(), selected_match.end(), row)
                    == selected_match.end()) {
                    selected_match.push_back(row);
                }
            }

            long out_row = 0;
            for (long src_row : selected_match) {
                acb_set(acb_mat_entry(rhs, out_row, 0),
                        match_rhs[src_row].raw());
                for (long c = 0; c < n_var; ++c) {
                    acb_set(acb_mat_entry(Asys, out_row, c),
                            match_rows[src_row][c].raw());
                }
                ++out_row;
            }
            for (std::size_t r = 0; r < pair_rows.size(); ++r) {
                acb_set(acb_mat_entry(rhs, out_row, 0), pair_rhs[r].raw());
                for (long c = 0; c < n_var; ++c) {
                    acb_set(acb_mat_entry(Asys, out_row, c),
                            pair_rows[r][c].raw());
                }
                ++out_row;
            }

            // Solve via rationalize -> exact rational Dixon -> arb conversion
            // (mirrors MMA DESolver.m:906 Solve, exact over rationals).
            acb_mat_t sol_mat;
            acb_mat_init(sol_mat, n_var, 1);

            int ratdig = working_pre();
            long rat_bits = decimal_digits_to_bits(ratdig) + 32;

            fmpq_mat_t Mblock, Bblock, Sblock;
            fmpq_mat_init(Mblock, 2 * n_var, 2 * n_var);
            fmpq_mat_init(Bblock, 2 * n_var, 1);
            fmpq_mat_init(Sblock, 2 * n_var, 1);

            bool rationalize_ok = true;
            {
                arb_t scale;
                arb_init(scale);
                arb_set_si(scale, 10);
                arb_pow_ui(scale, scale,
                           static_cast<unsigned long>(ratdig),
                           rat_bits);

                fmpz_t pow10;
                fmpz_init(pow10);
                fmpz_set_si(pow10, 10);
                fmpz_pow_ui(pow10, pow10,
                            static_cast<unsigned long>(ratdig));

                arb_t scaled, part;
                arb_init(scaled);
                arb_init(part);
                fmpz_t numer;
                fmpz_init(numer);

                auto rationalize_acb = [&](acb_srcptr src,
                                           fmpq_t out_re,
                                           fmpq_t out_im) -> bool {
                    arb_set(part, acb_realref(src));
                    arb_mul(scaled, part, scale, rat_bits);
                    if (!arb_is_finite(scaled)) return false;
                    arf_get_fmpz(numer, arb_midref(scaled), ARF_RND_NEAR);
                    fmpq_set_fmpz_frac(out_re, numer, pow10);

                    arb_set(part, acb_imagref(src));
                    arb_mul(scaled, part, scale, rat_bits);
                    if (!arb_is_finite(scaled)) return false;
                    arf_get_fmpz(numer, arb_midref(scaled), ARF_RND_NEAR);
                    fmpq_set_fmpz_frac(out_im, numer, pow10);
                    return true;
                };

                fmpq_t a_re, a_im;
                fmpq_init(a_re);
                fmpq_init(a_im);
                for (long r = 0; r < n_var && rationalize_ok; ++r) {
                    for (long c = 0; c < n_var; ++c) {
                        if (!rationalize_acb(acb_mat_entry(Asys, r, c),
                                             a_re, a_im)) {
                            rationalize_ok = false;
                            break;
                        }
                        fmpq_set(fmpq_mat_entry(Mblock, r, c), a_re);
                        fmpq_neg(fmpq_mat_entry(Mblock, r, c + n_var), a_im);
                        fmpq_set(fmpq_mat_entry(Mblock, r + n_var, c), a_im);
                        fmpq_set(fmpq_mat_entry(Mblock, r + n_var, c + n_var),
                                 a_re);
                    }
                    if (!rationalize_ok) break;
                    if (!rationalize_acb(acb_mat_entry(rhs, r, 0),
                                         a_re, a_im)) {
                        rationalize_ok = false;
                        break;
                    }
                    fmpq_set(fmpq_mat_entry(Bblock, r, 0), a_re);
                    fmpq_set(fmpq_mat_entry(Bblock, r + n_var, 0), a_im);
                }
                fmpq_clear(a_re);
                fmpq_clear(a_im);
                fmpz_clear(numer);
                fmpz_clear(pow10);
                arb_clear(scaled);
                arb_clear(part);
                arb_clear(scale);
            }

            int ok = 0;
            if (rationalize_ok) {
                ok = fmpq_mat_solve_dixon(Sblock, Mblock, Bblock);
                if (ok) {
                    acb_t cur;
                    acb_init(cur);
                    for (long i = 0; i < n_var; ++i) {
                        arb_set_fmpq(acb_realref(cur),
                                     fmpq_mat_entry(Sblock, i, 0),
                                     prec);
                        arb_set_fmpq(acb_imagref(cur),
                                     fmpq_mat_entry(Sblock, i + n_var, 0),
                                     prec);
                        acb_set(acb_mat_entry(sol_mat, i, 0), cur);
                    }
                    acb_clear(cur);
                }
            }
            fmpq_mat_clear(Mblock);
            fmpq_mat_clear(Bblock);
            fmpq_mat_clear(Sblock);
            if (!ok) {
                std::ostringstream msg;
                msg << "calcx00: rationalized boundary matching system "
                    << (rationalize_ok
                        ? "singular under exact rational solve"
                        : "could not be rationalized (interval too wide?)");
                acb_mat_clear(Asys); acb_mat_clear(rhs); acb_mat_clear(sol_mat);
                throw std::runtime_error(msg.str());
            }
            for (long i = 0; i < n_var; ++i) {
                acb_set(sol_v[i].raw(), acb_mat_entry(sol_mat, i, 0));
            }
            acb_mat_clear(Asys); acb_mat_clear(rhs); acb_mat_clear(sol_mat);
        }

        // Compute f0[block, beh[m].mu, p, n] = kk[m, p, n] . variables[m] + bb[m, p, n]
        for (long m = 0; m < (long)beh.size(); ++m) {
            std::vector<AcbValue> var_m(Nblk);
            auto it = std::find(essentialset.begin(), essentialset.end(), m);
            if (it != essentialset.end()) {
                long mi = static_cast<long>(it - essentialset.begin());
                for (long k_ = 0; k_ < (long)Nblk; ++k_) {
                    acb_set(var_m[k_].raw(), sol_v[mi * Nblk + k_].raw());
                }
            }

            const long logk = beh[m].log_power;
            std::vector<std::vector<std::vector<AcbValue>>> per_integral(Nblk);
            for (long i = 0; i < (long)Nblk; ++i) {
                per_integral[i].resize(logk + 1);
                for (long p = 0; p <= logk; ++p) {
                    per_integral[i][p].resize(xorder + 1);
                }
            }

            for (long p = 0; p <= logk; ++p) {
                for (long n = 0; n <= xorder; ++n) {
                    auto& kk_pn = chain.kk_at(m, p, n);
                    auto& bb_pn = chain.bb_at(m, p, n);
                    if (kk_pn.empty()) kk_pn.resize(Nblk * Nblk);
                    if (bb_pn.empty()) bb_pn.resize(Nblk);
                    std::vector<AcbValue> kv = matvec(kk_pn, var_m, Nblk, prec);
                    std::vector<AcbValue> coeff_vec(Nblk);
                    for (long i = 0; i < (long)Nblk; ++i) {
                        acb_add(coeff_vec[i].raw(), kv[i].raw(), bb_pn[i].raw(), prec);
                    }
                    if (rot != nullptr) {
                        coeff_vec = matvec(rot->u, coeff_vec, Nblk, prec);
                    }
                    for (long i = 0; i < (long)Nblk; ++i) {
                        acb_set(per_integral[i][p][n].raw(), coeff_vec[i].raw());
                    }
                }
            }

            for (long i = 0; i < (long)Nblk; ++i) {
                AsyTerm term;
                term.mu = beh[m].mu.clone();
                term.exp.resize(logk + 1);
                for (long p = 0; p <= logk; ++p) {
                    term.exp[p] = clone_vec(per_integral[i][p]);
                }
                result[block[i]].push_back(std::move(term));
            }
        }
    }

    return result;
}

std::vector<AsyExpansion>
calcx00(const std::vector<BlockEquation>&     nheq,
        const std::vector<BlockEquationNum>&  nheqn,
        const std::vector<AcbValue>&          bc,
        acb_srcptr                            x0,
        const AsymptoticBehaviorList&         behavior,
        long                                  prec) {
    return calcx00_impl(nheq, nheqn, bc, x0, behavior, nullptr, prec);
}

// ===========================================================================
//  Learn / CalcZero / PickZero
// ===========================================================================

BlockBehavior learn_from_rule_s(const AsyExpansion& asy) {
    int test = test_x_order();
    int learn = learn_x_order();
    long upto = std::min<long>(test, learn) + 1;
    int chop = chop_pre();

    BlockBehavior out;
    for (const auto& term : asy) {
        long k = static_cast<long>(term.exp.size());
        while (k >= 1) {
            bool any = false;
            const auto& row = term.exp[k - 1];
            for (long n = 0; n < (long)row.size() && n < upto; ++n) {
                if (!acb_is_chop_zero(row[n].raw(), chop)) { any = true; break; }
            }
            if (!any) --k;
            else break;
        }
        long log_power = k - 1;
        if (log_power >= 0) {
            BehaviorEntry e;
            e.mu = term.mu.clone();
            e.log_power = log_power;
            out.push_back(std::move(e));
        }
    }
    return out;
}

AsymptoticBehaviorList
learn_from_rule_s_all(const std::vector<AsyExpansion>& asy_list,
                      const std::vector<std::vector<std::size_t>>& blocks) {
    AsymptoticBehaviorList out(blocks.size());
    long prec = working_prec_bits();

    for (std::size_t bi = 0; bi < blocks.size(); ++bi) {
        BlockBehavior joined;
        for (auto idx : blocks[bi]) {
            BlockBehavior contrib = learn_from_rule_s(asy_list[idx]);
            for (auto& e : contrib) joined.push_back(BehaviorEntry{e.mu.clone(), e.log_power});
        }
        BlockBehavior unioned;
        std::vector<bool> done(joined.size(), false);
        for (std::size_t i = 0; i < joined.size(); ++i) {
            if (done[i]) continue;
            BehaviorEntry combined;
            combined.mu = joined[i].mu.clone();
            combined.log_power = joined[i].log_power;
            done[i] = true;
            for (std::size_t j = i + 1; j < joined.size(); ++j) {
                if (done[j]) continue;
                AcbValue diff;
                acb_sub(diff.raw(), joined[i].mu.raw(), joined[j].mu.raw(), prec);
                if (acb_is_chop_zero(diff.raw(), chop_pre())) {
                    if (joined[j].log_power > combined.log_power)
                        combined.log_power = joined[j].log_power;
                    done[j] = true;
                }
            }
            unioned.push_back(std::move(combined));
        }
        out[bi] = std::move(unioned);
    }
    return out;
}

std::vector<AsyExpansion>
calc_zero(const RationalMatrix& de,
          const std::vector<AcbValue>& bc,
          acb_srcptr x0,
          long prec) {
    if (!silent_mode()) log_line("CalcZero: starting");

    NormalizationResult norm;
    std::vector<BlockJordanRotation> rotations;
    bool use_block_rotations = false;
    bool use_acb_shearing = false;
    try {
        norm = normalize_mat(de, prec);
    } catch (const std::exception& e) {
        const std::string msg = e.what();
        if (msg.find("algebraic Jordan rotation") == std::string::npos &&
            msg.find("algebraic leading residue") == std::string::npos &&
            msg.find("algebraic-number shearing") == std::string::npos) {
            throw;
        }
        if (!silent_mode()) {
            log_line("CalcZero: falling back to block-local algebraic Jordan rotations");
        }
        auto fallback = normalize_mat_for_calc_zero(de, prec);
        norm = std::move(fallback.base);
        rotations = std::move(fallback.rotations);
        use_block_rotations = true;
        use_acb_shearing =
            (msg.find("algebraic-number shearing") != std::string::npos);
        if (use_acb_shearing) {
            return calc_zero_acb_shearing_single_block(norm, bc, x0, prec);
        }
    }

    std::vector<AcbValue> bcT(de.rows());
    for (std::size_t i = 0; i < de.rows(); ++i) {
        AcbValue acc;
        for (std::size_t j = 0; j < de.cols(); ++j) {
            std::ostringstream where;
            where << "calc_zero: invT(x0) entry=(" << i << "," << j << ")";
            AcbValue tij = evaluate_with_context(norm.invT(i, j), x0,
                                                  prec, where.str());
            AcbValue tmp;
            acb_mul(tmp.raw(), tij.raw(), bc[j].raw(), prec);
            acb_add(acc.raw(), acc.raw(), tmp.raw(), prec);
        }
        bcT[i] = std::move(acc);
    }

    auto nheq  = nh_equations(norm.B, EquationMode::Singular);
    auto nheqn = nh_equations_num(nheq, prec);

    AsymptoticBehaviorList behavior = use_block_rotations
        ? asymptotic_behavior_for_calc_zero(norm.B, rotations, prec)
        : asymptotic_behavior(norm.B, prec);

    if (learn_x_order() >= 0) {
        numeric::GlobalScope guard;
        guard.expansion.x_order = learn_x_order();
        guard.global.silent_mode = true;
        guard.commit();
        auto trial = use_block_rotations
            ? calcx00_impl(nheq, nheqn, bcT, x0, behavior, &rotations, prec)
            : calcx00_impl(nheq, nheqn, bcT, x0, behavior, nullptr, prec);
        auto blocks = analyze_block(norm.B);
        behavior = learn_from_rule_s_all(trial, blocks);
    }

    auto raw = use_block_rotations
        ? calcx00_impl(nheq, nheqn, bcT, x0, behavior, &rotations, prec)
        : calcx00_impl(nheq, nheqn, bcT, x0, behavior, nullptr, prec);

    auto T_ps = to_power_series_matrix(norm.T, prec);

    auto out = ps_map_rule_set(T_ps, raw, prec);
    if (!silent_mode()) log_line("CalcZero: finished");
    return out;
}

AcbValue pick_zero_rule_s(const AsyExpansion& asy, long prec) {
    std::vector<long> integer_mus;
    std::vector<std::size_t> integer_idx;

    int chop = chop_pre();
    if (chop <= 0) chop = 20;
    const double eps_d = std::pow(10.0, -chop);

    for (std::size_t i = 0; i < asy.size(); ++i) {
        arf_t mid_im;
        arf_init(mid_im);
        arf_abs(mid_im, arb_midref(acb_imagref(asy[i].mu.raw())));
        bool im_zero = (arf_cmp_d(mid_im, eps_d) < 0);
        arf_clear(mid_im);
        if (!im_zero) continue;

        fmpz_t z;
        fmpz_init(z);
        arf_get_fmpz(z, arb_midref(acb_realref(asy[i].mu.raw())), ARF_RND_NEAR);

        arb_t resid;
        arb_init(resid);
        arb_set_fmpz(resid, z);
        arb_sub(resid, acb_realref(asy[i].mu.raw()), resid, prec);
        arf_t mr;
        arf_init(mr);
        arf_abs(mr, arb_midref(resid));
        bool re_int = (arf_cmp_d(mr, eps_d) < 0);
        arf_clear(mr); arb_clear(resid);

        if (re_int) {
            integer_mus.push_back(fmpz_get_si(z));
            integer_idx.push_back(i);
        }
        fmpz_clear(z);
    }

    if (integer_mus.empty()) return AcbValue();
    if (integer_mus.size() > 1) {
        throw std::runtime_error(
            "pick_zero_rule_s: multiple integer mu values found; "
            "the asymptotic regions have not been fully merged");
    }

    long key = integer_mus[0];
    std::size_t which = integer_idx[0];
    if (key > 0) return AcbValue();

    if (key == 0) {
        const auto& exp = asy[which].exp;
        if (exp.empty() || exp[0].empty()) return AcbValue();
        return exp[0][0].clone();
    }

    const auto& exp = asy[which].exp;
    if (exp.empty()) return AcbValue();
    long want = -key;
    if (want >= (long)exp[0].size()) return AcbValue();
    return exp[0][want].clone();
}

}  // namespace amflow::ode
