// SPDX-License-Identifier: MIT
// qft::vacuum — single-mass vacuum integrals (Layer 14f).
//
//
// Mirrors AMFlow.m lines 819-824.  Available (loop, prop) pairs:
//   (1,1) (2,3) (3,4) (3,5) (4,5).

#ifndef AMFLOW_QFT_VACUUM_HPP
#define AMFLOW_QFT_VACUUM_HPP

#include <flint/acb.h>

#include "amflow/numeric/acb_wrap.hpp"

namespace amflow::qft {

void vacuum(long loop_count, long prop_count,
            const acb_t eps, acb_t out, long prec);

inline numeric::AcbValue vacuum(long loop_count, long prop_count,
                                 const acb_t eps, long prec) {
    numeric::AcbValue v;
    vacuum(loop_count, prop_count, eps, v.raw(), prec);
    return v;
}

bool vacuum_known(long loop_count, long prop_count) noexcept;

}  // namespace amflow::qft

#endif  // AMFLOW_QFT_VACUUM_HPP
