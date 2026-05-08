// SPDX-License-Identifier: MIT
// numeric::acb_wrap — implementation.
//

#include "amflow/numeric/acb_wrap.hpp"

#include <cmath>
#include <ostream>
#include <string>

#include "amflow/numeric/log.hpp"

#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/arf.h>

namespace amflow::numeric {

// ---------------------------------------------------------------------------
//  Chop predicate.  See comment in acb_wrap.hpp / numeric.md §8c.
// ---------------------------------------------------------------------------

bool acb_is_chop_zero(acb_srcptr x, int chop_digits) {
    if (chop_digits <= 0) return false;

    arf_t mid_re, mid_im;
    arf_init(mid_re); arf_init(mid_im);
    arf_abs(mid_re, arb_midref(acb_realref(x)));
    arf_abs(mid_im, arb_midref(acb_imagref(x)));

    bool result;
    if (chop_digits <= 300) {
        // Double can represent 10^-300 with full mantissa precision; for
        // larger digit counts we fall through to the arb_t path below.
        const double eps = std::pow(10.0, -chop_digits);
        result = (arf_cmp_d(mid_re, eps) < 0) && (arf_cmp_d(mid_im, eps) < 0);
    } else {
        arb_t eps_arb;
        arb_init(eps_arb);
        arb_set_si(eps_arb, 10);
        arb_inv(eps_arb, eps_arb, decimal_digits_to_bits(chop_digits) + 32);
        arb_pow_ui(eps_arb, eps_arb, static_cast<unsigned long>(chop_digits),
                   decimal_digits_to_bits(chop_digits) + 32);

        arf_srcptr eps_mid = arb_midref(eps_arb);
        result = (arf_cmp(mid_re, eps_mid) < 0) && (arf_cmp(mid_im, eps_mid) < 0);

        arb_clear(eps_arb);
    }
    arf_clear(mid_re); arf_clear(mid_im);
    return result;
}

bool AcbValue::is_chop_zero(int chop_digits) const {
    return acb_is_chop_zero(handle_, chop_digits);
}

// ---------------------------------------------------------------------------
//  Pretty printing
// ---------------------------------------------------------------------------

std::string AcbValue::to_string(long display_digits) const {
    // FLINT 3 dropped acb_get_str.  Build "(re + im*I)" by formatting the
    // two arb halves independently.
    //
    // AMFLOW_DEBUG_PRECISE=1 sets ARB_STR_MORE|ARB_STR_NO_RADIUS, which
    // prints midpoint digits beyond the ball radius (diagnostic only —
    // untrusted digits).
    int arb_flags = 0;
    if (::amflow::numeric::log::trace_enabled("AMFLOW_DEBUG_PRECISE")) arb_flags = 1 | 2;
    char* re = arb_get_str(acb_realref(handle_), display_digits, arb_flags);
    char* im = arb_get_str(acb_imagref(handle_), display_digits, arb_flags);
    std::string out = "(";
    out += (re != nullptr) ? re : "<null>";
    out += " + ";
    out += (im != nullptr) ? im : "<null>";
    out += "*I)";
    if (re != nullptr) flint_free(re);
    if (im != nullptr) flint_free(im);
    return out;
}

std::ostream& operator<<(std::ostream& os, const AcbValue& v) {
    return os << v.to_string(20);
}

}  // namespace amflow::numeric
