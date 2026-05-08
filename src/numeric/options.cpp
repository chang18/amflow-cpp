// SPDX-License-Identifier: MIT
// numeric::Options + GlobalScope + logging — implementation.
//
// amflow::numeric.

#include "amflow/numeric/options.hpp"

#include <cmath>
#include <iostream>
#include <ostream>
#include <stdexcept>

#include <flint/fmpq.h>

namespace amflow::numeric {

namespace {

// Process-wide singletons via function-local statics (avoids static-init order).
GlobalOptions& mutable_global() {
    static GlobalOptions g{};
    return g;
}
ExpansionOptions& mutable_expansion() {
    static ExpansionOptions e{};
    return e;
}
RunningOptions& mutable_running() {
    static RunningOptions r{};
    return r;
}

constexpr double kLog2_10 = 3.3219280948873623478703194294894;

}  // namespace

// ---------------------------------------------------------------------------
GlobalOptions&    global_options()    { return mutable_global();    }
ExpansionOptions& expansion_options() { return mutable_expansion(); }
RunningOptions&   running_options()   { return mutable_running();   }

void set_global_options(const GlobalOptions& opt)       { mutable_global()    = opt; }
void set_expansion_options(const ExpansionOptions& opt) { mutable_expansion() = opt; }
void set_running_options(const RunningOptions& opt)     { mutable_running()   = opt; }

void set_default_options() {
    mutable_global()    = GlobalOptions{};
    mutable_expansion() = ExpansionOptions{};
    mutable_running()   = RunningOptions{};
}

// ---------------------------------------------------------------------------
long decimal_digits_to_bits(int digits) {
    if (digits <= 0) return 53;
    long bits = static_cast<long>(std::ceil(static_cast<double>(digits) * kLog2_10));
    if (bits < 53) bits = 53;
    return bits;
}

long working_prec_bits() {
    return decimal_digits_to_bits(working_pre());
}

// ---------------------------------------------------------------------------
void d0_eps_shift_fmpq(fmpq_t out) {
    const std::string& s = global_options().d0;

    fmpq_t d;
    fmpq_init(d);
    if (fmpq_set_str(d, s.c_str(), 10) != 0) {
        fmpq_clear(d);
        throw std::runtime_error(
            "d0_eps_shift_fmpq: failed to parse D0 rational from '" + s + "'");
    }

    // out = (4 - D0) / 2
    fmpq_t four;
    fmpq_init(four);
    fmpq_set_si(four, 4, 1);
    fmpq_sub(out, four, d);

    fmpq_t half;
    fmpq_init(half);
    fmpq_set_si(half, 1, 2);
    fmpq_mul(out, out, half);

    fmpq_clear(half);
    fmpq_clear(four);
    fmpq_clear(d);
}

// ---------------------------------------------------------------------------
GlobalScope::GlobalScope()
    : global(global_options()),
      expansion(expansion_options()),
      running(running_options()),
      saved_global_(global_options()),
      saved_expansion_(expansion_options()),
      saved_running_(running_options()) {}

GlobalScope::~GlobalScope() {
    if (committed_) {
        mutable_global()    = saved_global_;
        mutable_expansion() = saved_expansion_;
        mutable_running()   = saved_running_;
    }
}

void GlobalScope::commit() {
    if (committed_) return;
    committed_         = true;
    mutable_global()    = global;
    mutable_expansion() = expansion;
    mutable_running()   = running;
}

// ---------------------------------------------------------------------------
std::ostream& log_stream() {
    return std::clog;
}

void log_line(const std::string& msg) {
    if (silent_mode()) return;
    log_stream() << msg << '\n';
}

}  // namespace amflow::numeric
