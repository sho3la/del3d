// del3d - exact arithmetic on non-overlapping floating-point expansions.
//
// A value is represented as a sum of doubles ("components") whose magnitudes
// are strictly increasing and which do not overlap in significance: no two
// components have significand bits at the same binary exponent. The sum is
// therefore exact - it loses nothing that the individual doubles hold - and the
// sign of the whole value is simply the sign of the last, largest, component.
//
// Every operation below is built from error-free transformations: pairs of
// IEEE double operations that return not only the rounded result but also the
// exact round-off it discarded. No extended-precision type, no multiprecision
// integers, no allocation beyond the component vector.
//
// This is the exact fallback stage of the predicates in predicates.cpp and
// predicates2.cpp; the filtered first stage answers the vast majority of calls
// and never reaches this file.
#ifndef DEL3D_EXPANSION_H
#define DEL3D_EXPANSION_H

#include "del3d/predicates.h"

#include <cstddef>
#include <vector>

namespace del3d {
namespace detail {

// ---- error-free transformations ------------------------------------------

/// Splits a + b into x + err with x = fl(a+b) and err the exact round-off, so
/// that x + err == a + b exactly. Makes no assumption about the relative
/// magnitudes of a and b.
inline void two_sum(double a, double b, double& x, double& err) {
    x = a + b;
    const double bv = x - a;
    const double av = x - bv;
    err = (a - av) + (b - bv);
}

/// Splits a - b into x + err exactly. Written as two_sum(a, -b): negation is
/// exact in IEEE-754, so this loses nothing, and it avoids the sign slips the
/// hand-expanded difference form invites.
inline void two_diff(double a, double b, double& x, double& err) {
    two_sum(a, -b, x, err);
}

/// Splits a into two halves hi + lo == a with disjoint significands, each
/// holding at most 26 significant bits, so that products of halves are exact.
inline void split(double a, double& hi, double& lo) {
    const double c = 134217729.0 * a;          // 2^27 + 1
    const double abig = c - a;
    hi = c - abig;
    lo = a - hi;
}

/// Splits a * b into x + err with x = fl(a*b), exactly, by multiplying the
/// four half-products above and subtracting them off the rounded result.
inline void two_product(double a, double b, double& x, double& err) {
    x = a * b;
    double ahi, alo, bhi, blo;
    split(a, ahi, alo);
    split(b, bhi, blo);
    const double e1 = x - (ahi * bhi);
    const double e2 = e1 - (alo * bhi);
    const double e3 = e2 - (ahi * blo);
    err = (alo * blo) - e3;
}

// ---- expansions ----------------------------------------------------------

/// An exact value: the sum of its components, ordered by increasing magnitude
/// and mutually non-overlapping. Zero components are dropped as they appear, so
/// an empty expansion is exactly zero and is_zero() is an exact test.
class Expansion {
public:
    Expansion() = default;

    explicit Expansion(double v) { if (v != 0.0) c_.push_back(v); }

    /// The exact difference a - b of two doubles, as a two-component value.
    static Expansion from_diff(double a, double b) {
        double x, e;
        two_diff(a, b, x, e);
        Expansion r;
        if (e != 0.0) r.c_.push_back(e);
        if (x != 0.0) r.c_.push_back(x);
        return r;
    }

    bool is_zero() const { return c_.empty(); }

    /// The nearest double to the exact value. Since the components are
    /// non-overlapping and ordered by increasing magnitude, summing them in
    /// that order - smallest first - cannot cancel, so one pass suffices.
    double to_double() const {
        double s = 0.0;
        for (double v : c_) s += v;
        return s;
    }

    /// Exact sign: the largest component dominates all the others combined.
    Sign sign() const {
        if (c_.empty()) return ZERO;
        const double top = c_.back();            // largest component
        return top > 0.0 ? POSITIVE : NEGATIVE;
    }

    /// Exact sum: absorb the other value's components one at a time, keeping
    /// the result normalised at every step.
    Expansion operator+(const Expansion& o) const {
        Expansion r = *this;
        for (double b : o.c_) r.grow(b);
        return r;
    }

    /// Exact difference. Negating every component is exact and preserves both
    /// the ordering and the non-overlap property, so it needs no renormalising.
    Expansion operator-(const Expansion& o) const {
        Expansion neg;
        neg.c_.reserve(o.c_.size());
        for (double b : o.c_) neg.c_.push_back(-b);
        return *this + neg;
    }

    /// Exact product: distribute over the other value's components, scaling
    /// this whole expansion by one double at a time and accumulating.
    Expansion operator*(const Expansion& o) const {
        Expansion acc;
        for (double b : o.c_) acc = acc + scaled(b);
        return acc;
    }

private:
    /// Adds a single double, re-normalising: run a two_sum chain through the
    /// components, emitting each round-off as a new small component and
    /// carrying the running sum upwards.
    void grow(double b) {
        if (b == 0.0) return;
        std::vector<double> out;
        out.reserve(c_.size() + 1);
        double q = b;
        for (double a : c_) {
            double sum, err;
            two_sum(q, a, sum, err);
            if (err != 0.0) out.push_back(err);
            q = sum;
        }
        if (q != 0.0) out.push_back(q);
        c_.swap(out);
    }

    /// Multiplies by a single double: each component contributes an exact
    /// two_product, whose two halves are merged into the running expansion.
    /// The result has at most twice as many components as the input.
    Expansion scaled(double b) const {
        Expansion r;
        if (b == 0.0 || c_.empty()) return r;
        r.c_.reserve(2 * c_.size());
        double q, hh;
        two_product(c_[0], b, q, hh);
        if (hh != 0.0) r.c_.push_back(hh);
        for (std::size_t i = 1; i < c_.size(); ++i) {
            double pi, pe;
            two_product(c_[i], b, pi, pe);
            double s, e;
            two_sum(q, pe, s, e);
            if (e != 0.0) r.c_.push_back(e);
            double t, te;
            two_sum(pi, s, t, te);
            if (te != 0.0) r.c_.push_back(te);
            q = t;
        }
        if (q != 0.0) r.c_.push_back(q);
        return r;
    }

    std::vector<double> c_;
};

} // namespace detail
} // namespace del3d

#endif // DEL3D_EXPANSION_H
