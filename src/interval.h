// del3d - directed-rounding interval arithmetic on doubles.
//
// Why a library that computes exactly also needs an *inexact* number type: the
// dual (circumcentre) values consumed downstream are produced by a lazy exact
// pipeline, which evaluates an expression in interval arithmetic first and only
// falls back on exact arithmetic when the interval is too wide. Its answer is
// therefore the midpoint of an interval, not the exact value, whenever the
// interval is already accurate to 1e-5 relative. On a sliver the circumcentre
// is enormous and the interval correspondingly wide in absolute terms, so that
// midpoint can differ from the exact circumcentre by ~1e-9 relative - and those
// far-away duals are exactly the ones that decide a nearest-Voronoi-vertex
// query. Reproducing such a pipeline bit for bit therefore means reproducing
// its *approximation*, not just the exact answer; this file is the arithmetic
// that circumcenter_lazy() in triangulation.cpp needs to do that.
//
// The whole file is only valid while the FPU rounds toward +infinity; see
// RoundingUpward below.
#ifndef DEL3D_INTERVAL_H
#define DEL3D_INTERVAL_H

#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>

namespace del3d {
namespace detail {

/// Scoped switch of the FPU rounding mode to "toward +infinity", restoring the
/// previous mode on destruction. The interval operators below store the lower
/// bound negated, so that a single upward-rounding operation widens *both*
/// ends outwards; with any other rounding mode they would produce intervals
/// that do not actually contain the true result.
class RoundingUpward {
public:
    RoundingUpward() : saved_(std::fegetround()) { std::fesetround(FE_UPWARD); }
    ~RoundingUpward() { std::fesetround(saved_); }
    RoundingUpward(const RoundingUpward&) = delete;
    RoundingUpward& operator=(const RoundingUpward&) = delete;
private:
    int saved_;
};

/// Forces a value through memory. Without this the optimiser is free to
/// reassociate the surrounding arithmetic, keep it in a wider register, or fold
/// it at compile time under the default rounding mode - any of which silently
/// invalidates the directed-rounding bounds.
inline double opacify(double x) {
    volatile double v = x;
    return v;
}

/**
 * A closed interval [inf, sup] guaranteed to contain the exact value.
 *
 * The lower bound is stored negated (ni_ = -inf) so that every bound - upper
 * and lower alike - is produced by one round-toward-+infinity operation, which
 * always widens the interval and never narrows it.
 */
class Interval {
public:
    Interval() : ni_(0.0), s_(0.0) {}
    /// A degenerate interval containing exactly one representable value.
    explicit Interval(double d) : ni_(-d), s_(d) {}
    /// Takes the true bounds (inf, sup); the lower one is negated on the way in.
    Interval(double i, double s) : ni_(-i), s_(s) {}

    double inf() const { return -ni_; }
    double sup() const { return s_; }
    bool is_point() const { return sup() == inf(); }

    /// The whole real line: what division by an interval straddling zero
    /// returns, since no finite bound is valid there.
    static Interval largest() {
        const double inf = std::numeric_limits<double>::infinity();
        return Interval(-inf, inf);
    }

    friend Interval operator+(const Interval& a, const Interval& b) {
        return Interval(-opacify(opacify(-a.inf()) + opacify(-b.inf())),
                        opacify(opacify(a.sup()) + opacify(b.sup())));
    }

    friend Interval operator-(const Interval& a, const Interval& b) {
        return Interval(-opacify(opacify(-a.inf()) + opacify(b.sup())),
                        opacify(opacify(a.sup()) + opacify(-b.inf())));
    }

    /// Interval product. Branching on the sign of each operand picks the pair
    /// of endpoint products that bound the result, so only two multiplications
    /// are needed unless both intervals straddle zero, in which case all four
    /// corner products are compared.
    friend Interval operator*(const Interval& a, const Interval& b) {
        if (a.inf() >= 0.0) {                              // a >= 0
            double aa = a.inf(), bb = a.sup();
            if (bb <= 0.0) return Interval(0.0);
            if (b.inf() < 0.0) {
                aa = bb;
                if (b.sup() < 0.0) bb = a.inf();
            }
            const double r = (b.sup() == 0) ? 0.0 : mul(bb, b.sup());
            return Interval(-mul(aa, -b.inf()), r);
        }
        if (a.sup() <= 0.0) {                              // a <= 0
            double aa = a.sup(), bb = a.inf();
            if (b.inf() < 0.0) {
                aa = bb;
                if (b.sup() <= 0.0) bb = a.sup();
            } else if (b.sup() <= 0) {
                return Interval(0.0);
            }
            return Interval(-mul(-bb, b.sup()), mul(-aa, -b.inf()));
        }
        // 0 in a
        if (b.inf() >= 0.0) {                              // b >= 0
            if (b.sup() <= 0.0) return Interval(0.0);
            return Interval(-mul(-a.inf(), b.sup()), mul(a.sup(), b.sup()));
        }
        if (b.sup() <= 0.0) {                              // b <= 0
            return Interval(-mul(a.sup(), -b.inf()), mul(-a.inf(), -b.inf()));
        }
        // 0 in a and 0 in b: the extremes can come from any corner.
        const double tmp1 = mul(-a.inf(), b.sup());
        const double tmp2 = mul(a.sup(), -b.inf());
        const double tmp3 = mul(-a.inf(), -b.inf());
        const double tmp4 = mul(a.sup(), b.sup());
        return Interval(-(std::max)(tmp1, tmp2), (std::max)(tmp3, tmp4));
    }

    /// Interval quotient, by the same sign analysis as the product. If the
    /// divisor contains zero the quotient is unbounded, hence largest().
    friend Interval operator/(const Interval& a, const Interval& b) {
        if (b.inf() > 0.0) {                               // b > 0
            double aa = b.sup(), bb = b.inf();
            if (a.inf() < 0.0) {
                aa = bb;
                if (a.sup() < 0.0) bb = b.sup();
            }
            return Interval(-div(-a.inf(), aa), div(a.sup(), bb));
        }
        if (b.sup() < 0.0) {                               // b < 0
            double aa = b.sup(), bb = b.inf();
            if (a.inf() < 0.0) {
                bb = aa;
                if (a.sup() < 0.0) aa = b.inf();
            }
            return Interval(-div(a.sup(), -aa), div(a.inf(), bb));
        }
        return largest();                                  // 0 in b
    }

    /// Square. Tighter than x * x, which cannot know that both factors are the
    /// same number: for an interval straddling zero the lower bound is 0, not
    /// the negative corner product. The distinction is observable, so the
    /// expression being reproduced must be squared with this and not with `*`.
    friend Interval square(const Interval& d) {
        if (d.inf() >= 0.0) return Interval(-mul(-d.inf(), d.inf()), mul(d.sup(), d.sup()));
        if (d.sup() <= 0.0) return Interval(-mul(d.sup(), -d.sup()), mul(-d.inf(), -d.inf()));
        const double m = (std::max)(-d.inf(), d.sup());
        return Interval(0.0, mul(m, m));
    }

    /// The representative double for this interval: its midpoint.
    double to_double() const { return (sup() + inf()) * 0.5; }

    double magnitude() const { return (std::max)(std::fabs(inf()), std::fabs(sup())); }
    double width() const { return sup() - inf(); }
    double radius() const { return width() / 2; }

    /// True iff the interval pins the value down to the given relative
    /// precision - the test that decides whether the midpoint is good enough
    /// or the exact fallback is required.
    bool has_smaller_relative_precision(double prec) const {
        return magnitude() == 0 || radius() < prec * magnitude();
    }

private:
    static double mul(double a, double b) { return opacify(opacify(a) * opacify(b)); }
    static double div(double a, double b) { return opacify(opacify(a) / opacify(b)); }

    double ni_;   // -inf
    double s_;    //  sup
};

} // namespace detail
} // namespace del3d

#endif // DEL3D_INTERVAL_H
