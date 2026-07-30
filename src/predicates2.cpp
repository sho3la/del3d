// del3d - the coplanar predicates. See include/del3d/predicates2.h.
//
// Same two-stage scheme as predicates.cpp: a filtered floating-point evaluation
// with an error bound, falling back on exact expansion arithmetic when the sign
// is in doubt.
//
// The projection choice, the sign conventions and the perturbation rule are all
// spelled out explicitly rather than derived per call site. Any self-consistent
// choice would give *a* valid triangulation, but only one fixed choice gives a
// reproducible one, so these details are part of the specification.
#include "del3d/predicates2.h"
#include "expansion.h"

#include <algorithm>
#include <cmath>

namespace del3d {
namespace {

using detail::Expansion;

constexpr double kEps = 1.1102230246251565e-16;
// Error bound for the 2x2 orientation determinant.
constexpr double kOrient2dBound = (3.0 + 16.0 * kEps) * kEps;
// The in-circle test below is a 4x4 determinant of the same shape as the
// in-sphere one, so it reuses that bound.
constexpr double kInCircleBound = (16.0 + 224.0 * kEps) * kEps;

inline Sign sign_of(double d) { return d > 0.0 ? POSITIVE : (d < 0.0 ? NEGATIVE : ZERO); }
/// Composes two orientation signs (used to re-sign an answer against a
/// reference triple whose own orientation may be negative).
inline Sign mul(Sign a, Sign b) { return Sign(int(a) * int(b)); }

/// |q-p, r-p| in exact expansion arithmetic.
Sign orient2d_exact(double px, double py, double qx, double qy, double rx, double ry) {
    const Expansion qpx = Expansion::from_diff(qx, px);
    const Expansion qpy = Expansion::from_diff(qy, py);
    const Expansion rpx = Expansion::from_diff(rx, px);
    const Expansion rpy = Expansion::from_diff(ry, py);
    return (qpx * rpy - qpy * rpx).sign();
}

/**
 * The in-circle determinant for three coplanar points and a query point, in
 * exact expansion arithmetic.
 *
 * Working in the plane without projecting: translate so that t is the origin,
 * lift p, q, r onto the paraboloid (offset plus squared length), and add a
 * fourth row v = pq x pr, the plane's normal, lifted by |v|^2. That normal row
 * confines the 4x4 determinant to the plane, making it the in-circle test of
 * the three points' circumcircle rather than a sphere test.
 *
 * Two conventions below are load-bearing and must not be "simplified":
 *   * the rows are ordered (pt, rt, qt, v) - r before q;
 *   * each 2x2 minor is formed as m01 = a10*a01 - a00*a11, i.e. the opposite
 *     operand order to the usual a00*a11 - a10*a01, which negates the whole
 *     determinant.
 * Both are absorbed into the sign convention of the public function.
 */
Sign in_circle_exact(const double* p, const double* q, const double* r, const double* t) {
    auto E = [](double a, double b) { return Expansion::from_diff(a, b); };

    const Expansion ptx = E(p[0], t[0]), pty = E(p[1], t[1]), ptz = E(p[2], t[2]);
    const Expansion qtx = E(q[0], t[0]), qty = E(q[1], t[1]), qtz = E(q[2], t[2]);
    const Expansion rtx = E(r[0], t[0]), rty = E(r[1], t[1]), rtz = E(r[2], t[2]);

    const Expansion pt2 = ptx * ptx + pty * pty + ptz * ptz;
    const Expansion qt2 = qtx * qtx + qty * qty + qtz * qtz;
    const Expansion rt2 = rtx * rtx + rty * rty + rtz * rtz;

    const Expansion pqx = E(q[0], p[0]), pqy = E(q[1], p[1]), pqz = E(q[2], p[2]);
    const Expansion prx = E(r[0], p[0]), pry = E(r[1], p[1]), prz = E(r[2], p[2]);

    // v = pq x pr, the plane normal.
    const Expansion vx = pqy * prz - pqz * pry;
    const Expansion vy = pqz * prx - pqx * prz;
    const Expansion vz = pqx * pry - pqy * prx;
    const Expansion v2 = vx * vx + vy * vy + vz * vz;

    const Expansion a00 = ptx, a01 = pty, a02 = ptz, a03 = pt2;
    const Expansion a10 = rtx, a11 = rty, a12 = rtz, a13 = rt2;
    const Expansion a20 = qtx, a21 = qty, a22 = qtz, a23 = qt2;
    const Expansion a30 = vx,  a31 = vy,  a32 = vz,  a33 = v2;

    // Laplace expansion by the first two columns: 2x2 minors, then 3x3, then
    // the full determinant along the last column.
    const Expansion m01 = a10 * a01 - a00 * a11;
    const Expansion m02 = a20 * a01 - a00 * a21;
    const Expansion m03 = a30 * a01 - a00 * a31;
    const Expansion m12 = a20 * a11 - a10 * a21;
    const Expansion m13 = a30 * a11 - a10 * a31;
    const Expansion m23 = a30 * a21 - a20 * a31;

    const Expansion m012 = m12 * a02 - m02 * a12 + m01 * a22;
    const Expansion m013 = m13 * a02 - m03 * a12 + m01 * a32;
    const Expansion m023 = m23 * a02 - m03 * a22 + m02 * a32;
    const Expansion m123 = m23 * a12 - m13 * a22 + m12 * a32;

    const Expansion det = m123 * a03 - m023 * a13 + m013 * a23 - m012 * a33;
    return det.sign();
}

inline bool lexico_less(const double* a, const double* b) {
    if (a[0] != b[0]) return a[0] < b[0];
    if (a[1] != b[1]) return a[1] < b[1];
    return a[2] < b[2];
}

} // namespace

Sign orient2d(double px, double py, double qx, double qy, double rx, double ry) {
    const double qpx = qx - px, qpy = qy - py;
    const double rpx = rx - px, rpy = ry - py;
    const double a = qpx * rpy, b = qpy * rpx;
    const double det = a - b;
    const double permanent = std::fabs(a) + std::fabs(b);
    const double errbound = kOrient2dBound * permanent;
    if (det > errbound || -det > errbound) return sign_of(det);
    return orient2d_exact(px, py, qx, qy, rx, ry);
}

Sign coplanar_orientation(const double* p, const double* q, const double* r) {
    // Try the xy projection, then yz, then xz, and take the first in which the
    // triple is not collinear. A degenerate projection returns ZERO and simply
    // falls through; since the triple is non-collinear in 3D (precondition), at
    // least one projection is decisive.
    Sign oxy = orient2d(p[0], p[1], q[0], q[1], r[0], r[1]);
    if (oxy != ZERO) return oxy;
    Sign oyz = orient2d(p[1], p[2], q[1], q[2], r[1], r[2]);
    if (oyz != ZERO) return oyz;
    return orient2d(p[0], p[2], q[0], q[2], r[0], r[2]);
}

Sign coplanar_orientation(const double* p, const double* q,
                          const double* r, const double* s) {
    // The reference triple (p,q,r) selects the projection; s is then measured
    // in that same projection and the two signs are multiplied, so the answer
    // is "s on the same side of line (p,q) as r" regardless of which projection
    // was picked or how it flips handedness.
    Sign oxy = orient2d(p[0], p[1], q[0], q[1], r[0], r[1]);
    if (oxy != ZERO) return mul(oxy, orient2d(p[0], p[1], q[0], q[1], s[0], s[1]));
    Sign oyz = orient2d(p[1], p[2], q[1], q[2], r[1], r[2]);
    if (oyz != ZERO) return mul(oyz, orient2d(p[1], p[2], q[1], q[2], s[1], s[2]));
    Sign oxz = orient2d(p[0], p[2], q[0], q[2], r[0], r[2]);
    return mul(oxz, orient2d(p[0], p[2], q[0], q[2], s[0], s[2]));
}

bool collinear(const double* p, const double* q, const double* r) {
    // Collinear in 3D means all three axis projections are degenerate at once.
    // This cannot be phrased as coplanar_orientation(p,q,r) == ZERO: that
    // predicate has non-collinearity as its precondition and returns the sign
    // of whichever projection is decisive, so it never reports collinearity.
    if (orient2d(r[0], r[1], p[0], p[1], q[0], q[1]) != ZERO) return false;
    if (orient2d(r[0], r[2], p[0], p[2], q[0], q[2]) != ZERO) return false;
    return orient2d(r[1], r[2], p[1], p[2], q[1], q[2]) == ZERO;
}

bool equal(const double* p, const double* q) {
    return p[0] == q[0] && p[1] == q[1] && p[2] == q[2];
}

CollinearPosition collinear_position(const double* p, const double* q, const double* r) {
    // Assuming p, q, r collinear, the ordering along the line is decided by any
    // coordinate in which p and r differ, so compare on the first such one.
    if (equal(q, p)) return CP_SOURCE;
    if (equal(q, r)) return CP_TARGET;
    int k = 0;
    if (p[0] != r[0]) k = 0;
    else if (p[1] != r[1]) k = 1;
    else k = 2;
    const double a = p[k], b = q[k], c = r[k];
    if (a < c) {                 // the line runs in the increasing direction
        if (b < a) return CP_BEFORE;
        if (b > c) return CP_AFTER;
        return CP_MIDDLE;
    }
    if (b > a) return CP_BEFORE; // ... and here in the decreasing one
    if (b < c) return CP_AFTER;
    return CP_MIDDLE;
}

Sign coplanar_in_circle_unperturbed(const double* p0, const double* p1,
                                    const double* p2, const double* p) {
    // Filtered evaluation of exactly the determinant in_circle_exact() builds,
    // including its row order and minor convention.
    const double ptx = p0[0] - p[0], pty = p0[1] - p[1], ptz = p0[2] - p[2];
    const double qtx = p1[0] - p[0], qty = p1[1] - p[1], qtz = p1[2] - p[2];
    const double rtx = p2[0] - p[0], rty = p2[1] - p[1], rtz = p2[2] - p[2];
    const double pt2 = ptx * ptx + pty * pty + ptz * ptz;
    const double qt2 = qtx * qtx + qty * qty + qtz * qtz;
    const double rt2 = rtx * rtx + rty * rty + rtz * rtz;

    const double pqx = p1[0] - p0[0], pqy = p1[1] - p0[1], pqz = p1[2] - p0[2];
    const double prx = p2[0] - p0[0], pry = p2[1] - p0[1], prz = p2[2] - p0[2];
    const double vx = pqy * prz - pqz * pry;
    const double vy = pqz * prx - pqx * prz;
    const double vz = pqx * pry - pqy * prx;
    const double v2 = vx * vx + vy * vy + vz * vz;

    const double a00 = ptx, a01 = pty, a02 = ptz, a03 = pt2;
    const double a10 = rtx, a11 = rty, a12 = rtz, a13 = rt2;
    const double a20 = qtx, a21 = qty, a22 = qtz, a23 = qt2;
    const double a30 = vx,  a31 = vy,  a32 = vz,  a33 = v2;

    const double m01 = a10 * a01 - a00 * a11;
    const double m02 = a20 * a01 - a00 * a21;
    const double m03 = a30 * a01 - a00 * a31;
    const double m12 = a20 * a11 - a10 * a21;
    const double m13 = a30 * a11 - a10 * a31;
    const double m23 = a30 * a21 - a20 * a31;
    const double m012 = m12 * a02 - m02 * a12 + m01 * a22;
    const double m013 = m13 * a02 - m03 * a12 + m01 * a32;
    const double m023 = m23 * a02 - m03 * a22 + m02 * a32;
    const double m123 = m23 * a12 - m13 * a22 + m12 * a32;
    const double det = m123 * a03 - m023 * a13 + m013 * a23 - m012 * a33;

    const double rows[4][4] = { {a00, a01, a02, a03},
                                {a10, a11, a12, a13},
                                {a20, a21, a22, a23},
                                {a30, a31, a32, a33} };
    // Conservative magnitude bound: the product of the rows' absolute sums
    // dominates every term of the expanded determinant, so it plays the role
    // the permanent plays in the 3D predicates.
    double mag = 0.0;
    for (int i = 0; i < 4; ++i) {
        double row = 0.0;
        for (int j = 0; j < 4; ++j) row += std::fabs(rows[i][j]);
        mag = (mag == 0.0) ? row : mag * row;
    }
    const double errbound = kInCircleBound * mag;
    if (det > errbound || -det > errbound) return sign_of(det);

    return in_circle_exact(p0, p1, p2, p);
}

Sign coplanar_in_circle(const double* p0, const double* p1,
                        const double* p2, const double* p) {
    const Sign bs = coplanar_in_circle_unperturbed(p0, p1, p2, p);
    if (bs != ZERO) return bs;

    // Exactly cocircular: the 2D form of the perturbation in predicates.cpp.
    // Walk the four points from the lexicographically largest downwards; the
    // first one whose replacement by p gives a decisive orientation settles it.
    // `local` re-signs that orientation against the reference triple, which may
    // itself be negatively oriented.
    const double* pts[4] = { p0, p1, p2, p };
    std::sort(pts, pts + 4, lexico_less);

    const Sign local = coplanar_orientation(p0, p1, p2);

    for (int i = 3; i > 0; --i) {
        if (pts[i] == p) return NEGATIVE;   // perturbing the query point lifts it out
        Sign o;
        if (pts[i] == p2 && (o = coplanar_orientation(p0, p1, p)) != ZERO) return mul(o, local);
        if (pts[i] == p1 && (o = coplanar_orientation(p0, p, p2)) != ZERO) return mul(o, local);
        if (pts[i] == p0 && (o = coplanar_orientation(p, p1, p2)) != ZERO) return mul(o, local);
    }
    return Sign(-int(local));
}

} // namespace del3d
