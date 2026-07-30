// del3d - exact 3D geometric predicates. See include/del3d/predicates.h.
//
// Each predicate is evaluated twice at most. The first pass computes the
// determinant in plain double arithmetic and, alongside it, the "permanent" -
// the same expression with every term replaced by its absolute value - which
// bounds how much round-off the evaluation can have accumulated. If the
// determinant exceeds that bound in magnitude its sign is certain and the
// function returns. Otherwise the same determinant is recomputed in exact
// expansion arithmetic (src/expansion.h), which is slower by an order of
// magnitude but is reached only for near-degenerate configurations.
#include "del3d/predicates.h"
#include "expansion.h"

#include <algorithm>
#include <cmath>

namespace del3d {
namespace {

using detail::Expansion;

// 2^-53: the round-off unit of IEEE double.
constexpr double kEps = 1.1102230246251565e-16;
// First-stage error bounds, as multiples of the permanent. The constants come
// from bounding the round-off of each product and sum in the expressions below.
constexpr double kOrient3dBound  = (7.0 + 56.0 * kEps) * kEps;
constexpr double kInSphereBound  = (16.0 + 224.0 * kEps) * kEps;

// Sign convention. The determinants evaluated here are the "difference from the
// last point" forms: orient3d computes |a-d, b-d, c-d| and in_sphere works
// relative to e. The interface, however, is specified in terms of |b-a, c-a,
// d-a|, which is an odd permutation of the former and therefore its negation;
// the in-sphere sign inherits the same flip, because the side of the sphere
// reverses with the orientation of (a,b,c,d). The public functions apply that
// negation at the end, so the sign returned is the one the perturbation rule
// and the triangulation assume.
inline Sign flip(Sign s) { return Sign(-int(s)); }

inline Sign sign_of(double d) { return d > 0.0 ? POSITIVE : (d < 0.0 ? NEGATIVE : ZERO); }

// ---- orient3d ------------------------------------------------------------

/// The 3x3 determinant |a-d, b-d, c-d| in exact expansion arithmetic. Every
/// coordinate difference is exact (from_diff of two doubles), and expansion
/// products and sums are exact, so the sign of the result is the true sign.
Sign orient3d_exact(const double* a, const double* b, const double* c, const double* d) {
    const Expansion adx = Expansion::from_diff(a[0], d[0]);
    const Expansion ady = Expansion::from_diff(a[1], d[1]);
    const Expansion adz = Expansion::from_diff(a[2], d[2]);
    const Expansion bdx = Expansion::from_diff(b[0], d[0]);
    const Expansion bdy = Expansion::from_diff(b[1], d[1]);
    const Expansion bdz = Expansion::from_diff(b[2], d[2]);
    const Expansion cdx = Expansion::from_diff(c[0], d[0]);
    const Expansion cdy = Expansion::from_diff(c[1], d[1]);
    const Expansion cdz = Expansion::from_diff(c[2], d[2]);

    // det = adz*(bdx*cdy - cdx*bdy) + bdz*(cdx*ady - adx*cdy) + cdz*(adx*bdy - bdx*ady)
    const Expansion t1 = bdx * cdy - cdx * bdy;
    const Expansion t2 = cdx * ady - adx * cdy;
    const Expansion t3 = adx * bdy - bdx * ady;
    const Expansion det = adz * t1 + bdz * t2 + cdz * t3;
    return det.sign();
}

// ---- in_sphere -----------------------------------------------------------

/**
 * The in-sphere determinant in exact expansion arithmetic.
 *
 * Translating the configuration so that e is the origin reduces the 5x5
 * determinant to a 4x4 one in the lifted coordinates: each remaining point
 * contributes its offset from e plus the squared length of that offset (its
 * height on the paraboloid). Expanding along the lifted column gives the
 * 2x2 minors ab..bd and then the four 3x3 cofactors abc, bcd, cda, dab.
 */
Sign in_sphere_exact(const double* a, const double* b, const double* c,
                     const double* d, const double* e) {
    const Expansion aex = Expansion::from_diff(a[0], e[0]);
    const Expansion aey = Expansion::from_diff(a[1], e[1]);
    const Expansion aez = Expansion::from_diff(a[2], e[2]);
    const Expansion bex = Expansion::from_diff(b[0], e[0]);
    const Expansion bey = Expansion::from_diff(b[1], e[1]);
    const Expansion bez = Expansion::from_diff(b[2], e[2]);
    const Expansion cex = Expansion::from_diff(c[0], e[0]);
    const Expansion cey = Expansion::from_diff(c[1], e[1]);
    const Expansion cez = Expansion::from_diff(c[2], e[2]);
    const Expansion dex = Expansion::from_diff(d[0], e[0]);
    const Expansion dey = Expansion::from_diff(d[1], e[1]);
    const Expansion dez = Expansion::from_diff(d[2], e[2]);

    // Paraboloid lifts: the squared distance of each point from e.
    const Expansion alift = aex * aex + aey * aey + aez * aez;
    const Expansion blift = bex * bex + bey * bey + bez * bez;
    const Expansion clift = cex * cex + cey * cey + cez * cez;
    const Expansion dlift = dex * dex + dey * dey + dez * dez;

    // 2x2 minors in the xy plane.
    const Expansion ab = aex * bey - bex * aey;
    const Expansion bc = bex * cey - cex * bey;
    const Expansion cd = cex * dey - dex * cey;
    const Expansion da = dex * aey - aex * dey;
    const Expansion ac = aex * cey - cex * aey;
    const Expansion bd = bex * dey - dex * bey;

    // 3x3 cofactors, each omitting one of the four points.
    const Expansion abc = aez * bc - bez * ac + cez * ab;
    const Expansion bcd = bez * cd - cez * bd + dez * bc;
    const Expansion cda = cez * da + dez * ac + aez * cd;
    const Expansion dab = dez * ab + aez * bd + bez * da;

    const Expansion det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd);
    return det.sign();
}

/// Lexicographic order on (x, y, z). Used to make the perturbation below a
/// function of the point set alone, independent of the argument order.
inline bool lexico_less(const double* p, const double* q) {
    if (p[0] != q[0]) return p[0] < q[0];
    if (p[1] != q[1]) return p[1] < q[1];
    return p[2] < q[2];
}

} // namespace

Sign orient3d(const double* a, const double* b, const double* c, const double* d) {
    const double adx = a[0] - d[0], bdx = b[0] - d[0], cdx = c[0] - d[0];
    const double ady = a[1] - d[1], bdy = b[1] - d[1], cdy = c[1] - d[1];
    const double adz = a[2] - d[2], bdz = b[2] - d[2], cdz = c[2] - d[2];

    const double bdxcdy = bdx * cdy, cdxbdy = cdx * bdy;
    const double cdxady = cdx * ady, adxcdy = adx * cdy;
    const double adxbdy = adx * bdy, bdxady = bdx * ady;

    const double det = adz * (bdxcdy - cdxbdy)
                     + bdz * (cdxady - adxcdy)
                     + cdz * (adxbdy - bdxady);

    // The same expression with absolute values: an upper bound on the size of
    // the terms, hence on the round-off that could have been accumulated.
    const double permanent = (std::fabs(bdxcdy) + std::fabs(cdxbdy)) * std::fabs(adz)
                           + (std::fabs(cdxady) + std::fabs(adxcdy)) * std::fabs(bdz)
                           + (std::fabs(adxbdy) + std::fabs(bdxady)) * std::fabs(cdz);
    const double errbound = kOrient3dBound * permanent;
    if (det > errbound || -det > errbound) return flip(sign_of(det));

    return flip(orient3d_exact(a, b, c, d));
}

Sign in_sphere_unperturbed(const double* a, const double* b, const double* c,
                           const double* d, const double* e) {
    // The filtered counterpart of in_sphere_exact: identical expression, plain
    // doubles, plus the permanent that bounds its error.
    const double aex = a[0] - e[0], aey = a[1] - e[1], aez = a[2] - e[2];
    const double bex = b[0] - e[0], bey = b[1] - e[1], bez = b[2] - e[2];
    const double cex = c[0] - e[0], cey = c[1] - e[1], cez = c[2] - e[2];
    const double dex = d[0] - e[0], dey = d[1] - e[1], dez = d[2] - e[2];

    const double ab = aex * bey - bex * aey;
    const double bc = bex * cey - cex * bey;
    const double cd = cex * dey - dex * cey;
    const double da = dex * aey - aex * dey;
    const double ac = aex * cey - cex * aey;
    const double bd = bex * dey - dex * bey;

    const double abc = aez * bc - bez * ac + cez * ab;
    const double bcd = bez * cd - cez * bd + dez * bc;
    const double cda = cez * da + dez * ac + aez * cd;
    const double dab = dez * ab + aez * bd + bez * da;

    const double alift = aex * aex + aey * aey + aez * aez;
    const double blift = bex * bex + bey * bey + bez * bez;
    const double clift = cex * cex + cey * cey + cez * cez;
    const double dlift = dex * dex + dey * dey + dez * dez;

    const double det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd);

    const double aezplus = std::fabs(aez), bezplus = std::fabs(bez);
    const double cezplus = std::fabs(cez), dezplus = std::fabs(dez);
    const double aexbeyplus = std::fabs(aex * bey), bexaeyplus = std::fabs(bex * aey);
    const double bexceyplus = std::fabs(bex * cey), cexbeyplus = std::fabs(cex * bey);
    const double cexdeyplus = std::fabs(cex * dey), dexceyplus = std::fabs(dex * cey);
    const double dexaeyplus = std::fabs(dex * aey), aexdeyplus = std::fabs(aex * dey);
    const double aexceyplus = std::fabs(aex * cey), cexaeyplus = std::fabs(cex * aey);
    const double bexdeyplus = std::fabs(bex * dey), dexbeyplus = std::fabs(dex * bey);

    const double permanent =
        ((cexdeyplus + dexceyplus) * bezplus +
         (dexbeyplus + bexdeyplus) * cezplus +
         (bexceyplus + cexbeyplus) * dezplus) * alift +
        ((dexaeyplus + aexdeyplus) * cezplus +
         (aexceyplus + cexaeyplus) * dezplus +
         (cexdeyplus + dexceyplus) * aezplus) * blift +
        ((aexbeyplus + bexaeyplus) * dezplus +
         (bexdeyplus + dexbeyplus) * aezplus +
         (dexaeyplus + aexdeyplus) * bezplus) * clift +
        ((bexceyplus + cexbeyplus) * aezplus +
         (cexaeyplus + aexceyplus) * bezplus +
         (aexbeyplus + bexaeyplus) * cezplus) * dlift;

    const double errbound = kInSphereBound * permanent;
    if (det > errbound || -det > errbound) return flip(sign_of(det));

    return flip(in_sphere_exact(a, b, c, d, e));
}

Sign in_sphere(const double* a, const double* b, const double* c,
               const double* d, const double* e) {
    const Sign s = in_sphere_unperturbed(a, b, c, d, e);
    if (s != ZERO) return s;

    // Exactly cospherical: decide by symbolic perturbation (see the header).
    // Sorting by (x,y,z) makes the outcome depend only on the point set, and
    // the walk from the lexicographically largest point downwards implements
    // "displace the largest point first, by the largest infinitesimal".
    const double* p[5] = { a, b, c, d, e };
    std::sort(p, p + 5, lexico_less);

    for (int i = 4; i > 2; --i) {
        // Perturbing the query point itself lifts it off the sphere outwards.
        if (p[i] == e) return NEGATIVE;
        Sign o;
        // Otherwise: replace the perturbed point by e in the oriented
        // tetrahedron and take that orientation, if it is decisive.
        if (p[i] == d && (o = orient3d(a, b, c, e)) != ZERO) return o;
        if (p[i] == c && (o = orient3d(a, b, e, d)) != ZERO) return o;
        if (p[i] == b && (o = orient3d(a, e, c, d)) != ZERO) return o;
        if (p[i] == a && (o = orient3d(e, b, c, d)) != ZERO) return o;
    }
    // Two iterations always decide when a,b,c,d are affinely independent, which
    // the precondition guarantees; this is unreachable.
    return NEGATIVE;
}

} // namespace del3d
