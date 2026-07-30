// del3d - exact geometric predicates in 3D.
//
// Both predicates are *filtered*: they first evaluate the determinant in plain
// double arithmetic together with a bound on the accumulated round-off error,
// and only when the computed value is smaller than that bound - i.e. when the
// sign is not certain - do they re-evaluate the same determinant in exact
// expansion arithmetic (src/expansion.h). The exact stage uses nothing but IEEE
// double operations: no extended precision types, no arbitrary-precision
// integers.
//
// in_sphere() additionally decides the degenerate (exactly cospherical) case by
// symbolic perturbation instead of reporting ZERO, so the predicate is a total
// order on configurations. This is what makes the resulting triangulation
// well defined - and reproducible - on cospherical input rather than merely
// self-consistent. See in_sphere() for the rule.
#ifndef DEL3D_PREDICATES_H
#define DEL3D_PREDICATES_H

namespace del3d {

enum Sign { NEGATIVE = -1, ZERO = 0, POSITIVE = 1 };

/// Sign of the determinant |b-a, c-a, d-a|: POSITIVE iff d lies on the positive
/// side of the plane through (a,b,c) oriented by that triple, ZERO iff the four
/// points are exactly coplanar. Exact for any finite double input.
Sign orient3d(const double* a, const double* b, const double* c, const double* d);

/**
 * POSITIVE iff e lies strictly inside the sphere circumscribing a,b,c,d.
 * \pre orient3d(a,b,c,d) == POSITIVE
 *
 * Never returns ZERO. When the five points are exactly cospherical the sign is
 * decided by a symbolic perturbation, i.e. by asking what the sign would be if
 * each point in turn were displaced by an infinitesimal amount, largest point
 * first in lexicographic order:
 *
 *     sort {a,b,c,d,e} lexicographically by (x,y,z)
 *     for i = 4 down to 3:                       // two iterations suffice
 *         if points[i] is e            -> NEGATIVE
 *         if points[i] is d and orient3d(a,b,c,e) != 0 -> that sign
 *         if points[i] is c and orient3d(a,b,e,d) != 0 -> that sign
 *         if points[i] is b and orient3d(a,e,c,d) != 0 -> that sign
 *         if points[i] is a and orient3d(e,b,c,d) != 0 -> that sign
 *
 * A Delaunay triangulation is uniquely determined by its in-sphere predicate,
 * so once this function is fixed the set of tetrahedra is fixed too, whatever
 * insertion order or cavity algorithm produces it.
 */
Sign in_sphere(const double* a, const double* b, const double* c,
               const double* d, const double* e);

/// The unperturbed predicate: returns ZERO when the five points are exactly
/// cospherical. Exposed for testing; the triangulation always uses in_sphere().
Sign in_sphere_unperturbed(const double* a, const double* b, const double* c,
                           const double* d, const double* e);

} // namespace del3d

#endif // DEL3D_PREDICATES_H
