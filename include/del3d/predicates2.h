// del3d - the coplanar (2D-embedded-in-3D) predicates.
//
// The incremental build needs these in two situations: while the triangulation
// is still degenerate - every point inserted so far collinear or coplanar, so
// there is no tetrahedron to test - and whenever a new point falls exactly on a
// facet of the convex hull, where the 3D in-sphere test degenerates and the
// decision has to be made inside the plane of that facet.
//
// Like the 3D pair in predicates.h these are filtered-exact, and
// coplanar_in_circle() resolves its degenerate (cocircular) case by symbolic
// perturbation so that it never returns ZERO.
//
// A note on projections: the 2D predicates work on an axis-aligned projection
// of the plane. Which projection is chosen is not observable in the final
// triangulation as long as one fixed rule is applied everywhere, so the rule is
// stated explicitly on each function below and must not be varied.
#ifndef DEL3D_PREDICATES2_H
#define DEL3D_PREDICATES2_H

#include "del3d/predicates.h"

namespace del3d {

/// Sign of the 2x2 determinant |q-p, r-p|: the orientation of the planar
/// triple (p,q,r), ZERO iff they are collinear. Exact.
Sign orient2d(double px, double py, double qx, double qy, double rx, double ry);

/**
 * Orientation of three coplanar 3D points, evaluated in the first coordinate
 * projection - xy, then yz, then xz - in which they are not collinear.
 * \pre p, q, r are not collinear (otherwise every projection is degenerate).
 */
Sign coplanar_orientation(const double* p, const double* q, const double* r);

/// Orientation of s with respect to the directed line (p,q), signed
/// consistently with the reference triple (p,q,r): the projection is selected
/// by (p,q,r) and s is then measured in that same projection, so the two
/// answers can be compared and multiplied.
Sign coplanar_orientation(const double* p, const double* q,
                          const double* r, const double* s);

/// True iff the three points are exactly collinear. Defined for every input,
/// unlike coplanar_orientation(), which has non-collinearity as a precondition.
bool collinear(const double* p, const double* q, const double* r);

/// True iff the two points have bitwise-equal coordinates.
bool equal(const double* p, const double* q);

/// Where q lies on the line through p and r, assuming the three are collinear:
/// before p, at p, strictly between, at r, or after r.
enum CollinearPosition { CP_BEFORE, CP_SOURCE, CP_MIDDLE, CP_TARGET, CP_AFTER };
CollinearPosition collinear_position(const double* p, const double* q, const double* r);

/**
 * POSITIVE iff p lies strictly inside the circle through the coplanar points
 * p0, p1, p2.
 * \pre coplanar_orientation(p0,p1,p2) != ZERO
 *
 * Never returns ZERO: the exactly-cocircular case is decided by the 2D
 * counterpart of the perturbation in predicates.h - sort the four points
 * lexicographically, then walk down from the largest, evaluating the
 * orientation of the triple with that point replaced by p.
 */
Sign coplanar_in_circle(const double* p0, const double* p1,
                        const double* p2, const double* p);

/// The unperturbed version: ZERO when the four points are exactly cocircular.
Sign coplanar_in_circle_unperturbed(const double* p0, const double* p1,
                                    const double* p2, const double* p);

} // namespace del3d

#endif // DEL3D_PREDICATES2_H
