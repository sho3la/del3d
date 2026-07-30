// del3d - 3D Delaunay triangulation. Public interface.
#ifndef DEL3D_DELAUNAY_H
#define DEL3D_DELAUNAY_H

#include "del3d/spatial_sort.h"   // Point3

#include <array>
#include <cstddef>
#include <vector>

namespace del3d {

/**
 * The circumcentre of a tetrahedron - the centre of the unique sphere through
 * its four vertices, i.e. the vertex of the Voronoi diagram dual to the cell.
 *
 * Numerator and denominator are both formed in exact expansion arithmetic and
 * only the final quotient is rounded, so the result is the correctly rounded
 * value of the exact expression. This is not interchangeable with evaluating
 * the same formula in double or long double: on a sliver (a nearly flat
 * tetrahedron) the denominator 2 a.(b x c) cancels catastrophically, and since
 * a sliver's circumcentre is legitimately far away, an inexact evaluation can
 * be wrong by many orders of magnitude. Those distant duals are exactly the
 * ones that dominate a nearest-Voronoi-vertex search, so the difference is not
 * cosmetic.
 *
 * Returns p0 if the four points are coplanar (no circumcentre exists).
 */
Point3 circumcenter(const Point3& p0, const Point3& p1,
                    const Point3& p2, const Point3& p3);

/**
 * The circumcentre as produced by a lazy exact-arithmetic evaluation of the
 * same formula, which is deliberately *not* the exact value.
 *
 * The lazy scheme evaluates the expression in interval arithmetic first and
 * returns the midpoint of the resulting interval as soon as that interval is
 * accurate to 1e-5 relative, falling back on exact arithmetic only when it is
 * not. On a sliver the returned value can therefore differ from the exact
 * circumcentre by ~1e-9 relative. This function reproduces that behaviour
 * exactly: it evaluates the formula in the interval arithmetic of
 * src/interval.h, applies the same accept-or-refine decision per coordinate,
 * and calls circumcenter() above for the coordinates that need refining.
 *
 * Use this when the goal is bit-compatibility with such a pipeline; use
 * circumcenter() when the goal is the most accurate answer.
 */
Point3 circumcenter_lazy(const Point3& p, const Point3& q,
                         const Point3& r, const Point3& s);

/**
 * The Delaunay triangulation of a point set.
 *
 * Construction is incremental: points are inserted in the spatial-sort order of
 * spatial_sort.h, each insertion locating the point by a walk, collecting the
 * cells whose circumsphere contains it (the conflict region) and re-starring
 * that region from the new vertex. All predicates are exact and all degenerate
 * cases are resolved by symbolic perturbation, so the output is a deterministic
 * function of the input point sequence - the same tetrahedra *and the same cell
 * order* on every platform, including on cospherical and cocircular input.
 */
class Delaunay {
public:
    /// A finite tetrahedron, as indices into the point set passed to build().
    using Tet = std::array<int, 4>;

    /// Builds the triangulation. Points equal to an earlier point are ignored.
    void build(const std::vector<Point3>& points);

    /// The finite cells - the tetrahedra - in cell-allocation order.
    const std::vector<Tet>& finite_cells() const { return cells_; }

    /// Every cell in allocation order, infinite ones included, with the
    /// infinite vertex reported as -1. The infinite cells are the ones sharing
    /// a facet with the convex hull; finite_cells() is this list filtered.
    /// Exposed mainly for tests, which use it to inspect the slot layout that
    /// the finite-cell order is a view of.
    const std::vector<Tet>& all_cells() const { return all_cells_; }

    /// For each input point, the indices (into finite_cells()) of the cells
    /// incident to it. Empty for points that duplicated an earlier one.
    const std::vector<std::vector<int>>& incident_cells() const { return incident_; }

    /// Number of distinct vertices actually inserted, duplicates excluded.
    std::size_t number_of_vertices() const { return n_vertices_; }

private:
    std::vector<Tet>              cells_;
    std::vector<Tet>              all_cells_;
    std::vector<std::vector<int>> incident_;
    std::size_t                   n_vertices_ = 0;
};

} // namespace del3d

#endif // DEL3D_DELAUNAY_H
