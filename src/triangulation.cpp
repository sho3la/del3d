// del3d - incremental Delaunay construction and the dual (circumcentre).
//
// The build inserts points one at a time in spatial-sort order. Each insertion:
//
//   1. locates the point - a walk from a hint cell towards the cell containing
//      it, reporting whether it landed in a cell, on a facet, on an edge, on an
//      existing vertex, outside the convex hull, or outside the affine hull;
//   2. collects the conflict region - the connected set of cells whose
//      circumsphere contains the point, found by a flood fill from the located
//      cell (for an infinite cell, "circumsphere contains" degenerates to a
//      half-space test);
//   3. re-stars that region from the new vertex (tds.h), which is what restores
//      the Delaunay property locally, and deletes the old cells.
//
// Points that arrive before the configuration is full-dimensional are handled
// by a separate ladder: dimension -2 -> -1 -> 0 -> 1 -> 2 -> 3, raising the
// dimension whenever a point falls outside the current affine hull.
//
// The order in which cells are created and destroyed is observable in the
// output (see compact_container.h), so the traversals here are fixed by
// specification: where an arbitrary-but-consistent choice exists, the choice is
// written out rather than left to chance.
#include "del3d/delaunay.h"
#include "del3d/predicates.h"
#include "del3d/predicates2.h"
#include "del3d/spatial_sort.h"

#include "expansion.h"
#include "interval.h"
#include "rand48.h"
#include "tds.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace del3d {
namespace {

using detail::Cell;
using detail::Index;
using detail::kNull;
using detail::Mark;
using detail::Tds;
using detail::Utils;
using detail::Vertex;

/// What a point location returned: the dimension of the face the point landed
/// on, or that it lies beyond the convex hull / outside the current affine hull.
enum LocateType { LT_VERTEX = 0, LT_EDGE, LT_FACET, LT_CELL,
                  LT_OUTSIDE_CONVEX_HULL, LT_OUTSIDE_AFFINE_HULL };

/// Owns the data structure while the triangulation is being built.
class Builder {
public:
    explicit Builder(const std::vector<Point3>& pts) : pts_(pts) {
        // Creating the infinite vertex takes the dimension from -2 to -1.
        tds_.infinite = insert_increase_dimension(kNull);
    }

    /// Inserts every point in spatial-sort order, carrying the previous
    /// vertex forward as the hint for the next point location.
    void build() {
        std::vector<int> order(pts_.size());
        for (std::size_t i = 0; i < order.size(); ++i) order[i] = int(i);
        spatial_sort(pts_, order);

        Index hint = kNull;
        int step = 0;
        for (int idx : order) {
#ifdef DEL3D_TRACE
            std::fprintf(stdout, "BEGIN step %d (point %d)\n", step + 1, idx);
#endif
            hint = insert(idx, hint);
            // The point index is recorded after every insertion, including when
            // the point duplicated an earlier one and insert() returned the
            // existing vertex - so a repeated point ends up carrying the LAST
            // input index that mapped to it, not the first.
            if (hint != kNull) tds_.vertices[hint].info = idx;
            ++step;
#ifdef DEL3D_TRACE
            std::fprintf(stdout, "step %d (point %d) dim %d:\n", step, idx, tds_.dimension);
            std::vector<Index> slot_of;
            tds_.cells.for_each([&](Index c) { slot_of.push_back(c); });
            auto slot_index = [&](Index c) {
                for (std::size_t k = 0; k < slot_of.size(); ++k)
                    if (slot_of[k] == c) return int(k);
                return -9;
            };
            int pos = 0;
            tds_.cells.for_each([&](Index c) {
                std::fprintf(stdout, "   slot %3d :", pos++);
                for (int i = 0; i < 4; ++i) {
                    const Index vv = tds_.cells[c].v[i];
                    std::fprintf(stdout, " %3d",
                                 (vv == tds_.infinite || vv == kNull)
                                     ? -1 : tds_.vertices[vv].info);
                }
                std::fprintf(stdout, "   n=[");
                for (int i = 0; i < 4; ++i)
                    std::fprintf(stdout, "%d ", slot_index(tds_.cells[c].n[i]));
                std::fprintf(stdout, "]\n");
            });
            if (hint != kNull)
                std::fprintf(stdout, "   >>> hint vertex %d ->cell() = position %d\n",
                             tds_.vertices[hint].info, slot_index(tds_.vertices[hint].cell));
#endif
#ifdef DEL3D_VALIDATE
            const std::string err = tds_.validate();
            if (!err.empty()) {
                std::fprintf(stderr,
                    "del3d: TDS invariant broken after insert #%d (point %d, dim %d): %s\n",
                    step, idx, tds_.dimension, err.c_str());
                std::abort();
            }
#endif
        }
    }

    Tds& tds() { return tds_; }

private:
    const std::vector<Point3>& pts_;
    Tds tds_;
    // Scratch buffers reused across insertions. Rebuilding them per insertion
    // cost ~10^5 short-lived allocations on a 35k-point input; reusing them
    // changes no result, only where the memory comes from.
    std::vector<Index> scratch_conflict_;
    std::vector<std::pair<Index, int>> scratch_facets_;
    std::vector<Index> scratch_stack_;

    /// The coordinates of a vertex, as three consecutive doubles. Point3 is a
    /// plain struct of three doubles, so its address is a valid double[3].
    const double* pt(Index v) const {
        return reinterpret_cast<const double*>(&pts_[tds_.vertices[v].info]);
    }
    const double* pt_of(int idx) const {
        return reinterpret_cast<const double*>(&pts_[idx]);
    }

    // ---- raising the dimension ------------------------------------------
    /**
     * Adds a vertex that lies outside the current affine hull, rebuilding the
     * structure one dimension higher. `star` is the vertex the new cells are
     * coned from (the infinite vertex, except for the very first call).
     *
     * One case per starting dimension:
     *   -2 -> -1  create the infinite vertex and its single dummy cell;
     *   -1 ->  0  the first finite point: two cells facing each other;
     *    0 ->  1  a segment: three cells forming the 1D chain plus its two
     *             infinite ends;
     *    1 ->  2  cone every edge of the 1D chain to the new point, walking the
     *             chain from `star` all the way round;
     *    2 ->  3  give every existing triangle the new point as its fourth
     *             vertex, and mirror the ones not incident to `star` into new
     *             cells on the other side, then link the mirrored layer up.
     * Returns the new vertex.
     */
    Index insert_increase_dimension(Index star) {
        Index v = tds_.create_vertex();
        const int dim = tds_.dimension;
        tds_.dimension = dim + 1;

        switch (dim) {
        case -2: {                       // the very first vertex (the infinite one)
            Index c = tds_.create_cell(v, kNull, kNull, kNull);
            tds_.vertices[v].cell = c;
            break;
        }
        case -1: {                       // second vertex, i.e. the first finite one
            Index d = tds_.create_cell(v, kNull, kNull, kNull);
            tds_.vertices[v].cell = d;
            tds_.set_adjacency(d, 0, tds_.vertices[star].cell, 0);
            break;
        }
        case 0: {
            Index c = tds_.vertices[star].cell;
            Index d = tds_.cells[c].n[0];
            tds_.cells[c].v[1] = tds_.cells[d].v[0];
            tds_.cells[d].v[1] = v;
            tds_.cells[d].n[1] = c;
            Index e = tds_.create_cell(v, star, kNull, kNull);
            tds_.set_adjacency(e, 0, c, 1);
            tds_.set_adjacency(e, 1, d, 0);
            tds_.vertices[v].cell = d;
            break;
        }
        case 1: {
            Index c = tds_.vertices[star].cell;
            const int i = tds_.cells[c].index_of_vertex(star);   // 0 or 1
            const int j = (i == 0) ? 1 : 0;
            Index d = tds_.cells[c].n[j];

            tds_.cells[c].v[2] = v;

            // Walk the chain of 1D cells from c back round to d, giving each
            // one v as its third vertex and creating its mirror image coned
            // from `star` on the other side.
            Index e = tds_.cells[c].n[i];
            Index cnew = c;
            Index enew = kNull;

            while (e != d) {
                enew = tds_.create_cell();
                tds_.cells[enew].v[i] = tds_.cells[e].v[j];
                tds_.cells[enew].v[j] = tds_.cells[e].v[i];
                tds_.cells[enew].v[2] = star;
                tds_.set_adjacency(enew, i, cnew, j);
                tds_.set_adjacency(enew, 2, e, 2);
                tds_.cells[e].v[2] = v;
                e = tds_.cells[e].n[i];
                cnew = enew;
            }
            tds_.cells[d].v[2] = v;
            tds_.set_adjacency(enew, j, d, 2);

            c = tds_.vertices[star].cell;
            tds_.cells[c].n[2] = tds_.cells[tds_.cells[c].n[i]].n[2];
            tds_.cells[c].n[j] = d;
            tds_.vertices[v].cell = d;
            break;
        }
        case 2: {
            std::vector<Index> new_cells;
            new_cells.reserve(16);

            bool first = true;
            // Snapshot first: the loop creates cells, and iterating the live
            // container while doing so would visit them too.
            std::vector<Index> snapshot;
            tds_.cells.for_each([&](Index c) { snapshot.push_back(c); });

            for (Index it : snapshot) {
                if (first) { tds_.vertices[v].cell = it; first = false; }
                // Cells made inside this loop are recognised by their still
                // unset neighbour 0 and skipped.
                if (tds_.cells[it].n[0] == kNull) continue;
                tds_.cells[it].n[3] = kNull;
                tds_.cells[it].v[3] = v;
                if (!tds_.cells[it].has_vertex(star)) {
                    // Mirror image on the other side of the old plane; vertices
                    // 1 and 2 are swapped to keep the orientation positive.
                    Index cnew = tds_.create_cell(tds_.cells[it].v[0], tds_.cells[it].v[2],
                                                  tds_.cells[it].v[1], star);
                    tds_.set_adjacency(cnew, 3, it, 3);
                    tds_.cells[cnew].n[0] = kNull;
                    new_cells.push_back(cnew);
                }
            }

            // Link the mirrored layer: each new cell's neighbours are the
            // mirrors of the corresponding neighbours of the cell it mirrors.
            for (Index nc : new_cells) {
                Index n = tds_.cells[nc].n[3];        // the cell opposite star
                for (int i = 0; i < 3; ++i) {
                    const int j = (i == 0) ? 0 : 3 - i;   // undo the 1/2 swap
                    Index c = tds_.cells[tds_.cells[n].n[i]].n[3];
                    if (c != kNull) tds_.cells[nc].n[j] = c;
                    else tds_.set_adjacency(nc, j, tds_.cells[n].n[i], 3);
                }
            }
            break;
        }
        default: break;
        }
        return v;
    }

    /// Flips the orientation of the whole triangulation by swapping vertices 0
    /// and 1 (and their neighbours) in every cell. Used when a dimension raise
    /// would otherwise leave all cells negatively oriented.
    void reorient() {
        tds_.cells.for_each([&](Index c) {
            std::swap(tds_.cells[c].v[0], tds_.cells[c].v[1]);
            std::swap(tds_.cells[c].n[0], tds_.cells[c].n[1]);
        });
    }

    Index infinite_cell() const {
        return tds_.vertices[tds_.infinite].cell;
    }

    // ---- point location --------------------------------------------------
    /**
     * The orientation test in plain double precision: the sign of
     * det(q-p, r-p, s-p), with no filtering and no exact fallback.
     *
     * Used only by the pre-walk below, where a wrong sign near zero costs an
     * extra step but not correctness. The grouping of the products is fixed:
     * because this is unfiltered floating point, a different association
     * rounds differently, flips the sign of a near-zero determinant and sends
     * the walk into a different cell - which is observable downstream.
     */
    static Sign inexact_orientation(const double* p, const double* q,
                                    const double* r, const double* s) {
        const double a00 = q[0] - p[0], a01 = q[1] - p[1], a02 = q[2] - p[2];
        const double a10 = r[0] - p[0], a11 = r[1] - p[1], a12 = r[2] - p[2];
        const double a20 = s[0] - p[0], a21 = s[1] - p[1], a22 = s[2] - p[2];
        const double m01 = a00 * a11 - a10 * a01;
        const double m02 = a00 * a21 - a20 * a01;
        const double m12 = a10 * a21 - a20 * a11;
        const double det = m01 * a22 - m02 * a12 + m12 * a02;
        if (det > 0) return POSITIVE;
        if (det < 0) return NEGATIVE;
        return ZERO;
    }

    /**
     * Cheap pre-walk ("structural filtering"): the same walk as exact_locate()
     * but with inexact orientation tests and a step budget, used to get close
     * to the target cell before the exact walk takes over.
     *
     * It is part of the specification, not an optional speed-up. The cell this
     * returns is the cell the exact walk starts from; a different start means a
     * different walk, a different cell for the conflict search to start from, a
     * different boundary-facet order and ultimately a different cell order.
     *
     * Deliberately non-stochastic - it scans facets 0..3 in order - so unlike
     * the exact walk it draws nothing from the random engine.
     */
    Index inexact_locate(const double* t, Index start, int n_of_turns) const {
        if (tds_.dimension < 3) return start;

        if (start == kNull) start = infinite_cell();
        // Start from a finite cell: orientation tests against the infinite
        // vertex are meaningless.
        if (tds_.cells[start].has_vertex(tds_.infinite))
            start = tds_.cells[start].n[tds_.cells[start].index_of_vertex(tds_.infinite)];

        Index previous = kNull;
        Index c = start;

    try_next_cell:
        --n_of_turns;
        {
            const double* pts[4] = { pt(tds_.cells[c].v[0]), pt(tds_.cells[c].v[1]),
                                     pt(tds_.cells[c].v[2]), pt(tds_.cells[c].v[3]) };
            for (int i = 0; i != 4; ++i) {
                const Index next = tds_.cells[c].n[i];
                if (previous == next) continue;   // never step straight back

                // Replace vertex i by the query point: a negative orientation
                // means the point is on the far side of facet i.
                const double* backup = pts[i];
                pts[i] = t;
                if (inexact_orientation(pts[0], pts[1], pts[2], pts[3]) != NEGATIVE) {
                    pts[i] = backup;
                    continue;
                }
                // Crossing that facet leaves the convex hull; hand over here.
                if (tds_.cells[next].has_vertex(tds_.infinite)) return next;

                previous = c;
                c = next;
                if (n_of_turns) goto try_next_cell;
            }
        }
        return c;
    }

    /// Locates a point: pre-walk in double precision, then the exact walk.
    Index locate(const double* p, LocateType& lt, int& li, int& lj, Index start) {
        const Index ch = inexact_locate(p, start, kMaxVisitedCells);
        return exact_locate(p, lt, li, lj, ch);
    }

    /// Step budget for the pre-walk.
    static constexpr int kMaxVisitedCells = 2500;

    /**
     * The exact point-location walk.
     *
     * At dimension 3 this is a remembering stochastic walk: from the current
     * cell, test the four facets in a random cyclic order (skipping the one
     * just crossed - that is the "remembering" part, which prevents immediate
     * backtracking) and step through the first facet the query point lies
     * beyond. When no facet is crossed the point is inside the current cell,
     * and the number of orientation tests that came out exactly ZERO says which
     * face it landed on: none = interior, one = a facet, two = an edge, three =
     * a vertex. The indices of that face are returned in li and lj.
     *
     * Lower dimensions are delegated to locate_2/1/0.
     */
    Index exact_locate(const double* p, LocateType& lt, int& li, int& lj, Index start) {
        if (tds_.dimension >= 1) {
            if (start == kNull) start = infinite_cell();
            if (tds_.cells[start].has_vertex(tds_.infinite))
                start = tds_.cells[start].n[tds_.cells[start].index_of_vertex(tds_.infinite)];
        }
        // A fresh engine per call, so the walk's random sequence is a function
        // of the walk alone and not of how many points were inserted before.
        detail::Rand48 rng;

        if (tds_.dimension == 2) return locate_2(p, lt, li, lj, start, rng);
        if (tds_.dimension == 1) return locate_1(p, lt, li, lj, start);
        if (tds_.dimension == 0) return locate_0(p, lt, li);
        if (tds_.dimension == -1) { lt = LT_OUTSIDE_AFFINE_HULL; return kNull; }

        assert(tds_.dimension == 3);
        Index previous = kNull;
        Index c = start;
        Sign o[4] = { ZERO, ZERO, ZERO, ZERO };

#ifdef DEL3D_WALK
        std::fprintf(stdout, "   walk start cell verts:");
        for (int q = 0; q < 4; ++q)
            std::fprintf(stdout, " %g", tds_.cells[c].v[q] == tds_.infinite
                                            ? -1.0 : pt(tds_.cells[c].v[q])[0]);
        std::fprintf(stdout, "\n");
#endif
        bool try_next_cell = true;
        while (try_next_cell) {
            try_next_cell = false;
            const double* pts[4] = { pt(tds_.cells[c].v[0]), pt(tds_.cells[c].v[1]),
                                     pt(tds_.cells[c].v[2]), pt(tds_.cells[c].v[3]) };
            // Random starting facet: this is what makes the walk terminate on
            // adversarial configurations instead of cycling.
            int i = die4(rng);
#ifdef DEL3D_WALK
            std::fprintf(stdout, "   walk at cell %d die4=%d\n", c, i);
#endif
            for (int j = 0; !try_next_cell && j != 4; ++j, i = (i + 1) & 3) {
                Index next = tds_.cells[c].n[i];
                if (previous == next) { o[i] = POSITIVE; continue; }

                const double* backup = pts[i];
                pts[i] = p;
                o[i] = orient3d(pts[0], pts[1], pts[2], pts[3]);
                if (o[i] != NEGATIVE) {
                    pts[i] = backup;
                } else {
                    if (tds_.cells[next].has_vertex(tds_.infinite)) {
                        li = tds_.cells[next].index_of_vertex(tds_.infinite);
                        lt = LT_OUTSIDE_CONVEX_HULL;
                        return next;
                    }
                    previous = c;
                    c = next;
                    try_next_cell = true;
                }
            }
        }

        // Degenerate orientations count the dimensions the point is pinned in.
        const int sum = (o[0] == ZERO) + (o[1] == ZERO) + (o[2] == ZERO) + (o[3] == ZERO);
        switch (sum) {
        case 0: lt = LT_CELL; break;
        case 1:
            lt = LT_FACET;
            li = (o[0] == ZERO) ? 0 : (o[1] == ZERO) ? 1 : (o[2] == ZERO) ? 2 : 3;
            break;
        case 2:
            // The edge is the intersection of the two zero facets, named by the
            // two indices whose orientation was NOT zero.
            lt = LT_EDGE;
            li = (o[0] != ZERO) ? 0 : (o[1] != ZERO) ? 1 : 2;
            lj = (o[li + 1] != ZERO) ? li + 1 : (o[li + 2] != ZERO) ? li + 2 : li + 3;
            break;
        default:
            lt = LT_VERTEX;
            li = (o[0] != ZERO) ? 0 : (o[1] != ZERO) ? 1 : (o[2] != ZERO) ? 2 : 3;
            break;
        }
        return c;
    }

    /// Point location at dimension 2: the same walk inside the plane, using the
    /// coplanar orientation of each of the triangle's three edges.
    Index locate_2(const double* p, LocateType& lt, int& li, int& lj,
                   Index start, detail::Rand48& rng) {
        Index c = start;
        // First: is the point in the plane of the triangulation at all? If not,
        // inserting it raises the dimension to 3.
        if (orient3d(pt(tds_.cells[c].v[0]), pt(tds_.cells[c].v[1]),
                     pt(tds_.cells[c].v[2]), p) != ZERO) {
            lt = LT_OUTSIDE_AFFINE_HULL;
            li = 3;                       // dimension 2 has a single facet, 3
            return c;
        }
        for (;;) {
            if (tds_.cells[c].has_vertex(tds_.infinite)) {
                // Reaching an infinite triangle means the point is beyond the
                // hull; li, lj name the hull edge it is beyond.
                const int inf = tds_.cells[c].index_of_vertex(tds_.infinite);
                lt = LT_OUTSIDE_CONVEX_HULL;
                li = Utils::cw(inf);
                lj = Utils::ccw(inf);
                return c;
            }
            // Test the three edges from a random starting one; a negative
            // orientation means the point is across that edge.
            int i = die3(rng);
            const double* p0 = pt(tds_.cells[c].v[i]);
            const double* p1 = pt(tds_.cells[c].v[Utils::ccw(i)]);
            const double* p2 = pt(tds_.cells[c].v[Utils::cw(i)]);
            Sign o[3];
            o[0] = coplanar_orientation(p0, p1, p);
            if (o[0] == NEGATIVE) { c = tds_.cells[c].n[Utils::cw(i)]; continue; }
            o[1] = coplanar_orientation(p1, p2, p);
            if (o[1] == NEGATIVE) { c = tds_.cells[c].n[i]; continue; }
            o[2] = coplanar_orientation(p2, p0, p);
            if (o[2] == NEGATIVE) { c = tds_.cells[c].n[Utils::ccw(i)]; continue; }

            // Inside this triangle; zero orientations pin it to an edge/vertex.
            const int sum = (o[0] == ZERO) + (o[1] == ZERO) + (o[2] == ZERO);
            switch (sum) {
            case 0: lt = LT_FACET; li = 3; break;
            case 1:
                lt = LT_EDGE;
                li = (o[0] == ZERO) ? i : (o[1] == ZERO) ? Utils::ccw(i) : Utils::cw(i);
                lj = Utils::ccw(li);
                break;
            default:
                lt = LT_VERTEX;
                li = (o[0] != ZERO) ? Utils::cw(i) : (o[1] != ZERO) ? i : Utils::ccw(i);
                break;
            }
            return c;
        }
    }

    /// Point location at dimension 1: walk along the chain of collinear
    /// segments until the point is between the current segment's endpoints.
    Index locate_1(const double* p, LocateType& lt, int& li, int& lj, Index start) {
        Index c = start;
        if (!collinear(p, pt(tds_.cells[c].v[0]), pt(tds_.cells[c].v[1]))) {
            lt = LT_OUTSIDE_AFFINE_HULL;   // off the line: the dimension rises
            return c;
        }
        for (;;) {
            if (tds_.cells[c].has_vertex(tds_.infinite)) {
                lt = LT_OUTSIDE_CONVEX_HULL;   // past one end of the chain
                return c;
            }
            switch (collinear_position(pt(tds_.cells[c].v[0]), p, pt(tds_.cells[c].v[1]))) {
            case CP_AFTER:  c = tds_.cells[c].n[0]; continue;
            case CP_BEFORE: c = tds_.cells[c].n[1]; continue;
            case CP_MIDDLE: lt = LT_EDGE;   li = 0; lj = 1; return c;
            case CP_SOURCE: lt = LT_VERTEX; li = 0; return c;
            default:        lt = LT_VERTEX; li = 1; return c;
            }
        }
    }

    /// Point location at dimension 0: there is one finite vertex, so the point
    /// either coincides with it or lies outside the affine hull.
    Index locate_0(const double* p, LocateType& lt, int& li) {
        Index found = kNull;
        tds_.vertices.for_each([&](Index v) {
            if (found == kNull && v != tds_.infinite) found = v;
        });
        if (!equal(p, pt(found))) lt = LT_OUTSIDE_AFFINE_HULL;
        else { lt = LT_VERTEX; li = 0; }
        return tds_.vertices[found].cell;
    }

    /// A value in [0,2], used to pick the walk's starting edge at dimension 2.
    /// The engine's range is far larger than 3, so the bias of the modulo
    /// reduction is negligible; the reduction itself is fixed by specification
    /// because it determines which facet the walk tries first.
    static int die3(detail::Rand48& rng) {
        const std::uint32_t val = rng() - detail::Rand48::min();
        return int(val % 3u);
    }

    /// The same for the four facets of a cell at dimension 3.
    static int die4(detail::Rand48& rng) {
        const std::uint32_t val = rng() - detail::Rand48::min();
        return int(val % 4u);
    }

    // ---- the Delaunay conflict test --------------------------------------
    /**
     * True iff the point lies strictly inside cell c's circumsphere, i.e. iff
     * inserting it must destroy c.
     *
     * For a finite cell this is the in-sphere predicate. For an infinite cell -
     * one incident to the hull - the "circumsphere" is the open half-space
     * beyond the corresponding hull facet, so the test degenerates to an
     * orientation. The table gives, for each position of the infinite vertex,
     * the three finite vertices in the order that makes that orientation the
     * conflict answer directly. If the point is exactly coplanar with the hull
     * facet the half-space test is inconclusive and the decision is made inside
     * the facet's plane, by the perturbed in-circle predicate.
     */
    bool in_conflict(Index c, const double* p) const {
        const Cell& cell = tds_.cells[c];
        for (int i = 0; i < 4; ++i) {
            if (cell.v[i] != tds_.infinite) continue;
            static const int kOrder[4][3] = { {2,1,3}, {2,3,0}, {1,0,3}, {0,1,2} };
            const double* a = pt(cell.v[kOrder[i][0]]);
            const double* b = pt(cell.v[kOrder[i][1]]);
            const double* d = pt(cell.v[kOrder[i][2]]);
            const Sign s = orient3d(a, b, d, p);
            if (s != ZERO) return s == POSITIVE;
            return coplanar_in_circle(a, b, d, p) == POSITIVE;
        }
        return in_sphere(pt(cell.v[0]), pt(cell.v[1]), pt(cell.v[2]), pt(cell.v[3]), p)
               == POSITIVE;
    }

    /**
     * The conflict test at dimension 2: strictly inside the circumcircle of a
     * finite triangle, or - for an infinite triangle - strictly beyond its hull
     * edge. A point exactly on the line of that edge conflicts only if it lies
     * strictly between the edge's two endpoints.
     */
    bool in_conflict_2(Index c, const double* p) const {
        const Cell& cell = tds_.cells[c];
        if (!cell.has_vertex(tds_.infinite)) {
            // In a valid triangulation the triple (0,1,2) is positively oriented.
            return coplanar_in_circle(pt(cell.v[0]), pt(cell.v[1]), pt(cell.v[2]), p)
                   == POSITIVE;
        }
        // v1, v2 are the finite vertices, ordered so that (v1, v2, infinite) is
        // positively oriented.
        const int i3 = cell.index_of_vertex(tds_.infinite);
        const double* v1 = pt(cell.v[Utils::ccw(i3)]);
        const double* v2 = pt(cell.v[Utils::cw(i3)]);
        const Sign o = coplanar_orientation(v1, v2, p);
        if (o != ZERO) return o == POSITIVE;
        return collinear_position(v1, p, v2) == CP_MIDDLE;
    }

    // ---- the conflict region ---------------------------------------------
    /**
     * Flood fill from cell d over all cells in conflict with p, at dimension 3
     * or 2.
     *
     * The region is connected (it is star-shaped about p), so a depth-first
     * walk over neighbours finds all of it. Every visited cell is marked, so
     * each is tested exactly once: InConflict cells go into `cells_out` and
     * have their own neighbours explored; cells that fail the test are marked
     * OnBoundary and the facet leading to them is appended to `facets_out` as
     * the (cell, local facet index) pair *on the inside*.
     *
     * Both outputs matter beyond their contents: the star rebuild either cones
     * from the last facet of the list or walks the whole list in order, so the
     * order the facets are discovered in is part of the result.
     */
    void find_conflicts(Index d, const double* p,
                        std::vector<Index>& cells_out,
                        std::vector<std::pair<Index, int>>& facets_out) {
#ifdef DEL3D_VALIDATE
        {   // Every mark must be Clear on entry: this function is the only one
            // that sets them and the star rebuild is what clears them again.
            int dirty = 0;
            tds_.cells.for_each([&](Index c) {
                if (tds_.cells[c].mark != Mark::Clear) {
                    if (dirty < 6)
                        std::fprintf(stderr,
                            "del3d: cell %d already marked %d on entry to find_conflicts\n",
                            c, int(tds_.cells[c].mark));
                    ++dirty;
                }
            });
            if (dirty) std::fprintf(stderr, "del3d: %d cells dirty on entry\n", dirty);
        }
#endif
        // dimension+1 neighbours per cell, so this serves dimensions 3 and 2.
        const int nnb = tds_.dimension + 1;
        std::vector<Index>& stack = scratch_stack_;
        stack.clear();
        stack.push_back(d);
        tds_.cells[d].mark = Mark::InConflict;
        cells_out.push_back(d);

        do {
            const Index c = stack.back();
            stack.pop_back();
            for (int i = 0; i < nnb; ++i) {
                const Index test = tds_.cells[c].n[i];
                if (tds_.cells[test].mark == Mark::InConflict) continue;
                if (tds_.cells[test].mark == Mark::Clear) {
                    if (tds_.dimension == 3 ? in_conflict(test, p)
                                            : in_conflict_2(test, p)) {
                        stack.push_back(test);
                        tds_.cells[test].mark = Mark::InConflict;
                        cells_out.push_back(test);
                        continue;
                    }
                    tds_.cells[test].mark = Mark::OnBoundary;
                }
                facets_out.push_back({ c, i });
            }
        } while (!stack.empty());
    }

    // ---- insertion -------------------------------------------------------
    /**
     * Inserts point `idx`, starting the search from the cell cached on vertex
     * `hint`. Returns the vertex holding the point - the existing one if the
     * point duplicates an earlier vertex, in which case nothing is modified.
     *
     * Dispatch is on the current dimension. At 3 and 2 the point goes through
     * the conflict-region path with the matching conflict tester; below that
     * the location type decides directly, since there is no circumsphere to
     * test against.
     */
    Index insert(int idx, Index hint) {
        const double* p = pt_of(idx);
#ifdef DEL3D_VALIDATE
        {   // Distinguish "corrupted by this insertion" from "already corrupt".
            const std::string e0 = tds_.validate();
            if (!e0.empty())
                std::fprintf(stderr, "del3d: ALREADY broken on ENTRY to insert of point %d: "
                                     "%s\n", idx, e0.c_str());
        }
#endif

        LocateType lt = LT_OUTSIDE_AFFINE_HULL;
        int li = -1, lj = -1;
        // The hint is a *vertex*; the walk needs a cell, so take the incident
        // cell cached on it. (Passing a vertex index where a cell index is
        // expected would start the walk in an unrelated - possibly freed - cell
        // and the conflict region would then contain dead cells.)
        const Index start = (hint == kNull) ? infinite_cell() : tds_.vertices[hint].cell;
        Index c = locate(p, lt, li, lj, start);

        Index v;
        if (tds_.dimension >= 2) {
            if (tds_.dimension == 2 && lt == LT_OUTSIDE_AFFINE_HULL) {
                v = insert_outside_affine_hull(p);      // this point makes it 3D
            } else if (lt == LT_VERTEX) {
                return tds_.cells[c].v[li];             // duplicate point
            } else {
                std::vector<Index>& conflict = scratch_conflict_;
                std::vector<std::pair<Index, int>>& facets = scratch_facets_;
                conflict.clear();
                facets.clear();
                find_conflicts(c, p, conflict, facets);
#ifdef DEL3D_TRACE
                std::fprintf(stdout, "   [locate lt=%d cell=%d]  region %zu, facets:",
                             int(lt), c, conflict.size());
                for (const auto& f : facets) std::fprintf(stdout, " (%d,%d)", f.first, f.second);
                std::fprintf(stdout, "\n");
#endif
                v = tds_.create_vertex();
                // The two star rebuilds create cells in different orders, so
                // which one runs is fixed by the size of the hole rather than
                // chosen for speed. See Tds::insert_in_small_hole.
                if (tds_.dimension == 3 && detail::Tds::is_small_hole(facets.size()))
                    tds_.insert_in_small_hole(conflict, facets, v);
                else
                    tds_.insert_in_hole(conflict, facets.back().first,
                                        facets.back().second, v);
            }
        } else {
            // Dimension <= 1: no cells to be in conflict with.
            switch (lt) {
            case LT_VERTEX:
                return tds_.cells[c].v[li];             // duplicate point
            case LT_EDGE:
                v = tds_.insert_in_edge(c, li, lj);
                break;
            case LT_OUTSIDE_CONVEX_HULL:
                // Beyond the end of the 1D chain: the point splits the infinite
                // edge that closes it off.
                v = tds_.insert_in_edge(c, 0, 1);
                break;
            default:
                v = insert_outside_affine_hull(p);
                break;
            }
        }
        return v;
    }

    /**
     * Inserts a point that lies outside the current affine hull, raising the
     * dimension by one.
     *
     * The rebuild cones the existing structure to the new point on one fixed
     * side. If the new point is actually on the other side, every cell comes
     * out negatively oriented, so the orientation of the new point against the
     * old hull is tested first and the whole structure is flipped afterwards
     * when needed.
     */
    Index insert_outside_affine_hull(const double* p) {
        bool reorient_needed = false;
        if (tds_.dimension == 1) {
            Index c = infinite_cell();
            Index n = tds_.cells[c].n[tds_.cells[c].index_of_vertex(tds_.infinite)];
            const Sign o = coplanar_orientation(pt(tds_.cells[n].v[0]),
                                                pt(tds_.cells[n].v[1]), p);
            reorient_needed = (o == NEGATIVE);
        }
        if (tds_.dimension == 2) {
            Index c = infinite_cell();
            Index n = tds_.cells[c].n[tds_.cells[c].index_of_vertex(tds_.infinite)];
            const Sign o = orient3d(pt(tds_.cells[n].v[0]), pt(tds_.cells[n].v[1]),
                                    pt(tds_.cells[n].v[2]), p);
            reorient_needed = (o == NEGATIVE);
        }
        Index v = insert_increase_dimension(tds_.infinite);
        if (reorient_needed) reorient();
        return v;
    }
};

} // namespace

Point3 circumcenter(const Point3& p0, const Point3& p1,
                    const Point3& p2, const Point3& p3) {
    using detail::Expansion;

    // Translate p0 to the origin: a = p1-p0, b = p2-p0, c = p3-p0. Each
    // difference of two doubles is exact, so nothing is lost here.
    const Expansion ax = Expansion::from_diff(p1.x, p0.x);
    const Expansion ay = Expansion::from_diff(p1.y, p0.y);
    const Expansion az = Expansion::from_diff(p1.z, p0.z);
    const Expansion bx = Expansion::from_diff(p2.x, p0.x);
    const Expansion by = Expansion::from_diff(p2.y, p0.y);
    const Expansion bz = Expansion::from_diff(p2.z, p0.z);
    const Expansion cx = Expansion::from_diff(p3.x, p0.x);
    const Expansion cy = Expansion::from_diff(p3.y, p0.y);
    const Expansion cz = Expansion::from_diff(p3.z, p0.z);

    // b x c, c x a, a x b.
    const Expansion bcx = by * cz - bz * cy;
    const Expansion bcy = bz * cx - bx * cz;
    const Expansion bcz = bx * cy - by * cx;
    const Expansion cax = cy * az - cz * ay;
    const Expansion cay = cz * ax - cx * az;
    const Expansion caz = cx * ay - cy * ax;
    const Expansion abx = ay * bz - az * by;
    const Expansion aby = az * bx - ax * bz;
    const Expansion abz = ax * by - ay * bx;

    // Squared edge lengths from p0.
    const Expansion a2 = ax * ax + ay * ay + az * az;
    const Expansion b2 = bx * bx + by * by + bz * bz;
    const Expansion c2 = cx * cx + cy * cy + cz * cz;

    // The standard closed form:
    //   num = |a|^2 (b x c) + |b|^2 (c x a) + |c|^2 (a x b),  den = 2 a.(b x c)
    // den is twice the signed volume form, so it vanishes exactly when the four
    // points are coplanar.
    const Expansion nx = a2 * bcx + b2 * cax + c2 * abx;
    const Expansion ny = a2 * bcy + b2 * cay + c2 * aby;
    const Expansion nz = a2 * bcz + b2 * caz + c2 * abz;
    const Expansion den = (ax * bcx + ay * bcy + az * bcz) * Expansion(2.0);

    if (den.is_zero()) return p0;          // coplanar: no circumcentre
    // Everything above was exact; this is the single rounding step.
    const double d = den.to_double();
    return { p0.x + nx.to_double() / d,
             p0.y + ny.to_double() / d,
             p0.z + nz.to_double() / d };
}

namespace {

/// A 3x3 determinant over intervals. The association of the products is fixed:
/// with interval arithmetic a different grouping yields a different (still
/// valid but differently sized) interval, and the width decides whether the
/// midpoint is accepted or the exact fallback runs - so it changes the result,
/// not just its rounding.
detail::Interval det3(const detail::Interval& a00, const detail::Interval& a01,
                      const detail::Interval& a02, const detail::Interval& a10,
                      const detail::Interval& a11, const detail::Interval& a12,
                      const detail::Interval& a20, const detail::Interval& a21,
                      const detail::Interval& a22) {
    const detail::Interval m01 = a00 * a11 - a10 * a01;
    const detail::Interval m02 = a00 * a21 - a20 * a01;
    const detail::Interval m12 = a10 * a21 - a20 * a11;
    return m01 * a22 - m02 * a12 + m12 * a02;
}

} // namespace

Point3 circumcenter_lazy(const Point3& p, const Point3& q,
                         const Point3& r, const Point3& s) {
    using detail::Interval;

    Point3 out;
    // Which coordinates were not pinned down well enough by the interval pass
    // and need the exact one. The exact circumcentre costs roughly 15x an
    // interval evaluation, and the interval suffices for the great majority of
    // cells, so it is computed lazily below rather than up front.
    bool need_exact[3] = { false, false, false };
    {
        detail::RoundingUpward guard;

        // The same closed form as circumcenter(), evaluated in intervals with
        // p as the origin. The inputs are exact doubles, so they start as
        // degenerate (single-value) intervals.
        const Interval px(p.x), py(p.y), pz(p.z);
        const Interval qpx = Interval(q.x) - px, qpy = Interval(q.y) - py, qpz = Interval(q.z) - pz;
        const Interval qp2 = square(qpx) + square(qpy) + square(qpz);
        const Interval rpx = Interval(r.x) - px, rpy = Interval(r.y) - py, rpz = Interval(r.z) - pz;
        const Interval rp2 = square(rpx) + square(rpy) + square(rpz);
        const Interval spx = Interval(s.x) - px, spy = Interval(s.y) - py, spz = Interval(s.z) - pz;
        const Interval sp2 = square(spx) + square(spy) + square(spz);

        // Cramer's rule: each numerator is the determinant with one coordinate
        // column replaced by the squared lengths; the denominator is the
        // volume determinant.
        const Interval num_x = det3(qpy, qpz, qp2, rpy, rpz, rp2, spy, spz, sp2);
        const Interval num_y = det3(qpx, qpz, qp2, rpx, rpz, rp2, spx, spz, sp2);
        const Interval num_z = det3(qpx, qpy, qp2, rpx, rpy, rp2, spx, spy, sp2);
        const Interval den   = det3(qpx, qpy, qpz, rpx, rpy, rpz, spx, spy, spz);

        const Interval inv = Interval(1.0) / (Interval(2.0) * den);
        const Interval cx = px + num_x * inv;
        const Interval cy = py - num_y * inv;   // the y minor comes out negated
        const Interval cz = pz + num_z * inv;

        // Per coordinate: take the bound if the interval collapsed to a single
        // double, else its midpoint if it is already accurate to 1e-5 relative,
        // else defer to exact arithmetic.
        const double kPrec = 0.00001;
        const Interval iv[3] = { cx, cy, cz };
        double* co[3] = { &out.x, &out.y, &out.z };
        for (int i = 0; i < 3; ++i) {
            if (iv[i].is_point())                                *co[i] = iv[i].inf();
            else if (iv[i].has_smaller_relative_precision(kPrec)) *co[i] = iv[i].to_double();
            else                                                 need_exact[i] = true;
        }
    }
    if (need_exact[0] || need_exact[1] || need_exact[2]) {
        const Point3 ex = circumcenter(p, q, r, s);
        if (need_exact[0]) out.x = ex.x;
        if (need_exact[1]) out.y = ex.y;
        if (need_exact[2]) out.z = ex.z;
    }
    return out;
}

void Delaunay::build(const std::vector<Point3>& points) {
    cells_.clear();
    all_cells_.clear();
    incident_.assign(points.size(), {});
    n_vertices_ = 0;
    if (points.empty()) return;

    Builder builder(points);
    builder.build();
    detail::Tds& tds = builder.tds();

    // The vertex count is meaningful whatever the dimension, so it is taken
    // before the early return below.
    tds.vertices.for_each([&](detail::Index v) {
        if (v != tds.infinite) ++n_vertices_;
    });

    if (tds.dimension < 3) return;   // all points collinear or coplanar: no cells

    // Walk the cells in container order - which is the order the contract fixes
    // - reporting the infinite vertex as -1 and filtering the finite cells out
    // into cells_ while recording each vertex's incident cells.
    tds.cells.for_each([&](detail::Index c) {
        Tet a;
        for (int i = 0; i < 4; ++i) {
            const detail::Index vv = tds.cells[c].v[i];
            a[i] = (vv == tds.infinite) ? -1 : tds.vertices[vv].info;
        }
        all_cells_.push_back(a);
        if (tds.cells[c].has_vertex(tds.infinite)) return;
        Tet t;
        for (int i = 0; i < 4; ++i) t[i] = tds.vertices[tds.cells[c].v[i]].info;
        const int id = int(cells_.size());
        cells_.push_back(t);
        for (int i = 0; i < 4; ++i) incident_[t[i]].push_back(id);
    });
}

} // namespace del3d
