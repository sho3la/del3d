// del3d - the tetrahedral data structure.
//
// A cell-and-neighbour representation of a triangulation of the whole sphere:
// every cell has four vertices and four neighbours, the neighbour at index i
// being the one opposite vertex i. The convex hull is closed off by a single
// "infinite" vertex, so that hull facets belong to infinite cells and every
// facet has exactly two incident cells - which removes all boundary special
// cases from the walk and from the conflict-region search.
//
// `dimension` tracks the affine dimension of the point set inserted so far:
// -2 empty, -1 only the infinite vertex, then 0, 1, 2, 3 as points are added.
// Below dimension 3 the cells are degenerate (a cell uses dimension+1 of its
// four slots) and the insertion routines below have a case per dimension.
//
// Scope is limited to what an incremental Delaunay build needs: creation,
// adjacency, the three local splits (in a cell, in a facet, in an edge) and the
// star-hole rebuild. No removal, no I/O.
//
// Element order is part of the library's contract (see compact_container.h), so
// cells and vertices come from order-preserving containers and every routine
// that creates several cells does so in a fixed, documented sequence.
#ifndef DEL3D_TDS_H
#define DEL3D_TDS_H

#include "compact_container.h"

#include <array>
#include <cassert>
#include <cstdint>
#include <utility>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace del3d {
namespace detail {

using Index = int;
constexpr Index kNull = -1;

/// Static index tables over a cell's four local vertex slots.
struct Utils {
    /// tab_next_around_edge[i][j]: the local index of the facet to leave
    /// through when turning around the oriented edge (i,j). The diagonal is
    /// never queried and holds the sentinel 5.
    static constexpr int tab_next_around_edge[4][4] = {
        {5, 2, 3, 1},
        {3, 5, 0, 2},
        {1, 3, 5, 0},
        {2, 0, 1, 5} };
    /// tab_vertex_triple_index[i][k]: the k-th vertex of facet i, ordered so
    /// that the triple is positively oriented as seen from outside the cell.
    static constexpr int tab_vertex_triple_index[4][3] = {
        {1, 3, 2},
        {0, 2, 3},
        {0, 3, 1},
        {0, 1, 2} };
    static constexpr int ccw_map[3] = {1, 2, 0};
    static constexpr int cw_map[3]  = {2, 0, 1};

    /// Next / previous index in the cyclic order (0,1,2), used at dimension 2
    /// where a cell is a triangle.
    static int ccw(int i) { return ccw_map[i]; }
    static int cw(int i)  { return cw_map[i]; }
    /// Which neighbour to step to when rotating around the oriented edge (i,j).
    static int next_around_edge(int i, int j) { return tab_next_around_edge[i][j]; }
    static int vertex_triple_index(int i, int j) { return tab_vertex_triple_index[i][j]; }
};

struct Vertex {
    Index cell = kNull;      ///< one arbitrary incident cell; the entry point
                             ///< for any traversal of this vertex's star
    int   info = -1;         ///< the caller's point index; -1 for the infinite vertex
};

/// Per-cell mark used while collecting the conflict region of a new point.
/// InConflict: the new point lies inside this cell's circumsphere, so the cell
/// will be destroyed. OnBoundary: tested and rejected, remembered so it is not
/// tested twice. Marks must be back to Clear before the next insertion.
enum class Mark : unsigned char { Clear, InConflict, OnBoundary };

struct Cell {
    std::array<Index, 4> v{ {kNull, kNull, kNull, kNull} };  ///< vertices
    std::array<Index, 4> n{ {kNull, kNull, kNull, kNull} };  ///< n[i] is opposite v[i]
    Mark mark = Mark::Clear;

    /// Local slot holding the given vertex. Undefined if it is not present.
    int index_of_vertex(Index vh) const {
        for (int i = 0; i < 4; ++i) if (v[i] == vh) return i;
        assert(false); return 0;
    }
    /// Local slot whose neighbour is the given cell, i.e. the shared facet.
    int index_of_neighbor(Index ch) const {
        for (int i = 0; i < 4; ++i) if (n[i] == ch) return i;
        assert(false); return 0;
    }
    bool has_vertex(Index vh) const {
        return v[0] == vh || v[1] == vh || v[2] == vh || v[3] == vh;
    }
};

class Tds {
public:
    CompactContainer<Vertex> vertices;
    CompactContainer<Cell>   cells;

    int dimension = -2;      ///< -2 empty, then -1, 0, 1, 2, 3
    Index infinite = kNull;  ///< the vertex closing off the convex hull

    // ---- element creation (order matters: see compact_container.h) --------

    Index create_vertex() { return vertices.create(); }

    Index create_cell() { return cells.create(); }

    /// Creates a cell with the given vertices and no neighbours yet.
    Index create_cell(Index v0, Index v1, Index v2, Index v3) {
        const Index c = cells.create();
        cells[c].v = { {v0, v1, v2, v3} };
        return c;
    }

    void delete_cell(Index c) { cells.destroy(c); }

    // ---- connectivity ----------------------------------------------------

    /// Records the two halves of one adjacency: c0's i0-th neighbour becomes
    /// c1 and c1's i1-th becomes c0.
    void set_adjacency(Index c0, int i0, Index c1, int i1) {
        cells[c0].n[i0] = c1;
        cells[c1].n[i1] = c0;
    }

    bool is_infinite(Index vh) const { return vh == infinite; }
    bool is_infinite_cell(Index ch) const { return cells[ch].has_vertex(infinite); }

    /**
     * Rebuilds the star of a new vertex over a conflict region (dimension 3).
     *
     * The region is the set of cells marked InConflict; its boundary is the set
     * of facets with an InConflict cell on one side and an unmarked cell on the
     * other. One new cell is created per boundary facet, joining that facet to
     * the new vertex `v`, and the new cells are stitched to each other and to
     * the cells outside the region.
     *
     * Recursive form: `c` is a cell in conflict and `li` the local index of a
     * boundary facet of it (so cells[c].n[li] lies outside). The new cell is a
     * copy of c with vertex li replaced by v. Each of its other three
     * neighbours is found by rotating around the corresponding edge until the
     * region is left; if the cell that should sit there does not exist yet, it
     * is created by recursing.
     *
     * Cells are created in this traversal's order, which is part of the output
     * contract, so the traversal must not be reordered.
     *
     * `prev_ind2` is the facet the caller will connect itself, skipped here.
     * Returns the new cell.
     */
    Index create_star_3(Index v, Index c, int li, int prev_ind2 = -1) {
        assert(dimension == 3);
        assert(cells[c].mark == Mark::InConflict);
        assert(cells[cells[c].n[li]].mark != Mark::InConflict);

        const Index cnew = create_cell(cells[c].v[0], cells[c].v[1],
                                       cells[c].v[2], cells[c].v[3]);
        cells[cnew].v[li] = v;
        const Index c_li = cells[c].n[li];
        set_adjacency(cnew, li, c_li, cells[c_li].index_of_neighbor(c));

        // Find the other three neighbours of cnew.
        for (int ii = 0; ii < 4; ++ii) {
            if (ii == prev_ind2 || cells[cnew].n[ii] != kNull) continue;
            vertices[cells[cnew].v[ii]].cell = cnew;

            // The two vertices shared with the neighbour across facet ii,
            // oriented so that (ii, vj1, vj2, li) is positively oriented.
            const Index vj1 = cells[c].v[Utils::next_around_edge(ii, li)];
            const Index vj2 = cells[c].v[Utils::next_around_edge(li, ii)];
            Index cur = c;
            int   zz  = ii;
            Index n   = cells[cur].n[zz];
            // Turn around the oriented edge (vj1, vj2) until leaving the region.
            while (cells[n].mark == Mark::InConflict) {
                assert(n != c);
                cur = n;
                zz  = Utils::next_around_edge(cells[n].index_of_vertex(vj1),
                                              cells[n].index_of_vertex(vj2));
                n   = cells[cur].n[zz];
            }
            // n is now outside the region and cur is the last cell inside it,
            // so (cur, zz) is another boundary facet.
            cells[n].mark = Mark::Clear;      // this facet is now handled

            const int jj1 = cells[n].index_of_vertex(vj1);
            const int jj2 = cells[n].index_of_vertex(vj2);
            const Index vvv = cells[n].v[Utils::next_around_edge(jj1, jj2)];
            Index nnn = cells[n].n[Utils::next_around_edge(jj2, jj1)];
            const int zzz = cells[nnn].index_of_vertex(vvv);
            if (nnn == cur) {
                // The neighbour relation points back into the region: the cell
                // that belongs there has not been created yet, so create it.
                nnn = create_star_3(v, nnn, zz, zzz);
            }
            set_adjacency(nnn, zzz, cnew, ii);
        }
#ifdef DEL3D_VALIDATE
        // On return every neighbour of the new cell must be set, except the one
        // the caller is about to set (prev_ind2).
        for (int ii = 0; ii < 4; ++ii) {
            if (ii == prev_ind2) continue;
            if (cells[cnew].n[ii] == kNull) {
                std::fprintf(stderr,
                    "del3d: create_star_3 left cnew=%d neighbour %d unset "
                    "(c=%d li=%d prev_ind2=%d)\n", cnew, ii, c, li, prev_ind2);
                std::abort();
            }
        }
#endif
        return cnew;
    }

    /**
     * Visits every cell incident to the oriented edge (v[i], v[j]) of cell c,
     * exactly once, by repeatedly stepping to the neighbour across
     * next_around_edge(). The cells around an interior edge form a cycle, so
     * the loop terminates when it returns to c.
     */
    template <class F>
    void cells_around_edge(Index c, int i, int j, F&& f) const {
        const Index vs = cells[c].v[i], vt = cells[c].v[j];
        Index cur = c;
        do {
            f(cur);
            cur = cells[cur].n[Utils::next_around_edge(cells[cur].index_of_vertex(vs),
                                                       cells[cur].index_of_vertex(vt))];
        } while (cur != c);
    }

    /**
     * Splits cell c into four by a new vertex in its interior (dimension 3).
     *
     * c itself is reused as one of the four - its vertex 0 is overwritten by
     * the new vertex - and three cells are created, in the order c3, c2, c1
     * (the one opposite the old vertex 3 first). Each new cell inherits one of
     * c's outward neighbours; the three internal facets are stitched pairwise.
     * Returns the new vertex.
     */
    Index insert_in_cell(Index c) {
        assert(dimension == 3);
        const Index v = create_vertex();

        const Index v0 = cells[c].v[0], v1 = cells[c].v[1];
        const Index v2 = cells[c].v[2], v3 = cells[c].v[3];
        const Index n1 = cells[c].n[1], n2 = cells[c].n[2], n3 = cells[c].n[3];

        const Index c3 = create_cell(v0, v1, v2, v);
        const Index c2 = create_cell(v0, v1, v,  v3);
        const Index c1 = create_cell(v0, v,  v2, v3);

        // Each new cell shares the facet opposite v0 with the reused cell c.
        set_adjacency(c3, 0, c, 3);
        set_adjacency(c2, 0, c, 2);
        set_adjacency(c1, 0, c, 1);

        // ... and one facet with each of the other two new cells.
        set_adjacency(c2, 3, c3, 2);
        set_adjacency(c1, 3, c3, 1);
        set_adjacency(c1, 2, c2, 1);

        // The outward facets keep their old neighbours.
        set_adjacency(n1, cells[n1].index_of_neighbor(c), c1, 1);
        set_adjacency(n2, cells[n2].index_of_neighbor(c), c2, 2);
        set_adjacency(n3, cells[n3].index_of_neighbor(c), c3, 3);

        cells[c].v[0] = v;
        vertices[v0].cell = c1;
        vertices[v].cell = c;
        return v;
    }

    /**
     * Splits the facet (c, i) by a new vertex lying exactly on it.
     *
     * At dimension 3 the facet is shared by c and its neighbour d, and each of
     * the two cells is split into three: two new cells plus the original,
     * reused. At dimension 2 the "facet" is the whole cell (a triangle), which
     * is split into three triangles. Returns the new vertex.
     */
    Index insert_in_facet(Index c, int i) {
        const Index v = create_vertex();

        if (dimension == 3) {
            // Order i, i1, i2, i3 so that the quadruple is positively oriented;
            // the new cells then replace vertices in that order.
            int i1, i2, i3;
            if ((i & 1) == 0) { i1 = (i + 1) & 3; i2 = (i + 2) & 3; i3 = 6 - i - i1 - i2; }
            else              { i1 = (i + 1) & 3; i2 = (i + 3) & 3; i3 = 6 - i - i1 - i2; }
            const Index vi = cells[c].v[i],  v1 = cells[c].v[i1];
            const Index v2 = cells[c].v[i2], v3 = cells[c].v[i3];

            // This side of the facet: two new cells, then c is reused.
            Index nc = cells[c].n[i1];
            const Index cnew1 = create_cell(vi, v, v2, v3);
            set_adjacency(cnew1, 1, nc, cells[nc].index_of_neighbor(c));
            set_adjacency(cnew1, 3, c, i1);
            vertices[v3].cell = cnew1;

            nc = cells[c].n[i2];
            const Index cnew2 = create_cell(vi, v1, v, v3);
            set_adjacency(cnew2, 2, nc, cells[nc].index_of_neighbor(c));
            set_adjacency(cnew2, 3, c, i2);
            set_adjacency(cnew1, 2, cnew2, 1);

            cells[c].v[i3] = v;

            // The other side of the facet, split the same way.
            const Index d = cells[c].n[i];
            const int j  = cells[d].index_of_neighbor(c);
            const int j1 = cells[d].index_of_vertex(v1);
            const int j2 = cells[d].index_of_vertex(v2);
            const int j3 = 6 - j - j1 - j2;

            Index nd = cells[d].n[j1];
            const Index dnew1 = create_cell(cells[d].v[j], v, v3, v2);
            set_adjacency(dnew1, 1, nd, cells[nd].index_of_neighbor(d));
            set_adjacency(dnew1, 2, d, j1);
            set_adjacency(dnew1, 0, cnew1, 0);   // across the split facet

            nd = cells[d].n[j2];
            const Index dnew2 = create_cell(cells[d].v[j], v1, v3, v);
            set_adjacency(dnew2, 3, nd, cells[nd].index_of_neighbor(d));
            set_adjacency(dnew2, 2, d, j2);
            set_adjacency(dnew2, 0, cnew2, 0);
            set_adjacency(dnew1, 3, dnew2, 1);

            cells[d].v[j3] = v;
            vertices[v].cell = d;
        } else {                                  // dimension 2
            Index n = cells[c].n[2];
            const Index cnew = create_cell(cells[c].v[0], cells[c].v[1], v, kNull);
            set_adjacency(cnew, 2, n, cells[n].index_of_neighbor(c));
            set_adjacency(cnew, 0, c, 2);
            vertices[cells[c].v[0]].cell = cnew;

            n = cells[c].n[1];
            const Index dnew = create_cell(cells[c].v[0], v, cells[c].v[2], kNull);
            set_adjacency(dnew, 1, n, cells[n].index_of_neighbor(c));
            set_adjacency(dnew, 0, c, 1);
            set_adjacency(dnew, 2, cnew, 1);

            cells[c].v[0] = v;
            vertices[v].cell = c;
        }
        return v;
    }

    /**
     * Splits the edge (v[i], v[j]) of cell c by a new vertex lying on it.
     *
     * At dimension 3 every cell around that edge has to be split, so the ring
     * of incident cells is collected, marked as the conflict region and
     * re-starred. At dimension 2 the edge has two incident triangles, each
     * split in two. At dimension 1 the triangulation is a chain of edges and
     * the split is local. Returns the new vertex.
     */
    Index insert_in_edge(Index c, int i, int j) {
        if (dimension == 3) {
            std::vector<Index> conflict;
            conflict.reserve(32);
            cells_around_edge(c, i, j, [&](Index cc) {
                conflict.push_back(cc);
                cells[cc].mark = Mark::InConflict;
            });
            const Index v = create_vertex();
            return insert_in_hole(conflict, c, i, v);
        }
        if (dimension == 2) {
            const Index v = create_vertex();
            const int k = 3 - i - j;                  // the third vertex of the triangle
            const Index d = cells[c].n[k];            // the triangle across the edge
            const int kd = cells[d].index_of_neighbor(c);
            const int id = cells[d].index_of_vertex(cells[c].v[i]);
            const int jd = cells[d].index_of_vertex(cells[c].v[j]);

            // Each old triangle keeps the half towards j; the new one takes the
            // half towards i.
            const Index cnew = create_cell();
            cells[cnew].v[i] = cells[c].v[i];
            vertices[cells[c].v[i]].cell = cnew;
            cells[cnew].v[j] = v;
            cells[cnew].v[k] = cells[c].v[k];
            cells[c].v[i] = v;

            const Index dnew = create_cell();
            cells[dnew].v[id] = cells[d].v[id];
            cells[dnew].v[jd] = v;
            cells[dnew].v[kd] = cells[d].v[kd];
            cells[d].v[id] = v;

            Index nj = cells[c].n[j];
            set_adjacency(cnew, i, c, j);
            set_adjacency(cnew, j, nj, cells[nj].index_of_neighbor(c));

            nj = cells[d].n[jd];
            set_adjacency(dnew, id, d, jd);
            set_adjacency(dnew, jd, nj, cells[nj].index_of_neighbor(d));

            set_adjacency(cnew, k, dnew, kd);         // across the new edge
            vertices[v].cell = cnew;
            return v;
        }
        // dimension 1: split the segment c into c + cnew.
        const Index v = create_vertex();
        const Index cnew = create_cell(v, cells[c].v[1], kNull, kNull);
        vertices[cells[c].v[1]].cell = cnew;
        cells[c].v[1] = v;
        set_adjacency(cnew, 0, cells[c].n[0], 1);
        set_adjacency(cnew, 1, c, 0);
        vertices[v].cell = cnew;
        return v;
    }

    /**
     * The dimension-2 counterpart of create_star_3: rebuilds the star of a new
     * vertex over a conflict region of triangles.
     *
     * The boundary of a planar region is a single cycle of edges, so instead of
     * recursing this walks that cycle: starting from the boundary edge (c, li)
     * it turns around the current boundary vertex until it leaves the region,
     * emits one triangle (v, v1, next), links it to the previous one, and
     * advances - closing the ring between the last and first triangles at the
     * end. Returns the last cell created.
     */
    Index create_star_2(Index v, Index c, int li) {
        assert(dimension == 2);
        Index cnew = kNull;
        int i1 = Utils::ccw(li);                  // v, i1, i2 positively oriented
        Index bound = c;
        Index v1 = cells[c].v[i1];
        const int ind = cells[cells[c].n[li]].index_of_neighbor(c);
        Index pnew = kNull;
        do {
            Index cur = bound;
            // Turn around v1 until the boundary of the region is reached.
            while (cells[cells[cur].n[Utils::cw(i1)]].mark == Mark::InConflict) {
                cur = cells[cur].n[Utils::cw(i1)];
                i1 = cells[cur].index_of_vertex(v1);
            }
            cells[cells[cur].n[Utils::cw(i1)]].mark = Mark::Clear;

            cnew = create_cell(v, v1, cells[cur].v[Utils::ccw(i1)], kNull);
            const Index nb = cells[cur].n[Utils::cw(i1)];
            set_adjacency(cnew, 0, nb, cells[nb].index_of_neighbor(cur));
            cells[cnew].n[1] = kNull;             // set when the next one is made
            cells[cnew].n[2] = pnew;
            vertices[v1].cell = cnew;
            if (pnew != kNull) cells[pnew].n[1] = cnew;

            bound = cur;
            i1 = Utils::ccw(i1);
            v1 = cells[bound].v[i1];
            pnew = cnew;
        } while (v1 != cells[c].v[Utils::ccw(li)]);   // back at the start
        // Close the ring between the last and the first created cell.
        const Index first = cells[cells[c].n[li]].n[ind];
        set_adjacency(cnew, 1, first, 2);
        return cnew;
    }

    /**
     * Structural self-check, used by the tests and the DEL3D_VALIDATE builds.
     * Returns a description of the first broken invariant, or an empty string.
     *
     * Checks that every live cell has non-null, distinct vertices; that every
     * neighbour index refers to a live cell and points back; that no conflict
     * mark survived the last insertion; and that every vertex's cached incident
     * cell is alive and really contains it. Running this after each insertion
     * localises the insertion that corrupts the structure, instead of the
     * failure surfacing later as an out-of-range index or a walk that never
     * terminates.
     */
    std::string validate() const {
        std::string err;
        auto fail = [&](const std::string& s) { if (err.empty()) err = s; };
        const int nv = dimension + 1;          // used vertex/neighbour slots per cell
        if (dimension < 1) return err;

        cells.for_each([&](Index c) {
            const Cell& cell = cells[c];
            for (int i = 0; i < nv; ++i) {
                if (cell.v[i] == kNull) { fail("cell " + std::to_string(c) +
                                               " has a null vertex " + std::to_string(i)); return; }
                for (int j = i + 1; j < nv; ++j)
                    if (cell.v[i] == cell.v[j])
                        fail("cell " + std::to_string(c) + " repeats vertex " +
                             std::to_string(cell.v[i]));
            }
            for (int i = 0; i < nv; ++i) {
                const Index n = cell.n[i];
                if (n == kNull) { fail("cell " + std::to_string(c) +
                                       " has a null neighbour " + std::to_string(i)); return; }
                if (n < 0 || n >= Index(cells.capacity()) || !cells.is_used(n)) {
                    fail("cell " + std::to_string(c) + " neighbour " + std::to_string(i) +
                         " -> dead cell " + std::to_string(n)); return;
                }
                bool back = false;
                for (int k = 0; k < nv; ++k) if (cells[n].n[k] == c) back = true;
                if (!back)
                    fail("cell " + std::to_string(c) + " neighbour " + std::to_string(i) +
                         " -> " + std::to_string(n) + " which does not point back");
            }
        });
        // Marks are set by the conflict search and cleared by the star rebuild,
        // so between insertions every one of them must be Clear again.
        cells.for_each([&](Index c) {
            if (cells[c].mark != Mark::Clear)
                fail("cell " + std::to_string(c) + " left marked (" +
                     std::to_string(int(cells[c].mark)) + ") after the insertion");
        });
        vertices.for_each([&](Index v) {
            const Index c = vertices[v].cell;
            if (c == kNull || !cells.is_used(c))
                fail("vertex " + std::to_string(v) + " has a dead cell");
            else if (!cells[c].has_vertex(v))
                fail("vertex " + std::to_string(v) + " -> cell " + std::to_string(c) +
                     " which does not contain it");
        });
        return err;
    }

    /**
     * Lookup table mapping a directed vertex pair (u,v) to the new cell that
     * carries it, used to stitch a freshly created star together.
     *
     * Open addressing with linear probing over a buffer that survives between
     * insertions; "clearing" it is a stamp bump rather than a deallocation.
     * This replaces a node-based map rebuilt inside every insertion, which on a
     * 35k-point input cost a few million short-lived allocations.
     *
     * Iteration order differs from an ordered map's, which is harmless here:
     * the second pass only needs to pair each (u,v) with its opposite (v,u),
     * and the adjacency it then records is symmetric.
     */
    class PairMap {
    public:
        /// Which new cell (by position in the facet list) holds the pair, and
        /// the local slot in it opposite that pair's edge.
        struct Val { int facet; int idx; };

        /// Prepares the map for about `nentries` insertions, growing to keep
        /// the load factor near 25% and invalidating all previous entries.
        void reset(std::size_t nentries) {
            std::size_t cap = 64;
            while (cap < nentries * 4) cap <<= 1;
            if (slots_.size() < cap) { slots_.assign(cap, Slot{}); stamp_ = 0; }
            mask_ = slots_.size() - 1;
            if (++stamp_ == 0) {                    // stamp wrapped: really clear
                slots_.assign(slots_.size(), Slot{});
                stamp_ = 1;
            }
        }

        void set(Index u, Index v, int facet, int idx) {
            std::size_t i = hash(u, v) & mask_;
            for (;;) {
                Slot& s = slots_[i];
                // A slot from an earlier stamp counts as empty.
                if (s.stamp != stamp_) { s = Slot{ u, v, { facet, idx }, stamp_ }; return; }
                if (s.u == u && s.v == v) { s.val = { facet, idx }; return; }
                i = (i + 1) & mask_;
            }
        }

        bool get(Index u, Index v, Val& out) const {
            std::size_t i = hash(u, v) & mask_;
            for (;;) {
                const Slot& s = slots_[i];
                if (s.stamp != stamp_) return false;
                if (s.u == u && s.v == v) { out = s.val; return true; }
                i = (i + 1) & mask_;
            }
        }

        template <class F>
        void for_each(F f) const {
            for (const Slot& s : slots_)
                if (s.stamp == stamp_) f(s.u, s.v, s.val);
        }

    private:
        struct Slot { Index u = 0, v = 0; Val val{ 0, 0 }; unsigned stamp = 0; };
        /// Asymmetric mix, so that (u,v) and (v,u) land in different slots.
        static std::size_t hash(Index u, Index v) {
            return std::size_t(std::uint32_t(u) * 0x9E3779B1u) ^
                   (std::size_t(std::uint32_t(v) * 0x85EBCA6Bu) << 1);
        }
        std::vector<Slot> slots_;
        std::size_t mask_ = 0;
        unsigned stamp_ = 0;
    };

    PairMap pair_scratch_;                 // reused by insert_in_small_hole
    std::vector<Index> new_cells_scratch_; // ditto

    /// Conflict regions with at most this many boundary facets - nearly all of
    /// them in practice - are re-starred by insert_in_small_hole below.
    static constexpr std::size_t kMaxSmallHoleFacets = 128;
    static bool is_small_hole(std::size_t nfacets) { return nfacets <= kMaxSmallHoleFacets; }

    /// The same facet as seen from the cell on the other side.
    std::pair<Index, int> mirror_facet(Index c, int i) const {
        const Index n = cells[c].n[i];
        return { n, cells[n].index_of_neighbor(c) };
    }

    /**
     * Iterative star rebuild for a region with few boundary facets.
     *
     * Two passes. The first walks the boundary facet list in order, creating
     * one cell per facet, attaching it to the cell outside the region, and
     * registering its three directed vertex pairs in the pair map. The second
     * matches each directed pair (u,v) with its opposite (v,u) - the two new
     * cells sharing that edge - and records the adjacency between them. The
     * conflict cells are unmarked and deleted at the end.
     *
     * This is not merely a faster path than create_star_3: it creates cells in
     * *facet-list order*, whereas the recursive version creates them in the
     * order its traversal happens to reach them. Since cell order is part of
     * the contract, which of the two runs is itself part of the specification -
     * see is_small_hole().
     */
    template <class CellRange, class FacetRange>
    Index insert_in_small_hole(const CellRange& conflict, const FacetRange& facets,
                               Index nv) {
        std::vector<Index>& new_cells = new_cells_scratch_;
        new_cells.assign(facets.size(), kNull);
        PairMap& pair_map = pair_scratch_;
        pair_map.reset(3 * facets.size());

        for (std::size_t k = 0; k < facets.size(); ++k) {
            // Take the facet from the outside, where the cell survives.
            const std::pair<Index, int> f = mirror_facet(facets[k].first, facets[k].second);
            cells[f.first].mark = Mark::Clear;                  // was OnBoundary
            const Index u = cells[f.first].v[Utils::vertex_triple_index(f.second, 0)];
            const Index v = cells[f.first].v[Utils::vertex_triple_index(f.second, 1)];
            const Index w = cells[f.first].v[Utils::vertex_triple_index(f.second, 2)];
            // The outside cell survives, so it is a safe incident cell to cache.
            vertices[u].cell = f.first;
            vertices[v].cell = f.first;
            vertices[w].cell = f.first;

            // v, u, w - swapped relative to the outward triple, so that the new
            // cell is positively oriented as seen from the new vertex.
            const Index nc = create_cell(v, u, w, nv);
            new_cells[k] = nc;
            vertices[nv].cell = nc;
            cells[nc].n[3] = f.first;                 // facet 3 faces outwards
            cells[f.first].n[f.second] = nc;

            // Register the three edges of the facet, each with the local slot
            // opposite it (i.e. the third vertex) in the new cell.
            pair_map.set(u, v, int(k), cells[nc].index_of_vertex(w));
            pair_map.set(v, w, int(k), cells[nc].index_of_vertex(u));
            pair_map.set(w, u, int(k), cells[nc].index_of_vertex(v));
        }

        // Pair each directed edge with its opposite; u < v visits each edge
        // once. Both halves of the adjacency are written here.
        pair_map.for_each([&](Index u, Index v, const PairMap::Val& val) {
            if (!(u < v)) return;
            PairMap::Val other;
            if (!pair_map.get(v, u, other)) return;
            const Index a = new_cells[val.facet];
            const Index b = new_cells[other.facet];
            cells[a].n[val.idx] = b;
            cells[b].n[other.idx] = a;
        });

        for (Index c : conflict) cells[c].mark = Mark::Clear;
        for (Index c : conflict) delete_cell(c);
        return nv;
    }

    /**
     * Recursive star rebuild: star the hole from the boundary facet (begin, i),
     * then delete the conflict cells.
     *
     * The deletion happens *after* the star is complete, so the freed slots are
     * reused by the next insertion rather than by this one. That ordering is
     * what determines where subsequent cells land, so it must not be
     * interleaved.
     */
    template <class CellRange>
    Index insert_in_hole(const CellRange& conflict, Index begin, int i, Index newv) {
        const Index cnew = (dimension == 3) ? create_star_3(newv, begin, i)
                                            : create_star_2(newv, begin, i);
        vertices[newv].cell = cnew;
        for (Index ch : conflict) delete_cell(ch);
        return newv;
    }
};

} // namespace detail
} // namespace del3d

#endif // DEL3D_TDS_H
