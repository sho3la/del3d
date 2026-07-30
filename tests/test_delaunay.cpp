// del3d - end-to-end test of the triangulation.
//
// Rather than compare against another library, this checks the properties that
// define the output, each one against the exact-integer reference in
// tests/exact.h:
//
//   * every finite cell is positively oriented;
//   * the empty-sphere property holds globally - no input point lies strictly
//     inside the circumsphere of any cell;
//   * the cell complex is closed - every triangular facet is shared by exactly
//     two cells, counting the infinite ones;
//   * every distinct input point is a vertex, and duplicates are not;
//   * incident_cells() agrees with finite_cells();
//   * building twice gives byte-identical output, cell order included;
//   * and on a fixed set of inputs the cell sequence hashes to a pinned value,
//     which is what "the same tetrahedra and the same cell order on every
//     compiler and platform" means in practice.
//
// The cases are chosen for their degeneracies: a lattice is massively
// cospherical, and a plane plus one apex forces the whole low-dimension
// bootstrap (0 -> 1 -> 2 -> 3) before any tetrahedron exists.
#include "del3d/delaunay.h"

#include "exact.h"
#include "pointsets.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <map>
#include <random>
#include <set>
#include <vector>

static int failures = 0;

static void fail(const char* what, const char* detail = nullptr) {
    ++failures;
    if (detail) std::printf("    FAIL: %s (%s)\n", what, detail);
    else        std::printf("    FAIL: %s\n", what);
}

static std::uint64_t hash_cells(const std::vector<del3d::Delaunay::Tet>& cells) {
    std::uint64_t h = 1469598103934665603ull;                 // FNV-1a, 64 bit
    auto mix = [&h](std::uint64_t v) {
        for (int b = 0; b < 8; ++b) {
            h ^= (v >> (8 * b)) & 0xffull;
            h *= 1099511628211ull;
        }
    };
    mix(cells.size());
    for (const del3d::Delaunay::Tet& t : cells)
        for (int v : t) mix(static_cast<std::uint64_t>(static_cast<std::int64_t>(v)));
    return h;
}

// Exhaustive empty-sphere check: every input point against every cell. Exact,
// and quadratic, so it is run in full only on the smaller cases; on the larger
// ones a fixed subset of cells is checked instead.
static void check_empty_spheres(const std::vector<del3d::Point3>& pts,
                                const std::vector<del3d::Delaunay::Tet>& cells,
                                std::size_t cell_budget, long long& on_sphere) {
    std::vector<const double*> p(pts.size());
    for (std::size_t i = 0; i < pts.size(); ++i) p[i] = &pts[i].x;

    const std::size_t stride = cells.size() > cell_budget
                             ? (cells.size() + cell_budget - 1) / cell_budget : 1;

    int reported = 0;
    for (std::size_t ci = 0; ci < cells.size(); ci += stride) {
        const del3d::Delaunay::Tet& t = cells[ci];
        const double* a = p[t[0]];
        const double* b = p[t[1]];
        const double* c = p[t[2]];
        const double* d = p[t[3]];
        if (exact::orient3d(a, b, c, d) != 1) {
            if (reported++ < 5) fail("cell is not positively oriented");
            continue;
        }
        for (std::size_t i = 0; i < p.size(); ++i) {
            if (int(i) == t[0] || int(i) == t[1] || int(i) == t[2] || int(i) == t[3]) continue;
            const int s = exact::in_sphere_raw(a, b, c, d, p[i]);
            if (s > 0) {
                if (reported++ < 5) fail("point strictly inside a circumsphere");
            } else if (s == 0) {
                ++on_sphere;      // cospherical: allowed, and the interesting case
            }
        }
    }
}

static void check_case(const char* what, const std::vector<del3d::Point3>& pts,
                       std::size_t cell_budget, std::uint64_t expect_hash) {
    del3d::Delaunay d;
    d.build(pts);
    const std::vector<del3d::Delaunay::Tet>& cells = d.finite_cells();
    const std::vector<del3d::Delaunay::Tet>& all   = d.all_cells();

    // --- the empty-sphere property, and cell orientation ---
    long long on_sphere = 0;
    check_empty_spheres(pts, cells, cell_budget, on_sphere);

    // --- the complex is closed: every facet is shared by exactly two cells ---
    {
        std::map<std::array<int, 3>, int> facets;
        for (const del3d::Delaunay::Tet& t : all)
            for (int skip = 0; skip < 4; ++skip) {
                std::array<int, 3> f{};
                for (int k = 0, w = 0; k < 4; ++k) if (k != skip) f[w++] = t[k];
                std::sort(f.begin(), f.end());
                ++facets[f];
            }
        std::size_t bad = 0;
        for (const auto& kv : facets) if (kv.second != 2) ++bad;
        if (bad) fail("facets not shared by exactly two cells");
    }

    // --- every distinct point is a vertex, and nothing else is ---
    {
        std::set<std::array<double, 3>> distinct;
        for (const del3d::Point3& q : pts) distinct.insert({q.x, q.y, q.z});

        std::set<int> used;
        for (const del3d::Delaunay::Tet& t : all)
            for (int v : t) if (v >= 0) used.insert(v);

        // A triangulation of fewer than two distinct points has no cells at all,
        // so only the non-degenerate case can be checked this way.
        if (!all.empty()) {
            if (used.size() != d.number_of_vertices())
                fail("vertices in the cell table differ from number_of_vertices()");
            std::set<std::array<double, 3>> used_pts;
            for (int v : used) used_pts.insert({pts[v].x, pts[v].y, pts[v].z});
            if (used_pts != distinct) fail("some distinct input point is not a vertex");
        }
        if (d.number_of_vertices() != distinct.size())
            fail("number_of_vertices() differs from the count of distinct points");
    }

    // --- incident_cells() agrees with finite_cells() ---
    {
        std::vector<std::vector<int>> want(pts.size());
        for (std::size_t ci = 0; ci < cells.size(); ++ci)
            for (int v : cells[ci])
                if (want[v].empty() || want[v].back() != int(ci)) want[v].push_back(int(ci));

        const std::vector<std::vector<int>>& got = d.incident_cells();
        if (got.size() != pts.size()) {
            fail("incident_cells() has the wrong length");
        } else {
            for (std::size_t i = 0; i < pts.size(); ++i) {
                std::vector<int> a = got[i], b = want[i];
                std::sort(a.begin(), a.end());
                std::sort(b.begin(), b.end());
                if (a != b) { fail("incident_cells() disagrees with finite_cells()"); break; }
            }
        }
    }

    // --- determinism ---
    del3d::Delaunay again;
    again.build(pts);
    if (again.finite_cells() != cells || again.all_cells() != all)
        fail("rebuilding the same points gave a different answer");

    // --- the pinned cell sequence ---
    const std::uint64_t h = hash_cells(cells);
    if (expect_hash != 0 && h != expect_hash)
        fail("cell sequence differs from the pinned value");

    std::printf("  %-24s %5zu pts  %6zu cells  %6lld cospherical  hash %016llx%s\n",
                what, pts.size(), cells.size(), on_sphere,
                static_cast<unsigned long long>(h),
                expect_hash == 0 ? "  (not pinned)" : "");
}

static std::vector<del3d::Point3> to_points(const std::vector<pointsets::Pt>& v) {
    std::vector<del3d::Point3> out;
    out.reserve(v.size());
    for (const pointsets::Pt& p : v) out.push_back({p[0], p[1], p[2]});
    return out;
}

int main() {
    std::printf("del3d::Delaunay - properties, determinism and the pinned cell order\n");

    // Degenerate input that produces no tetrahedron at all: these must not
    // crash and must report the right vertex count.
    {
        check_case("empty", {}, 64, 0);
        check_case("one point", {{0, 0, 0}}, 64, 0);
        check_case("two points", {{0, 0, 0}, {1, 0, 0}}, 64, 0);
        std::vector<del3d::Point3> dup(20, del3d::Point3{1.0, 2.0, 3.0});
        check_case("20 identical points", dup, 64, 0);
        std::vector<del3d::Point3> line;
        for (int i = 0; i < 10; ++i) line.push_back({double(i), 2.0 * i, 3.0 * i});
        check_case("collinear", line, 64, 0);
        std::vector<del3d::Point3> flat;
        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 5; ++j) flat.push_back({double(i), double(j), 0.0});
        check_case("coplanar", flat, 64, 0);
    }

    // Generic position. The hashes are the pinned cell sequences: any platform
    // that produces a different tetrahedron, or the same ones in a different
    // order, fails here.
    {
        const struct { int n; std::uint64_t hash; } gen[] = {
            {   5, 0x5a7aed0478362803ull },
            {   8, 0x0b63a88c5b428febull },
            {  20, 0xdc8940eebc8963d7ull },
            { 100, 0x92dcfb61db16f9c3ull },
            { 500, 0xab8f0ef0a819629aull },
        };
        for (const auto& g : gen) {
            char name[64];
            std::snprintf(name, sizeof name, "random n=%d", g.n);
            check_case(name, to_points(pointsets::random_cube(std::size_t(g.n), 4u)),
                       g.n <= 100 ? 4096u : 64u, g.hash);
        }
    }

    // Degenerate by construction.
    {
        check_case("lattice 5x5x5", to_points(pointsets::lattice(5, 5, 5)),
                   4096, 0xf0f09ebcc3974328ull);
        std::vector<del3d::Point3> plane;
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j) plane.push_back({double(i), double(j), 0.0});
        plane.push_back({2.0, 2.0, 3.0});
        check_case("plane + 1 apex", plane, 4096, 0x235f7d68eec125daull);
        check_case("box shell 5^3", to_points(pointsets::box_shell(5)),
                   256, 0x2839796a61e16c25ull);
    }

    // Quantised surface samples.
    check_case("cylinder (quantised)", to_points(pointsets::cylinder(300)),
               128, 0x9403338c5b1fde7eull);
    check_case("sphere (quantised)", to_points(pointsets::sphere(300)),
               128, 0x85b7cd8bf3bdc99full);
    check_case("torus (quantised)", to_points(pointsets::torus(300)),
               128, 0x67adf648ba1948a7ull);

    std::printf("\n%s (%d failures)\n",
                failures ? "FAILED" : "ALL PROPERTIES HOLD", failures);
    return failures ? 1 : 0;
}
