// del3d - coplanar predicate conformance test.
//
// Checks coplanar_orientation (both arities), collinear and coplanar_in_circle
// (raw *and* perturbed) against the independent exact-integer reference in
// tests/exact.h, on points that are exactly coplanar by construction - which is
// the only case these predicates are ever asked about.
//
// The reference derives the in-circle answer from an in-sphere test against an
// apex raised along the plane normal, where the library builds the equivalent
// 4x4 determinant directly; the two share no formula and no arithmetic.
//
// Four stages: answers known by hand, which also pin the sign conventions; an
// integer grid in the z = 0 plane; the same grid on the tilted plane z = x + y,
// so that no single coordinate projection is the degenerate one; and a ring of
// exactly cocircular points, which forces the perturbation branch on nearly
// every call.
#include "del3d/predicates2.h"

#include "exact.h"
#include "pointsets.h"

#include <array>
#include <cstdio>
#include <vector>

static int failures = 0;
static long long n_or3 = 0, n_or4 = 0, n_ic = 0, deg_or3 = 0, deg_ic = 0;

static void fail(const char* what) {
    if (++failures < 10) std::printf("  MISMATCH: %s\n", what);
}

static void expect(bool ok, const char* what) {
    if (!ok) { ++failures; std::printf("  WRONG: %s\n", what); }
}

static void check_triple(const double* p, const double* q, const double* r) {
    // collinear() is defined for every input; coplanar_orientation() is not -
    // it has non-collinearity as a precondition, so screen for that first.
    const bool ref_col = exact::collinear(p, q, r);
    if (del3d::collinear(p, q, r) != ref_col) fail("collinear");
    ++n_or3;
    if (ref_col) { ++deg_or3; return; }
    if (del3d::coplanar_orientation(p, q, r) != exact::coplanar_orientation(p, q, r))
        fail("coplanar_orientation(p,q,r)");
}

static void check_quad(const double* p, const double* q, const double* r, const double* s) {
    if (exact::collinear(p, q, r)) return;   // precondition of both
    ++n_or4;
    if (del3d::coplanar_orientation(p, q, r, s) != exact::coplanar_orientation(p, q, r, s))
        fail("coplanar_orientation(p,q,r,s)");
    ++n_ic;
    if (exact::coplanar_in_circle_raw(p, q, r, s) == 0) ++deg_ic;
    if (del3d::coplanar_in_circle_unperturbed(p, q, r, s) != exact::coplanar_in_circle_raw(p, q, r, s))
        fail("coplanar_in_circle (raw)");
    if (del3d::coplanar_in_circle(p, q, r, s) != exact::coplanar_in_circle(p, q, r, s))
        fail("coplanar_in_circle (perturbed)");
}

int main() {
    // ---- 0. answers that can be read off the geometry ----------------------
    {
        const double a[3] = {0, 0, 0}, b[3] = {1, 0, 0}, c[3] = {0, 1, 0};
        const double centre[3] = {0.5, 0.5, 0.0};    // the circumcentre: inside
        const double far[3]    = {9.0, 9.0, 0.0};    // plainly outside
        const double online[3] = {2.0, 0.0, 0.0};    // collinear with a, b

        expect(exact::collinear(a, b, online), "reference collinear");
        expect(!exact::collinear(a, b, c), "reference non-collinear");
        expect(exact::coplanar_orientation(a, b, c) == 1, "reference orientation ccw");
        expect(exact::coplanar_orientation(a, c, b) == -1, "reference orientation cw");
        expect(exact::coplanar_in_circle_raw(a, b, c, centre) == 1, "reference in-circle inside");
        expect(exact::coplanar_in_circle_raw(a, b, c, far) == -1, "reference in-circle outside");
        // ... and the same answers, in the opposite orientation and on a plane
        // that no coordinate projection leaves axis-aligned.
        expect(exact::coplanar_in_circle_raw(a, c, b, centre) == 1, "in-circle is orientation-free");
        const double ta[3] = {0, 0, 0}, tb[3] = {1, 0, 1}, tc[3] = {0, 1, 1};
        const double tin[3] = {0.5, 0.5, 1.0};
        expect(exact::coplanar_in_circle_raw(ta, tb, tc, tin) == 1, "reference in-circle, tilted");

        expect(del3d::collinear(a, b, online), "del3d collinear");
        expect(del3d::coplanar_in_circle(a, b, c, centre) == del3d::POSITIVE, "del3d in-circle inside");
        expect(del3d::coplanar_in_circle(a, b, c, far) == del3d::NEGATIVE, "del3d in-circle outside");
        std::printf("stage 0 (sign conventions): %s\n", failures ? "FAILED" : "ok");
    }

    // ---- 1. an exact integer grid in the z = 0 plane -----------------------
    {
        std::vector<std::array<double, 3>> g;
        for (int i = -3; i <= 3; ++i)
            for (int j = -3; j <= 3; ++j) g.push_back({double(i), double(j), 0.0});
        pointsets::Rng rng(11);
        for (int t = 0; t < 20000; ++t) {
            const double* a = g[rng.below(g.size())].data();
            const double* b = g[rng.below(g.size())].data();
            const double* c = g[rng.below(g.size())].data();
            const double* d = g[rng.below(g.size())].data();
            check_triple(a, b, c);
            check_quad(a, b, c, d);
        }
        std::printf("z=0 integer grid: %lld orient3 (%lld collinear), %lld in-circle (%lld cocircular)\n",
                    n_or3, deg_or3, n_ic, deg_ic);
    }

    // ---- 2. a tilted plane, so no coordinate projection is degenerate ------
    {
        const long long o0 = n_or3, i0 = n_ic, d0 = deg_or3, e0 = deg_ic;
        std::vector<std::array<double, 3>> g;
        for (int i = -3; i <= 3; ++i)
            for (int j = -3; j <= 3; ++j)
                g.push_back({double(i), double(j), double(i) + double(j)});  // z = x + y
        pointsets::Rng rng(13);
        for (int t = 0; t < 20000; ++t) {
            const double* a = g[rng.below(g.size())].data();
            const double* b = g[rng.below(g.size())].data();
            const double* c = g[rng.below(g.size())].data();
            const double* d = g[rng.below(g.size())].data();
            check_triple(a, b, c);
            check_quad(a, b, c, d);
        }
        std::printf("tilted plane z=x+y: %lld orient3 (%lld collinear), %lld in-circle (%lld cocircular)\n",
                    n_or3 - o0, deg_or3 - d0, n_ic - i0, deg_ic - e0);
    }

    // ---- 3. points on a circle: maximally cocircular -----------------------
    {
        const long long i0 = n_ic, e0 = deg_ic;
        std::vector<std::array<double, 3>> c;
        // Exact rational points on x^2+y^2 = 25.
        const int pts[][2] = { {5,0},{4,3},{3,4},{0,5},{-3,4},{-4,3},{-5,0},
                               {-4,-3},{-3,-4},{0,-5},{3,-4},{4,-3} };
        for (const auto& q : pts) c.push_back({double(q[0]), double(q[1]), 0.0});
        c.push_back({0.0, 0.0, 0.0});          // the centre, strictly inside
        c.push_back({9.0, 9.0, 0.0});          // strictly outside
        for (std::size_t a = 0; a < c.size(); ++a)
        for (std::size_t b = 0; b < c.size(); ++b)
        for (std::size_t d = 0; d < c.size(); ++d)
        for (std::size_t e = 0; e < c.size(); ++e) {
            if (a == b || a == d || b == d) continue;
            check_quad(c[a].data(), c[b].data(), c[d].data(), c[e].data());
        }
        std::printf("cocircular ring: %lld in-circle (%lld exactly cocircular)\n",
                    n_ic - i0, deg_ic - e0);
    }

    std::printf("\n%s  (%lld orient3 + %lld orient4 + %lld in-circle, %d failures)\n",
                failures ? "FAILED" : "ALL COPLANAR PREDICATES MATCH THE REFERENCE",
                n_or3, n_or4, n_ic, failures);
    return failures ? 1 : 0;
}
