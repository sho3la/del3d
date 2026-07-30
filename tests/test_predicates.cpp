// del3d - 3D predicate conformance test.
//
// Checks del3d::orient3d and in_sphere - raw and perturbed - against the
// independent exact-integer reference in tests/exact.h, on
//   (a) configurations whose answer is known by hand, which also pin the sign
//       conventions of the reference itself,
//   (b) hand-made configurations that are exactly degenerate by construction,
//       and
//   (c) random 5-tuples of quantised surface samples, which is where the
//       cospherical cases that matter in practice actually occur.
//
// The counts of exactly degenerate cases printed at the end are the point of
// the test: those are the configurations where an unfiltered predicate silently
// returns the wrong sign.
#include "del3d/predicates.h"

#include "exact.h"
#include "pointsets.h"

#include <array>
#include <cstdio>
#include <vector>

static int failures = 0;
static long long checked_o = 0, checked_s = 0, degenerate_o = 0, degenerate_s = 0;

static void check_orient(const double* a, const double* b, const double* c, const double* d) {
    const int mine = del3d::orient3d(a, b, c, d);
    const int ref  = exact::orient3d(a, b, c, d);
    ++checked_o;
    if (ref == 0) ++degenerate_o;
    if (mine != ref) {
        if (++failures < 10)
            std::printf("  ORIENT MISMATCH del3d=%d ref=%d\n", mine, ref);
    }
}

static void check_sphere(const double* a, const double* b, const double* c,
                         const double* d, const double* e) {
    // The predicate's precondition.
    if (del3d::orient3d(a, b, c, d) != del3d::POSITIVE) return;
    ++checked_s;
    if (exact::in_sphere_raw(a, b, c, d, e) == 0) ++degenerate_s;

    if (del3d::in_sphere_unperturbed(a, b, c, d, e) != exact::in_sphere_raw(a, b, c, d, e)) {
        if (++failures < 10) std::printf("  INSPHERE(raw) MISMATCH\n");
    }
    if (del3d::in_sphere(a, b, c, d, e) != exact::in_sphere(a, b, c, d, e)) {
        if (++failures < 10) std::printf("  INSPHERE(perturbed) MISMATCH\n");
    }
}

static void expect(bool ok, const char* what) {
    if (!ok) { ++failures; std::printf("  WRONG: %s\n", what); }
}

int main() {
    // ---- 0. the reference's own sign conventions ---------------------------
    // Both implementations could in principle agree on a convention that is not
    // the documented one, so anchor it on answers that can be read off the
    // geometry: the unit corner tetrahedron is positively oriented, its
    // circumsphere is centred at (0.5,0.5,0.5), and the two query points below
    // are plainly inside and outside it.
    {
        const double o[3] = {0, 0, 0}, x[3] = {1, 0, 0}, y[3] = {0, 1, 0}, z[3] = {0, 0, 1};
        const double inside[3] = {0.25, 0.25, 0.25};
        const double outside[3] = {5.0, 5.0, 5.0};
        expect(exact::orient3d(o, x, y, z) == 1, "reference orient3d of the unit corner");
        expect(exact::orient3d(o, y, x, z) == -1, "reference orient3d under a swap");
        expect(exact::in_sphere_raw(o, x, y, z, inside) == 1, "reference in_sphere inside");
        expect(exact::in_sphere_raw(o, x, y, z, outside) == -1, "reference in_sphere outside");
        expect(del3d::orient3d(o, x, y, z) == del3d::POSITIVE, "del3d orient3d of the unit corner");
        expect(del3d::in_sphere(o, x, y, z, inside) == del3d::POSITIVE, "del3d in_sphere inside");
        expect(del3d::in_sphere(o, x, y, z, outside) == del3d::NEGATIVE, "del3d in_sphere outside");
        std::printf("stage 0 (sign conventions): %s\n", failures ? "FAILED" : "ok");
    }

    // ---- 1. exactly-degenerate configurations by construction --------------
    {
        const double o[3] = {0, 0, 0}, x[3] = {1, 0, 0}, y[3] = {0, 1, 0};
        const double coplanar[3] = {1, 1, 0};                 // in the z=0 plane
        check_orient(o, x, y, coplanar);                      // must be ZERO
        const double up[3] = {0, 0, 1};
        check_orient(o, x, y, up);

        // Four corners of a unit square on a circle, plus a fifth cocircular
        // point: the classic cospherical degeneracy.
        const double a[3] = {1, 0, 0}, b[3] = {0, 1, 0}, c[3] = {-1, 0, 0};
        const double d[3] = {0, 0, 1}, e[3] = {0, -1, 0};
        check_sphere(a, b, c, d, e);
        check_sphere(a, b, d, c, e);
        // Lattice points: many exact ties.
        for (int i = -2; i <= 2; ++i)
        for (int j = -2; j <= 2; ++j)
        for (int k = -2; k <= 2; ++k) {
            const double p[3] = {double(i), double(j), double(k)};
            const double q[3] = {double(j), double(k), double(i)};
            check_orient(o, x, y, p);
            check_sphere(a, b, d, p, q);
        }
        std::printf("stage 1 (constructed degeneracies): %lld orient, %lld insphere, "
                    "%lld/%lld exactly degenerate\n",
                    checked_o, checked_s, degenerate_o, degenerate_s);
    }

    // ---- 2. quantised surface samples and lattices -------------------------
    struct Case { const char* name; std::vector<pointsets::Pt> pts; };
    const Case cases[] = {
        // A coarse quantisation - eighths - so that ties are the common case
        // rather than the rare one.
        { "cylinder (quantised)", pointsets::cylinder(400, 8) },
        { "sphere (quantised)",   pointsets::sphere(400, 8) },
        { "torus (quantised)",    pointsets::torus(400, 8) },
        { "box shell 6^3",        pointsets::box_shell(6) },
    };
    for (const Case& c : cases) {
        const long long o0 = checked_o, s0 = checked_s, d0 = degenerate_o, e0 = degenerate_s;
        pointsets::Rng rng(12345);
        for (int t = 0; t < 20000; ++t) {
            const double* p[5];
            for (int k = 0; k < 5; ++k) p[k] = c.pts[rng.below(c.pts.size())].data();
            check_orient(p[0], p[1], p[2], p[3]);
            check_sphere(p[0], p[1], p[2], p[3], p[4]);
        }
        std::printf("  %-22s %4zu pts  %6lld orient (%5lld exactly coplanar), "
                    "%6lld insphere (%5lld exactly cospherical)\n",
                    c.name, c.pts.size(), checked_o - o0, degenerate_o - d0,
                    checked_s - s0, degenerate_s - e0);
    }

    std::printf("\n%s  (%lld orient + %lld insphere checks, %d failures)\n",
                failures ? "FAILED" : "ALL PREDICATES MATCH THE REFERENCE",
                checked_o, checked_s, failures);
    return failures ? 1 : 0;
}
