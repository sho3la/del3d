// del3d - spatial sort test.
//
// The insertion order is part of the library's output contract: it is what
// makes the cell order reproducible. So what is checked here is not that the
// order is "good" but that it is exactly what it is supposed to be -
//
//   * it is a permutation of the input indices, whatever the input;
//   * it is a function of the input alone, so sorting the same points twice
//     gives the same answer;
//   * on a fixed set of inputs it hashes to a pinned value, which is how a
//     platform whose answer differs gets caught;
//
// - plus one property check, that the order really is spatially local, since a
// permutation could satisfy everything above and still be useless.
//
// The sizes below straddle the algorithm's thresholds (the Hilbert recursion
// limit of 8 and the multiscale threshold of 64), and the grid and duplicate
// sets produce masses of exactly-equal coordinates, which is where a median
// selection with unspecified order among equal elements would diverge.
#include "del3d/spatial_sort.h"

#include "pointsets.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <numeric>
#include <random>
#include <vector>

static int failures = 0;

static void fail(const char* what, const char* detail) {
    ++failures;
    std::printf("    FAIL: %s (%s)\n", what, detail);
}

static std::vector<int> sorted_order(const std::vector<del3d::Point3>& pts) {
    std::vector<int> order(pts.size());
    std::iota(order.begin(), order.end(), 0);
    del3d::spatial_sort(pts, order);
    return order;
}

static std::uint64_t hash_order(const std::vector<int>& order) {
    std::uint64_t h = 1469598103934665603ull;                 // FNV-1a, 64 bit
    for (std::uint64_t v : std::vector<std::uint64_t>(order.begin(), order.end()))
        for (int b = 0; b < 4; ++b) { h ^= (v >> (8 * b)) & 0xffull; h *= 1099511628211ull; }
    return h;
}

/// Mean distance between consecutive points in a given order.
static double mean_step(const std::vector<del3d::Point3>& pts, const std::vector<int>& order) {
    if (order.size() < 2) return 0.0;
    double sum = 0.0;
    for (std::size_t i = 1; i < order.size(); ++i) {
        const del3d::Point3& a = pts[order[i - 1]];
        const del3d::Point3& b = pts[order[i]];
        const double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
        sum += std::sqrt(dx * dx + dy * dy + dz * dz);
    }
    return sum / double(order.size() - 1);
}

static void check(const char* what, const std::vector<del3d::Point3>& pts,
                  std::uint64_t expect_hash, bool check_locality) {
    const std::vector<int> order = sorted_order(pts);

    // A permutation of 0..n-1.
    {
        std::vector<int> seen = order;
        std::sort(seen.begin(), seen.end());
        std::vector<int> want(pts.size());
        std::iota(want.begin(), want.end(), 0);
        if (seen != want) fail(what, "not a permutation of the input indices");
    }

    // A function of the points alone.
    if (sorted_order(pts) != order) fail(what, "not reproducible");

    // Actually local: consecutive points must be much closer together than a
    // random pair, or the sort is not doing its job. The bound is deliberately
    // loose, and only applied where the sample is large enough for the mean to
    // mean anything - this is a sanity check, not a quality metric.
    double ratio = 0.0;
    if (check_locality && pts.size() >= 200) {
        std::vector<int> rnd(pts.size());
        std::iota(rnd.begin(), rnd.end(), 0);
        std::mt19937 rng(5);
        std::shuffle(rnd.begin(), rnd.end(), rng);
        const double base = mean_step(pts, rnd);
        ratio = base > 0.0 ? mean_step(pts, order) / base : 0.0;
        if (base > 0.0 && ratio > 0.5) fail(what, "the order is not spatially local");
    }

    const std::uint64_t h = hash_order(order);
    if (expect_hash != 0 && h != expect_hash) fail(what, "differs from the pinned order");

    std::printf("  %-24s %6zu pts  hash %016llx  locality %.3f%s\n",
                what, pts.size(), static_cast<unsigned long long>(h), ratio,
                expect_hash == 0 ? "  (not pinned)" : "");
}

static std::vector<del3d::Point3> to_points(const std::vector<pointsets::Pt>& v) {
    std::vector<del3d::Point3> out;
    out.reserve(v.size());
    for (const pointsets::Pt& p : v) out.push_back({p[0], p[1], p[2]});
    return out;
}

int main() {
    std::printf("del3d::spatial_sort\n");

    const struct { int n; std::uint64_t hash; } gen[] = {
        {    0, 0x14650fb0739d0383ull },
        {    1, 0x315446a086a23133ull },
        {    2, 0x29034675a49f07c2ull },
        {    5, 0x744ddb73fb935d27ull },
        {    8, 0x7ec129e73709cb53ull },
        {    9, 0x580cd772e82efaabull },
        {   63, 0x18cf12b455d4e1acull },
        {   64, 0x4fde96ee0b2e8be3ull },
        {   65, 0xafa28526c6130f43ull },
        {  200, 0x23d15ee1f09fef23ull },
        { 1000, 0x8b7191bde7a7449bull },
    };
    for (const auto& g : gen) {
        char name[64];
        std::snprintf(name, sizeof name, "random n=%d", g.n);
        check(name, to_points(pointsets::random_cube(std::size_t(g.n), 7u)), g.hash, true);
    }

    check("lattice 12x12x6", to_points(pointsets::lattice(12, 12, 6)),
          0xdb779fbe6eb82833ull, true);
    check("box shell 8^3", to_points(pointsets::box_shell(8)),
          0x1667fe1e40586b6aull, true);
    check("300 identical points",
          std::vector<del3d::Point3>(300, del3d::Point3{1.0, 2.0, 3.0}),
          0xc9325f52533ad94full, false);

    check("cylinder (quantised)", to_points(pointsets::cylinder(2000)),
          0x0b3b773499c64173ull, true);
    check("sphere (quantised)", to_points(pointsets::sphere(2000)),
          0x93fb3414e085d46bull, true);
    check("torus (quantised)", to_points(pointsets::torus(2000)),
          0x0e117e42cb5f457full, true);

    std::printf("\n%s (%d failures)\n",
                failures ? "FAILED" : "THE ORDER IS THE PINNED ONE", failures);
    return failures ? 1 : 0;
}
