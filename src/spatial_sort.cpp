// del3d - the spatial sort. See include/del3d/spatial_sort.h.
#include "del3d/spatial_sort.h"
#include "rand48.h"

#include <algorithm>
#include <cstddef>

namespace del3d {
namespace {

using detail::Rand48;
using detail::ShuffleRng;

using Iter = std::vector<int>::iterator;

/// Ranges of at most 8 points are left in whatever order they arrived in.
constexpr std::ptrdiff_t kHilbertLimit      = 8;
/// Below 64 points the multiscale pass stops peeling off a coarse subsample.
constexpr std::ptrdiff_t kMultiscaleThresh  = 64;
/// Fraction of a range peeled off as the coarser scale.
constexpr double         kMultiscaleRatio   = 0.125;

/// Compares two point indices on one coordinate axis; `up` reverses the
/// comparison, which is how the Hilbert traversal flips direction per octant.
struct Cmp {
    const std::vector<Point3>* pts;
    int axis;
    bool up;
    double coord(int i) const {
        const Point3& p = (*pts)[i];
        return axis == 0 ? p.x : (axis == 1 ? p.y : p.z);
    }
    bool operator()(int a, int b) const {
        return up ? (coord(b) < coord(a)) : (coord(a) < coord(b));
    }
};

/// Lomuto partition around the element at `pivot_it`, returning the pivot's
/// final position. `right` addresses the LAST element, not one past it.
Iter partition_deterministic(Iter left, Iter right, Iter pivot_it, Cmp& cmp) {
    std::iter_swap(pivot_it, right);
    Iter result = left;
    for (Iter it = left; it != right; ++it) {
        if (cmp(*it, *right)) {
            std::iter_swap(result, it);
            ++result;
        }
    }
    std::iter_swap(right, result);
    return result;
}

/// Quickselect: places the element that belongs at `nth` there, with smaller
/// elements before it, always pivoting on the middle element.
///
/// std::nth_element is NOT interchangeable here. It leaves the relative order
/// of elements comparing equal unspecified, and on input with repeated
/// coordinates - grids, extruded meshes, quantised scans - that difference
/// changes the permutation, hence the insertion order, hence the cell order.
void nth_element_deterministic(Iter left, Iter nth, Iter right, Cmp& cmp) {
    if (left == right) return;
    --right;                       // now addresses the last element
    if (left == right) return;
    for (;;) {
        Iter pivot_it = left + ((right - left) / 2);
        Iter new_pivot = partition_deterministic(left, right, pivot_it, cmp);
        if (new_pivot == nth) return;
        if (nth < new_pivot) right = new_pivot - 1;
        else                 left  = new_pivot + 1;
    }
}

/// Splits a range at its median along one axis and returns the split point.
Iter hilbert_split(Iter begin, Iter end, Cmp cmp) {
    if (begin >= end) return begin;
    Iter middle = begin + (end - begin) / 2;
    nth_element_deterministic(begin, middle, end, cmp);
    return middle;
}

/**
 * One level of the Hilbert-curve sort.
 *
 * Splits the range into eight octants by three successive median cuts - first
 * on axis x, then on y within each half, then on z within each quarter - and
 * recurses into each octant with the axis roles rotated and the directions
 * flipped so that consecutive octants are visited end to end. That is what
 * makes the resulting order follow a Hilbert curve, and therefore makes
 * consecutive points spatially close, which is what the incremental
 * construction's point-location walk exploits.
 *
 * `x` names the current primary axis and upx/upy/upz the traversal direction
 * along each; the recursive calls permute them per octant.
 */
void recursive_sort(const std::vector<Point3>& pts, Iter begin, Iter end,
                    int x, bool upx, bool upy, bool upz) {
    const int y = (x + 1) % 3, z = (x + 2) % 3;
    if (end - begin <= kHilbertLimit) return;

    Iter m0 = begin, m8 = end;
    Iter m4 = hilbert_split(m0, m8, Cmp{&pts, x,  upx});
    Iter m2 = hilbert_split(m0, m4, Cmp{&pts, y,  upy});
    Iter m1 = hilbert_split(m0, m2, Cmp{&pts, z,  upz});
    Iter m3 = hilbert_split(m2, m4, Cmp{&pts, z, !upz});
    Iter m6 = hilbert_split(m4, m8, Cmp{&pts, y, !upy});
    Iter m5 = hilbert_split(m4, m6, Cmp{&pts, z,  upz});
    Iter m7 = hilbert_split(m6, m8, Cmp{&pts, z, !upz});

    recursive_sort(pts, m0, m1, z,  upz,  upx,  upy);
    recursive_sort(pts, m1, m2, y,  upy,  upz,  upx);
    recursive_sort(pts, m2, m3, y,  upy,  upz,  upx);
    recursive_sort(pts, m3, m4, x,  upx, !upy, !upz);
    recursive_sort(pts, m4, m5, x,  upx, !upy, !upz);
    recursive_sort(pts, m5, m6, y, !upy,  upz, !upx);
    recursive_sort(pts, m6, m7, y, !upy,  upz, !upx);
    recursive_sort(pts, m7, m8, z, !upz, !upx,  upy);
}

void hilbert_sort(const std::vector<Point3>& pts, Iter begin, Iter end) {
    recursive_sort(pts, begin, end, 0, false, false, false);
}

/**
 * Multiscale pass: peel the leading 12.5% off the range, order that recursively
 * first, and Hilbert-sort the remainder after it. The triangulation is then
 * built from a coarse subsample outwards, which keeps the point-location walk
 * short - inserting a Hilbert-sorted sequence alone would make every walk
 * traverse a long thin chain of recently created cells.
 */
void multiscale_sort(const std::vector<Point3>& pts, Iter begin, Iter end) {
    Iter middle = begin;
    if (end - begin >= kMultiscaleThresh) {
        middle = begin + static_cast<std::ptrdiff_t>(double(end - begin) * kMultiscaleRatio);
        multiscale_sort(pts, begin, middle);
    }
    hilbert_sort(pts, middle, end);
}

} // namespace

void spatial_sort(const std::vector<Point3>& points, std::vector<int>& order) {
    if (order.size() < 2) return;
    // The shuffle removes any structure in the caller's ordering (a mesh's
    // vertex list is usually already spatially coherent, which would defeat the
    // multiscale pass) and does so from a fixed seed, so the result is still a
    // deterministic function of the input.
    Rand48 rng;
    ShuffleRng shuffle_rng(rng);
    detail::random_shuffle(order.begin(), order.end(), shuffle_rng);
    multiscale_sort(points, order.begin(), order.end());
}

} // namespace del3d
