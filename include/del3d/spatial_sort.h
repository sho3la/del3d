// del3d - the insertion order used by the incremental construction.
//
// This is not purely a performance detail. On degenerate input the Delaunay
// triangulation is not unique, and even when the set of tetrahedra is unique
// the order in which cells are *allocated* depends on the order points are
// inserted. del3d treats both as part of its output contract, so the order is
// fixed by an algorithm that is fully specified here rather than by anything
// the standard library leaves unspecified.
//
// The sequence is:
//
//   1. shuffle the indices with a fixed-seed 48-bit linear congruential
//      generator (src/rand48.h), using a deterministic Fisher-Yates loop;
//   2. multiscale pass: while the range holds at least 64 elements, split off
//      the leading 12.5% and sort that recursively first, so that the
//      triangulation grows from a coarse subsample outwards;
//   3. Hilbert sort each range: recursive median splits along a Hilbert curve
//      traversal, stopping at ranges of 8 or fewer points.
//
// The median split uses an explicitly written selection routine rather than
// std::nth_element, whose relative order among equal elements is unspecified;
// that difference is observable on any input with repeated coordinates.
#ifndef DEL3D_SPATIAL_SORT_H
#define DEL3D_SPATIAL_SORT_H

#include <cstddef>
#include <vector>

namespace del3d {

struct Point3 { double x, y, z; };

/**
 * Reorders `order` - indices into `points` - into the spatial-sort sequence.
 * `order` is normally 0..n-1 on entry. `points` itself is not modified.
 */
void spatial_sort(const std::vector<Point3>& points, std::vector<int>& order);

} // namespace del3d

#endif // DEL3D_SPATIAL_SORT_H
