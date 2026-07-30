# del3d

3D Delaunay triangulation of a point set. One job, no dependencies: standard
C++17 and nothing else, tests included.

```cpp
#include <del3d/delaunay.h>

std::vector<del3d::Point3> pts = /* ... */;

del3d::Delaunay d;
d.build(pts);

for (const del3d::Delaunay::Tet& t : d.finite_cells()) {
    // t[0..3] are indices into pts
    del3d::Point3 dual = del3d::circumcenter(pts[t[0]], pts[t[1]],
                                             pts[t[2]], pts[t[3]]);
}
```

## What it guarantees

* **Exact predicates.** Orientation and in-sphere are evaluated with a filtered
  floating-point stage backed by exact expansion arithmetic, so the sign is
  always the true sign, whatever the input coordinates.
* **Total predicates.** Cospherical and cocircular configurations are resolved
  by symbolic perturbation rather than reported as ties, so the triangulation is
  well defined on degenerate input - grids, extruded geometry, quantised scans -
  and not merely self-consistent.
* **Determinism.** The output is a function of the input point sequence alone:
  the same tetrahedra *and the same cell order*, on every compiler and platform.
  Nothing whose behaviour the standard leaves unspecified is on the path from
  points to cells.

That last point is the unusual one and it shapes the whole design. Callers
frequently break ties by "whichever incident cell comes first" - the Voronoi
pole of a vertex, for instance, is picked from its incident cells - so cell
order is observable, and therefore part of the contract rather than an
implementation detail.

## Exactness requires strict floating point

`src/expansion.h` recovers the round-off of each operation by recomputing it
from the rounded result. That is exact only if every operation is correctly
rounded and nothing reassociates it, contracts it into a fused multiply-add, or
evaluates it in a wider register. `src/interval.h` additionally changes the FPU
rounding mode at run time.

`CMakeLists.txt` therefore sets `/fp:strict` (MSVC) or `-ffp-contract=off
-frounding-math` (GCC/Clang, plus SSE2 on 32-bit x86). **Do not build del3d with
`-ffast-math`.**

## Building

```sh
cmake -S . -B build
cmake --build build
```

Produces a static library `del3d` exporting `include/`. `DEL3D_BUILD_TESTS=ON`
adds the tests.

Two optional diagnostics, compiled out by default:

| define | effect |
|---|---|
| `DEL3D_VALIDATE` | runs `Tds::validate()` after every insertion and aborts on the first broken invariant, naming the insertion that broke it |
| `DEL3D_TRACE` / `DEL3D_WALK` | dumps the cell table after each insertion / traces the location walk |

## Tests

The tests need nothing beyond the standard library either - no data files, no
third-party package, nothing to install.

The predicates are checked against `tests/exact.h`, an arbitrary-precision
*integer* reference written for the purpose. It shares no code with the library:
every finite double is exactly `m * 2^e`, so each predicate's coordinates are
converted to integers on a common scale - which changes the determinant only by
a positive factor, and so not its sign - and the determinant is then evaluated
with no rounding and no filtering anywhere in the decision path. It is slow and
obvious, which is what a reference should be.

| test | what it checks |
|---|---|
| `test_container` | the allocation order the rest of the design depends on: block sizes, front-to-back handout within a block, LIFO reuse, memory-order iteration |
| `test_predicates` | `orient3d` and `in_sphere`, raw and perturbed, against the exact reference, on constructed degeneracies and on random 5-tuples of quantised surface samples |
| `test_predicates2` | the coplanar predicates, on an axis-aligned integer grid, a tilted grid (so no single projection is the degenerate one), and a ring of exactly cocircular points |
| `test_spatial_sort` | the order is a permutation, is reproducible, is spatially local, and hashes to a pinned value |
| `test_delaunay` | orientation, the global empty-sphere property, closure of the cell complex, vertex coverage, `incident_cells()`, determinism, and the pinned cell sequence |

`test_delaunay` verifies the output rather than comparing it: for each cell it
asks the exact reference whether any input point lies strictly inside its
circumsphere, checks that every triangular facet is shared by exactly two cells
once the infinite ones are counted, and checks that every distinct input point
is a vertex. Determinism is pinned by hashing the cell sequence and comparing
against a stored value, so a platform that produces the same tetrahedra in a
different order still fails. `test_spatial_sort` pins its permutation the same
way. Both generate their point sets from `std::mt19937` with the mapping onto
doubles written out explicitly, since the standard's distribution adaptors are
not specified bit for bit and pinned values would otherwise be meaningless.

Status: all five pass with zero failures.

The degenerate counts the tests print are the point of them: in one run
`test_predicates` covered 1 281 exactly coplanar orientations and 315 exactly
cospherical in-sphere tests on the lattice-quantised box shell alone, and
`test_predicates2` covered 18 324 exactly cocircular configurations. Those are
the cases where an unfiltered predicate silently returns the wrong sign and the
triangulation quietly stops being Delaunay.

The input the tests lean on is deliberately quantised - snapped to a power-of-two
grid - rather than uniformly random. Random points are in general position and
never reach the exact stage at all; CAD models, lattices and anything exported
from a scanner are full of exactly coincident planes and spheres, and that is
what these predicates exist for.