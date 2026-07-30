// del3d viewer - a small mesh reader.
//
// Only the viewer needs this. del3d triangulates a point set and has no opinion
// about where the points came from, so nothing in the library or the tests
// reads a file; keeping the reader here is what lets that stay true.
//
// The formats are the ones an exported mesh actually arrives in. Faces are read
// where the format has them - not because the triangulation needs them, but so
// the input surface can be drawn alongside its tetrahedralisation.
#ifndef DEL3D_VIEWER_MESH_IO_H
#define DEL3D_VIEWER_MESH_IO_H

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace meshio {

struct Mesh {
    std::vector<std::array<double, 3>> vertices;
    std::vector<std::array<int, 3>>    faces;    // empty for a plain point set
    bool empty() const { return vertices.empty(); }
};

/**
 * Loads a mesh, choosing the reader by file extension:
 *
 *   .obj          Wavefront OBJ, ASCII. Polygons are triangulated as a fan;
 *                 v/vt/vn face syntax and negative indices are handled.
 *   .off          OFF / COFF / NOFF.
 *   .ply          ASCII and binary-little-endian PLY.
 *   .stl          Binary and ASCII STL. Welded, since STL repeats every corner.
 *   .xyz .txt .pts  One "x y z" per line, no faces.
 *
 * On failure returns an empty mesh and sets `error`.
 */
Mesh load(const std::string& path, std::string& error);

/// Merges bitwise-identical vertices and rewrites the faces to match. Exact
/// equality is the right test here: the point is to undo a format that stored
/// the same coordinate twice, not to repair a mesh.
void weld(Mesh& m);

/**
 * Removes vertices whose coordinates are not finite, and any face touching one,
 * remapping the rest. Returns how many went.
 *
 * Exporters really do write NaN and inf. del3d requires finite input, and a
 * single one of them would also take a renderer's bounding box with it, so they
 * have to go - but visibly, which is why the count comes back rather than being
 * swallowed.
 */
std::size_t drop_non_finite(Mesh& m);

/// The formats above, as a human-readable list for the UI.
const char* supported_formats();

} // namespace meshio

#endif // DEL3D_VIEWER_MESH_IO_H
