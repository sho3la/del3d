// del3d - a viewer for the triangulation, built on Polyscope.
//
// Not part of the library and not built by default: Polyscope pulls in a window
// system, OpenGL, ImGui and GLM, none of which del3d itself has any use for.
// Turn it on with -DDEL3D_BUILD_VIEWER=ON.
//
// It shows the tetrahedra, the input points, and the Voronoi vertices dual to
// the cells, and colours the cells by two measures that make the degenerate
// output visible: circumradius, and the radius-edge ratio, which is the usual
// way to spot slivers.
//
// A tetrahedralisation is opaque from outside - the interesting cells are the
// ones you cannot see - so there are three separate ways to look inside, and
// they answer different questions:
//
//   * edges on the volume mesh: the cell boundaries on the visible hull only;
//   * the "tet edges" curve network: every edge of every cell, interior ones
//     included, which with the volume mesh switched off is the true wireframe;
//   * a slice plane, optionally cutting through cells rather than dropping
//     whole ones, which shows the interior as a solid cross-section.
//
// Transparency is the fourth, and combines with all of them.
//
// Run with no arguments for the built-in point sets, or pass a mesh - OBJ, OFF,
// PLY, STL or a plain point list, see mesh_io.h - which can also be loaded from
// the UI once the window is open. Only the vertices reach the triangulation;
// the faces are drawn beside it so the input surface and the tetrahedra of its
// vertices can be compared. `--check` builds and registers everything against
// Polyscope's mock backend and exits without opening a window, which is how the
// viewer can be smoke-tested on a machine with no display.
#include <del3d/delaunay.h>

#include "mesh_io.h"

#include <polyscope/curve_network.h>
#include <polyscope/point_cloud.h>
#include <polyscope/polyscope.h>
#include <polyscope/slice_plane.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/view.h>
#include <polyscope/volume_mesh.h>

#include <imgui.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <random>
#include <tuple>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace {

using Pt = std::array<double, 3>;

/// Above this many cells the interior wireframe is built on demand rather than
/// with every rebuild.
constexpr std::size_t kWireframeAutoLimit = 150000;

/// Voronoi vertices further than this many model diagonals from the input are
/// left out of the point cloud. They are real - a sliver's dual belongs far
/// away - but drawing them is useless and lets one cell set the scene's scale.
constexpr double kDualRadiusLimit = 10.0;

// --------------------------------------------------------------------------
// Input.
// --------------------------------------------------------------------------

/// Rounds to a multiple of 1/q, q a power of two, so the result is exact - which
/// is what makes the quantised sets genuinely, not nearly, degenerate.
double snap(double v, int q) { return std::floor(v * q + 0.5) / q; }

std::vector<del3d::Point3> make_points(int kind, int n, unsigned seed) {
    std::vector<del3d::Point3> pts;
    std::mt19937 rng(seed);
    auto u01 = [&rng]() { return double(rng()) * (1.0 / 4294967296.0); };
    const double two_pi = 6.283185307179586;

    switch (kind) {
        case 0:                                   // random in a cube
            for (int i = 0; i < n; ++i)
                pts.push_back({2.0 * u01() - 1.0, 2.0 * u01() - 1.0, 2.0 * u01() - 1.0});
            break;
        case 1: {                                 // integer lattice
            const int m = std::max(2, int(std::lround(std::cbrt(double(n)))));
            for (int i = 0; i < m; ++i)
                for (int j = 0; j < m; ++j)
                    for (int k = 0; k < m; ++k)
                        pts.push_back({double(i) / m, double(j) / m, double(k) / m});
            break;
        }
        case 2:                                   // quantised sphere
            for (int i = 0; i < n; ++i) {
                const double z = 2.0 * u01() - 1.0, r = std::sqrt(1.0 - z * z);
                const double t = two_pi * u01();
                pts.push_back({snap(r * std::cos(t), 16), snap(r * std::sin(t), 16),
                               snap(z, 16)});
            }
            break;
        default:                                  // quantised torus
            for (int i = 0; i < n; ++i) {
                const double a = two_pi * u01(), b = two_pi * u01();
                const double R = 1.0 + 0.35 * std::cos(b);
                pts.push_back({snap(R * std::cos(a), 16), snap(R * std::sin(a), 16),
                               snap(0.35 * std::sin(b), 16)});
            }
            break;
    }
    return pts;
}

// --------------------------------------------------------------------------
// State, and the registration that follows a rebuild.
// --------------------------------------------------------------------------

struct State {
    std::vector<del3d::Point3> points;
    del3d::Delaunay            tri;
    int          kind = 0;
    int          count = 2000;
    unsigned     seed = 1;

    // An imported mesh. Its faces are kept only so the input surface can be
    // drawn next to the tetrahedralisation of its vertices.
    bool         from_file = false;
    std::string  file;
    char         path_field[512] = {0};
    meshio::Mesh mesh;
    std::string  load_error;
    bool         show_surface = true;
    std::size_t  dropped_vertices = 0;   // non-finite coordinates in the file
    std::size_t  duals_dropped = 0;      // circumcentres too far away to draw

    // How the interior is being looked at. Held here rather than read back off
    // the structures because rebuild() destroys and recreates them.
    bool  edges = true;          // cell edges on the volume mesh
    float edge_width = 1.0f;
    float transparency = 1.0f;   // 1 = opaque
    bool  solid = true;          // the volume mesh itself
    bool  wireframe = false;     // the all-edges curve network
    bool  cut_cells = true;      // slice planes cut cells open, not drop them
    polyscope::SlicePlane* plane = nullptr;
};

State g;

/// Every distinct edge of every cell - the interior ones too, which is what
/// separates this from the volume mesh's own edge display.
std::vector<std::array<int, 2>> tet_edges(const std::vector<del3d::Delaunay::Tet>& cells) {
    std::set<std::array<int, 2>> uniq;
    for (const del3d::Delaunay::Tet& t : cells)
        for (int i = 0; i < 4; ++i)
            for (int j = i + 1; j < 4; ++j)
                uniq.insert(t[i] < t[j] ? std::array<int, 2>{t[i], t[j]}
                                        : std::array<int, 2>{t[j], t[i]});
    return std::vector<std::array<int, 2>>(uniq.begin(), uniq.end());
}

/// Pushes the interior-viewing settings onto whatever is currently registered.
void apply_style() {
    if (polyscope::hasVolumeMesh("triangulation")) {
        polyscope::VolumeMesh* vm = polyscope::getVolumeMesh("triangulation");
        vm->setEnabled(g.solid);
        vm->setEdgeWidth(g.edges ? double(g.edge_width) : 0.0);
        vm->setTransparency(g.transparency);
        // With whole-element culling off, a slice plane cuts the tetrahedra
        // themselves and you see the cross-section; with it on, cells in front
        // of the plane are dropped entirely and you see whole neighbours.
        vm->setCullWholeElements(!g.cut_cells);
    }
    if (polyscope::hasCurveNetwork("tet edges"))
        polyscope::getCurveNetwork("tet edges")->setEnabled(g.wireframe);
    if (polyscope::hasSurfaceMesh("input surface"))
        polyscope::getSurfaceMesh("input surface")->setEnabled(g.show_surface);
}

double edge_min(const del3d::Point3& a, const del3d::Point3& b,
                const del3d::Point3& c, const del3d::Point3& d) {
    const del3d::Point3* p[4] = { &a, &b, &c, &d };
    double best = 0.0;
    bool first = true;
    for (int i = 0; i < 4; ++i)
        for (int j = i + 1; j < 4; ++j) {
            const double dx = p[i]->x - p[j]->x;
            const double dy = p[i]->y - p[j]->y;
            const double dz = p[i]->z - p[j]->z;
            const double len = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (first || len < best) { best = len; first = false; }
        }
    return best;
}

/**
 * Fixes the scene's extents to the input points.
 *
 * Polyscope otherwise derives them from everything registered, and the derived
 * structures here are not bounded by the input: one hull sliver's circumcentre
 * can sit 1e12 model-widths away, at which point the length scale is set by
 * that dual, the model becomes a sub-pixel speck, and the home view degenerates
 * into "lookAt() yielded an invalid view". The input points are the only honest
 * definition of how big the scene is, so they are what gets used.
 */
void set_scene_extents(const Pt& lo, const Pt& hi, double diag, bool reset_view) {
    polyscope::options::automaticallyComputeSceneExtents = false;

    // A single point, or every point identical, has no extent to speak of.
    const double scale = (diag > 0.0 && std::isfinite(diag)) ? diag : 1.0;
    polyscope::state::lengthScale = float(scale);

    Pt c{ 0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]), 0.5 * (lo[2] + hi[2]) };
    for (int k = 0; k < 3; ++k) if (!std::isfinite(c[k])) c[k] = 0.0;
    const float h = float(0.5 * scale);
    polyscope::state::boundingBox = std::make_tuple(
        glm::vec3{float(c[0]) - h, float(c[1]) - h, float(c[2]) - h},
        glm::vec3{float(c[0]) + h, float(c[1]) + h, float(c[2]) + h});

    // Only on a genuine change of input: resetting while someone is dragging a
    // slider would yank the camera out from under them.
    if (reset_view) polyscope::view::resetCameraToHomeView();
}

/// Loads `path` into the state. Leaves the current input untouched on failure,
/// so a mistyped path in the UI does not blank the view.
bool load_file(const std::string& path) {
    std::string err;
    meshio::Mesh m = meshio::load(path, err);
    if (m.empty()) {
        g.load_error = err.empty() ? ("nothing readable in " + path) : err;
        return false;
    }
    // Exporters do emit NaN and inf. del3d requires finite coordinates, and one
    // of them would take the whole scene's bounding box with it, so they go -
    // reported, not silently.
    g.dropped_vertices = meshio::drop_non_finite(m);
    if (m.empty()) {
        g.load_error = "every vertex in " + path + " was non-finite";
        return false;
    }

    g.load_error.clear();
    g.mesh = std::move(m);
    g.file = path;
    g.from_file = true;
    std::snprintf(g.path_field, sizeof g.path_field, "%s", path.c_str());
    return true;
}

void rebuild(bool reset_view) {
    if (g.from_file) {
        // Only the vertices go into the triangulation - del3d takes a point
        // set, and the faces are along for the ride so the surface can be drawn.
        g.points.clear();
        g.points.reserve(g.mesh.vertices.size());
        for (const Pt& v : g.mesh.vertices) g.points.push_back({v[0], v[1], v[2]});
    } else {
        g.points = make_points(g.kind, g.count, g.seed);
    }

    g.tri.build(g.points);
    const std::vector<del3d::Delaunay::Tet>& cells = g.tri.finite_cells();

    std::vector<Pt> verts;
    verts.reserve(g.points.size());
    for (const del3d::Point3& p : g.points) verts.push_back({p.x, p.y, p.z});

    // The extents of the input, which is the only thing that should decide
    // where the camera goes - see set_scene_extents().
    Pt lo{0, 0, 0}, hi{0, 0, 0};
    for (std::size_t i = 0; i < verts.size(); ++i)
        for (int k = 0; k < 3; ++k) {
            if (i == 0) { lo[k] = hi[k] = verts[i][k]; }
            else { lo[k] = std::min(lo[k], verts[i][k]); hi[k] = std::max(hi[k], verts[i][k]); }
        }
    const double diag = std::sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) +
                                  (hi[1] - lo[1]) * (hi[1] - lo[1]) +
                                  (hi[2] - lo[2]) * (hi[2] - lo[2]));

    std::vector<Pt> duals;
    std::vector<double> radius, radius_edge;
    duals.reserve(cells.size());
    radius.reserve(cells.size());
    radius_edge.reserve(cells.size());
    g.duals_dropped = 0;

    for (const del3d::Delaunay::Tet& t : cells) {
        const del3d::Point3& a = g.points[t[0]];
        const del3d::Point3  c = del3d::circumcenter(a, g.points[t[1]],
                                                     g.points[t[2]], g.points[t[3]]);
        const double dx = c.x - a.x, dy = c.y - a.y, dz = c.z - a.z;
        const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
        const double e = edge_min(a, g.points[t[1]], g.points[t[2]], g.points[t[3]]);
        radius.push_back(r);
        radius_edge.push_back(e > 0.0 ? r / e : 0.0);

        // A sliver's circumcentre is legitimately far away - that is the whole
        // reason circumcenter() is evaluated exactly - and a surface mesh is
        // full of slivers, so the dual of a hull cell can sit many orders of
        // magnitude outside the model, or overflow outright. Drawing those is
        // pointless and lets one cell dictate the scale of the entire scene, so
        // the cloud keeps only the duals that are near the input.
        const bool finite = std::isfinite(c.x) && std::isfinite(c.y) && std::isfinite(c.z);
        if (finite && (diag <= 0.0 || r <= kDualRadiusLimit * diag))
            duals.push_back({c.x, c.y, c.z});
        else
            ++g.duals_dropped;
    }

    polyscope::removeAllStructures();

    polyscope::PointCloud* pc = polyscope::registerPointCloud("input points", verts);
    pc->setPointRadius(0.003);

    if (g.from_file && !g.mesh.faces.empty())
        polyscope::registerSurfaceMesh("input surface", g.mesh.vertices, g.mesh.faces);

    if (!cells.empty()) {
        polyscope::VolumeMesh* vm = polyscope::registerTetMesh("triangulation", verts, cells);
        vm->addCellScalarQuantity("circumradius", radius);
        // The sliver measure: a well-shaped tetrahedron sits near 0.7, and the
        // value grows without bound as the cell flattens.
        vm->addCellScalarQuantity("radius / shortest edge", radius_edge)->setEnabled(true);

        // The interior wireframe, disabled by default: it is a lot of lines,
        // and only worth drawing once the solid cells are turned down or off.
        // An imported mesh can be large enough that collecting the edges is
        // itself the slow part, so past a threshold it is built only when
        // actually asked for.
        if (g.wireframe || cells.size() <= kWireframeAutoLimit) {
            polyscope::CurveNetwork* cn =
                polyscope::registerCurveNetwork("tet edges", verts, tet_edges(cells));
            cn->setRadius(0.0007);
            cn->setColor(glm::vec3{0.15f, 0.15f, 0.15f});
        }

        if (!duals.empty()) {
            polyscope::PointCloud* dual = polyscope::registerPointCloud("Voronoi vertices", duals);
            dual->setPointRadius(0.002);
            dual->setEnabled(false);
        }
    }

    set_scene_extents(lo, hi, diag, reset_view);
    apply_style();
}

void ui() {
    ImGui::PushItemWidth(160);

    ImGui::TextUnformatted("del3d");
    ImGui::Text("%zu points, %zu vertices, %zu cells",
                g.points.size(), g.tri.number_of_vertices(),
                g.tri.finite_cells().size());
    if (g.duals_dropped)
        ImGui::TextDisabled("%zu Voronoi vertices lie off-scene (slivers)", g.duals_dropped);
    ImGui::Separator();

    // ---- import ------------------------------------------------------------
    ImGui::TextUnformatted("Import");
    ImGui::PushItemWidth(300);
    const bool entered = ImGui::InputText("##path", g.path_field, sizeof g.path_field,
                                          ImGuiInputTextFlags_EnterReturnsTrue);
    ImGui::PopItemWidth();
    ImGui::SameLine();
    if (ImGui::Button("load") || entered) {
        if (load_file(g.path_field)) rebuild(true);
    }
    ImGui::TextDisabled("%s", meshio::supported_formats());

    if (!g.load_error.empty()) ImGui::TextWrapped("could not load: %s", g.load_error.c_str());

    if (g.from_file) {
        ImGui::Text("loaded: %s", g.file.c_str());
        ImGui::Text("%zu mesh vertices, %zu triangles",
                    g.mesh.vertices.size(), g.mesh.faces.size());
        if (g.dropped_vertices)
            ImGui::TextWrapped("%zu non-finite vertices were dropped",
                               g.dropped_vertices);
        if (!g.mesh.faces.empty()) {
            if (ImGui::Checkbox("show the input surface", &g.show_surface)) apply_style();
        }
        if (ImGui::Button("back to a generated set")) { g.from_file = false; rebuild(true); }
    }

    ImGui::Separator();

    // ---- generated point sets ----------------------------------------------
    if (!g.from_file) {
        const char* kinds[] = { "random cube", "lattice", "quantised sphere", "quantised torus" };
        bool dirty = false;
        dirty |= ImGui::Combo("point set", &g.kind, kinds, IM_ARRAYSIZE(kinds));
        dirty |= ImGui::SliderInt("count", &g.count, 8, 20000);
        int seed = int(g.seed);
        if (ImGui::SliderInt("seed", &seed, 1, 100)) { g.seed = unsigned(seed); dirty = true; }
        if (dirty || ImGui::Button("rebuild")) rebuild(true);
        ImGui::TextWrapped("The lattice and quantised sets are massively cospherical: "
                           "that is the case the perturbed predicates exist for.");
    }

    ImGui::Separator();
    ImGui::TextUnformatted("Looking inside");

    bool style = false;
    style |= ImGui::Checkbox("solid cells", &g.solid);
    ImGui::SameLine();
    style |= ImGui::Checkbox("cell edges", &g.edges);
    if (g.edges) style |= ImGui::SliderFloat("edge width", &g.edge_width, 0.1f, 4.0f);

    if (ImGui::Checkbox("wireframe (all edges, interior included)", &g.wireframe)) {
        // On a big triangulation the edge list was skipped at rebuild time, so
        // asking for it now is what builds it.
        if (g.wireframe && !polyscope::hasCurveNetwork("tet edges")) rebuild(false);
        else style = true;
    }
    style |= ImGui::SliderFloat("transparency", &g.transparency, 0.05f, 1.0f);

    // Turning transparency on has a cost, so the mode is only raised when it is
    // actually being used.
    if (g.transparency < 1.0f)
        polyscope::options::transparencyMode = polyscope::TransparencyMode::Pretty;

    bool has_plane = (g.plane != nullptr);
    if (ImGui::Checkbox("slice plane", &has_plane)) {
        if (has_plane) {
            // Visible by default: the plane's gizmo is how it gets dragged, and
            // an invisible one looks like the checkbox did nothing.
            g.plane = polyscope::addSceneSlicePlane(true);
        } else {
            polyscope::removeLastSceneSlicePlane();
            g.plane = nullptr;
        }
        style = true;
    }
    if (g.plane) {
        ImGui::SameLine();
        style |= ImGui::Checkbox("cut cells open", &g.cut_cells);
    }

    if (ImGui::Button("wireframe only")) {
        g.solid = false; g.wireframe = true; g.transparency = 1.0f; style = true;
    }
    ImGui::SameLine();
    if (ImGui::Button("solid")) {
        g.solid = true; g.wireframe = false; g.transparency = 1.0f; style = true;
    }
    ImGui::SameLine();
    if (ImGui::Button("x-ray")) {
        g.solid = true; g.wireframe = false; g.edges = true; g.transparency = 0.35f;
        style = true;
    }

    if (style) apply_style();

    ImGui::PopItemWidth();
}

} // namespace

int main(int argc, char** argv) {
    bool check = false;
    std::string path;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--check") check = true;
        else path = arg;
    }

    // Loaded before the window opens, so a bad path is a message on the
    // terminal rather than an empty view.
    if (!path.empty() && !load_file(path)) {
        std::fprintf(stderr, "%s\n", g.load_error.c_str());
        return 1;
    }

    polyscope::init(check ? "openGL_mock" : "");
    polyscope::state::userCallback = ui;

    rebuild(true);
    if (g.points.empty()) {
        std::fprintf(stderr, "no points to show\n");
        return 1;
    }

    if (check) {
        std::printf("ok: %zu points, %zu vertices, %zu cells registered\n",
                    g.points.size(), g.tri.number_of_vertices(),
                    g.tri.finite_cells().size());
        return 0;
    }

    polyscope::show();
    return 0;
}
