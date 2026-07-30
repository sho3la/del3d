// del3d - the smallest useful program.
//
// Builds the Delaunay triangulation of a point set and reads three things back
// out of it: the tetrahedra, the Voronoi vertices dual to them, and the cells
// incident to one input point.
//
// Run it with no arguments for a built-in point set, or pass a text file with
// one "x y z" per line.
#include <del3d/delaunay.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

static std::vector<del3d::Point3> read_xyz(const char* path) {
    std::vector<del3d::Point3> pts;
    std::ifstream f(path);
    double x, y, z;
    while (f >> x >> y >> z) pts.push_back({x, y, z});
    return pts;
}

static std::vector<del3d::Point3> demo_points() {
    // A lattice: cospherical everywhere, so it is a fair demonstration that the
    // degenerate cases are handled rather than avoided.
    std::vector<del3d::Point3> pts;
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            for (int k = 0; k < 6; ++k)
                pts.push_back({double(i), double(j), double(k)});
    return pts;
}

int main(int argc, char** argv) {
    const std::vector<del3d::Point3> pts =
        (argc > 1) ? read_xyz(argv[1]) : demo_points();

    if (pts.empty()) {
        std::fprintf(stderr, "no points\n");
        return 1;
    }

    // The whole API, in three lines.
    del3d::Delaunay d;
    d.build(pts);
    const std::vector<del3d::Delaunay::Tet>& cells = d.finite_cells();

    std::printf("%zu points, %zu distinct vertices, %zu tetrahedra\n",
                pts.size(), d.number_of_vertices(), cells.size());

    // Each cell's dual is the centre of its circumsphere: the Voronoi diagram
    // falls out of the triangulation this way.
    double rmin = 0.0, rmax = 0.0;
    for (std::size_t i = 0; i < cells.size(); ++i) {
        const del3d::Delaunay::Tet& t = cells[i];
        const del3d::Point3 c = del3d::circumcenter(pts[t[0]], pts[t[1]],
                                                    pts[t[2]], pts[t[3]]);
        const del3d::Point3& a = pts[t[0]];
        const double dx = c.x - a.x, dy = c.y - a.y, dz = c.z - a.z;
        const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (i == 0 || r < rmin) rmin = r;
        if (i == 0 || r > rmax) rmax = r;
    }
    std::printf("circumradius: min %g, max %g\n", rmin, rmax);

    // Cells are reported in a fixed order, and so is each vertex's list of
    // incident cells - the same on every platform, which is what makes it safe
    // to pick "the first incident cell" and get a reproducible answer.
    const std::vector<int>& around0 = d.incident_cells()[0];
    std::printf("point 0 is in %zu cells; the first is (%d %d %d %d)\n",
                around0.size(),
                around0.empty() ? -1 : cells[around0[0]][0],
                around0.empty() ? -1 : cells[around0[0]][1],
                around0.empty() ? -1 : cells[around0[0]][2],
                around0.empty() ? -1 : cells[around0[0]][3]);

    return 0;
}
