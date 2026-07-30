// del3d tests - the point sets the tests run on.
//
// The interesting input for an exact-predicate library is not random points -
// random points are in general position and never exercise the exact stage.
// What matters is input whose coordinates repeat and align: CAD models, voxel
// and lattice data, and anything quantised, which is most scanned or exported
// geometry. Every generator here therefore snaps its samples to a power-of-two
// grid, which is exactly representable in binary floating point and so produces
// genuinely - not nearly - coplanar and cospherical configurations, in bulk.
//
// All of them are deterministic functions of their arguments, so a failure is
// always reproducible.
#ifndef DEL3D_TEST_POINTSETS_H
#define DEL3D_TEST_POINTSETS_H

#include <array>
#include <cmath>
#include <cstddef>
#include <random>
#include <vector>

namespace pointsets {

using Pt = std::array<double, 3>;

constexpr double kTwoPi = 6.283185307179586;

/**
 * A deterministic generator with a fully specified mapping onto doubles.
 *
 * std::mt19937 is specified bit for bit by the standard, but the distribution
 * adaptors are not: uniform_real_distribution may consume a different number of
 * words, or map them differently, in each standard library. The tests compare
 * results against pinned values, so the mapping is written out here instead.
 */
class Rng {
public:
    explicit Rng(unsigned seed) : g_(seed) {}
    /// Uniform in [0,1) at 32-bit resolution.
    double u01() { return double(g_()) * (1.0 / 4294967296.0); }
    double range(double lo, double hi) { return lo + (hi - lo) * u01(); }
    std::size_t below(std::size_t n) { return std::size_t(g_() % n); }
private:
    std::mt19937 g_;
};

/// Rounds to a multiple of 1/q, with q a power of two so the result is exact.
inline double snap(double v, int q) { return std::floor(v * q + 0.5) / q; }

/// n samples of the lateral surface of a cylinder, quantised.
inline std::vector<Pt> cylinder(std::size_t n, int q = 32, unsigned seed = 1) {
    std::vector<Pt> out;
    out.reserve(n);
    Rng rng(seed);
    for (std::size_t i = 0; i < n; ++i) {
        const double a = kTwoPi * rng.u01();
        out.push_back({ snap(2.0 * std::cos(a), q),
                        snap(2.0 * std::sin(a), q),
                        snap(rng.range(-1.5, 1.5), q) });
    }
    return out;
}

/// n samples of a sphere, quantised.
inline std::vector<Pt> sphere(std::size_t n, int q = 32, unsigned seed = 2) {
    std::vector<Pt> out;
    out.reserve(n);
    Rng rng(seed);
    for (std::size_t i = 0; i < n; ++i) {
        const double z = rng.range(-1.0, 1.0);
        const double r = std::sqrt(1.0 - z * z);
        const double t = kTwoPi * rng.u01();
        out.push_back({ snap(2.0 * r * std::cos(t), q),
                        snap(2.0 * r * std::sin(t), q),
                        snap(2.0 * z, q) });
    }
    return out;
}

/// n samples of a torus, quantised.
inline std::vector<Pt> torus(std::size_t n, int q = 32, unsigned seed = 3) {
    std::vector<Pt> out;
    out.reserve(n);
    Rng rng(seed);
    for (std::size_t i = 0; i < n; ++i) {
        const double u = kTwoPi * rng.u01(), v = kTwoPi * rng.u01();
        const double R = 2.0 + 0.75 * std::cos(v);
        out.push_back({ snap(R * std::cos(u), q),
                        snap(R * std::sin(u), q),
                        snap(0.75 * std::sin(v), q) });
    }
    return out;
}

/// The surface of an axis-aligned box on an integer lattice: flat faces and
/// sharp edges, so a large fraction of every subset is exactly coplanar.
inline std::vector<Pt> box_shell(int m) {
    std::vector<Pt> out;
    for (int i = 0; i <= m; ++i)
        for (int j = 0; j <= m; ++j)
            for (int k = 0; k <= m; ++k)
                if (i == 0 || i == m || j == 0 || j == m || k == 0 || k == m)
                    out.push_back({ double(i), double(j), double(k) });
    return out;
}

/// A full integer lattice: the standard worst case for an in-sphere predicate.
inline std::vector<Pt> lattice(int nx, int ny, int nz) {
    std::vector<Pt> out;
    out.reserve(std::size_t(nx) * ny * nz);
    for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j)
            for (int k = 0; k < nz; ++k)
                out.push_back({ double(i), double(j), double(k) });
    return out;
}

/// Uniform random points in a cube - generic position, no degeneracies.
inline std::vector<Pt> random_cube(std::size_t n, unsigned seed = 4) {
    std::vector<Pt> out;
    out.reserve(n);
    Rng rng(seed);
    for (std::size_t i = 0; i < n; ++i)
        out.push_back({ rng.range(-1.0, 1.0), rng.range(-1.0, 1.0), rng.range(-1.0, 1.0) });
    return out;
}

} // namespace pointsets

#endif // DEL3D_TEST_POINTSETS_H
