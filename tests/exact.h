// del3d tests - an independent exact reference for the predicates.
//
// The library evaluates its determinants with a floating-point filter backed by
// error-free transformations (src/expansion.h). To check that machinery, the
// tests need a second implementation that shares none of it. This header is
// that: arbitrary-precision *integer* arithmetic, no filtering, no floating
// point in the decision path at all.
//
// Turning double coordinates into integers is exact and costs nothing in
// generality. Every finite double is m * 2^e with m a 53-bit integer, so given
// the coordinates of one predicate call we take the smallest e among them and
// write each coordinate as an integer times that common 2^e. Scaling every
// coordinate of a predicate by the same positive factor multiplies its
// determinant by a positive power of that factor, so the sign - which is the
// whole answer - is untouched.
//
// Nothing here is meant to be fast; it is meant to be obviously right.
#ifndef DEL3D_TEST_EXACT_H
#define DEL3D_TEST_EXACT_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace exact {

// ---------------------------------------------------------------------------
// Arbitrary-precision signed integer: sign plus a little-endian base-2^32
// magnitude. Only what the determinants below need - add, subtract, multiply,
// shift, sign.
// ---------------------------------------------------------------------------
class Int {
public:
    Int() = default;

    Int(long long v) {
        unsigned long long m;
        if (v > 0)      { sg_ = 1;  m = static_cast<unsigned long long>(v); }
        else if (v < 0) { sg_ = -1; m = 0ull - static_cast<unsigned long long>(v); }
        else            return;
        d_.push_back(static_cast<std::uint32_t>(m & 0xffffffffu));
        d_.push_back(static_cast<std::uint32_t>(m >> 32));
        trim();
    }

    int sign() const { return sg_; }

    Int operator-() const { Int r = *this; r.sg_ = -r.sg_; return r; }

    /// Multiplication by 2^bits.
    Int operator<<(unsigned bits) const {
        if (sg_ == 0) return Int();
        const unsigned words = bits / 32, rest = bits % 32;
        Int r;
        r.sg_ = sg_;
        r.d_.assign(words, 0u);
        std::uint32_t carry = 0;
        for (std::uint32_t w : d_) {
            const std::uint64_t v = (static_cast<std::uint64_t>(w) << rest) | carry;
            r.d_.push_back(static_cast<std::uint32_t>(v & 0xffffffffu));
            carry = static_cast<std::uint32_t>(v >> 32);
        }
        if (carry) r.d_.push_back(carry);
        r.trim();
        return r;
    }

    friend Int operator+(const Int& a, const Int& b) {
        if (a.sg_ == 0) return b;
        if (b.sg_ == 0) return a;
        Int r;
        if (a.sg_ == b.sg_) {
            r.sg_ = a.sg_;
            r.d_ = add_mag(a.d_, b.d_);
        } else {
            const int c = cmp_mag(a.d_, b.d_);
            if (c == 0) return Int();
            if (c > 0) { r.sg_ = a.sg_; r.d_ = sub_mag(a.d_, b.d_); }
            else       { r.sg_ = b.sg_; r.d_ = sub_mag(b.d_, a.d_); }
        }
        r.trim();
        return r;
    }

    friend Int operator-(const Int& a, const Int& b) { return a + (-b); }

    friend Int operator*(const Int& a, const Int& b) {
        if (a.sg_ == 0 || b.sg_ == 0) return Int();
        Int r;
        r.sg_ = a.sg_ * b.sg_;
        r.d_.assign(a.d_.size() + b.d_.size(), 0u);
        for (std::size_t i = 0; i < a.d_.size(); ++i) {
            std::uint64_t carry = 0;
            for (std::size_t j = 0; j < b.d_.size(); ++j) {
                const std::uint64_t cur = r.d_[i + j]
                    + static_cast<std::uint64_t>(a.d_[i]) * b.d_[j] + carry;
                r.d_[i + j] = static_cast<std::uint32_t>(cur & 0xffffffffu);
                carry = cur >> 32;
            }
            std::size_t k = i + b.d_.size();
            while (carry) {
                const std::uint64_t cur = r.d_[k] + carry;
                r.d_[k] = static_cast<std::uint32_t>(cur & 0xffffffffu);
                carry = cur >> 32;
                ++k;
            }
        }
        r.trim();
        return r;
    }

private:
    int sg_ = 0;
    std::vector<std::uint32_t> d_;

    void trim() {
        while (!d_.empty() && d_.back() == 0u) d_.pop_back();
        if (d_.empty()) sg_ = 0;
    }

    static int cmp_mag(const std::vector<std::uint32_t>& a,
                       const std::vector<std::uint32_t>& b) {
        if (a.size() != b.size()) return a.size() < b.size() ? -1 : 1;
        for (std::size_t i = a.size(); i-- > 0; )
            if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
        return 0;
    }

    static std::vector<std::uint32_t> add_mag(const std::vector<std::uint32_t>& a,
                                              const std::vector<std::uint32_t>& b) {
        std::vector<std::uint32_t> r;
        r.reserve(std::max(a.size(), b.size()) + 1);
        std::uint64_t carry = 0;
        for (std::size_t i = 0; i < std::max(a.size(), b.size()); ++i) {
            const std::uint64_t s = carry + (i < a.size() ? a[i] : 0u)
                                          + (i < b.size() ? b[i] : 0u);
            r.push_back(static_cast<std::uint32_t>(s & 0xffffffffu));
            carry = s >> 32;
        }
        if (carry) r.push_back(static_cast<std::uint32_t>(carry));
        return r;
    }

    /// Precondition: a >= b.
    static std::vector<std::uint32_t> sub_mag(const std::vector<std::uint32_t>& a,
                                              const std::vector<std::uint32_t>& b) {
        std::vector<std::uint32_t> r;
        r.reserve(a.size());
        std::int64_t borrow = 0;
        for (std::size_t i = 0; i < a.size(); ++i) {
            std::int64_t s = static_cast<std::int64_t>(a[i]) - borrow
                           - (i < b.size() ? static_cast<std::int64_t>(b[i]) : 0);
            if (s < 0) { s += (std::int64_t(1) << 32); borrow = 1; } else borrow = 0;
            r.push_back(static_cast<std::uint32_t>(s));
        }
        return r;
    }
};

// ---------------------------------------------------------------------------
// Exact conversion of a group of coordinates to integers on a common scale.
// ---------------------------------------------------------------------------

/// Writes a finite double as m * 2^e with m an integer of at most 53 bits.
inline void decompose(double x, long long& m, int& e) {
    if (x == 0.0) { m = 0; e = 0; return; }
    int bexp = 0;
    const double f = std::frexp(x, &bexp);        // x = f * 2^bexp, 0.5 <= |f| < 1
    m = static_cast<long long>(std::ldexp(f, 53));  // exact: f * 2^53 is an integer
    e = bexp - 53;
}

/// The common scale for one predicate call. Feed it every coordinate that the
/// determinant will touch, then use it to convert them.
class Frame {
public:
    void add(const double* p) {
        for (int k = 0; k < 3; ++k) {
            long long m; int e;
            decompose(p[k], m, e);
            if (m == 0) continue;                 // zero is exact on any scale
            if (!seen_ || e < emin_) { emin_ = e; seen_ = true; }
        }
    }
    void add(double x) {
        long long m; int e;
        decompose(x, m, e);
        if (m == 0) return;
        if (!seen_ || e < emin_) { emin_ = e; seen_ = true; }
    }

    Int operator()(double x) const {
        long long m; int e;
        decompose(x, m, e);
        if (m == 0) return Int();
        return Int(m) << static_cast<unsigned>(e - emin_);
    }

private:
    int  emin_ = 0;
    bool seen_ = false;
};

// ---------------------------------------------------------------------------
// Determinants.
// ---------------------------------------------------------------------------
struct V3 { Int x, y, z; };

inline V3 sub(const V3& a, const V3& b) { return { a.x - b.x, a.y - b.y, a.z - b.z }; }
inline Int norm2(const V3& a) { return a.x * a.x + a.y * a.y + a.z * a.z; }

inline V3 cross(const V3& a, const V3& b) {
    return { a.y * b.z - a.z * b.y,
             a.z * b.x - a.x * b.z,
             a.x * b.y - a.y * b.x };
}

inline Int det3(const V3& a, const V3& b, const V3& c) {
    return a.x * (b.y * c.z - b.z * c.y)
         - a.y * (b.x * c.z - b.z * c.x)
         + a.z * (b.x * c.y - b.y * c.x);
}

/// | ax ay az |a|^2 |
/// | bx by bz |b|^2 |   expanded along the last column.
/// | cx cy cz |c|^2 |
/// | dx dy dz |d|^2 |
inline Int det4_lifted(const V3& a, const V3& b, const V3& c, const V3& d) {
    return   norm2(d) * det3(a, b, c)
           - norm2(c) * det3(a, b, d)
           + norm2(b) * det3(a, c, d)
           - norm2(a) * det3(b, c, d);
}

// ---------------------------------------------------------------------------
// The reference predicates. Signs follow the conventions documented in
// include/del3d/predicates.h and predicates2.h.
// ---------------------------------------------------------------------------

/// Sign of |b-a, c-a, d-a|.
inline int orient3d(const double* a, const double* b, const double* c, const double* d) {
    Frame f;
    f.add(a); f.add(b); f.add(c); f.add(d);
    const V3 A{f(a[0]), f(a[1]), f(a[2])}, B{f(b[0]), f(b[1]), f(b[2])};
    const V3 C{f(c[0]), f(c[1]), f(c[2])}, D{f(d[0]), f(d[1]), f(d[2])};
    return det3(sub(B, A), sub(C, A), sub(D, A)).sign();
}

/// POSITIVE iff e is strictly inside the sphere through a,b,c,d, ZERO iff the
/// five points are exactly cospherical.
/// \pre orient3d(a,b,c,d) > 0
inline int in_sphere_raw(const double* a, const double* b, const double* c,
                         const double* d, const double* e) {
    Frame f;
    f.add(a); f.add(b); f.add(c); f.add(d); f.add(e);
    const V3 E{f(e[0]), f(e[1]), f(e[2])};
    const V3 A = sub(V3{f(a[0]), f(a[1]), f(a[2])}, E);
    const V3 B = sub(V3{f(b[0]), f(b[1]), f(b[2])}, E);
    const V3 C = sub(V3{f(c[0]), f(c[1]), f(c[2])}, E);
    const V3 D = sub(V3{f(d[0]), f(d[1]), f(d[2])}, E);
    // The lifted determinant is negative when e is inside a positively oriented
    // (in this library's convention) tetrahedron - see tests/test_predicates.cpp,
    // which pins that on a configuration whose answer is known by hand.
    return -det4_lifted(A, B, C, D).sign();
}

/// The cospherical perturbation rule of in_sphere(), written out here against
/// the reference orientation so the library's version is compared with a second
/// implementation rather than with itself.
inline int in_sphere(const double* a, const double* b, const double* c,
                     const double* d, const double* e) {
    const int s = in_sphere_raw(a, b, c, d, e);
    if (s != 0) return s;

    const double* p[5] = { a, b, c, d, e };
    std::sort(p, p + 5, [](const double* x, const double* y) {
        if (x[0] != y[0]) return x[0] < y[0];
        if (x[1] != y[1]) return x[1] < y[1];
        return x[2] < y[2];
    });
    for (int i = 4; i > 2; --i) {
        if (p[i] == e) return -1;
        int o;
        if (p[i] == d && (o = orient3d(a, b, c, e)) != 0) return o;
        if (p[i] == c && (o = orient3d(a, b, e, d)) != 0) return o;
        if (p[i] == b && (o = orient3d(a, e, c, d)) != 0) return o;
        if (p[i] == a && (o = orient3d(e, b, c, d)) != 0) return o;
    }
    return -1;
}

/// Sign of the 2x2 determinant |q-p, r-p|.
inline int orient2d(double px, double py, double qx, double qy, double rx, double ry) {
    Frame f;
    f.add(px); f.add(py); f.add(qx); f.add(qy); f.add(rx); f.add(ry);
    const Int PX = f(px), PY = f(py);
    const Int QX = f(qx) - PX, QY = f(qy) - PY;
    const Int RX = f(rx) - PX, RY = f(ry) - PY;
    return (QX * RY - QY * RX).sign();
}

inline bool collinear(const double* p, const double* q, const double* r) {
    return orient2d(p[0], p[1], q[0], q[1], r[0], r[1]) == 0
        && orient2d(p[1], p[2], q[1], q[2], r[1], r[2]) == 0
        && orient2d(p[0], p[2], q[0], q[2], r[0], r[2]) == 0;
}

/// Orientation in the first of the xy, yz, xz projections that is not
/// degenerate. \pre p, q, r not collinear.
inline int coplanar_orientation(const double* p, const double* q, const double* r) {
    int o = orient2d(p[0], p[1], q[0], q[1], r[0], r[1]);
    if (o != 0) return o;
    o = orient2d(p[1], p[2], q[1], q[2], r[1], r[2]);
    if (o != 0) return o;
    return orient2d(p[0], p[2], q[0], q[2], r[0], r[2]);
}

/// s against the directed line (p,q), signed consistently with (p,q,r).
inline int coplanar_orientation(const double* p, const double* q,
                                const double* r, const double* s) {
    int o = orient2d(p[0], p[1], q[0], q[1], r[0], r[1]);
    if (o != 0) return o * orient2d(p[0], p[1], q[0], q[1], s[0], s[1]);
    o = orient2d(p[1], p[2], q[1], q[2], r[1], r[2]);
    if (o != 0) return o * orient2d(p[1], p[2], q[1], q[2], s[1], s[2]);
    o = orient2d(p[0], p[2], q[0], q[2], r[0], r[2]);
    return o * orient2d(p[0], p[2], q[0], q[2], s[0], s[2]);
}

/**
 * POSITIVE iff p is strictly inside the circle through the coplanar points
 * p0, p1, p2; ZERO iff the four are exactly cocircular.
 *
 * Derived independently of the library's formulation: every sphere through
 * p0, p1, p2 meets their plane in exactly their circumcircle, so the in-circle
 * question for a coplanar p is the in-sphere question for p against the
 * tetrahedron (p0, p1, p2, apex) for any apex off the plane. Taking
 * apex = p0 + (p1-p0) x (p2-p0) makes that tetrahedron positively oriented by
 * construction - its determinant is the squared length of the normal - so no
 * orientation correction is needed.
 */
inline int coplanar_in_circle_raw(const double* p0, const double* p1,
                                  const double* p2, const double* p) {
    Frame f;
    f.add(p0); f.add(p1); f.add(p2); f.add(p);
    const V3 P{f(p[0]), f(p[1]), f(p[2])};
    const V3 A = sub(V3{f(p0[0]), f(p0[1]), f(p0[2])}, P);
    const V3 B = sub(V3{f(p1[0]), f(p1[1]), f(p1[2])}, P);
    const V3 C = sub(V3{f(p2[0]), f(p2[1]), f(p2[2])}, P);

    const V3 n = cross(sub(B, A), sub(C, A));
    const V3 D = { A.x + n.x, A.y + n.y, A.z + n.z };   // the apex, p-relative

    return -det4_lifted(A, B, C, D).sign();
}

/// The cocircular perturbation rule of coplanar_in_circle().
inline int coplanar_in_circle(const double* p0, const double* p1,
                              const double* p2, const double* p) {
    const int bs = coplanar_in_circle_raw(p0, p1, p2, p);
    if (bs != 0) return bs;

    const double* pts[4] = { p0, p1, p2, p };
    std::sort(pts, pts + 4, [](const double* a, const double* b) {
        if (a[0] != b[0]) return a[0] < b[0];
        if (a[1] != b[1]) return a[1] < b[1];
        return a[2] < b[2];
    });
    const int local = coplanar_orientation(p0, p1, p2);
    for (int i = 3; i > 0; --i) {
        if (pts[i] == p) return -1;
        int o;
        if (pts[i] == p2 && (o = coplanar_orientation(p0, p1, p)) != 0) return o * local;
        if (pts[i] == p1 && (o = coplanar_orientation(p0, p, p2)) != 0) return o * local;
        if (pts[i] == p0 && (o = coplanar_orientation(p, p1, p2)) != 0) return o * local;
    }
    return -local;
}

} // namespace exact

#endif // DEL3D_TEST_EXACT_H
