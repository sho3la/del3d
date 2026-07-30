// del3d viewer - mesh reading. See mesh_io.h.
#include "mesh_io.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <istream>
#include <map>
#include <sstream>
#include <utility>

namespace meshio {
namespace {

std::string lowercase(std::string s) {
    for (char& c : s) c = char(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

std::string extension_of(const std::string& path) {
    const std::size_t dot = path.find_last_of('.');
    const std::size_t sep = path.find_last_of("/\\");
    if (dot == std::string::npos || (sep != std::string::npos && dot < sep)) return "";
    return lowercase(path.substr(dot));
}

void add_polygon(Mesh& m, const std::vector<int>& poly) {
    // Fan triangulation. Fine for the convex quads and n-gons that exporters
    // emit; a non-convex n-gon would need ear clipping, and is not worth it
    // here since the faces are only ever drawn, never triangulated against.
    for (std::size_t k = 2; k < poly.size(); ++k)
        m.faces.push_back({poly[0], poly[k - 1], poly[k]});
}

// ---------------------------------------------------------------------------
// OBJ
// ---------------------------------------------------------------------------
Mesh load_obj(std::istream& in, std::string& error) {
    Mesh m;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ls(line);
        std::string tag;
        ls >> tag;
        if (tag == "v") {
            double x, y, z;
            if (ls >> x >> y >> z) m.vertices.push_back({x, y, z});
        } else if (tag == "f") {
            std::vector<int> poly;
            std::string tok;
            while (ls >> tok) {
                // "v", "v/vt", "v//vn" and "v/vt/vn" all start with the vertex.
                const std::size_t slash = tok.find('/');
                if (slash != std::string::npos) tok.resize(slash);
                if (tok.empty()) continue;
                long idx = std::strtol(tok.c_str(), nullptr, 10);
                if (idx < 0) idx = long(m.vertices.size()) + idx;   // relative
                else         idx -= 1;                              // 1-based
                if (idx >= 0 && idx < long(m.vertices.size())) poly.push_back(int(idx));
            }
            if (poly.size() >= 3) add_polygon(m, poly);
        }
    }
    if (m.vertices.empty()) error = "no vertices in the OBJ";
    return m;
}

// ---------------------------------------------------------------------------
// OFF
// ---------------------------------------------------------------------------
Mesh load_off(std::istream& in, std::string& error) {
    Mesh m;

    // Skip comments and blank lines, and tolerate the counts sitting on the
    // same line as the magic word, which some writers do.
    auto next_tokens = [&in](std::istringstream& out) -> bool {
        std::string line;
        while (std::getline(in, line)) {
            const std::size_t hash = line.find('#');
            if (hash != std::string::npos) line.resize(hash);
            if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;
            out.clear();
            out.str(line);
            return true;
        }
        return false;
    };

    std::istringstream ls;
    if (!next_tokens(ls)) { error = "empty OFF"; return m; }
    std::string magic;
    ls >> magic;
    if (lowercase(magic).find("off") == std::string::npos) {
        error = "not an OFF file";
        return m;
    }

    long nv = 0, nf = 0, ne = 0;
    if (!(ls >> nv >> nf >> ne)) {                 // counts on the next line
        if (!next_tokens(ls) || !(ls >> nv >> nf >> ne)) {
            error = "malformed OFF header";
            return m;
        }
    }

    m.vertices.reserve(std::size_t(std::max(0L, nv)));
    for (long i = 0; i < nv; ++i) {
        double x, y, z;
        if (!(in >> x >> y >> z)) { error = "OFF ended inside the vertex list"; return m; }
        m.vertices.push_back({x, y, z});
        std::string rest;                          // per-vertex colour or normal
        std::getline(in, rest);
    }
    for (long i = 0; i < nf; ++i) {
        long n = 0;
        if (!(in >> n)) break;
        std::vector<int> poly;
        poly.reserve(std::size_t(std::max(0L, n)));
        for (long k = 0; k < n; ++k) {
            long idx = 0;
            if (!(in >> idx)) break;
            if (idx >= 0 && idx < long(m.vertices.size())) poly.push_back(int(idx));
        }
        if (poly.size() >= 3) add_polygon(m, poly);
        std::string rest;                          // per-face colour
        std::getline(in, rest);
    }
    return m;
}

// ---------------------------------------------------------------------------
// PLY
// ---------------------------------------------------------------------------
struct PlyProperty {
    std::string name;
    std::string type;        // scalar type, or the item type of a list
    bool        is_list = false;
    std::string count_type;
};

struct PlyElement {
    std::string              name;
    long                     count = 0;
    std::vector<PlyProperty> props;
};

int type_size(const std::string& t) {
    if (t == "char" || t == "uchar" || t == "int8" || t == "uint8")     return 1;
    if (t == "short" || t == "ushort" || t == "int16" || t == "uint16") return 2;
    if (t == "int" || t == "uint" || t == "int32" || t == "uint32" ||
        t == "float" || t == "float32")                                 return 4;
    if (t == "double" || t == "float64")                                return 8;
    return 0;
}

bool is_float_type(const std::string& t) {
    return t == "float" || t == "float32" || t == "double" || t == "float64";
}

bool is_signed_type(const std::string& t) {
    return t == "char" || t == "int8" || t == "short" || t == "int16" ||
           t == "int"  || t == "int32";
}

/// Reads one little-endian value and widens it to double, which holds every
/// PLY integer type exactly.
bool read_binary_value(std::istream& in, const std::string& type, double& out) {
    const int n = type_size(type);
    if (n == 0) return false;
    unsigned char buf[8];
    in.read(reinterpret_cast<char*>(buf), n);
    if (!in) return false;

    std::uint64_t raw = 0;
    for (int i = n - 1; i >= 0; --i) raw = (raw << 8) | buf[i];   // little-endian

    if (is_float_type(type)) {
        if (n == 4) {
            float f;
            const std::uint32_t r32 = std::uint32_t(raw);
            std::memcpy(&f, &r32, 4);
            out = double(f);
        } else {
            double d;
            std::memcpy(&d, &raw, 8);
            out = d;
        }
    } else if (is_signed_type(type)) {
        const int shift = 64 - 8 * n;
        out = double(std::int64_t(raw << shift) >> shift);        // sign-extend
    } else {
        out = double(raw);
    }
    return true;
}

Mesh load_ply(std::ifstream& in, std::string& error) {
    Mesh m;

    std::string line;
    if (!std::getline(in, line) || lowercase(line).substr(0, 3) != "ply") {
        error = "not a PLY file";
        return m;
    }

    std::string format = "ascii";
    std::vector<PlyElement> elements;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        std::istringstream ls(line);
        std::string tag;
        ls >> tag;
        if (tag == "comment" || tag == "obj_info") continue;
        if (tag == "format") { ls >> format; continue; }
        if (tag == "element") {
            PlyElement e;
            ls >> e.name >> e.count;
            elements.push_back(e);
            continue;
        }
        if (tag == "property") {
            if (elements.empty()) continue;
            PlyProperty p;
            std::string t;
            ls >> t;
            if (t == "list") {
                p.is_list = true;
                ls >> p.count_type >> p.type >> p.name;
            } else {
                p.type = t;
                ls >> p.name;
            }
            elements.back().props.push_back(p);
            continue;
        }
        if (tag == "end_header") break;
    }

    if (format != "ascii" && format != "binary_little_endian") {
        error = "unsupported PLY format '" + format + "'"
                " (ascii and binary_little_endian are read)";
        return m;
    }
    const bool ascii = (format == "ascii");

    for (const PlyElement& e : elements) {
        const bool is_vertex = (e.name == "vertex");
        const bool is_face   = (e.name == "face");

        int ix = -1, iy = -1, iz = -1;
        for (std::size_t i = 0; i < e.props.size(); ++i) {
            if (e.props[i].name == "x") ix = int(i);
            if (e.props[i].name == "y") iy = int(i);
            if (e.props[i].name == "z") iz = int(i);
        }
        if (is_vertex && (ix < 0 || iy < 0 || iz < 0)) {
            error = "PLY vertex element has no x/y/z";
            return m;
        }

        for (long n = 0; n < e.count; ++n) {
            std::vector<double> scalars(e.props.size(), 0.0);
            std::vector<int>    poly;

            std::istringstream as;
            if (ascii) {
                if (!std::getline(in, line)) { error = "PLY ended early"; return m; }
                as.str(line);
                as.clear();
            }

            for (std::size_t p = 0; p < e.props.size(); ++p) {
                const PlyProperty& prop = e.props[p];
                if (prop.is_list) {
                    double cnt = 0;
                    if (ascii) { if (!(as >> cnt)) break; }
                    else if (!read_binary_value(in, prop.count_type, cnt)) {
                        error = "PLY ended inside a list"; return m;
                    }
                    for (long k = 0; k < long(cnt); ++k) {
                        double v = 0;
                        if (ascii) { if (!(as >> v)) break; }
                        else if (!read_binary_value(in, prop.type, v)) {
                            error = "PLY ended inside a list"; return m;
                        }
                        if (is_face) poly.push_back(int(v));
                    }
                } else {
                    double v = 0;
                    if (ascii) { if (!(as >> v)) break; }
                    else if (!read_binary_value(in, prop.type, v)) {
                        error = "PLY ended inside an element"; return m;
                    }
                    scalars[p] = v;
                }
            }

            if (is_vertex) m.vertices.push_back({scalars[ix], scalars[iy], scalars[iz]});
            if (is_face) {
                std::vector<int> clean;
                for (int idx : poly)
                    if (idx >= 0 && idx < int(m.vertices.size())) clean.push_back(idx);
                if (clean.size() >= 3) add_polygon(m, clean);
            }
        }
    }

    if (m.vertices.empty() && error.empty()) error = "no vertices in the PLY";
    return m;
}

// ---------------------------------------------------------------------------
// STL
// ---------------------------------------------------------------------------
Mesh load_stl_binary(std::ifstream& in, std::string& error) {
    Mesh m;
    in.clear();
    in.seekg(80, std::ios::beg);
    std::uint32_t n = 0;
    in.read(reinterpret_cast<char*>(&n), 4);
    if (!in) { error = "truncated STL header"; return m; }

    m.vertices.reserve(std::size_t(n) * 3);
    m.faces.reserve(n);
    for (std::uint32_t i = 0; i < n; ++i) {
        float buf[12];
        in.read(reinterpret_cast<char*>(buf), 48);
        std::uint16_t attr = 0;
        in.read(reinterpret_cast<char*>(&attr), 2);
        if (!in) break;                       // tolerate a truncated tail
        const int base = int(m.vertices.size());
        for (int k = 1; k <= 3; ++k)
            m.vertices.push_back({double(buf[3 * k]), double(buf[3 * k + 1]),
                                  double(buf[3 * k + 2])});
        m.faces.push_back({base, base + 1, base + 2});
    }
    if (m.vertices.empty()) error = "no triangles in the STL";
    return m;
}

Mesh load_stl_ascii(std::istream& in, std::string& error) {
    Mesh m;
    std::string line;
    std::vector<int> tri;
    while (std::getline(in, line)) {
        std::istringstream ls(line);
        std::string tag;
        ls >> tag;
        if (lowercase(tag) != "vertex") continue;
        double x, y, z;
        if (!(ls >> x >> y >> z)) continue;
        tri.push_back(int(m.vertices.size()));
        m.vertices.push_back({x, y, z});
        if (tri.size() == 3) {
            m.faces.push_back({tri[0], tri[1], tri[2]});
            tri.clear();
        }
    }
    if (m.vertices.empty()) error = "no vertices in the STL";
    return m;
}

/// An STL is binary unless it both starts with "solid" and holds a "facet"
/// keyword: the word "solid" alone is not decisive, since a binary file's
/// 80-byte header may happen to begin with it.
bool stl_is_binary(std::ifstream& in) {
    in.clear();
    in.seekg(0, std::ios::end);
    const std::streamoff size = in.tellg();
    in.seekg(0, std::ios::beg);

    char head[512] = {0};
    in.read(head, sizeof head - 1);
    const std::streamsize got = in.gcount();
    in.clear();
    in.seekg(0, std::ios::beg);

    const std::string text(head, std::size_t(got));
    if (lowercase(text.substr(0, 5)) != "solid") return true;
    if (lowercase(text).find("facet") != std::string::npos) return false;

    // Ambiguous: fall back on the exact size a binary file must have.
    if (size >= 84) {
        in.seekg(80, std::ios::beg);
        std::uint32_t n = 0;
        in.read(reinterpret_cast<char*>(&n), 4);
        in.clear();
        in.seekg(0, std::ios::beg);
        if (std::streamoff(84) + std::streamoff(n) * 50 == size) return true;
    }
    return false;
}

// ---------------------------------------------------------------------------
// Plain point lists
// ---------------------------------------------------------------------------
Mesh load_xyz(std::istream& in, std::string& error) {
    Mesh m;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ls(line);
        double x, y, z;
        if (ls >> x >> y >> z) m.vertices.push_back({x, y, z});
    }
    if (m.vertices.empty()) error = "no points in the file";
    return m;
}

} // namespace

void weld(Mesh& m) {
    std::map<std::array<double, 3>, int> seen;
    std::vector<int> remap(m.vertices.size());
    std::vector<std::array<double, 3>> unique;
    unique.reserve(m.vertices.size());

    for (std::size_t i = 0; i < m.vertices.size(); ++i) {
        auto it = seen.find(m.vertices[i]);
        if (it == seen.end()) {
            const int id = int(unique.size());
            seen.emplace(m.vertices[i], id);
            unique.push_back(m.vertices[i]);
            remap[i] = id;
        } else {
            remap[i] = it->second;
        }
    }

    if (unique.size() == m.vertices.size()) return;

    for (std::array<int, 3>& f : m.faces)
        for (int& v : f) v = remap[std::size_t(v)];

    // A face whose corners collapsed onto each other is no longer a triangle.
    m.faces.erase(std::remove_if(m.faces.begin(), m.faces.end(),
                                 [](const std::array<int, 3>& f) {
                                     return f[0] == f[1] || f[1] == f[2] || f[0] == f[2];
                                 }),
                  m.faces.end());

    m.vertices.swap(unique);
}

std::size_t drop_non_finite(Mesh& m) {
    std::vector<int> remap(m.vertices.size(), -1);
    std::vector<std::array<double, 3>> kept;
    kept.reserve(m.vertices.size());

    for (std::size_t i = 0; i < m.vertices.size(); ++i) {
        const std::array<double, 3>& v = m.vertices[i];
        if (std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2])) {
            remap[i] = int(kept.size());
            kept.push_back(v);
        }
    }

    const std::size_t dropped = m.vertices.size() - kept.size();
    if (dropped == 0) return 0;

    std::vector<std::array<int, 3>> faces;
    faces.reserve(m.faces.size());
    for (const std::array<int, 3>& f : m.faces) {
        const int a = remap[std::size_t(f[0])];
        const int b = remap[std::size_t(f[1])];
        const int c = remap[std::size_t(f[2])];
        if (a >= 0 && b >= 0 && c >= 0) faces.push_back({a, b, c});
    }

    m.vertices.swap(kept);
    m.faces.swap(faces);
    return dropped;
}

const char* supported_formats() { return "obj, off, ply, stl, xyz"; }

Mesh load(const std::string& path, std::string& error) {
    error.clear();

    std::ifstream in(path, std::ios::binary);
    if (!in) { error = "cannot open " + path; return Mesh(); }

    const std::string ext = extension_of(path);
    Mesh m;
    if (ext == ".obj") {
        m = load_obj(in, error);
    } else if (ext == ".off") {
        m = load_off(in, error);
    } else if (ext == ".ply") {
        m = load_ply(in, error);
    } else if (ext == ".stl") {
        m = stl_is_binary(in) ? load_stl_binary(in, error) : load_stl_ascii(in, error);
        weld(m);                       // STL stores every corner separately
    } else if (ext == ".xyz" || ext == ".txt" || ext == ".pts" || ext.empty()) {
        m = load_xyz(in, error);
    } else {
        error = "unknown extension '" + ext + "' (" + supported_formats() + ")";
    }
    return m;
}

} // namespace meshio
