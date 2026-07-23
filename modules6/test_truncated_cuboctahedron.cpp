// test_truncated_cuboctahedron.cpp
//
// Headless (no GL) runtime verification that Cube::getFacetsTruncatedCuboctahedron
// builds a correct truncated cuboctahedron from a unit cube's 8 corners:
//   SOLID  : 92 triangles, 48 vertices, 138 mesh-edges, V-E+F = 2,
//            every mesh edge shared by exactly 2 (closed manifold),
//            72 polyhedron edges (the non-coplanar-adjacent ones), all equal
//            length = 2r/(1+2√2), all triangle normals outward.
//   HOLLOW : 288 triangles, every mesh edge shared by 1 or 2 (closed outer shell
//            with 26 hole boundaries), no edge shared by > 2.
//
// Mirrors the offline verification in truncated_octahedron_mesh.md §9 (which used
// a throwaway /tmp script); this one stays in the repo so the check is repeatable.
//
// Build:  g++ -std=c++17 test_truncated_cuboctahedron.cpp -o test_truncated_cuboctahedron -lm
// Run:    ./test_truncated_cuboctahedron

#include <vector>
#include <array>
#include <set>
#include <map>
#include <cmath>
#include <cstdio>
#include <string>
#include <iostream>
#include "Vector3D.cpp"
#include "Quaternion.cpp"
#include "Facet.cpp"
#include "FacetBox.hpp"
#include "Cube.hpp"

using namespace std;

static bool approx(double a, double b, double eps = 1e-6) { return fabs(a - b) < eps; }

static string key(const Vector3D& v) {
    char b[96];
    snprintf(b, sizeof(b), "%.6f,%.6f,%.6f", v.x(), v.y(), v.z());
    return string(b);
}

int main() {
    int fails = 0;
    auto check = [&](bool cond, const string& m) {
        cout << (cond ? "PASS " : "FAIL ") << m << "\n";
        if (!cond) ++fails;
    };

    const double r = 1.0;
    Cube cube(r, Vector3D{0, 0, 0});            // unit cube, non-subdivided
    const double R = 1.0 + 2.0 * sqrt(2.0);    // canonical cube half-side carrying octagons
    const double expected_edge = 2.0 * r / R;  // = 2/(1+2√2) ≈ 0.5224078

    // ---------------- SOLID ----------------
    FacetBox fb = cube.getFacetsTruncatedCuboctahedron(false, 0.5);
    const size_t T = fb.size();
    check(T == 92, "solid: 92 triangles (got " + to_string(T) + ")");

    vector<Vector3D> verts;
    map<string, int> vidx;
    auto idx = [&](const Vector3D& v) -> int {
        string k = key(v);
        auto it = vidx.find(k);
        if (it != vidx.end()) return it->second;
        int i = (int)verts.size(); verts.push_back(v); vidx[k] = i; return i;
    };

    vector<array<int, 3>> tris;
    map<pair<int, int>, vector<int>> edgeTris;     // edge -> triangle indices
    for (size_t i = 0; i < fb.size(); ++i) {
        const Facet& f = fb[i];
        int a = idx(f[0]), b = idx(f[1]), c = idx(f[2]);
        tris.push_back({a, b, c});
        for (auto e : vector<pair<int, int>>{{a, b}, {b, c}, {c, a}}) {
            int lo = min(e.first, e.second), hi = max(e.first, e.second);
            edgeTris[{lo, hi}].push_back((int)i);
        }
    }
    const int V = (int)verts.size();
    const int E = (int)edgeTris.size();
    check(V == 48, "solid: 48 unique vertices (got " + to_string(V) + ")");
    check(E == 138, "solid: 138 mesh edges (got " + to_string(E) + ")");
    check(V - E + (int)T == 2, "solid: Euler V-E+F = 2 (got " + to_string(V - E + (int)T) + ")");
    int shared2 = 0, over2 = 0;
    for (auto& kv : edgeTris) {
        if (kv.second.size() == 2) ++shared2;
        else if (kv.second.size() > 2) ++over2;
    }
    check(shared2 == E && over2 == 0,
          "solid: every mesh edge shared by exactly 2 (" + to_string(shared2) + "/" + to_string(E) + ")");

    // polyhedron edges = mesh edges whose two adjacent triangles are NOT coplanar
    // (fan diagonals are coplanar within one face -> parallel normals; poly-edges
    // join two distinct faces -> non-parallel normals). All 72 must be equal.
    int polyEdges = 0; bool allEqual = true;
    double lmin = 1e9, lmax = 0.0;
    for (auto& kv : edgeTris) {
        if (kv.second.size() != 2) continue;
        Vector3D n1 = fb[kv.second[0]].getNormal();
        Vector3D n2 = fb[kv.second[1]].getNormal();
        double clen = abs(n1 % n2);                 // |cross(n1,n2)|: 0 -> coplanar
        if (clen < 1e-6) continue;                   // fan diagonal
        ++polyEdges;
        double L = abs(verts[kv.first.first] - verts[kv.first.second]);
        lmin = min(lmin, L); lmax = max(lmax, L);
        if (!approx(L, expected_edge, 1e-6)) allEqual = false;
    }
    check(polyEdges == 72, "solid: 72 polyhedron edges (got " + to_string(polyEdges) + ")");
    check(allEqual, "solid: all 72 polyhedron edges equal expected 2r/(1+2√2)");
    check(approx(lmin, expected_edge) && approx(lmax, expected_edge),
          "solid: edge lengths min/max match expected (min=" + to_string(lmin) + " max=" + to_string(lmax) + ")");

    // normals outward (cube centered at origin -> outward = normal . tri_centroid > 0)
    bool outOk = true;
    for (size_t i = 0; i < fb.size(); ++i) {
        const Facet& f = fb[i];
        Vector3D n = f.getNormal();
        Vector3D tc = (f[0] + f[1] + f[2]) / 3.0;
        if (n * tc <= 1e-9) outOk = false;
    }
    check(outOk, "solid: all triangle normals outward (normal . tri_centroid > 0)");

    // ---------------- HOLLOW ----------------
    FacetBox fbh = cube.getFacetsTruncatedCuboctahedron(true, 0.5);
    const size_t Th = fbh.size();
    check(Th == 288, "hollow: 288 triangles (got " + to_string(Th) + ")");

    map<pair<int, int>, int> hEdgeCount;
    verts.clear(); vidx.clear();
    auto hidx = [&](const Vector3D& v) -> int {
        string k = key(v);
        auto it = vidx.find(k);
        if (it != vidx.end()) return it->second;
        int i = (int)verts.size(); verts.push_back(v); vidx[k] = i; return i;
    };
    for (size_t i = 0; i < fbh.size(); ++i) {
        const Facet& f = fbh[i];
        int a = hidx(f[0]), b = hidx(f[1]), c = hidx(f[2]);
        for (auto e : vector<pair<int, int>>{{a, b}, {b, c}, {c, a}}) {
            int lo = min(e.first, e.second), hi = max(e.first, e.second);
            hEdgeCount[{lo, hi}]++;
        }
    }
    int b1 = 0, b2 = 0, bover = 0;
    for (auto& kv : hEdgeCount) {
        if (kv.second == 1) ++b1;
        else if (kv.second == 2) ++b2;
        else ++bover;
    }
    // 26 holes: 6 octagons*8 + 8 hexagons*6 + 12 squares*4 = 48+48+48 = 144 boundary edges
    check(b1 == 144, "hollow: 144 boundary edges (26 hole loops) (got " + to_string(b1) + ")");
    check(bover == 0, "hollow: no mesh edge shared by > 2 (clean manifold with boundary)");
    check(b1 + b2 == (int)hEdgeCount.size(), "hollow: all edges shared by 1 or 2");

    cout << "\n";
    if (fails == 0) { cout << "ALL CHECKS PASSED\n"; return 0; }
    cout << "*** " << fails << " CHECK(S) FAILED ***\n";
    return 1;
}