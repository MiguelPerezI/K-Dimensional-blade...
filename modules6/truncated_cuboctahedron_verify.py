#!/usr/bin/env python3
"""
truncated_cuboctahedron_verify.py

Generates the truncated cuboctahedron (great rhombicuboctahedron) from its
canonical coordinates, builds the 26 faces (6 octagons + 8 hexagons + 12
squares) with CCW-outward winding, verifies the topology/regularity, and emits
the four `static constexpr` tables that Cube.hpp uses for the
pushTruncatedCuboctahedronFacets surface method:

    tc_vert_norm_[48][3]   -- 48 normalized positions (perms of (±a, ±b, ±1)),
                              a = 1/(1+2√2), b = (1+1+√2)/(1+2√2); the parent
                              cube is half-side 1 in normalized space.
    tc_oct_faces_[6][8]     -- octagon  vertex indices, CCW-outward
    tc_hex_faces_[8][6]     -- hexagon  vertex indices, CCW-outward
    tc_sq_faces_[12][4]     -- square   vertex indices, CCW-outward

Canonical model (edge length 2, centered at origin):
    48 vertices = all permutations of (±1, ±(1+√2), ±(1+2√2)),  R = 1+2√2.
    Octagons  on planes  x = ±R,  y = ±R,  z = ±R           (6 faces)
    Hexagons  on planes  ±x ± y ± z = 3+3√2                 (8 faces)
    Squares   on planes  ±x ± y, ±x ± z, ±y ± z = 2+3√2     (12 faces)

The C++ rebuilds each vertex by trilinear interpolation of the parent cube's
8 corners at tc_vert_norm_[g]; the weights are constant per vertex, so each
surface vertex is a fixed affine combination of the 8 corners -> it deforms
with the lattice under sphere inversion, exactly like the existing surface
methods. Each vertex is computed once and referenced from all three of its
faces (octagon+hexagon+square, vertex config 4.6.8) -> watertight under
deformation.

Usage:  python3 truncated_cuboctahedron_verify.py
        (prints PASS/FAIL for every check, then the C++ tables)
"""

import math
import itertools

S2 = math.sqrt(2.0)
R  = 1.0 + 2.0 * S2          # cube half-side carrying the octagons
HEX_C = 3.0 + 3.0 * S2       # ±x±y±z plane constant for hexagons
SQ_C  = 2.0 + 3.0 * S2       # ±x±y (etc.) plane constant for squares
A = 1.0 / R                  # normalized "corner-cut" coord  ~0.26120
B = (1.0 + S2) / R           # normalized "edge-cut" coord    ~0.63060
EPS = 1e-9

# ---------------------------------------------------------------------------
# 1. Generate the 48 canonical vertices (perms of (±1, ±(1+√2), ±(1+2√2))).
#    Index order: for each of the 6 permutations of (1, s2, R) (with R the big
#    one), all 8 sign flips.  Deterministic -> stable table output.
# ---------------------------------------------------------------------------
VALS = (1.0, 1.0 + S2, R)    # the three distinct magnitudes

V = []                       # list of (x,y,z)
coord_to_idx = {}
for perm in itertools.permutations(VALS):       # 6 permutations
    for signs in itertools.product((+1, -1), repeat=3):   # 8 sign flips
        p = (signs[0] * perm[0], signs[1] * perm[1], signs[2] * perm[2])
        # round to avoid float-key misses
        key = (round(p[0], 9), round(p[1], 9), round(p[2], 9))
        if key not in coord_to_idx:
            coord_to_idx[key] = len(V)
            V.append(p)
assert len(V) == 48, f"expected 48 vertices, got {len(V)}"


def idx_of(p):
    return coord_to_idx[(round(p[0], 9), round(p[1], 9), round(p[2], 9))]


def vsub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def vadd(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def vscale(a, s):
    return (a[0] * s, a[1] * s, a[2] * s)


def vdot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def vcross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def vlen(a):
    return math.sqrt(vdot(a, a))


def vnorm(a):
    L = vlen(a)
    assert L > EPS, "normalize zero vector"
    return (a[0] / L, a[1] / L, a[2] / L)


def centroid(verts):
    c = (0.0, 0.0, 0.0)
    for p in verts:
        c = vadd(c, p)
    return vscale(c, 1.0 / len(verts))


def order_ccw_outward(verts):
    """Order coplanar verts CCW as seen from outside (polyhedron centered at
    origin -> outward normal = direction of the face centroid)."""
    c = centroid(verts)
    n = vnorm(c)                              # outward normal
    # in-plane basis u, v with u x v = n  -> atan2(p.v, p.u) CCW = outward
    axis = (1.0, 0.0, 0.0)
    if abs(vdot(n, axis)) > 0.9:
        axis = (0.0, 1.0, 0.0)
    u = vnorm(vsub(axis, vscale(n, vdot(n, axis))))   # axis - (axis.n)n
    v = vcross(n, u)                                  # n x u ; |u|=1 -> |v|=1
    def ang(p):
        d = vsub(p, c)
        return math.atan2(vdot(d, v), vdot(d, u))
    return sorted(verts, key=ang)


# ---------------------------------------------------------------------------
# 2. Build the 26 faces by plane filtering.
# ---------------------------------------------------------------------------
oct_faces = []   # list of list of vertex coords (8 each)
hex_faces = []
sq_faces = []

# octagons: coordinate == ±R  (order matches octagonal_faces_: z+,z-,x+,x-,y+,y-)
for (axis, sign) in [(2, +1), (2, -1), (0, +1), (0, -1), (1, +1), (1, -1)]:
    verts = [p for p in V if abs(p[axis] - sign * R) < 1e-6]
    assert len(verts) == 8, f"octagon {axis},{sign} got {len(verts)}"
    oct_faces.append([idx_of(q) for q in order_ccw_outward(verts)])

# hexagons: ±x ± y ± z = 3+3√2  (one per sign triple)
for sx in (+1, -1):
    for sy in (+1, -1):
        for sz in (+1, -1):
            verts = [p for p in V
                     if abs(sx * p[0] + sy * p[1] + sz * p[2] - HEX_C) < 1e-6]
            assert len(verts) == 6, f"hexagon {sx,sy,sz} got {len(verts)}"
            hex_faces.append([idx_of(q) for q in order_ccw_outward(verts)])

# squares: ±x±y, ±x±z, ±y±z = 2+3√2  (3 axis pairs x 4 sign combos = 12)
for (p_axis, q_axis) in [(0, 1), (0, 2), (1, 2)]:
    for sp in (+1, -1):
        for sq in (+1, -1):
            verts = [pt for pt in V
                     if abs(sp * pt[p_axis] + sq * pt[q_axis] - SQ_C) < 1e-6]
            assert len(verts) == 4, f"square {p_axis},{q_axis},{sp},{sq} got {len(verts)}"
            sq_faces.append([idx_of(q) for q in order_ccw_outward(verts)])

# ---------------------------------------------------------------------------
# 3. Verification.
# ---------------------------------------------------------------------------
def fmt(x):
    return f"{x:.10f}"

fails = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        fails.append(msg)


check(len(V) == 48, "48 vertices")
check(len(oct_faces) == 6 and all(len(f) == 8 for f in oct_faces), "6 octagons x 8 verts")
check(len(hex_faces) == 8 and all(len(f) == 6 for f in hex_faces), "8 hexagons x 6 verts")
check(len(sq_faces) == 12 and all(len(f) == 4 for f in sq_faces), "12 squares x 4 verts")
Fpoly = 26

# triangles (fan each face) and mesh topology
tris = []   # (i,j,k) global vertex indices
for f in oct_faces:
    for t in range(1, len(f) - 1):
        tris.append((f[0], f[t], f[t + 1]))
for f in hex_faces:
    for t in range(1, len(f) - 1):
        tris.append((f[0], f[t], f[t + 1]))
for f in sq_faces:
    for t in range(1, len(f) - 1):
        tris.append((f[0], f[t], f[t + 1]))
T = len(tris)
check(T == 92, f"92 solid triangles (got {T})")

# edges: every triangle edge must be shared by exactly 2 (closed manifold)
edge_count = {}
for (a, b, c) in tris:
    for (p, q) in [(a, b), (b, c), (c, a)]:
        e = (min(p, q), max(p, q))
        edge_count[e] = edge_count.get(e, 0) + 1
E = len(edge_count)
shared2 = sum(1 for v in edge_count.values() if v == 2)
check(E == 138, f"138 mesh edges (got {E})")
check(shared2 == E, f"every mesh edge shared by exactly 2 ({shared2}/{E})")
check(48 - E + T == 2, f"Euler V-E+F = 2 (got {48 - E + T})")

# polyhedron edges = face boundaries (consecutive verts, no fan diagonals)
poly_edges = set()
for f in oct_faces + hex_faces + sq_faces:
    n = len(f)
    for t in range(n):
        a, b = f[t], f[(t + 1) % n]
        poly_edges.add((min(a, b), max(a, b)))
check(len(poly_edges) == 72, f"72 polyhedron edges (got {len(poly_edges)})")

# all 72 polyhedron edges equal length (= 2)
edge_lens = []
for (a, b) in poly_edges:
    edge_lens.append(vlen(vsub(V[a], V[b])))
lmin, lmax = min(edge_lens), max(edge_lens)
check(abs(lmin - 2.0) < 1e-9 and abs(lmax - 2.0) < 1e-9,
      f"all 72 polyhedron edges length 2 (min={lmin:.6f} max={lmax:.6f})")

# faces regular: all boundary edges of a face equal, all verts equidistant from
# the face centroid (regular polygon).
def face_regular(face, name):
    pts = [V[i] for i in face]
    c = centroid(pts)
    rad = [vlen(vsub(p, c)) for p in pts]
    sides = [vlen(vsub(pts[t], pts[(t + 1) % len(face)]))
             for t in range(len(face))]
    ok_r = max(rad) - min(rad) < 1e-9
    ok_s = max(sides) - min(sides) < 1e-9
    check(ok_r, f"{name}: vertices equidistant from centroid (regular polygon)")
    check(ok_s, f"{name}: all boundary edges equal")
    return ok_r and ok_s

for fi, f in enumerate(oct_faces):
    face_regular(f, f"octagon {fi}")
for fi, f in enumerate(hex_faces):
    face_regular(f, f"hexagon {fi}")
for fi, f in enumerate(sq_faces):
    face_regular(f, f"square {fi}")

# normals outward: each triangle normal . (tri_centroid - origin) > 0
out_ok = True
for (a, b, c) in tris:
    pa, pb, pc = V[a], V[b], V[c]
    n = vcross(vsub(pb, pa), vsub(pc, pa))
    tc = vscale(vadd(vadd(pa, pb), pc), 1.0 / 3.0)
    if vdot(n, tc) <= EPS:
        out_ok = False
check(out_ok, "all triangle normals outward (normal . tri_centroid > 0)")

print()
print(f"a = {A:.10f}  (= 1/(1+2√2) = (2√2-1)/7)")
print(f"b = {B:.10f}  (= (1+√2)/(1+2√2) = (3+√2)/7)")
print(f"R = {R:.10f}  HEX_C={HEX_C:.6f}  SQ_C={SQ_C:.6f}")
print()
if fails:
    print(f"*** {len(fails)} CHECK(S) FAILED ***")
    for m in fails:
        print("   - " + m)
else:
    print("ALL CHECKS PASSED")
print()

# ---------------------------------------------------------------------------
# 4. Emit the C++ constexpr tables.
# ---------------------------------------------------------------------------
def emit_verts():
    lines = ["    static constexpr double tc_vert_norm_[48][3] = {"]
    for i, p in enumerate(V):
        nx, ny, nz = p[0] / R, p[1] / R, p[2] / R     # normalize to half-side 1
        sep = "," if i < len(V) - 1 else ""
        lines.append(f"        {{ {fmt(nx)}, {fmt(ny)}, {fmt(nz)} }}{sep}  // {i}")
    lines.append("    };")
    return "\n".join(lines)


def emit_int_table(name, faces, per_row):
    lines = [f"    static constexpr int {name}[{len(faces)}][{per_row}] = {{"]
    for i, f in enumerate(faces):
        sep = "," if i < len(faces) - 1 else ""
        body = ", ".join(f"{idx:2d}" for idx in f)
        lines.append(f"        {{ {body} }}{sep}  // face {i}")
    lines.append("    };")
    return "\n".join(lines)


print("// ---- paste into Cube.hpp (private, with the other surface-method tables) ----")
print()
print(emit_verts())
print()
print(emit_int_table("tc_oct_faces_", oct_faces, 8))
print()
print(emit_int_table("tc_hex_faces_", hex_faces, 6))
print()
print(emit_int_table("tc_sq_faces_", sq_faces, 4))
print()

exit_code = 0 if not fails else 1
exit(exit_code)