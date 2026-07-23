# Truncated cuboctahedron mesh (a fourth `Cube` surface method)

This documents the **truncated cuboctahedron** (great rhombicuboctahedron) surface
method added to `Cube.hpp` as a fourth way to build a cube's surface, parallel to
the plain cube, the octagonal (truncated-cube) mesh and the truncated-octahedron
mesh. It is the cantitruncated (omnitruncated) cube, `t_{0,1,2}{4,3}`.

Reference: <https://en.wikipedia.org/wiki/Truncated_cuboctahedron>

## 1. The solid

Archimedean solid, vertex configuration **4.6.8**:

| property | value |
|---|---|
| faces | **26** = 12 squares + 8 hexagons + 6 octagons |
| vertices | 48 |
| edges | 72 |
| Euler V−E+F | 2 |
| vertex figure | 4.6.8 (one square + one hexagon + one octagon meet at each vertex) |
| dihedral 4-6 | arccos(−√6/3) ≈ 144°44′ |
| dihedral 4-8 | arccos(−1/√2) = 135° |
| dihedral 6-8 | arccos(−√3/3) ≈ 125°16′ |

Cantitruncation of the cube: each **cube face → octagon**, each **cube edge →
square**, each **cube corner → hexagon**.

## 2. Canonical model and the cube-derived construction

Canonical coordinates (edge length 2, centred at origin): the 48 vertices are all
permutations of

```
(±1, ±(1+√2), ±(1+2√2))
```

The octagon faces lie on the planes `x = ±R`, `y = ±R`, `z = ±R` with
`R = 1 + 2√2` (the cube half-side whose faces carry the octagons); the hexagon
faces on `±x ± y ± z = 3 + 3√2`; the squares on `±x ± y`, `±x ± z`, `±y ± z =
2 + 3√2`.

Normalising by `R` so the parent cube is half-side 1, the 48 vertices are
**permutations of `(±a, ±b, ±1)`** with

```
a = 1/(1+2√2) = (2√2 − 1)/7 ≈ 0.2612038750   (the "corner-cut" coordinate)
b = (1+√2)/(1+2√2) = (3+√2)/7 ≈ 0.6306019375  (the "edge-cut" coordinate)
```

and `±1` selects which cube face the vertex lies on.

**Inversion compatibility (the design rationale).** Each of the 48 surface
vertices is computed **once** by trilinear interpolation of the cube's 8 corners
at its normalised position:

```
w_i = (1/8) · (1 + s_i.x · x)(1 + s_i.y · y)(1 + s_i.z · z)
vert = Σ_i w_i · corner_i
```

(`s_i` = sign pattern of corner `i` in the standard `initVertices` order 0..7 —
see `tc_corner_sign_` in `Cube.hpp`.) The trilinear weights are constant per
vertex, so every surface vertex is a **fixed affine combination of the 8 cube
corners** — exactly like the plain cube, octagonal and truncated-octahedron
methods. Consequently it deforms with the lattice under sphere inversion and
needs no extra storage. Computing each vertex once and referencing it from all
three of its faces (octagon + hexagon + square, vertex config 4.6.8) makes the
mesh **watertight under deformation** (no per-face cracks).

## 3. Tables and generator (`Cube.hpp`)

The construction is the same "static constexpr table + `static void push*Facets`"
pattern as the other three methods:

| artifact | name | size |
|---|---|---|
| corner sign patterns | `tc_corner_sign_[8][3]` | 8×3 |
| normalised vertex positions | `tc_vert_norm_[48][3]` | 48×3 (perms of (±a,±b,±1)) |
| octagon faces (CCW-outward) | `tc_oct_faces_[6][8]` | 6×8 |
| hexagon faces (CCW-outward) | `tc_hex_faces_[8][6]` | 8×6 |
| square faces (CCW-outward) | `tc_sq_faces_[12][4]` | 12×4 |
| hollow inset default | `default_tc_inset_` | 0.5 |
| generator | `pushTruncatedCuboctahedronFacets(fb, v, hollow, inset)` | — |

The four tables are **generated and verified** by
`truncated_cuboctahedron_verify.py` (see §9) and pasted into `Cube.hpp`; they are
correct by construction, so no runtime winding check is needed and the result does
not depend on `cell.center` (which goes stale under inversion).

The public getter family mirrors the truncated-octahedron family (no morph `s` —
the truncated cuboctahedron is regular only at the fixed `a/b`; the knobs are
`hollow` and `inset`):

- `buildFacetsTruncatedCuboctahedron`
- `getFacetsTruncatedCuboctahedron`
- `getCheckerboardFacetsTruncatedCuboctahedron`
- `getPlaneFacetsTruncatedCuboctahedron`
- `getSubCellFacetsTruncatedCuboctahedron`
- `refreshTriangulationTruncatedCuboctahedron`
- `writeSTL_s_truncated_cuboctahedron` (ASCII STL, modes `full` / `checkerboard` /
  `plane_xy` / `plane_xz` / `plane_yz`)

## 4. Triangle counts

Fan-triangulating each face from vertex 0 (CCW-outward):

| face type | per face | × count | triangles |
|---|---|---|---|
| octagon (8-gon) | 6 fan tris | 6 | 36 |
| hexagon (6-gon) | 4 fan tris | 8 | 32 |
| square (4-gon) | 2 fan tris | 12 | 24 |
| **solid total** | | | **92** |

| property | solid | hollow |
|---|---|---|
| triangles | 92 | 288 (6 oct×16 + 8 hex×12 + 12 sq×8) |
| vertices | 48 | 144 (48 outer + 96 inner: 6×8 + 8×6 + 12×4) |
| mesh edges | 138 (72 poly-edges + 66 fan diagonals) | — |
| Euler V−E+F | 2 | 2 − 26 = −24 (26 hole boundaries) |
| boundary edges | 0 (closed) | 144 (26 holes: 6×8 + 8×6 + 12×4) |

Hollow (mirrors the truncated-octahedron hollow convention): each face is inset
toward its own centroid by `inset` (a homothety, which preserves winding and
planarity) and the inner face is skipped as a hole.

## 5. Winding convention

Every face is emitted **CCW-outward**, fixed by the `tc_*_faces_` tables (computed
directly from the canonical 3-D coordinates by the verifier, which orders each
face's vertices counter-clockwise as seen from outside). No runtime dot-product
check is needed and the result does not depend on `cell.center`. The fan
triangulation (vertex 0 → every other consecutive pair) preserves the outward
winding; the hollow homothety toward each face centroid preserves it too.

## 9. Verified geometry (computational check)

The geometry was verified before/while implementing by
`truncated_cuboctahedron_verify.py` (mirrors the truncated-octahedron workflow in
`truncated_octahedron_mesh.md` §9), then confirmed at runtime by the in-repo
headless test `test_truncated_cuboctahedron.cpp`.

**`truncated_cuboctahedron_verify.py`** (generates + verifies, then emits the
`tc_vert_norm_` / `tc_oct_faces_` / `tc_hex_faces_` / `tc_sq_faces_` tables):

- ✅ 48 vertices, 6 octagons × 8, 8 hexagons × 6, 12 squares × 4 (26 faces)
- ✅ 92 solid triangles, 138 mesh edges, V−E+F = 2
- ✅ every mesh edge shared by exactly 2 triangles (closed manifold)
- ✅ 72 polyhedron edges (the non-coplanar-adjacent mesh edges)
- ✅ all 72 polyhedron edges equal length (= 2 in canonical, = 2r/(1+2√2) for
  cube half-side r)
- ✅ every face regular (vertices equidistant from face centroid; all boundary
  edges of each face equal)
- ✅ all triangle normals outward (normal · tri_centroid > 0)
- ✅ emits the four `static constexpr` tables pasted into `Cube.hpp`

**`test_truncated_cuboctahedron.cpp`** (runtime, headless, no GL):
`g++ -std=c++17 test_truncated_cuboctahedron.cpp -o test_truncated_cuboctahedron -lm && ./test_truncated_cuboctahedron`

- ✅ solid: 92 triangles, 48 vertices, 138 mesh edges, Euler = 2
- ✅ solid: every mesh edge shared by exactly 2 (closed manifold)
- ✅ solid: 72 polyhedron edges, all equal to 2r/(1+2√2) (≈0.5224 at r=1)
- ✅ solid: all triangle normals outward (normal · tri_centroid > 0)
- ✅ hollow: 288 triangles, 144 boundary edges (26 hole loops), no edge shared by >2

**STL facet counts** (`writeSTL_s_truncated_cuboctahedron`):

- ✅ solid: 92 facets
- ✅ hollow: 288 facets

## Build

```
g++ -std=c++17 main19_truncated_cuboctahedron.cpp -o main19_truncated_cuboctahedron -lGL -lGLU -lglut -lm
g++ -std=c++17 test_truncated_cuboctahedron.cpp   -o test_truncated_cuboctahedron   -lm
python3 truncated_cuboctahedron_verify.py
```

`main19_truncated_cuboctahedron` draws one truncated cuboctahedron (hollow/inset
knobs, STL on `c`); `test_truncated_cuboctahedron` is the headless runtime check;
`truncated_cuboctahedron_verify.py` generates + verifies the tables.