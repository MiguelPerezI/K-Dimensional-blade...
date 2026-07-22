# The Truncated Cube (Octagonal Mesh) — `Cube.hpp`

A complete map of how the **truncated cube** is defined, meshed, and wired into the
render paths inside `Cube.hpp`. The code calls it the *octagonal mesh*; geometrically
it is the Archimedean solid produced by cutting off the 8 corners of a cube — each of
the 6 square faces becomes an **octagon**, and the 8 cut corners become **8 triangles**.

> Companion to `cube_lineage_main15_stl_capture.md` (the plain 12-triangle cube) and
> `inversions_main15_stl_capture.md` (the sphere-inversion deformation). This doc zooms
> in on the *second* cube-surface method that sits alongside the plain one in `Cube.hpp`.
> Every entry cites `file:line` so you can jump straight to the code.

---

## 0. The one-sentence definition

The truncated cube is a **second, parallel method for building a cube's surface** that
lives entirely in `Cube.hpp`. It is built **on the fly from the same 8 corner vertices**
the plain 12-triangle path uses, so it deforms identically under the inversion animation
and needs **no extra storage**. Three `static constexpr` tables fix the topology; one
generator function (`pushOctagonalCubeFacets`) emits the triangles; two render paths
(CPU `Facet` and GPU VBO/IBO) consume it.

---

## 1. Where it is defined in `Cube.hpp`

The octagonal mesh is split across three regions of the file. There is no `Cube.cpp`
— `Cube.hpp` is header-only, so everything below is inline/`static`.

| Region | Lines | Contents |
|---|---|---|
| **GPU / VBO path** (public) | `Cube.hpp:1482-1661` | `fillOctagonalVertexLattice` (1511) and `fillCheckerboardIndicesOctagonal` (1602) — flat float/index buffers for the shader path |
| **CPU / Facet path** (public) | `Cube.hpp:1819-2003` | The `*Octagonal` getter/writer family: `buildFacetsOctagonal`, `getFacetsOctagonal`, `getCheckerboardFacetsOctagonal`, `getPlaneFacetsOctagonal`, `getSubCellFacetsOctagonal`, `refreshTriangulationOctagonal`, `writeSTL_s_octagonal` |
| **Tables + generator** (private) | `Cube.hpp:2036-2144` | `octagonal_faces_`, `corner_neighbors_`, `default_trunc_`, `default_inner_scale_`, and `pushOctagonalCubeFacets` |

The section banner that introduces it:

```
/* === OCTAGONAL (TRUNCATED-CUBE) MESH — a second method for building a cube ===
 * Each face becomes an octagon (outer ring + inner ring + center fan = 24
 * triangles) and 8 corner triangles close the solid => 152 triangles/cube.
 * Derived on the fly from the same 8 corner vertices the plain 12-triangle
 * path uses, so the slow-inversion animation is identical. ...
 */                                                          // Cube.hpp:1819
```

---

## 2. The substrate: standard cube vertex ordering (0-7)

Both cube-surface methods (plain and octagonal) read the **same 8 corners**, so the
tables only make sense against this fixed ordering, established in `initVertices`
(`Cube.hpp:2159`):

```
       3----------2            0: (x-, y-, z+)  front-bottom-left
      /|         /|            1: (x+, y-, z+)  front-bottom-right
     / |        / |            2: (x+, y+, z+)  front-top-right
    7----------6  |     Z      3: (x-, y+, z+)  front-top-left
    |  |       |  |     |      4: (x-, y-, z-)  back-bottom-left
    |  0-------|--1     |      5: (x+, y-, z-)  back-bottom-right
    | /        | /      |      6: (x+, y+, z-)  back-top-right
    |/         |/       +-----Y 7: (x-, y+, z-)  back-top-left
    4----------5       /
                      X
```

`verts_` (the main cube) and every `SubCell.vertices` use this exact order
(`Cube.hpp:2165-2172` and `Cube.hpp:2228-2237`). The octagonal tables index into it.

---

## 3. The three defining tables (private, `Cube.hpp:2036-2076`)

### 3a. `octagonal_faces_[6][4]` — the 4 corners of each face, CCW-outward

```cpp
static constexpr int octagonal_faces_[6][4] = {              // Cube.hpp:2047
    {0, 1, 2, 3},  // Front  (z = +radius, outward +z)
    {4, 7, 6, 5},  // Back   (z = -radius, outward -z)
    {1, 5, 6, 2},  // Right  (x = +radius, outward +x)
    {0, 3, 7, 4},  // Left   (x = -radius, outward -x)
    {3, 2, 6, 7},  // Top    (y = +radius, outward +y)
    {0, 4, 5, 1}   // Bottom (y = -radius, outward -y)
};
```

Each row is the four corner indices of one face in **counter-clockwise order when
viewed from outside the cube** — the same winding convention as the plain
`cube_triangles_[12][3]` table (`Cube.hpp:2021`). This ordering is what makes every
generated triangle's normal point outward with no runtime dot-product check.

### 3b. `corner_neighbors_[8][3]` — each corner's 3 edge-neighbours, outward winding

```cpp
static constexpr int corner_neighbors_[8][3] = {             // Cube.hpp:2060
    {1, 3, 4},  // 0: (-,-,+)
    {0, 5, 2},  // 1: (+,-,+)
    {3, 1, 6},  // 2: (+,+,+)   <-- NOT {1,3,6}; that winding is inward
    {2, 7, 0},  // 3: (-,+,+)
    {5, 0, 7},  // 4: (-,-,-)
    {4, 6, 1},  // 5: (+,-,-)
    {7, 2, 5},  // 6: (+,+,-)
    {6, 4, 3}   // 7: (-,+,-)
};
```

For corner `c` with position `Vc`, its corner triangle is built from three **cut
points** along the three edges leaving `Vc`:

```
Pa = Vc + t·(v[n[0]] - Vc)
Pb = Vc + t·(v[n[1]] - Vc)
Pc = Vc + t·(v[n[2]] - Vc)
triangle = (Pa, Pb, Pc)                                         // Cube.hpp:2137-2143
```

The neighbour **order** `n[0], n[1], n[2]` is chosen so that `(Pa, Pb, Pc)` has an
**outward** normal (along `Vc − cubeCenter`). This is hardcoded into the table, with an
explicit warning on corner 2 that `{1,3,6}` would be inward. Fixing the winding in the
table means:

- no runtime normal-direction check is needed, and
- the result does **not** depend on `cell.center`, which goes stale (wrong) once the
  inversion animation deforms the lattice.

### 3c. `default_trunc_` and `default_inner_scale_` — the regular-octagon constants

```cpp
static constexpr double default_trunc_       = 0.2928932188134524;  // = 1/(2+sqrt(2))   Cube.hpp:2075
static constexpr double default_inner_scale_ = 0.5;                                     // Cube.hpp:2076
```

**`default_trunc_`** is the corner-cut fraction `t ∈ (0, 0.5)`. The specific value
`1/(2+√2)` is the one that makes the octagon **regular**: it equalizes the 4 long outer
edges (the un-cut middles of the original square edges, length `2r·(1−2t)`) with the 4
short cut-corner edges (the new edges across each cut, length `2r·t·√2`):

```
2r·(1 − 2t) = 2r·t·√2   →   t = 1 / (2 + √2) ≈ 0.2928932188
```

It is hardcoded as a literal because `std::sqrt` is not `constexpr` until C++26 and
`M_SQRT2` is not guaranteed under strict `-std=c++17`. Change this to tune how deeply
the corners are cut (smaller = barely chamfered, 0.5 = cuts meet at edge midpoints).

**`default_inner_scale_`** shrinks the inner octagon ring toward the face centroid by a
homothety (a uniform scale about the centroid), so the inner octagon stays regular and
concentric with the outer one. `0.5` puts the inner ring halfway from the centroid to
the outer ring. Change this to tune the "frame" thickness of each face.

---

## 4. The geometry generator — `pushOctagonalCubeFacets` (`Cube.hpp:2096-2144`)

The single `static` function that turns 8 corners into 152 triangles. Everything else
on the CPU side is a loop over subcells that calls this; everything on the GPU side is
a reimplementation that emits the same triangles as indices.

```cpp
static void pushOctagonalCubeFacets(FacetBox& fb,
                                    const std::array<Vector3D, 8>& v,
                                    double trunc,
                                    double innerScale,
                                    bool hollow);                // Cube.hpp:2096
```

### 4a. Per face (6 faces × 24 triangles = 144)

For each face `f` with corners `P[0..3] = v[octagonal_faces_[f][0..3]]`:

**Step 1 — face centroid:**
```
C = (P0 + P1 + P2 + P3) / 4                                    // Cube.hpp:2104
```

**Step 2 — outer octagon ring** in CCW-outward order
`[F0, B1, F1, B2, F2, B3, F3, B0]` (`Cube.hpp:2109-2114`):
```
F[i]   = P[i]   + t·(P[i+1] − P[i])   // forward along edge i → i+1   (cut near P[i])
B[i+1] = P[i+1] + t·(P[i]   − P[i+1]) // backward along edge i+1 → i  (cut near P[i+1])
```

Geometrically, each original corner `P[i]` is replaced by two cut points — `F[i]` on
the edge to `P[i+1]` and `B[i]` on the edge to `P[i−1]` — and the corner triangle
(§4c) fills the gap between them. Between `F[i]` and `B[i+1]` lies a **long** octagon
edge (the un-cut middle of the original square edge); between `B[i+1]` and `F[i+1]`
lies a **short** cut-corner edge. With `t = 1/(2+√2)` these alternate equal lengths →
regular octagon.

```
        B3 ————— F3          (outer octagon, one face)
       /          \           F[i] = cut near corner P[i], toward P[i+1]
      F2           B0         B[i] = cut near corner P[i], toward P[i-1]
      |    ● C    |           C    = face centroid (center-fan apex)
      B2           F0         long edges:  F0-B1, F1-B2, F2-B3, F3-B0
       \          /           short edges: B0-F1, B1-F2, B2-F3, B3-F0
        F1 ————— B1
```

**Step 3 — inner ring** = homothety of the outer ring toward `C` (`Cube.hpp:2116-2118`):
```
inner[i] = C + innerScale·(outer[i] − C)
```
A uniform scale about `C`, so the inner octagon is regular and concentric.

**Step 4 — 16 annular triangles** between the outer and inner rings (`Cube.hpp:2121-2125`):
```
for i in 0..8 (j = i+1 mod 8):
    push(outer[i], outer[j], inner[j])    // quad (outer[i],outer[j],inner[j],inner[i])
    push(outer[i], inner[j], inner[i])    // split into 2 triangles, outward
```
8 octagon edges × 2 triangles = **16 annular triangles per face**.

**Step 5 — 8 center-fan triangles** filling the inner octagon (`Cube.hpp:2129-2134`):
```
if (!hollow)
    for i in 0..8 (j = i+1 mod 8):
        push(C, inner[i], inner[j])       // fan from centroid to each inner edge
```
**Skipped when `hollow = true`** — leaving an octagonal hole through each face so the
truncated cube is hollow (you can see through the face centers).

### 4b. Per corner (8 corner triangles)

```cpp
for (int c = 0; c < 8; ++c) {                                  // Cube.hpp:2137-2143
    const int* n = corner_neighbors_[c];
    Vector3D Vc = v[c];
    fb.push(Vc + trunc * (v[n[0]] - Vc),
            Vc + trunc * (v[n[1]] - Vc),
            Vc + trunc * (v[n[2]] - Vc));
}
```

One triangle per truncated corner, closing the solid. Winding is fixed by
`corner_neighbors_` (§3b) so the normal points outward along `Vc − cubeCenter`.

### 4c. Topology totals

| Mode | Per face | × 6 faces | + 8 corners | **Per cube** |
|---|---|---|---|---|
| **Solid** (`hollow=false`) | 16 annular + 8 fan = 24 | 144 | 8 | **152 triangles** |
| **Hollow** (`hollow=true`) | 16 annular + 0 fan = 16 | 96 | 8 | **104 triangles** |

Normals are computed **eagerly** by the `Facet` constructor from the current (possibly
deformed) vertex positions (`Facet.cpp:27-33`, via cross product + `unit()`), so the
mesh shades correctly even after the inversion animation bends it.

---

## 5. The two render paths

Both paths produce **byte-for-byte identical geometry** — same iteration order, same
selection predicate `(i%x==0 && k%z==0) || j%y==0`, same active-cell filter, same
winding. The CPU path constructs `Facet` objects (each does a cross product + `sqrt`);
the GPU path emits raw floats and indices so the shader computes normals from
`dFdx/dFdy`. The CPU cost is the only difference.

### 5a. CPU / Facet path (public, `Cube.hpp:1819-2003`)

Each getter returns a **fresh `FacetBox` by value** (matching `getCheckerboardFacets`'
contract, not the `const&` of `getFacets`):

| Method | Line | Builds from |
|---|---|---|
| `buildFacetsOctagonal(hollow)` | `Cube.hpp:1837` | the 8 main `verts_` → `facets_` (basic cube only) |
| `getFacetsOctagonal(hollow)` | `Cube.hpp:1851` | `verts_` if not subdivided, else all active subcells |
| `getCheckerboardFacetsOctagonal(x,y,z,hollow)` | `Cube.hpp:1878` | checkerboard-selected active subcells |
| `getPlaneFacetsOctagonal(axis,layer,hollow)` | `Cube.hpp:1894` | one 2D plane of subcells |
| `getSubCellFacetsOctagonal(x,y,z,hollow)` | `Cube.hpp:1907` | a single subcell (152 tris) |
| `refreshTriangulationOctagonal(hollow)` | `Cube.hpp:1922` | all active subcells → `facets_` (**~1 GB at n=37; not on live path**) |
| `writeSTL_s_octagonal(...)` | `Cube.hpp:1954` | same modes as `writeSTL_s` (full/checkerboard/plane) → ASCII STL |

The checkerboard/plane/octagonal getters reuse the **same subcell selection** as the
plain 12-triangle getters (`getCheckerboardSubcells`, `getPlane`); they only swap the
triangulation step — `pushOctagonalCubeFacets` instead of the `cube_triangles_` table.

### 5b. GPU / VBO path (public, `Cube.hpp:1482-1661`)

Two methods that mirror the plain `fillVertexLattice` / `fillCheckerboardIndices`
pair, but for the truncated mesh. They use plain `float` / `unsigned int` so `Cube.hpp`
stays GL-agnostic.

**`fillOctagonalVertexLattice(x,y,z, positions, hollow)`** (`Cube.hpp:1511`) — emits a
flat **126-vertex block per selected active cell** into `positions`:

```
Per-cell layout (126 verts = 378 floats):                     // Cube.hpp:1495-1500
  6 faces, 17 verts/face at offset f*17:
     outer[0..7]  @ f*17 + 0..7     (8 outer octagon verts)
     inner[0..7]  @ f*17 + 8..15    (8 inner octagon verts)
     centroid C   @ f*17 + 16       (1 face centroid)
  8 corners, 3 cut-points each at offset 102 + c*3:
     cp[0..2]     @ 102 + c*3       (the 3 corner-triangle verts)

  6*17 + 8*3 = 102 + 24 = 126 verts/cell
```

The vertex set is **independent of `hollow`** — `hollow` only changes which *triangles*
are emitted (the index buffer), not which *vertices* exist. So the VBO is allocated
once and never resized when toggling hollow; only the IBO rebuilds. (`(void)hollow;`
at `Cube.hpp:1579` makes this explicit.)

**`fillCheckerboardIndicesOctagonal(x,y,z, indices, hollow)`** (`Cube.hpp:1602`) — emits
the triangle index buffer. A running selected-cell ordinal `sel` sets
`base = sel * 126` so the indices line up with the vertex buffer produced by
`fillOctagonalVertexLattice` (same iteration order, same predicate, same active filter):

```
Per selected active cell:                                     // Cube.hpp:1626-1656
  6 faces, fb = base + f*17:
    16 annular indices   (8 edges × 2 tris)        (outer/inner at fb+0..15)
    if !hollow: 8 center-fan indices               (centroid at fb+16)
  8 corner triangles, cb = base + 102 + c*3:       (cut-points at cb+0..2)
    push(cb+0, cb+1, cb+2)

  solid:  6*(16+8) + 8 = 152 tris → 456 indices/cell
  hollow: 6*16       + 8 = 104 tris → 312 indices/cell
```

The winding matches `pushOctagonalCubeFacets` exactly (same `fb.push()` order), so the
GPU-rendered geometry is byte-for-byte the same as the CPU `getCheckerboardFacetsOctagonal`
path — only the per-frame `Facet` construction (cross + `sqrt` per triangle) is
eliminated. Divisors `< 1` are clamped to `1` to avoid modulo-by-zero (`Cube.hpp:1608-1610`).

### 5c. Equivalence at a glance

```
                        ┌─── pushOctagonalCubeFacets ───→ 152 Facets/cell  (CPU path)
  8 corner verts  ──────┤                                 (Cube.hpp:2096)
  (deformed)             │
                        └─── fillOctagonalVertexLattice  ──→ 126 floats/cell  (GPU path)
                             fillCheckerboardIndicesOctagonal → 456/312 ints/cell
                             (Cube.hpp:1511 / 1602)

  Same predicate, same winding, same active filter → identical triangles.
```

---

## 6. How the mains consume it

Both truncated-cube programs reuse the **same inversion pipeline** documented in
`inversions_main15_stl_capture.md` (`g_identityPositions → composeSigma →
subcells_[i][j][k].vertices[l]`); they only swap the surface method.

**`main15_slow_inversion_truncated_cube.cpp` (CPU):**
- Every frame in `Draw()`: `cube.getCheckerboardFacetsOctagonal(ii, kk, jj, g_hollow)`
  (`:407`) → immediate-mode OpenGL.
- Startup bake: `cube.writeSTL_s_octagonal(...)` (`:395`).
- `c`-key snapshot: `getCheckerboardFacetsOctagonal` again (`:1065, :1142`) so the STL
  matches the on-screen render triangle-for-triangle.

**`main15_slow_inversion_truncated_cube_gpu.cpp` (GPU):**
- Once: `fillCheckerboardIndicesOctagonal` builds the **static IBO** (`:829`); the IBO
  only rebuilds when `ii/jj/kk` or `--hollow` changes.
- Every frame: `fillOctagonalVertexLattice` refills the **dynamic VBO** with deformed
  vertices (`:859`) — the VBO is re-uploaded each frame because the geometry deforms.
- Startup bake: `writeSTL_s_octagonal` (`:566`).
- `c`-key snapshot still uses the CPU `getCheckerboardFacetsOctagonal` Facet path
  (`:1360, :1437`) so the snapshot matches the screen.

The `ii, kk, jj` argument order (note: not `ii, jj, kk`) is the same historical quirk
as the plain path — the predicate arguments are passed as `(ii, kk, jj)` in code.

---

## 7. Relationship to the plain 12-triangle path

| | Plain cube | Truncated cube (octagonal) |
|---|---|---|
| Triangulation table | `cube_triangles_[12][3]` (`Cube.hpp:2021`) | `octagonal_faces_[6][4]` + `corner_neighbors_[8][3]` (`Cube.hpp:2047, 2060`) |
| Triangles per cell | 12 | 152 (solid) / 104 (hollow) |
| Verts per cell (VBO) | 8 (`fillVertexLattice`) | 126 (`fillOctagonalVertexLattice`) |
| Index getter | `fillCheckerboardIndices` (`Cube.hpp:1454`) | `fillCheckerboardIndicesOctagonal` (`Cube.hpp:1602`) |
| Facet getter | `getCheckerboardFacets(x,y,z)` (`Cube.hpp:1388`) | `getCheckerboardFacetsOctagonal(x,y,z,hollow)` (`Cube.hpp:1878`) |
| STL writer | `writeSTL_s` (`Cube.hpp:614`) | `writeSTL_s_octagonal` (`Cube.hpp:1954`) |
| Source of corners | same 8 `SubCell.vertices` | same 8 `SubCell.vertices` |
| Selection predicate | `(i%x==0 && k%z==0) \|\| j%y==0` | identical |

The octagonal path is a **drop-in surface replacement**: it reads the same `subcells_`
grid, uses the same checkerboard/plane selection, and deforms under the same inversion
composition. Switching a main from plain to truncated is essentially swapping each
`getCheckerboardFacets` / `fillVertexLattice` / `fillCheckerboardIndices` / `writeSTL_s`
call for its `*Octagonal` counterpart and threading the `hollow` flag.

---

## 8. Tuning the look

Two `constexpr` in `Cube.hpp` control the appearance (`Cube.hpp:2075-2076`):

- **`default_trunc_`** (`1/(2+√2) ≈ 0.2929`) — corner-cut depth. The current value
  yields a *regular* octagon. Lower → corners barely chamfered (approaches a plain
  cube); higher → deeper cuts (at `0.5` the cuts meet at edge midpoints and the "face"
  collapses). Must stay in `(0, 0.5)`.
- **`default_inner_scale_`** (`0.5`) — inner ring size as a fraction of the outer ring
  relative to the face centroid. Lower → thicker annular "frame" and smaller center
  fan; higher → thinner frame, bigger fan. Must stay in `(0, 1)`. Only affects the
  solid look (the center fan exists only when `!hollow`).

`hollow` (a runtime flag, not a constant) drops the 8 center-fan triangles per face,
leaving an octagonal hole — 104 tris/cell instead of 152. Toggle it with `--hollow` on
the command line or the `O` key at runtime.

Because both constants are `constexpr` passed into `pushOctagonalCubeFacets` /
`fillOctagonalVertexLattice` as parameters, you can also override them per-call without
touching the defaults — none of the public API hardcodes them except by defaulting to
`default_trunc_` / `default_inner_scale_`.

---

## 9. Things to remember

1. **It is the same 8 corners.** The octagonal mesh is derived entirely from the
   subcell's 8 `Vector3D` vertices — no extra storage, no extra capture. Anything that
   deforms `subcells_` (the inversion animation) deforms the truncated cube for free.
2. **Winding is fixed in the tables, not computed.** `octagonal_faces_` (CCW-outward)
   and `corner_neighbors_` (outward corner triangles) make every normal point outward
   with no runtime dot-check. Corner 2's neighbour order is deliberately `{3,1,6}`,
   not `{1,3,6}` — the latter is inward (`Cube.hpp:2063`).
3. **`hollow` changes triangles, not vertices.** The 126-vertex VBO layout is identical
   solid or hollow; only the IBO (456 vs 312 indices/cell) changes. So toggling hollow
   rebuilds the IBO but not the VBO.
4. **CPU and GPU paths are byte-for-byte equivalent.** Same predicate, same iteration
   order, same winding. The GPU path exists to remove the per-frame `Facet`
   construction cost (cross + `sqrt` per triangle), not to change the geometry.
5. **`refreshTriangulationOctagonal` is a memory trap.** It stores 152·n³ `Facet`s
   (~1 GB at n=37). The live render paths read `subcells_` on the fly and never call it
   — it exists only for API parity with `refreshTriangulation` (`Cube.hpp:1916-1921`).
6. **`cell.center` is not trusted.** The corner-triangle winding relies on the table,
   not on `cell.center`, specifically because `cell.center` goes stale under inversion.
   Don't "simplify" the corner generation by computing a normal from `cell.center`.
7. **The regular-octagon constant is hardcoded** because `std::sqrt` isn't `constexpr`
   under `-std=c++17`. If you change `default_trunc_`, recompute the literal yourself;
   don't write `1.0/(2.0+std::sqrt(2.0))` — it won't compile as `constexpr`.

---

## 10. Quick-lookup: octagonal methods in `Cube.hpp`

| Method | Defined | Role |
|---|---|---|
| `fillOctagonalVertexLattice` | `Cube.hpp:1511` | GPU: 126 floats/selected cell into a flat VBO |
| `fillCheckerboardIndicesOctagonal` | `Cube.hpp:1602` | GPU: 456/312 indices/selected cell into IBO |
| `buildFacetsOctagonal` | `Cube.hpp:1837` | CPU: store 152 tris into `facets_` from `verts_` (basic cube) |
| `getFacetsOctagonal` | `Cube.hpp:1851` | CPU: fresh FacetBox, all active subcells (or `verts_`) |
| `getCheckerboardFacetsOctagonal` | `Cube.hpp:1878` | CPU: fresh FacetBox, checkerboard selection |
| `getPlaneFacetsOctagonal` | `Cube.hpp:1894` | CPU: fresh FacetBox, one plane |
| `getSubCellFacetsOctagonal` | `Cube.hpp:1907` | CPU: one subcell → 152 tris |
| `refreshTriangulationOctagonal` | `Cube.hpp:1922` | CPU: rebuild `facets_` from all subcells (**~1 GB at n=37; avoid**) |
| `writeSTL_s_octagonal` | `Cube.hpp:1954` | CPU: ASCII STL, same modes as `writeSTL_s` |
| `pushOctagonalCubeFacets` (private) | `Cube.hpp:2096` | the generator: 8 corners → 152/104 tris |
| `octagonal_faces_` (private) | `Cube.hpp:2047` | 6 faces × 4 corners, CCW-outward |
| `corner_neighbors_` (private) | `Cube.hpp:2060` | 8 corners × 3 neighbours, outward winding |
| `default_trunc_` (private) | `Cube.hpp:2075` | `1/(2+√2)` — regular-octagon cut fraction |
| `default_inner_scale_` (private) | `Cube.hpp:2076` | `0.5` — inner ring homothety toward centroid |
| `initVertices` (private) | `Cube.hpp:2159` | the standard 0-7 corner ordering both paths share |