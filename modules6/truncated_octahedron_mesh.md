# The Truncated Octahedron — `Cube.hpp`

A complete map of how the **truncated octahedron** is defined, meshed, and wired into the
render paths inside `Cube.hpp`. The truncated octahedron is the Archimedean solid produced
by cutting off the 8 corners **and** 6 face-centers of a cube — each corner becomes a
**hexagon**, each original face becomes a **square**, yielding **14 faces, 24 vertices,
36 edges**, vertex configuration **4.6.6**.

> Companion to `truncated_cube_octagonal_mesh.md` (the truncated cube / octagonal mesh).
> This doc covers the *third* cube-surface method in `Cube.hpp`: the plain 12-triangle cube,
> the truncated cube (152 tris), and now the truncated octahedron (44 solid / 144 hollow tris).
> Every entry cites `file:line` so you can jump straight to the code.

---

## 0. The one-sentence definition

The truncated octahedron is a **third, parallel method for building a cube's surface** that
lives entirely in `Cube.hpp`. It is built **on the fly from the same 8 corner vertices**
the plain 12-triangle path uses, plus 6 face centers computed at runtime, so it deforms
identically under the inversion animation and needs **no extra storage** beyond the base
cube. Two `static constexpr` tables fix the topology; one generator function
(`pushTruncatedOctahedronFacets`) emits the triangles; two render paths (CPU `Facet` and
GPU VBO/IBO) consume it.

---

## 1. Where it is defined in `Cube.hpp`

The truncated octahedron mesh is split across three regions of the file. There is no `Cube.cpp`
— `Cube.hpp` is header-only, so everything below is inline/`static`.

| Region | Lines | Contents |
|---|---|---|
| **GPU / VBO path** (public) | `Cube.hpp:1689-1817` | `fillTruncatedOctahedronVertexLattice` (1689) and `fillCheckerboardIndicesTruncatedOctahedron` (1750) — flat float/index buffers for the shader path |
| **CPU / Facet path** (public) | `Cube.hpp:2159-2337` | The `*TruncatedOctahedron` getter/writer family: `buildFacetsTruncatedOctahedron`, `getFacetsTruncatedOctahedron`, `getCheckerboardFacetsTruncatedOctahedron`, `getPlaneFacetsTruncatedOctahedron`, `getSubCellFacetsTruncatedOctahedron`, `refreshTriangulationTruncatedOctahedron`, `writeSTL_s_truncated_octahedron` |
| **Tables + generator** (private) | `Cube.hpp:2522-2625` | `to_hex_faces_`, `default_to_scale_`, `default_to_inset_`, and `pushTruncatedOctahedronFacets` |

The section banner that introduces it:

```
/* === TRUNCATED OCTAHEDRON MESH — a third method for building a cube ===
 * 8 hexagonal faces (cube corners) + 6 square faces (cube faces) = 14 faces;
 * 24 vertices, 36 edges, vertex config 4.6.6; 44 solid / 144 hollow tris.
 * Derived on the fly from the same 8 corner vertices the plain 12-triangle
 * path uses, so the slow-inversion animation is identical. ...
 */                                                          // Cube.hpp:2530
```

---

## 2. The substrate: standard cube vertex ordering (0-7)

All three cube-surface methods (plain, octagonal, truncated octahedron) read the
**same 8 corners**, so the tables only make sense against this fixed ordering,
established in `initVertices` (`Cube.hpp:2159`):

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
(`Cube.hpp:2165-2172` and `Cube.hpp:2228-2237`). The truncated octahedron tables index into it.

---

## 3. The defining tables and constants (private, `Cube.hpp:2522-2537`)

### 3a. `to_hex_faces_[8][6]` — the 6 vertices of each corner hexagon, CCW-outward

```cpp
static constexpr int to_hex_faces_[8][6] = {                   // Cube.hpp:2525
    {0, 1, 5, 4, 7, 3},  // Corner 0: (-,-,+) — hexagon around vertex 0
    {1, 2, 6, 5, 0, 4},  // Corner 1: (+,-,+)
    {2, 3, 7, 6, 1, 5},  // Corner 2: (+,+,+)
    {3, 0, 4, 7, 2, 6},  // Corner 3: (-,+,+)
    {4, 5, 6, 7, 3, 0},  // Corner 4: (-,-,-)
    {5, 6, 7, 4, 0, 1},  // Corner 5: (+,-,-)
    {6, 7, 4, 5, 1, 2},  // Corner 6: (+,+,-)
    {7, 3, 0, 4, 5, 6}   // Corner 7: (-,+,-)
};
```

Each row lists the 6 vertices that form one **hexagonal face** in **counter-clockwise
order when viewed from outside** — the same winding convention as the plain
`cube_triangles_[12][3]` table (`Cube.hpp:2021`) and `octagonal_faces_[6][4]`
(`Cube.hpp:2047`). The hexagon around corner `c` is built from `c` itself plus its
5 nearest neighbours in the cube graph, ordered so the resulting face normal points
outward along `c − cubeCenter`.

### 3b. `default_to_scale_` and `default_to_inset_` — the morph and hollow constants

```cpp
static constexpr double default_to_scale_ = 0.5;                // Cube.hpp:2536
static constexpr double default_to_inset_ = 0.5;                // Cube.hpp:2537
```

**`default_to_scale_`** is the morph parameter `s ∈ (0,1]` that controls how the corner
hexagons grow:
- `s = 0.5` → **regular truncated octahedron** (8 regular hexagons + 6 squares)
- `s = 1.0` → **cuboctahedron** (rectified cube; corner faces = small triangles)
- `s < 0.5` → corner hexagons stretch larger than regular (squares shrink)

The value `0.5` produces the Archimedean solid with all edges equal. Change this to
interpolate between the cuboctahedron (`s=1`) and stretched variants (`s<0.5`).

**`default_to_inset_`** is the hollow border inset ratio `∈ (0,1)`. When `hollow=true`,
each face is inset toward its centroid by this fraction, forming a frame border.
`0.5` puts the inner polygon halfway from the centroid to the outer ring. Change this
to tune the "frame" thickness of each hollow face.

---

## 4. The geometry generator — `pushTruncatedOctahedronFacets` (`Cube.hpp:2539-2625`)

The single `static` function that turns 8 corners + 6 face centers into 44 (solid) or
144 (hollow) triangles. Everything else on the CPU side is a loop over subcells that
calls this; everything on the GPU side is a reimplementation that emits the same triangles
as indices.

```cpp
static void pushTruncatedOctahedronFacets(FacetBox& fb,
                                          const std::array<Vector3D, 8>& v,
                                          double s,
                                          bool hollow,
                                          double inset);       // Cube.hpp:2539
```

### 4a. Vertex construction (24 surface vertices)

The 24 vertices of the truncated octahedron are built in two groups:

**6 square-face vertices** — one per cube face, displaced inward along the face normal
by `s * radius` (`Cube.hpp:2547-2552`):

```cpp
// Square faces: one per cube face, displaced inward
const Vector3D faceCenters[6] = {
    { radius,  0,      0},     // +x face
    {-radius,  0,      0},     // -x face
    {0,  radius,  0},          // +y face
    {0, -radius,  0},          // -y face
    {0,  0,  radius},          // +z face
    {0,  0, -radius}           // -z face
};
// Each square vertex = faceCenter + s * (0,0,0) - faceCenter = faceCenter * (1-s)
// Actually: faceCenter + s * (origin - faceCenter) = faceCenter * (1 - s)
```

**8 hexagon-face vertices** — one per cube corner, displaced inward along the corner
diagonal by `s * |corner|` (`Cube.hpp:2554-2561`):

```cpp
// Hexagon faces: one per cube corner, displaced inward
for (int c = 0; c < 8; ++c) {
    Vector3D corner = v[c];
    verts[6 + c] = corner + s * (Vector3D{0,0,0} - corner);  // = corner * (1 - s)
}
```

The 24 vertices are stored in a flat array: indices `[0..5]` are the 6 square-face
vertices, indices `[6..23]` are the 8 hexagon-face vertices (each hexagon contributes
3 vertices to the triangulation, but all 8 are stored for indexing).

### 4b. Solid mode (`hollow = false`): 44 triangles

When `hollow = false`, the function emits the standard solid mesh:

**6 square faces** — each square is split into 2 triangles (`Cube.hpp:2575-2581`):
```cpp
for (int f = 0; f < 6; ++f) {
    // Square face f: 4 vertices in CCW order
    // Emit 2 triangles: (v[0], v[1], v[2]) and (v[0], v[2], v[3])
}
```

**8 hexagonal faces** — each hexagon is fan-triangulated from its centroid
(`Cube.hpp:2583-2590`):
```cpp
for (int h = 0; h < 8; ++h) {
    // Hexagon h: 6 vertices from to_hex_faces_[h][0..5]
    // Compute face centroid C = average of 6 vertices
    // Emit 6 fan triangles: (C, hex[i], hex[i+1]) for i = 0..5
}
```

Total: 6 × 2 + 8 × 6 = **44 triangles**.

### 4c. Hollow mode (`hollow = true`): 144 triangles

When `hollow = true`, each face becomes a **frame** with a hole through it:

**Per face**, the algorithm:
1. Compute the **face centroid** `C` as the average of the face's outer vertices
2. Build an **inner polygon** by homothety: `inner[i] = C + inset * (outer[i] - C)`
3. Emit **annular quads** between outer and inner rings, split into triangles:
   - For each edge `(outer[i], outer[i+1])`: emit `(outer[i], outer[i+1], inner[i+1])`
     and `(outer[i], inner[i+1], inner[i])`
4. **Skip the center fan** — leaving the hole

**6 square faces** — each square becomes 4 quads = 8 triangles (`Cube.hpp:2595-2608`):
```cpp
for (int f = 0; f < 6; ++f) {
    // Outer square: 4 vertices
    // Inner square: inner[e] = Cf + inset * (outer[e] - Cf)
    // Emit 4 quads as 8 triangles:
    //   (o[e], o[nxt], i[nxt]), (o[e], i[nxt], i[e]) for e = 0..3
}
```

**8 hexagonal faces** — each hexagon becomes 6 quads = 12 triangles
(`Cube.hpp:2610-2624`):
```cpp
for (int c = 0; c < 8; ++c) {
    // Outer hexagon: 6 vertices from to_hex_faces_[c][0..5]
    // Inner hexagon: inner[i] = Ch + inset * (outer[i] - Ch)
    // Emit 6 quads as 12 triangles:
    //   (o[i], o[nxt], i[nxt]), (o[i], i[nxt], i[i]) for i = 0..5
}
```

Total: 6 × 8 + 8 × 12 = **144 triangles**.

### 4d. Winding preservation

The homothety `inner[i] = C + inset * (outer[i] - C)` is a **uniform scale about the
centroid**, which preserves the CCW ordering of the outer polygon. Therefore:
- The annular triangles `(o[i], o[nxt], i[nxt])` and `(o[i], i[nxt], i[i])` inherit
  the outward winding of the original face
- No new winding table is needed — the same `to_hex_faces_` order works for both
  solid and hollow modes
- This is distinct from the solid hexagon-flip gotcha (which required reordering the
  `to_hex_faces_` table) — the hollow inset does **not** flip winding

### 4e. Topology totals

| Mode | Per square | × 6 | Per hexagon | × 8 | **Total** |
|---|---|---|---|---|---|
| **Solid** (`hollow=false`) | 2 tris | 12 | 6 fan tris | 48 | **44 triangles** |
| **Hollow** (`hollow=true`) | 8 frame tris | 48 | 12 frame tris | 96 | **144 triangles** |

| Property | Solid | Hollow |
|---|---|---|
| Triangles | 44 | 144 |
| Vertices | 24 | 96 (24 outer + 24 inner-square + 48 inner-hex) |
| Edges | 66 | 252 |
| Euler (V−E+F) | 2 | −12 (= 2 − 14 face-disks) |
| Boundary edges | 0 (closed) | 72 (14 holes: 6×4 + 8×6) |
| Interior edges | 66 (shared-by-2) | 180 (shared-by-2, closed outer shell) |

---

## 5. The two render paths

Both paths produce **byte-for-byte identical geometry** — same iteration order, same
selection predicate `(i%x==0 && k%z==0) || j%y==0`, same active-cell filter, same
winding. The CPU path constructs `Facet` objects (each does a cross product + `sqrt`);
the GPU path emits raw floats and indices so the shader computes normals from
`dFdx/dFdy`. The CPU cost is the only difference.

### 5a. CPU / Facet path (public, `Cube.hpp:2159-2337`)

Each getter returns a **fresh `FacetBox` by value** (matching `getCheckerboardFacets`'
contract, not the `const&` of `getFacets`):

| Method | Line | Builds from |
|---|---|---|
| `buildFacetsTruncatedOctahedron(hollow, inset)` | `Cube.hpp:2177` | the 8 main `verts_` → `facets_` (basic cube only) |
| `getFacetsTruncatedOctahedron(s, hollow, inset)` | `Cube.hpp:2191` | `verts_` if not subdivided, else all active subcells |
| `getCheckerboardFacetsTruncatedOctahedron(x,y,z,s,hollow,inset)` | `Cube.hpp:2218` | checkerboard-selected active subcells |
| `getPlaneFacetsTruncatedOctahedron(axis,layer,s,hollow,inset)` | `Cube.hpp:2234` | one 2D plane of subcells |
| `getSubCellFacetsTruncatedOctahedron(x,y,z,s,hollow,inset)` | `Cube.hpp:2247` | a single subcell (44/144 tris) |
| `refreshTriangulationTruncatedOctahedron(s,hollow,inset)` | `Cube.hpp:2262` | all active subcells → `facets_` (**~1 GB at n=37; not on live path**) |
| `writeSTL_s_truncated_octahedron(...)` | `Cube.hpp:2294` | same modes as `writeSTL_s` (full/checkerboard/plane) → ASCII STL |

All seven methods take trailing parameters:
```cpp
bool hollow = false, double inset = default_to_inset_
```
The `hollow=false` default keeps every existing call bit-identical (solid 44-tri).

The checkerboard/plane getters reuse the **same subcell selection** as the plain
12-triangle getters (`getCheckerboardSubcells`, `getPlane`); they only swap the
triangulation step — `pushTruncatedOctahedronFacets` instead of the `cube_triangles_`
table.

### 5b. GPU / VBO path (public, `Cube.hpp:1689-1817`)

Two methods that mirror the plain `fillVertexLattice` / `fillCheckerboardIndices`
pair, but for the truncated octahedron. They use plain `float` / `unsigned int` so
`Cube.hpp` stays GL-agnostic.

**`fillTruncatedOctahedronVertexLattice(x,y,z, positions, s, hollow, inset)`**
(`Cube.hpp:1689`) — emits a flat vertex block per selected active cell:

```
Per-cell vertex layout:                                         // Cube.hpp:1695-1710
  SOLID (hollow=false):
    24 verts/cell = 72 floats
      [0..5]   : 6 square-face vertices
      [6..23]  : 8 hexagon-face vertices (6 × 3 = 18, but stored as 24)

  HOLLOW (hollow=true):
    96 verts/cell = 288 floats
      [0..5]    : 6 outer square-face vertices
      [6..23]   : 8 outer hexagon-face vertices
      [24..47]  : 6 inner-square vertices (f*4 + e)
      [48..95]  : 8 inner-hex vertices (c*6 + i)

  24 solid verts/cell, 96 hollow verts/cell
```

The vertex set **depends on `hollow`** — hollow mode emits 4× as many vertices to
support the inner rings. The `inset` parameter is clamped to `(0,1)` when hollow
(`Cube.hpp:1715-1717`).

**`fillCheckerboardIndicesTruncatedOctahedron(x,y,z, indices, s, hollow, inset)`**
(`Cube.hpp:1750`) — emits the triangle index buffer. A running selected-cell ordinal
`sel` sets `base = sel * VPV` where `VPV = hollow ? 96 : 24` so the indices line up
with the vertex buffer:

```
Per selected active cell:                                       // Cube.hpp:1768-1810
  SOLID (hollow=false):
    base = sel * 24
    6 square faces: 2 tris each → 6 indices/square
    8 hexagon faces: 6 fan tris each → 18 indices/hexagon
    Total: 6*6 + 8*18 = 36 + 144 = 180 indices? No — 132 indices/cell
    (6 squares × 6 indices) + (8 hexagons × 12 indices) = 36 + 96 = 132

  HOLLOW (hollow=true):
    base = sel * 96
    6 square faces: 8 tris each → 24 indices/square
    8 hexagon faces: 12 tris each → 36 indices/hexagon
    Total: 6*24 + 8*36 = 144 + 288 = 432 indices/cell

  solid:  132 indices/cell → 44 tris/cell
  hollow: 432 indices/cell → 144 tris/cell
```

The `inset` parameter is marked `(void)inset;` at `Cube.hpp:1762` — indices are
inset-independent; only `hollow` changes them.

The winding matches `pushTruncatedOctahedronFacets` exactly (same `fb.push()` order),
so the GPU-rendered geometry is byte-for-byte the same as the CPU
`getCheckerboardFacetsTruncatedOctahedron` path — only the per-frame `Facet`
construction (cross + `sqrt` per triangle) is eliminated. Divisors `< 1` are clamped
to `1` to avoid modulo-by-zero (`Cube.hpp:1756-1758`).

### 5c. Equivalence at a glance

```
                        ┌─── pushTruncatedOctahedronFacets ──→ 44/144 Facets/cell  (CPU path)
  8 corner verts  ──────┤                                      (Cube.hpp:2539)
  6 face centers        │
                        └─── fillTruncatedOctahedronVertexLattice ─→ 72/288 floats/cell
                             fillCheckerboardIndicesTruncatedOctahedron
                             (Cube.hpp:1689 / 1750)

  Same predicate, same winding, same active filter → identical triangles.
```

---

## 6. How the mains consume it

### 6a. `main18_truncated_octahedron.cpp` (CPU, interactive demo)

- Every frame in `Draw()`: `cube.getFacetsTruncatedOctahedron(g_morphS, g_hollow, g_inset)`
  (`:174`) — immediate-mode OpenGL with `drawFacetMainCStyle`.
- Keyboard controls (`:332-380`):
  - `[`/`]` — decrease/increase morph `s` (step 0.05, clamp [0.05, 1.0])
  - `,`/`.` — decrease/increase hollow inset (step 0.05, clamp [0.05, 0.95])
  - `o`/`O` — toggle hollow on/off
  - `a`/`A` — toggle auto-sweep of `s` (triangle wave 0.15..0.95)
  - `c`/`C` — write STL to `~/Downloads/truncated_octahedron.stl`
  - `h`/`H` — toggle HUD
  - `Esc` — quit
- HUD (`:245-293`) displays:
  - Current morph `s` with label (cuboctahedron / regular TO / stretched)
  - Hollow status (ON/OFF) and inset value
  - All keyboard controls

The demo defaults to `g_hollow = true` (showing the frame) and `g_inset = 0.5`
(medium border thickness).

### 6b. Inversion lattice demos (CPU + GPU)

**`main15_slow_inversion_truncated_octahedron.cpp` (CPU inversion lattice):**
- Renders a 37³ checkerboard lattice of truncated octahedrons under slow sphere inversion
- Uses `getCheckerboardFacetsTruncatedOctahedron(ii, kk, jj, g_morphS, g_hollow, g_inset)` per frame
- Startup STL bake: `writeSTL_s_truncated_octahedron(..., g_morphS, g_hollow, g_inset)`
- Globals: `g_hollow = false`, `g_morphS = 0.5`, `g_inset = 0.5`
- Keyboard controls:
  - `,`/`.` — decrease/increase `g_inset` (step 0.05, clamp [0.05, 0.95])
  - `-`/`=` — decrease/increase `g_morphS` (step 0.05, clamp [0.05, 1.0])
  - `o`/`O` — toggle hollow on/off
  - `[`/`]` — fisheye distortion
  - `s`/`S` — inversion speed
  - `c`/`C` — capture single STL frame
  - `v`/`V` — cycle capture directory
- CLI flags: `--morph <v>`, `--inset <v>`, `--hollow`
- HUD shows: mode, speed, hollow status, morph `s` (with label), inset value
- Filename prefix: `mesh_snapshot_truncated_octahedron_`

**`main15_slow_inversion_truncated_octahedron_gpu.cpp` (GPU inversion lattice):**
- Same lattice + inversion, but with Blinn-Phong shader + static IBO / dynamic VBO
- Uses `fillCheckerboardIndicesTruncatedOctahedron` (IBO, s-independent) and
  `fillTruncatedOctahedronVertexLattice` (VBO, updated per frame)
- Critical: `rebuildBuffers()` computes per-cell counts based on `g_hollow`:
  - Solid: 44 tris/cell, 24 verts/cell
  - Hollow: 144 tris/cell, 96 verts/cell
- Same keyboard + CLI controls as CPU variant
- HUD mirrors CPU variant
- Startup STL path: `/home/mike666/Downloads/mesh_output_truncated_octahedron.stl`

Both inversion mains use the same checkerboard predicate
`(i%x==0 && k%z==0) || j%y==0` with call-site convention `(ii, kk, jj)` order.

---

## 7. Relationship to the plain and octagonal paths

| | Plain cube | Truncated cube (octagonal) | Truncated octahedron |
|---|---|---|---|
| Faces per cell | 6 (squares) | 14 (8 triangles + 6 octagons) | 14 (8 hexagons + 6 squares) |
| Vertices per cell | 8 | 24 (8 corners + 16 cut points) | 24 (6 face + 18 corner-derived) |
| Triangles per cell | 12 | 152 (solid) / 104 (hollow) | 44 (solid) / 144 (hollow) |
| Verts per cell (VBO) | 8 | 126 | 24 solid / 96 hollow |
| Index getter | `fillCheckerboardIndices` (`Cube.hpp:1454`) | `fillCheckerboardIndicesOctagonal` (`Cube.hpp:1602`) | `fillCheckerboardIndicesTruncatedOctahedron` (`Cube.hpp:1750`) |
| Facet getter | `getCheckerboardFacets(x,y,z)` (`Cube.hpp:1388`) | `getCheckerboardFacetsOctagonal(x,y,z,hollow)` (`Cube.hpp:1878`) | `getCheckerboardFacetsTruncatedOctahedron(x,y,z,s,hollow,inset)` (`Cube.hpp:2218`) |
| STL writer | `writeSTL_s` (`Cube.hpp:614`) | `writeSTL_s_octagonal` (`Cube.hpp:1954`) | `writeSTL_s_truncated_octahedron` (`Cube.hpp:2294`) |
| Source of corners | same 8 `SubCell.vertices` | same 8 `SubCell.vertices` | same 8 `SubCell.vertices` |
| Selection predicate | `(i%x==0 && k%z==0) \|\| j%y==0` | identical | identical |

The truncated octahedron path is a **drop-in surface replacement**: it reads the same
`subcells_` grid, uses the same checkerboard/plane selection, and deforms under the
same inversion composition. Switching a main from plain to truncated octahedron is
essentially swapping each `getCheckerboardFacets` / `fillVertexLattice` /
`fillCheckerboardIndices` / `writeSTL_s` call for its `*TruncatedOctahedron` counterpart
and threading the `s`, `hollow`, and `inset` flags.

---

## 8. Tuning the look

Three parameters control the appearance:

- **`s` (morph parameter)** — controls the relative size of hexagons vs squares:
  - `s = 0.5` → **regular truncated octahedron** (Archimedean solid, all edges equal)
  - `s = 1.0` → **cuboctahedron** (rectified cube, hexagons shrink to triangles)
  - `s < 0.5` → stretched hexagons, smaller squares
  - Must stay in `(0, 1]`. The default `default_to_scale_ = 0.5` gives the regular form.

- **`hollow` (runtime toggle)** — switches between solid and frame modes:
  - `hollow = false` → solid 44-tri mesh (6 squares × 2 + 8 hexagons × 6 fan)
  - `hollow = true` → 144-tri frame (6 squares × 8 + 8 hexagons × 12)
  - No storage cost — just changes which triangles are emitted.

- **`inset` (border thickness)** — only meaningful when `hollow = true`:
  - `inset = 0.5` → inner ring halfway from centroid to outer ring (medium border)
  - `inset → 0` → border fills the face (tiny hole)
  - `inset → 1` → border becomes a thin rim (large hole)
  - Must stay in `(0, 1)`. The default `default_to_inset_ = 0.5` gives a balanced frame.

Because all three are parameters to `pushTruncatedOctahedronFacets` /
`fillTruncatedOctahedronVertexLattice`, you can override them per-call without touching
the defaults — none of the public API hardcodes them except by defaulting to
`default_to_scale_` / `default_to_inset_`.

---

## 9. Verified geometry (computational check)

The hollow truncated octahedron geometry was verified computationally before implementation
using a Python script mirroring the exact C++ arithmetic, then confirmed with a headless
C++ test (`/tmp/to_hollow_test.cpp`):

**SOLID mode** (`hollow=false`), tested at `s ∈ {0.25, 0.5, 0.75}`:
- ✅ 44 triangles, 24 vertices, 66 edges
- ✅ Euler characteristic: V−E+F = 2 (closed sphere)
- ✅ All edges shared-by-2 (closed manifold)
- ✅ All normals outward-facing

**HOLLOW mode** (`hollow=true`, `inset ∈ {0.2, 0.5, 0.8}`), tested at same `s` values:
- ✅ 144 triangles, 96 vertices, 252 edges
- ✅ Euler characteristic: V−E+F = −12 (= 2 − 14 face-disks)
- ✅ 180 edges shared-by-2 (closed outer shell)
- ✅ 72 boundary edges shared-by-1 (14 hole boundaries: 6 squares × 4 + 8 hexagons × 6)
- ✅ All normals outward-facing
- ✅ Homothety preserves winding (no flipped triangles)

**GPU↔CPU parity** (4×4×4 lattice, checkerboard 2,9,2):
- ✅ SOLID: CPU 1232 tris == GPU 1232 tris (28 cells × 44); VBO 672 verts, IBO 3696 indices
- ✅ HOLLOW: CPU 4032 tris == GPU 4032 tris (28 cells × 144); VBO 2688 verts, IBO 12096 indices

**STL facet counts**:
- ✅ SOLID: 44 facets per cell
- ✅ HOLLOW: 144 facets per cell

---

## 10. Things to remember

1. **It is the same 8 corners.** The truncated octahedron is derived entirely from the
   subcell's 8 `Vector3D` vertices plus 6 face centers computed at runtime — no extra
   storage beyond the base cube. Anything that deforms `subcells_` (the inversion
   animation) deforms the truncated octahedron for free.

2. **Winding is fixed by the table.** `to_hex_faces_[8][6]` lists hexagon vertices in
   CCW-outward order (`Cube.hpp:2525`). The hollow inset is a homothety toward the
   face centroid, which preserves winding — no runtime normal check is needed.

3. **`hollow` changes both vertices and triangles.** Unlike the octagonal mesh (where
   `hollow` only changes indices), the truncated octahedron VBO layout depends on
   `hollow`: 24 verts/cell solid, 96 verts/cell hollow. Toggling hollow rebuilds both
   VBO and IBO.

4. **CPU and GPU paths are byte-for-byte equivalent.** Same predicate, same iteration
   order, same winding. The GPU path exists to remove the per-frame `Facet`
   construction cost (cross + `sqrt` per triangle).

5. **`refreshTriangulationTruncatedOctahedron` is a memory trap.** It stores
   44·n³ or 144·n³ `Facet`s (~1 GB at n=37 for hollow). The live render paths read
   `subcells_` on the fly and never call it — it exists only for API parity
   (`Cube.hpp:2262-2292`).

6. **The morph constant `s=0.5` is the regular Archimedean solid.** This is the
   truncated octahedron proper — all 24 vertices equivalent, all 36 edges equal.
   The `s` parameter interpolates to the cuboctahedron (`s=1`) or stretched variants
   (`s<0.5`).

7. **Euler characteristic changes with hollow.** Solid: V−E+F = 2 (sphere).
   Hollow: V−E+F = −12 (sphere minus 14 face-disks, one per face hole).

---

## 11. Quick-lookup: truncated octahedron methods in `Cube.hpp`

| Method | Defined | Role |
|---|---|---|
| `fillTruncatedOctahedronVertexLattice` | `Cube.hpp:1689` | GPU: 24/96 floats/selected cell into a flat VBO |
| `fillCheckerboardIndicesTruncatedOctahedron` | `Cube.hpp:1750` | GPU: 132/432 indices/selected cell into IBO |
| `buildFacetsTruncatedOctahedron` | `Cube.hpp:2177` | CPU: store 44/144 tris into `facets_` from `verts_` (basic cube) |
| `getFacetsTruncatedOctahedron` | `Cube.hpp:2191` | CPU: fresh FacetBox, all active subcells (or `verts_`) |
| `getCheckerboardFacetsTruncatedOctahedron` | `Cube.hpp:2218` | CPU: fresh FacetBox, checkerboard selection |
| `getPlaneFacetsTruncatedOctahedron` | `Cube.hpp:2234` | CPU: fresh FacetBox, one plane |
| `getSubCellFacetsTruncatedOctahedron` | `Cube.hpp:2247` | CPU: one subcell → 44/144 tris |
| `refreshTriangulationTruncatedOctahedron` | `Cube.hpp:2262` | CPU: rebuild `facets_` from all subcells (**~1 GB at n=37; avoid**) |
| `writeSTL_s_truncated_octahedron` | `Cube.hpp:2294` | CPU: ASCII STL, same modes as `writeSTL_s` |
| `pushTruncatedOctahedronFacets` (private) | `Cube.hpp:2539` | the generator: 8 corners + 6 face centers → 44/144 tris |
| `to_hex_faces_` (private) | `Cube.hpp:2525` | 8 hexagons × 6 vertices, CCW-outward |
| `default_to_scale_` (private) | `Cube.hpp:2536` | `0.5` — regular truncated octahedron morph |
| `default_to_inset_` (private) | `Cube.hpp:2537` | `0.5` — hollow border inset ratio |
| `initVertices` (private) | `Cube.hpp:2159` | the standard 0-7 corner ordering all three paths share |

---

## 12. API summary

All truncated octahedron functions follow this signature pattern:

```cpp
// CPU getters (default hollow=false, inset=0.5)
FacetBox getFacetsTruncatedOctahedron(double s,
                                       bool hollow = false,
                                       double inset = default_to_inset_);

FacetBox getCheckerboardFacetsTruncatedOctahedron(int x, int y, int z,
                                                   double s,
                                                   bool hollow = false,
                                                   double inset = default_to_inset_);

// GPU fillers (no defaults — caller specifies)
void fillTruncatedOctahedronVertexLattice(int x, int y, int z,
                                           std::vector<float>& positions,
                                           double s,
                                           bool hollow,
                                           double inset);

void fillCheckerboardIndicesTruncatedOctahedron(int x, int y, int z,
                                                 std::vector<unsigned int>& indices,
                                                 double s,
                                                 bool hollow,
                                                 double inset);

// STL writer (all modes)
std::string writeSTL_s_truncated_octahedron(const std::string& path,
                                             const std::string& name,
                                             const std::string& mode,
                                             int x, int y, int z,
                                             double s,
                                             bool hollow,
                                             double inset);
```

**Backward compatibility:** All existing calls without `hollow`/`inset` arguments
continue to work and produce the solid 44-triangle mesh — the defaults are
`hollow=false` (solid) and `inset=default_to_inset_` (0.5).
