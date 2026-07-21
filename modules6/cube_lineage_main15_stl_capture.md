# Cube class lineage — `main15_stl_capture.cpp`

A reference for how `main15_stl_capture.cpp` **defines, creates, stores, and draws** its cube.
Use this as a map when working with the `Cube` class elsewhere — every entry cites
`file:line` so you can jump straight to the code.

> Scope: the live render path in this file. A legacy helper
> (`applySigmaTransformationToCube`, `main15_stl_capture.cpp:171`) is noted where
> relevant but is **not on the live path** — it's defined but never called here.

---

## 0. The one global cube

```cpp
int cube_dim = 2.0;                              // main15_stl_capture.cpp:152  (radius / half-side)
int N = 37;                                      //                         :153  (subdivision levels)
Cube cube(cube_dim, Vector3D{0,0.0001,0}, N);    //                         :154
FacetBox plane_subcells;                         //                         :155  (render scratch)
```

Includes that bring in the types (note `.cpp`s are included directly):

```cpp
#include "Vector3D.cpp"     // main15_stl_capture.cpp:6
#include "Quaternion.cpp"   //                         :8
#include "Facet.cpp"        //                         :9
#include "FacetBox.hpp"     //                         :10
#include "Cube.hpp"          //                         :12  (header-only — no Cube.cpp)
```

---

## 1. DEFINE — the type hierarchy

Everything bottoms out in three classes, one per include:

```
Vector3D   (Vector3D.hpp/.cpp)   — 3 doubles; atomic point + algebra
   ▲
   │ stored inside
Quaternion (Quaternion.hpp/.cpp) — 4-tuple; Cube/Facet use it as a point carrier (.V() → Vector3D)
   ▲
   │ stored inside
Facet      (Facet.hpp/.cpp)      — triangle: 3 Quaternion verts A,B,C + normal N
   ▲
   │ collected in
FacetBox   (FacetBox.hpp)        — std::vector<Facet> wrapper
   ▲
   │ stored inside
Cube       (Cube.hpp)            — verts_ + facets_ + subcells_ grid
```

### `Cube` internal storage (`Cube.hpp:1641-1648`)

| Member | Type | Role |
|---|---|---|
| `verts_` | `std::vector<Quaternion>` | 8 main corners (basic cube only) |
| `facets_` | `FacetBox` | triangulated faces (n³·12 after subdivision) — **not used by the live render path** |
| `subcells_` | `vector<vector<vector<SubCell>>>` | n³ grid — **the live data source for the animated render** |
| `subdivision_levels_` | `int` | n |
| `subdivision_built_` | `bool` | subdivision valid flag |
| `cube_triangles_[12][3]` | `static constexpr int` | lookup table: 8 ordered verts → 12 triangles. Every triangulation path funnels through this. |

### `Cube::SubCell` (`Cube.hpp:238-251`)

```cpp
struct SubCell {
    std::array<Vector3D, 8> vertices;  // 8 cube vertices, standard order 0-7
    Vector3D center;                   // subcell center in world space
    double   radius;                  // half side length (face→center, not corner)
    int i, j, k;                       // physical grid coords [0,n)
    bool active = true;                // skip in triangulation when false
};
```

This is what the animation **mutates** and what the render **reads**.

### File-scope storage outside the Cube

```cpp
static std::vector<float> g_identityPositions;  // main15_stl_capture.cpp:214
```

The **pristine identity lattice** — captured once, the source of truth the animation
recomposes from each frame (so deformation is never applied in-place → no frame-over-frame compounding).

---

## 2. CREATE — construction lineage

`Cube cube(cube_dim, Vector3D{…}, N)` → 3-arg constructor `Cube.hpp:299`:

```
Cube(radius, center, subdivisions)                 Cube.hpp:299
 ├─ initVertices(radius, center)                   Cube.hpp:1682
 │    └─ emplace_back 8 Quaternion corners into verts_   (standard ordering 0-7)
 ├─ buildFacets()                                  Cube.hpp:1706
 │    └─ for t in 0..12: Facet(verts_[cube_triangles_[t]]) → facets_.push()  Cube.hpp:1654
 └─ subdivide(37)                                   Cube.hpp:896
      ├─ subdivision_levels_ = 37
      ├─ facets_.clear()
      ├─ buildSubdividedCube(center, radius, 37)   Cube.hpp:1731
      │    ├─ step = 2·radius/37 ; subRadius = step/2
      │    ├─ subcells_.resize(37×37×37)
      │    └─ for i,j,k in [0,37):
      │         ├─ subCenter = center - radius + step·(i+0.5, j+0.5, k+0.5)
      │         ├─ subVerts[0..7] = 8 corners around subCenter
      │         ├─ subcells_[i][j][k] = SubCell(subVerts, subCenter, subRadius, i,j,k)
      │         └─ for t in 0..12: facets_.push(subVerts[tri…])   Cube.hpp:1654
      └─ subdivision_built_ = true
```

At construction the cube has: 8 `verts_`, 37³ `SubCell`s in `subcells_` (each with 8 live
`Vector3D` vertices), and 37³·12 `Facet`s in `facets_`. The animation later overwrites the
`subcells_` vertices; `facets_` is **not** used by the live render path (see §4).

---

## 3. STORE — how vertices get written (animation lineage)

Two writers feed `subcells_`; one reader-writer captures the identity lattice.

### Capture (once, in `Setup()` → `captureIdentityLattice()` `main:239`)

```
cube.fillVertexLattice(g_identityPositions)        Cube.hpp:1415
 └─ for i,j,k in [0,n): for l in 0..8:
        push cell.vertices[l].x/y/z as float        (n³·8·3 floats)
```

### Animate (every frame, `display()` → `updateAnimatedGeometry()` `main:249`, called at `main:763`)

```
updateAnimatedGeometry()
 ├─ currentSigmaParams(g_animTime)                 main:220
 │    ├─ cube.getSubcellCenter(1,1,1)              Cube.hpp:1593 → getSubCell → logicalToPhysical (956)
 │    └─ cube.getSubcellRadius(1,1,1)             Cube.hpp:1612 → cell.radius·√3
 ├─ for i,j,k in [0,n): for l in 0..8:
 │    p = g_identityPositions[idx…]
 │    cube.updateSubCellVertex(lx,ly,lz,l, composeSigma(p, sp))   Cube.hpp:1167
 │         ├─ getSubCellMutable(x,y,z)            Cube.hpp:1061 → logicalToPhysical → subcells_[i][j][k]
 │         └─ cell.vertices[vertex_idx] = new_pos   ← the actual store
 └─ composeSigma(p, sp)                           main:231  (two-pass sphere inversion, singular-center guarded)
      uses Vector3D operator-, operator* (dot), operator+, scalar operator*   Vector3D.cpp:66-118
      (mirrors the free sigma() at Vector3D.cpp:273, but never throws)
```

**Storage flow:** `g_identityPositions → composeSigma → subcells_[i][j][k].vertices[l]`,
recomputed every frame from the pristine copy.

> The live path writes `subcells_` directly and **never calls `refreshTriangulation()`**.
> That method (`Cube.hpp:1186`) only appears inside the dead helper
> `applySigmaTransformationToCube` (`main:171`).

---

## 4. DRAW — the render lineage (screen + STL)

The render reads `subcells_` live; it does **not** read `facets_`. One call site, two sinks.

```
Draw() / captureSTLSnapshot()
 └─ cube.getCheckerboardFacets(ii, kk, jj)        Cube.hpp:1388  (FacetBox return overload)
      └─ getCheckerboardSubcells(x,y,z)           Cube.hpp:1321
           predicate: (i%x==0 && k%z==0) || j%y==0   over subcells_[i][j][k]
      └─ for each active cell: for t in 0..12:
           checkerboard_facets.push(cell.vertices[tri…])   FacetBox.hpp:159
              └─ facets_.emplace_back(A,B,C) → Facet(Vector3D,…)  Facet.cpp:27
                   └─ normal N = unit(edge1 % edge2)            Facet.cpp:30-33
                        ├─ operator% (cross)                  Vector3D.cpp:113
                        └─ unit() → abs()/=                    Vector3D.cpp:148,127
```

The returned `FacetBox` is consumed by one of two sinks.

### To screen (`Draw()` `main:321`, every frame via `ProcessingProto()` `main:335` → `display()` `main:745`)

```
for i in plane_subcells:
   drawFacetMainCStyle(plane_subcells[i], i+500)    main:427
    ├─ f.getNormal()      → Facet.hpp:38 (N.V())
    ├─ f[0],f[1],f[2]     → Facet.cpp:40 (A.V/B.V/C.V)
    └─ OpenGL: glBegin(GL_TRIANGLES) + glNormal3f + glVertex3f
              glBegin(GL_LINE_LOOP)  …outline…
```

### To STL file — two flavors

- **Once at startup** (`Setup()` `main:316`):
  `cube.writeSTL_s(path, "MyCube", "checkerboard", ii, kk, jj)`
  → 6-arg overload `Cube.hpp:614` → `getCheckerboardFacets(x,y,z)` (same path as above) → ASCII STL loop.
- **On `c` key** (`captureSTLSnapshot()` `main:863` → `keyboard` `main:918`):
  `cube.getCheckerboardFacets(ii,kk,jj)` again, then the **local** writer `writeFacetBoxSTL()`
  `main:830` (a free function, not a Cube method) iterating `facets[i]`, `f.getNormal()`, `f[0..2]`.
  Reuses the identical `getCheckerboardFacets` call on the current animated `subcells_`, so the STL
  matches the OpenGL render triangle-for-triangle.

---

## 5. At a glance

```
                          ┌──────── Cube.hpp (header-only) ────────┐
   define:    Vector3D ─► Quaternion ─► Facet ─► FacetBox ─► Cube(verts_, facets_, subcells_, cube_triangles_[])
                          └────────────────────────────────────────┘

   create:    Cube(r,c,n) → initVertices → buildFacets → subdivide → buildSubdividedCube
                                                              (Cube.hpp:299→896→1731)

   store:     g_identityPositions ──fillVertexLattice──► (read once)   Cube.hpp:1415
                  │
                  ▼
             composeSigma ──updateSubCellVertex──► subcells_[i][j][k].vertices[l]   (every frame)
                                                              Cube.hpp:1167

   draw:      subcells_ ──getCheckerboardFacets──► FacetBox ──┬── drawFacetMainCStyle ──► OpenGL  (screen)
                          Cube.hpp:1388                        └── writeSTL_s / writeFacetBoxSTL ──► .stl (file)
```

---

## 6. Two things to remember

1. **The live render/snapshot path is `subcells_ → getCheckerboardFacets` and bypasses
   `facets_` entirely.** So `refreshTriangulation()` (the old way to refresh `facets_`
   after vertex edits) is dead code in this file; the animation just overwrites subcell
   vertices and the next `getCheckerboardFacets` reads them live.
2. **`sigma` exists in three forms:**
   - free `sigma()` — `Vector3D.cpp:273` (throws on coincidence with center)
   - `Facet::applySigma` / `Facet::sigma` — `Facet.cpp:175 / 189` (used by `FacetBox::applySigma`)
   - inlined `composeSigma` — `main:231` (singular-center guarded; what the animation actually uses)

---

## 7. Quick-lookup: Cube methods used by this file

| Method | Defined | Called from (main15_stl_capture.cpp) |
|---|---|---|
| `Cube(r, c, n)` constructor | `Cube.hpp:299` | `:154` (global) |
| `initVertices` (private) | `Cube.hpp:1682` | ctor |
| `buildFacets` (private) | `Cube.hpp:1706` | ctor |
| `subdivide` | `Cube.hpp:896` | ctor |
| `buildSubdividedCube` (private) | `Cube.hpp:1731` | `subdivide` |
| `cube_triangles_` table | `Cube.hpp:1654` | every triangulation path |
| `hasSubdivision` | `Cube.hpp:922` | guards in helpers |
| `getSubdivisionLevels` | `Cube.hpp:934` | `updateAnimatedGeometry` `:251` |
| `logicalToPhysical` | `Cube.hpp:956` | `getSubCell*`, `updateSubCellVertex` |
| `getSubCell` | `Cube.hpp:1010` | `applySigma…` `:183` (legacy) |
| `getSubCellMutable` | `Cube.hpp:1061` | `updateSubCellVertex` |
| `updateSubCellVertex` | `Cube.hpp:1167` | `updateAnimatedGeometry` `:261`, legacy `:185` |
| `getSubcellCenter` | `Cube.hpp:1593` | `currentSigmaParams` `:224` |
| `getSubcellRadius` | `Cube.hpp:1612` | `currentSigmaParams` `:225` |
| `fillVertexLattice` | `Cube.hpp:1415` | `captureIdentityLattice` `:240` |
| `getCheckerboardSubcells(x,y,z)` | `Cube.hpp:1321` | via `getCheckerboardFacets` |
| `getCheckerboardFacets(x,y,z)` | `Cube.hpp:1388` | `Draw()` `:328`, `captureSTLSnapshot` `:888` |
| `refreshTriangulation` | `Cube.hpp:1186` | legacy helper `:192` only — **dead on live path** |
| `writeSTL_s` (6-arg) | `Cube.hpp:614` | `Setup()` `:316` |
| `Facet::operator[]` | `Facet.cpp:40` | `drawFacetMainCStyle`, `writeFacetBoxSTL` |
| `Facet::getNormal` | `Facet.hpp:38` | `drawFacetMainCStyle`, `writeFacetBoxSTL` |
| `Facet(Vector3D,Vector3D,Vector3D)` | `Facet.cpp:27` | `FacetBox::push` |
| `FacetBox::push(3 verts)` | `FacetBox.hpp:159` | `getCheckerboardFacets` |
| `FacetBox::operator[]` / `size()` | `FacetBox.hpp:138 / 204` | `Draw()`, `writeFacetBoxSTL` |
| `sigma` (free) | `Vector3D.cpp:273` | via `Facet::applySigma` |
| `operator%` (cross) | `Vector3D.cpp:113` | `Facet` normal |
| `unit` | `Vector3D.cpp:148` | `Facet` normal |
| `operator*` (dot) | `Vector3D.cpp:107` | `composeSigma`, `sigma` |

---

## 8. Coordinate-system note (gotcha)

`getSubCell`, `getSubCellMutable`, `getSubcellCenter`, `getSubcellRadius`,
`updateSubCellVertex`, `setSubCellActive` take **logical** coords `(x,y,z)` where
`(0,0,0)` is the center cell — they internally call `logicalToPhysical` (`Cube.hpp:956`,
`center = n/2`, `i = center + x`). `getSubCellPhysical` / `getSubCellMutablePhysical`
take **physical** `[0,n)` indices directly.

`fillVertexLattice`, `getCheckerboardSubcells`, `getCheckerboardFacets`, and
`refreshTriangulation` iterate **physical** `[0,n)³` internally — no logical conversion.
The animation mirrors this: `updateAnimatedGeometry` (`main:249`) walks physical `[0,n)³`
and converts to logical only to call `updateSubCellVertex` (`lx = i - center`, etc.).