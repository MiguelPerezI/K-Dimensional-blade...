# Inversions — `main15_stl_capture.cpp`

A focused map of the **cube lattice** construction and the **animated sphere-inversion
composition** that deforms it each frame, plus recipes for extending the composition to
more inversion passes. Every entry cites `file:line` so you can jump straight to the code.

> Companion to `cube_lineage_main15_stl_capture.md`, which covers the broader
> DEFINE → CREATE → STORE → DRAW lineage. This doc zooms in on §3 (STORE) of that
> file — the inversion math itself — and how to extend it.

---

## 1. Where the cube lattice is built

Two lattices are built: the **live subdivided lattice** (inside `Cube`, what the render
reads and the animation mutates) and a **pristine flat copy** the animation recomposes
from each frame so deformation never compounds in place.

### 1a. Live subdivided lattice — global construction

```cpp
int cube_dim = 2.0;                              // main15_stl_capture.cpp:152  (half-side / radius)
int N = 37;                                      //                         :153  (subdivision levels)
Cube cube(cube_dim, Vector3D{0,0.0001,0}, N);    //                         :154
```

The 3-arg constructor `Cube.hpp:299` runs `initVertices` → `buildFacets` →
`subdivide(37)` → `buildSubdividedCube` (`Cube.hpp:1731`), which fills `subcells_` — a
37×37×37 grid where each `SubCell` holds 8 `Vector3D` vertices (`Cube.hpp:238-251`).

```
Cube(radius, center, 37)                          Cube.hpp:299
 └─ subdivide(37)                                  Cube.hpp:896
     └─ buildSubdividedCube(center, radius, 37)    Cube.hpp:1731
          ├─ step = 2·radius/37 ; subRadius = step/2
          ├─ subcells_.resize(37×37×37)
          └─ for i,j,k in [0,37):
               ├─ subCenter = center - radius + step·(i+0.5, j+0.5, k+0.5)
               └─ subcells_[i][j][k] = SubCell(subVerts[0..7], subCenter, subRadius, i,j,k)
```

### 1b. Pristine identity lattice — captured once

```cpp
static std::vector<float> g_identityPositions;   // main15_stl_capture.cpp:214  (flat n³·8·3 floats)

static void captureIdentityLattice() {           //                         :239
    cube.fillVertexLattice(g_identityPositions); // Cube.hpp:1415
}
```

Called from `Setup()` at `:314`. This flat vector is the **source of truth** the
animation reads each frame; the live `subcells_` vertices are overwritten, never
read-back-and-redeformed, so there is no frame-over-frame compounding.

```
g_identityPositions ──fillVertexLattice──► captured once
        │
        ▼
   read every frame in updateAnimatedGeometry, composed via composeSigma,
   pushed into subcells_[i][j][k].vertices[l]
```

---

## 2. How the inversions are applied

The pipeline runs every frame in `updateAnimatedGeometry()` (`main15_stl_capture.cpp:249`),
driven by a ~30 Hz timer (`animTimer` `:269`, which advances `g_animTime` and is dispatched
through `display()`):

```cpp
SigmaParams sp = currentSigmaParams(g_animTime);                       // :250
...
cube.updateSubCellVertex(lx, ly, lz, l, composeSigma(p, sp));         // :261
```

### 2a. The two inversion spheres — `currentSigmaParams(t)` `:220`

```cpp
static SigmaParams currentSigmaParams(double t) {
    SigmaParams s;
    s.c1 = Vector3D(0.0, 0.0, 0.1 * sin(t));   // pass 1: center oscillates in z, r = 0.5
    s.r1 = 0.5;
    s.c2 = cube.getSubcellCenter(1, 1, 1);     // pass 2: fixed at subcell (1,1,1)
    s.r2 = cube.getSubcellRadius(1, 1, 1);
    return s;
}
```

| Pass | Center | Radius | Behavior |
|---|---|---|---|
| 1 | `(0, 0, 0.1·sin t)` | `0.5` | oscillates along z with the animation clock |
| 2 | subcell (1,1,1) center | subcell (1,1,1) radius (`radius·√3`) | fixed in space |

### 2b. The math — `composeSigma(p, s)` `:231`

A **two-pass sphere inversion**. Each pass is the classic inversion
`a + (r² / |p−a|²) · (p−a)`, with a singular-center guard that leaves the point untouched
when it coincides with a center (the free `sigma()` in `Vector3D.cpp:273` *throws* in that
case; here we guard so the animation never crashes):

```cpp
static inline Vector3D composeSigma(const Vector3D& p, const SigmaParams& s) {
    Vector3D d1 = p - s.c1;  double ds1 = d1 * d1;
    Vector3D q  = (ds1 < 1e-12) ? p : s.c1 + (s.r1 * s.r1 / ds1) * d1;   // pass 1
    Vector3D d2 = q - s.c2;  double ds2 = d2 * d2;
    return (ds2 < 1e-12) ? q : s.c2 + (s.r2 * s.r2 / ds2) * d2;         // pass 2
}
```

The underlying single inversion (for reference):

```cpp
// Vector3D.cpp:273 — the free sigma() (throws on singular center; NOT used by the live path)
Vector3D sigma(const Vector3D& x, const Vector3D& a, double r) {
    Vector3D diff = x - a;
    double dist_squared = diff * diff;
    if (abs(dist_squared) < 1e-12)
        throw std::runtime_error("sigma: point coincides with sphere center");
    double factor = r * r / dist_squared;
    return a + factor * diff;
}
```

### 2c. Storage flow each frame — `updateAnimatedGeometry()` `:249`

```
updateAnimatedGeometry()
 ├─ sp = currentSigmaParams(g_animTime)              :250
 ├─ n = cube.getSubdivisionLevels(); center = n/2    :251-252
 ├─ for i,j,k in [0,n):  for l in 0..8:               :255-264
 │     p = g_identityPositions[idx…]
 │     cube.updateSubCellVertex(lx, ly, lz, l, composeSigma(p, sp))   Cube.hpp:1167
 │         └─ cell.vertices[vertex_idx] = new_pos     ← the actual store
 └─ idx += 3
```

**Net flow:** `g_identityPositions → composeSigma → subcells_[i][j][k].vertices[l]`,
recomputed every frame from the pristine copy. The render (`getCheckerboardFacets`,
`Draw()` `:328`) and the `c`-key STL snapshot (`captureSTLSnapshot()` `:863`) read
`subcells_` live, so whatever `composeSigma` produces shows up on screen and in the
snapshot with no further wiring.

> **Dead alternative.** `applySigmaTransformationToCube()` (`:171`) does a *single*
> `sigma()` pass and then calls `refreshTriangulation()`. It is **not on the live path**
> (never called here), and the live render bypasses `refreshTriangulation()` entirely
> (it reads `subcells_` directly). Do **not** extend that helper to add inversions —
> extend the `composeSigma` / `SigmaParams` pair described in §3.

---

## 3. How to add more inversions

The composition is hard-coded to 2 passes. **Three places must change in lockstep.**
Nothing downstream (`updateAnimatedGeometry`, the render, the STL snapshot) needs
editing — they all read whatever `composeSigma` returns.

### 3a. Recipe: add a 3rd inversion pass

**(1) `SigmaParams` struct** (`:212`) — add a field:

```cpp
struct SigmaParams { Vector3D c1; double r1; Vector3D c2; double r2;
                     Vector3D c3; double r3; };
```

**(2) `currentSigmaParams(t)`** (`:220`) — set the new center/radius (animate with `t`
if you want it to move, or anchor it to a subcell like pass 2):

```cpp
s.c3 = Vector3D(0.2 * cos(t), 0.0, 0.0);   // oscillates in x
s.r3 = 0.35;
// — or, fixed to a subcell: —
// s.c3 = cube.getSubcellCenter(1, -1, 1);
// s.r3 = cube.getSubcellRadius(1, -1, 1);
```

**(3) `composeSigma(p, s)`** (`:231`) — chain one more pass on the result:

```cpp
static inline Vector3D composeSigma(const Vector3D& p, const SigmaParams& s) {
    Vector3D d1 = p   - s.c1;  double ds1 = d1 * d1;
    Vector3D q  = (ds1 < 1e-12) ? p   : s.c1 + (s.r1 * s.r1 / ds1) * d1;   // pass 1
    Vector3D d2 = q   - s.c2;  double ds2 = d2 * d2;
    Vector3D r2v = (ds2 < 1e-12) ? q   : s.c2 + (s.r2 * s.r2 / ds2) * d2;  // pass 2
    Vector3D d3 = r2v - s.c3;  double ds3 = d3 * d3;
    return (ds3 < 1e-12) ? r2v : s.c3 + (s.r3 * s.r3 / ds3) * d3;          // pass 3
}
```

### 3b. Generalize to N inversions (cleaner if you'll add more than 3)

Swap the fixed struct for vectors and loop in `composeSigma`. The same downstream
call sites work unchanged.

```cpp
// main15_stl_capture.cpp:212  — replace the struct
struct SigmaParams {
    std::vector<Vector3D> c;   // inversion centers
    std::vector<double>   r;   // inversion radii
};

// :220 — push as many (center, radius) pairs as you want
static SigmaParams currentSigmaParams(double t) {
    SigmaParams s;
    s.c.push_back(Vector3D(0.0, 0.0, 0.1 * sin(t)));  s.r.push_back(0.5);
    s.c.push_back(cube.getSubcellCenter(1, 1, 1));    s.r.push_back(cube.getSubcellRadius(1, 1, 1));
    s.c.push_back(Vector3D(0.2 * cos(t), 0.0, 0.0));   s.r.push_back(0.35);
    // ... add more here
    return s;
}

// :231 — one guard per pass, applied in order
static inline Vector3D composeSigma(Vector3D p, const SigmaParams& s) {
    for (size_t i = 0; i < s.c.size(); ++i) {
        Vector3D d = p - s.c[i];  double ds = d * d;
        p = (ds < 1e-12) ? p : s.c[i] + (s.r[i] * s.r[i] / ds) * d;
    }
    return p;
}
```

### 3c. Things to keep in mind

- **Order matters.** Sphere inversions do not commute — `σ₂∘σ₁ ≠ σ₁∘σ₂` in general.
  The sequence in `currentSigmaParams` *is* the composition, so reordering the pushes
  changes the deformation.
- **The singular guard is per-pass.** Keep the `(ds < 1e-12) ? p : …` form on every
  pass; without it a vertex landing exactly on a center produces a NaN/Inf and the
  mesh explodes. The free `sigma()` (`Vector3D.cpp:273`) throws instead — fine for
  one-shot use, wrong for a 60 fps animation.
- **No retriangulation needed.** The live path reads `subcells_` directly via
  `getCheckerboardFacets`; adding passes does not require touching `refreshTriangulation`
  or the STL writers. New passes appear on screen and in `c`-key snapshots automatically.
- **Radius sign / magnitude is the lever.** Small radii near the lattice produce tight
  local folds; large radii produce gentle global bends. Negative radii compose an
  inversion with a reflection.

---

## 4. Quick reference

| Symbol | Location | Role |
|---|---|---|
| `cube` global | `main:154` | the live 37³ subdivided lattice |
| `g_identityPositions` | `main:214` | pristine flat copy, captured once |
| `captureIdentityLattice()` | `main:239` | one-time capture (called from `Setup()` `:314`) |
| `SigmaParams` struct | `main:212` | the set of inversion spheres |
| `currentSigmaParams(t)` | `main:220` | builds the spheres (animated w/ `t`) |
| `composeSigma(p, s)` | `main:231` | applies the passes — **extend here to add inversions** |
| `updateAnimatedGeometry()` | `main:249` | per-frame: read identity → compose → write subcells |
| `animTimer()` | `main:269` | ~30 Hz clock advancing `g_animTime` |
| `Draw()` | `main:321` | reads `subcells_` via `getCheckerboardFacets(ii,kk,jj)` |
| `captureSTLSnapshot()` | `main:863` | `c` key: same selection → STL file |
| free `sigma()` | `Vector3D.cpp:273` | single-pass reference (throws; not on live path) |
| `applySigmaTransformationToCube()` | `main:171` | **dead** single-pass helper — do not extend |