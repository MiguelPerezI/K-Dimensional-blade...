# 4D-Menger-Sponge
The fourth dimensional menger sponge, enjoy...

//LINUX commands
c++ main.cpp -o go -lGL -lGLU -lglut
./go

// To compile and run main9.cpp (Hyperbolic Dodecahedron):
c++ main9.cpp -o main9 -lGL -lGLU -lglut
./main9

// To compile and run main10.cpp (Cube + Dodecahedron Demo):
c++ main10.cpp -o main10 -lGL -lGLU -lglut
./main10

// To compile and run main11.cpp (Advanced Cube Subdivision Demo with NEW String-Based Coordinates):
g++ -std=c++17 main11.cpp -o main11 -lGL -lGLU -lglut
./main11



e, E    : move camera up and down
r, R    : rotate
v       : 4D rotation
i       : order 2
I       : order 3
P       : up one floor
p       : down one floor
f       : zoom out
F       : zoom In

---

## 🚀 GPU-Accelerated Build — `main15_gpu` (NVIDIA Tesla T4)

`main15_gpu.cpp` is a GPU-accelerated version of `main15.cpp`. It renders the same
checkerboard-reflected subdivided cube, with the same look and the same controls,
but on the **NVIDIA Tesla T4** instead of the CPU — which makes it dramatically
faster. The original `main15.cpp` is left untouched as a reference/baseline.

### Why the original was slow (root cause)

The original `main15` was slow for **two compounding reasons**, both verified on
this server:

1. **It never used the GPU.** Its OpenGL context rendered through Mesa's
   `llvmpipe` **software rasterizer** (CPU), not the installed Tesla T4. The T4
   was fully usable — the app simply never selected the NVIDIA GLX backend.
   Confirmed: querying the active renderer returned
   `VENDOR=Mesa, RENDERER=llvmpipe`; forcing `__GLX_VENDOR_LIBRARY_NAME=nvidia`
   switched the same program to `RENDERER=Tesla T4/PCIe/SSE2, OpenGL 4.6.0`.
2. **Pathological immediate-mode rendering.** Every frame the CPU rebuilt
   ~112,860 `Facet` objects (each does a cross product + `sqrt`) and issued
   **~225K `glBegin/glEnd` batches** (fill + outline, two batches per triangle).
   There were no VBOs, no display lists, no shaders. The 37³ cube geometry is
   built **once** in `Setup()` and never changes per frame; only the
   checkerboard *selection* (`ii/jj/kk`) changes, and only on key press.

`main15_gpu` fixes both:

- **Selects the NVIDIA GLX backend** via `setenv("__GLX_VENDOR_LIBRARY_NAME",
  "nvidia",1)` before `glutInit` (and via the `run_gpu.sh` launcher), so the T4
  rasterizes instead of llvmpipe.
- **Replaces the immediate-mode loop** with a static VBO (uploaded once), an
  index buffer rebuilt **only when `ii/jj/kk` change**, and a `#version 460 core`
  GLSL shader that computes per-triangle face normals on the GPU via
  `dFdx/dFdy`. The ~225K `glBegin/glEnd` batches collapse to **2
  `glDrawElements` calls** per frame (one fill pass, one outline pass).

### What it does (feature parity with `main15`)

- Builds one `Cube` of radius `2.0` subdivided `37×37×37` (`50,653` subcells).
- Applies two sphere-inversion (`sigma`) reflections in `Setup()` (one centered
  at the origin, one at subcell `(1,1,1)`) — identical to `main15`.
- Renders the **parameterized checkerboard selection** `(i%ii==0 && k%jj==0)
  || j%kk==0` of active subcells, each as 12 triangles, with gray lit fills and
  black triangle outlines (same look as `main15`).
- Same interactive camera (mouse rotate/pan/zoom, wheel zoom), same keyboard
  controls, same right-click menu, same fisheye "projection" toggle, same HUD.
- Writes an STL export of the selected mesh to
  `/home/mike666/Downloads/mesh_output_gpu.stl` (one-time, in `Setup()`).

### How to execute

```bash
cd /home/mike666/K-Dimensional-blade.../modules6

# 1. Compile (pure g++ + GL libs; no CUDA, no GLEW needed)
g++ -o main15_gpu main15_gpu.cpp -lGL -lGLU -lglut -lm

# 2. Run on the Tesla T4
./run_gpu.sh
#   run_gpu.sh is equivalent to:
#   DISPLAY=:0.0 __GLX_VENDOR_LIBRARY_NAME=nvidia ./main15_gpu
```

On startup the program prints which renderer it actually got, e.g.:

```
[main15_gpu] GL VENDOR   = NVIDIA Corporation
[main15_gpu] GL RENDERER = Tesla T4/PCIe/SSE2
[main15_gpu] GL VERSION  = 4.6.0 NVIDIA 580.173.02
[main15_gpu] VBO uploaded: 405224 vertices (4.64 MB)
[main15_gpu] IBO rebuilt: 112860 triangles (ii,jj,kk=4,5,9)
```

If you ever see `RENDERER = llvmpipe`, the GPU was not selected — run it through
`./run_gpu.sh` (or export `__GLX_VENDOR_LIBRARY_NAME=nvidia` yourself). The
in-program `setenv` means the binary self-selects the T4 even without the
launcher, but `run_gpu.sh` is the supported, belt-and-suspenders path.

> The original baseline still works for comparison:
> ```bash
> g++ -o main15 main15.cpp -lGL -lGLU -lglut -lm && ./main15   # still llvmpipe/CPU
> ```

### What to install before use

The program needs a working NVIDIA OpenGL stack and FreeGLUT. On this server
(Ubuntu 24.04, Tesla T4, driver 580.173.02, CUDA 13.0) everything is already
present. On a fresh machine, install:

```bash
# NVIDIA driver + userspace GL (provides libGLX_nvidia.so, the GL 4.6 backend)
sudo apt update
sudo apt install nvidia-driver-580     # or match your installed kernel module

# OpenGL / GLUT development headers and runtime libs
sudo apt install libgl1-mesa-dev libglu1-mesa-dev freeglut3-dev

# (Optional) only to *inspect* which GL renderer is active — not required to run
sudo apt install mesa-utils            # provides glxinfo
```

Notes:
- **No GLEW / GLFW / glad required.** Modern GL functions (VBOs, shaders,
  uniforms) are loaded with a tiny manual `glXGetProcAddressARB` loader bundled
  in `main15_gpu.cpp` (see *Code architecture* below). The build only needs
  `-lGL -lGLU -lglut -lm`.
- **No CUDA toolkit required.** The speedup comes from GPU-side OpenGL
  rasterization + shaders, not compute kernels. `nvcc` / CUDA are not used at
  compile time or runtime. (CUDA 13.0 happens to be installed here, but the GPU
  build does not depend on it.)
- **A display is required.** The app opens a real GLUT window. This server runs
  a real Xorg on `:0` (lightdm + sunshine), so `DISPLAY=:0.0` works. On a
  headless box you'd need `Xvfb`/`xvfb-run` or an EGL/Headless approach (out of
  scope for this build).

### Dependencies map

| Dependency | Provided by | Required? | Purpose |
|---|---|---|---|
| NVIDIA driver + `libGLX_nvidia.so` | `nvidia-driver-*` | **Yes** | Hardware GL 4.6 backend on the T4 |
| `libGL.so` (GLVND) | `libgl1` / `libglvnd0` | Yes | OpenGL 1.x symbols + `glXGetProcAddressARB` loader |
| `libGLU.so` | `libglu1-mesa-dev` | Yes | `gluPerspective` (fisheye projection) |
| FreeGLUT (`libglut.so`, `GL/glut.h`) | `freeglut3-dev` | Yes | Window, context, input callbacks |
| `libm` | glibc | Yes | `math.h` / `M_PI` |
| GLEW / GLFW / glad | — | **No** | Replaced by the manual loader |
| CUDA toolkit (`nvcc`, cudart) | `cuda-toolkit` | **No** | Not used by the GL render path |
| `mesa-utils` (`glxinfo`) | `mesa-utils` | Optional | Only to inspect the active renderer |
| C++17 compiler | `g++` (≥ 7) | Yes | Raw string literals, structured bindings |

### Code architecture

**Files added/changed for the GPU build:**

| File | Role |
|---|---|
| `main15_gpu.cpp` | New program. Copy of `main15.cpp` with the render path rewritten (see below). |
| `run_gpu.sh` | Launcher: sets `DISPLAY` + `__GLX_VENDOR_LIBRARY_NAME=nvidia`, execs `./main15_gpu`. |
| `Cube.hpp` | Two small **additive** `const` accessors added (no existing behavior changed): |

```cpp
void fillVertexLattice(std::vector<float>& positions) const;                    // n^3*8*3 floats
void fillCheckerboardIndices(int x,int y,int z, std::vector<unsigned int>& idx) const; // 12 tris/active cell
```

These iterate `subcells_[i][j][k]` and `cube_triangles_` exactly the way the
existing `getCheckerboardSubcells` / `getCheckerboardFacets` do (same `[0,n)`
indexing, same predicate, same triangulation table), so the GPU-rendered
geometry is identical to the CPU path — only the per-frame `Facet`
construction cost is removed. They use plain `float` / `unsigned int` so
`Cube.hpp` stays GL-agnostic.

**Render path (per frame) in `main15_gpu.cpp`:**

```
display()
 ├─ glClear
 ├─ set fixed-function MODELVIEW (glTranslate + glRotate)   ← same as main15
 ├─ drawAxes()                                              ← unchanged fixed-function
 ├─ drawGPU()
 │   ├─ read back PROJECTION & MODELVIEW with glGetFloatv
 │   ├─ uMVP = PROJECTION * MODELVIEW                        ← guarantees positional parity
 │   ├─ compute camera pos + light dir (inverse of view rotation)
 │   ├─ if selection dirty: rebuild IBO (only on i/I/j/J/k/K)
 │   ├─ glUseProgram; bind VBO/IBO; set uniforms
 │   ├─ fill pass:  glPolygonMode(FILL) + POLYGON_OFFSET_FILL; glDrawElements  (uMode=0)
 │   └─ outline pass: glPolygonMode(LINE); glDrawElements                    (uMode=1, black)
 ├─ drawHUD()                                               ← unchanged fixed-function
 └─ glutSwapBuffers
```

Key architectural decisions:

- **No external GL loader.** A ~25-line manual loader uses
  `glXGetProcAddressARB` (from `<GL/glx.h>`) to fetch the post-GL-1.1 function
  pointers (`glGenBuffers`, `glCreateShader`, `glUniformMatrix4fv`,
  `glVertexAttribPointer`, …). GL 1.x calls (`glDrawElements`, `glClear`,
  `glPolygonMode`, …) link directly from `-lGL`. Verified end-to-end on the T4.
- **Matrix parity trick.** The fixed-function matrix stack is kept exactly as
  in `main15` (`reshape` sets the projection with the fisheye distortion;
  `display` sets the modelview). The shader's `uMVP` is built by reading those
  matrices back with `glGetFloatv` and multiplying them. This guarantees the
  GPU-rendered triangles land at the same screen positions as the original
  fixed-function pipeline would have placed them, and it lets the axes/HUD
  keep using fixed-function (they're tiny and not worth porting).
- **Two-sided flat lighting in the fragment shader.** Normals are computed
  per-triangle from `normalize(cross(dFdx(vWorldPos), dFdy(vWorldPos)))` (zero
  CPU normal work), flipped to face the camera, and shaded with Blinn-Phong
  using the same ambient/diffuse/specular/shininess values as the original
  `initGL()` (`GL_LIGHT0` ambient 0.3, diffuse 0.7, specular 1.0, shininess 128,
  gray material `(200,200,200)`).
- **Index buffer rebuilt lazily.** `g_selectionDirty` is set by the `i/I/j/J/k/K`
  key handlers; the IBO is rebuilt once on the next frame (and on the first
  frame). Steady-state mouse rotate/pan/zoom never touches the IBO — it only
  re-uploads the `uMVP` uniform.

**Reused geometry classes (unchanged):** `Cube.hpp`, `Facet.hpp/cpp`,
`FacetBox.hpp`, `Vector3D.hpp/cpp`, `Quaternion.hpp/cpp`. The one-time `Setup()`
pipeline (`applySigmaTransformationToCube` → `refreshTriangulation` →
`writeSTL_s`) is identical to `main15`. `Dodecahedron`, `Torus`, `Vector4D`,
and `FacetBox.cpp` are linked in transitively but are **not on the per-frame
render path** (legacy).

### Performance numbers (this server)

| Metric | `main15` (CPU/llvmpipe) | `main15_gpu` (Tesla T4) |
|---|---|---|
| OpenGL renderer | Mesa `llvmpipe` (software) | `Tesla T4/PCIe/SSE2` (hardware) |
| Per-frame CPU facet construction | ~112,860 `Facet`s (cross+`sqrt` each) | **0** (normals in shader) |
| Per-frame draw calls | ~225,000 `glBegin/glEnd` batches | **2** `glDrawElements` |
| VBO upload | every frame (immediate mode) | **once** (405,224 verts, 4.64 MB) |
| IBO rebuild | n/a (rebuilt facets every frame) | **only on `ii/jj/kk` change** (112,860 tris) |
| Selection triangle count | ~112,860 | ~112,860 (identical geometry) |

### Troubleshooting

- **`RENDERER = llvmpipe`** — GPU not selected. Run via `./run_gpu.sh`, or
  `export __GLX_VENDOR_LIBRARY_NAME=nvidia`. The program also prints a warning
  in this case.
- **`freeglut: ERROR: Internal error <FBConfig with necessary capabilities not found>`**
  — the requested display mode couldn't be satisfied. Make sure you're using
  `GLUT_DEPTH` (not `GL_DEPTH`), that an X server is running on `$DISPLAY`, and
  that the NVIDIA driver is loaded (`nvidia-smi` works).
- **No window / `DISPLAY` not set** — on a headless box you need an X server
  (`Xvfb`/`xvfb-run`) or an EGL headless context. This build targets the
  on-server `:0` Xorg.
- **Link errors about `glGenBuffers` etc.** — those are loaded at runtime, not
  linked; if you see them, you've accidentally called a modern GL function
  without going through the loader. The build line is just
  `g++ -o main15_gpu main15_gpu.cpp -lGL -lGLU -lglut -lm`.
- **Shader compile error printed to stderr** — the program aborts if the
  `#version 460 core` program fails to compile/link. This means the active GL
  context is < 4.6 (e.g. an old software fallback); ensure the NVIDIA backend is
  selected.

---

## Cube Class Usage

The `Cube.hpp` header provides a complete cube implementation with triangulation and subdivision capabilities, featuring **full subdivision control** for real-time manipulation and selective rendering.

### Basic Usage

```cpp
#include "Cube.hpp"

// Create a basic cube
Cube basic_cube(1.0, Vector3D{0, 0, 0});  // radius=1.0, centered at origin

// Create a subdivided cube
Cube subdivided_cube(1.0, Vector3D{0, 0, 0}, 2);  // 2 subdivision levels
```

### Constructor Options

1. **Basic Cube**: `Cube(radius, center)`
   - Creates a cube with 12 triangular faces (2 triangles per face × 6 faces)
   - `radius` is the distance from center to face (half the side length)

2. **Subdivided Cube**: `Cube(radius, center, subdivisions)`
   - Creates a cube subdivided into `n³` smaller cubes (subcells)
   - Each subcell is triangulated using the regular 12 triangle mesh
   - Total triangles: `n³ × 12 triangles` (e.g., 2³ × 12 = 96 triangles)

---

## 🏗️ Advanced Subdivision Control

### Memory Layout and Storage

When a cube is subdivided with `n` levels, the subdivision creates a **3D grid** of subcells stored in memory as:

```cpp
std::vector<std::vector<std::vector<SubCell>>> subcells_[i][j][k]
```

**Storage Details:**
- **Grid Coordinates**: `(i,j,k)` where `i,j,k ∈ [0, n)`
- **Memory Layout**: Each `SubCell` contains:
  - 8 vertices (cube corners) as `std::array<Vector3D, 8>`
  - Center point, radius, grid coordinates, and active state
- **Vertex Order**: Standard cube vertex numbering (see diagram below)
- **Total Storage**: `n³` subcells, each with 8 vertices = `8n³` vertices total

### SubCell Structure

```cpp
struct SubCell {
    std::array<Vector3D, 8> vertices;  // 8 cube vertices in standard order
    Vector3D center;                   // Center point of subcell
    double radius;                     // Half side length
    int i, j, k;                      // Grid coordinates [0,n)
    bool active = true;               // Enable/disable for rendering
};
```

### Vertex Numbering Convention

Each subcell follows the standard cube vertex numbering:
```
       3----------2
      /|         /|
     / |        / |
    7----------6  |     Z
    |  |       |  |     |
    |  0-------|--1     |
    | /        | /      |
    |/         |/       +-----Y
    4----------5       /
                      X
```

### Core Access Methods

#### 1. Accessing Individual Subcells

```cpp
// Get read-only access to a specific subcell
const SubCell& cell = cube.getSubCell(i, j, k);
cout << "Center: " << cell.center << ", Active: " << cell.active << endl;

// Get mutable access for modifications
SubCell& cell = cube.getSubCellMutable(i, j, k);
cell.active = false;  // Disable this subcell
```

#### 2. Accessing Planes of Subcells

```cpp
// Get all subcells in a plane
// axis: 0=YZ plane, 1=XZ plane, 2=XY plane
// layer: which layer [0, n)

auto yz_plane = cube.getPlane(0, 1);  // YZ plane at X=1
auto xz_plane = cube.getPlane(1, 0);  // XZ plane at Y=0
auto xy_plane = cube.getPlane(2, 2);  // XY plane at Z=2

// Iterate through plane subcells
for (const SubCell& cell : xy_plane) {
    cout << "Cell (" << cell.i << "," << cell.j << "," << cell.k << ")" << endl;
}
```

#### 3. Real-Time Vertex Manipulation

```cpp
// Update a specific vertex of a subcell
cube.updateSubCellVertex(i, j, k, vertex_idx, new_position);

// Example: Move vertex 0 of subcell (1,1,1) upward
Vector3D new_pos = cube.getSubCell(1, 1, 1).vertices[0] + Vector3D(0, 0, 0.1);
cube.updateSubCellVertex(1, 1, 1, 0, new_pos);

// IMPORTANT: Refresh triangulation after vertex updates
cube.refreshTriangulation();
```

### Selective Rendering Methods

#### 1. Render Individual Subcells

```cpp
// Get triangles for a specific subcell
FacetBox subcell_faces = cube.getSubCellFacets(1, 1, 1);

// Render with specific color
for (size_t i = 0; i < subcell_faces.size(); ++i) {
    drawFacet(subcell_faces[i], 255, 0, 0, 0.8f);  // Red color
}
```

#### 2. Render Planes of Subcells

```cpp
// Method 1: Traditional numeric coordinates
// Get all triangles in a plane using axis/layer system
FacetBox plane_faces = cube.getPlaneFacets(2, 0);  // XY plane at Z=0

// Method 2: NEW - String-based coordinates (more intuitive)
// Get planes using string coordinate system
FacetBox xy_center = cube.getPlaneFacets("x", "y", "0");    // Center XY plane (Z=0)
FacetBox xz_center = cube.getPlaneFacets("x", "0", "z");    // Center XZ plane (Y=0)
FacetBox yz_center = cube.getPlaneFacets("0", "y", "z");    // Center YZ plane (X=0)

// Advanced usage - specific layers
FacetBox xy_top    = cube.getPlaneFacets("x", "y", "2");    // Upper XY plane (Z=+2)
FacetBox xy_bottom = cube.getPlaneFacets("x", "y", "-2");   // Lower XY plane (Z=-2)

// Render plane with gradient colors
for (size_t i = 0; i < plane_faces.size(); ++i) {
    float hue = float(i) / float(plane_faces.size()) * 360.0f;
    Color c = hsv2rgb(hue, 0.8f, 0.9f);
    int R = int(c.r * 255), G = int(c.g * 255), B = int(c.b * 255);
    drawFacet(plane_faces[i], R, G, B, 0.7f);
}
```

#### 3. Enable/Disable Subcells

```cpp
// Disable specific subcells for "hollow" effects
cube.setSubCellActive(1, 1, 1, false);  // Hide center subcell

// Create a "frame" by disabling interior cells
for (int i = 1; i < n-1; ++i) {
    for (int j = 1; j < n-1; ++j) {
        for (int k = 1; k < n-1; ++k) {
            cube.setSubCellActive(i, j, k, false);
        }
    }
}

// Refresh rendering after changes
cube.refreshTriangulation();
```

### Real-Time Animation Examples

#### 1. Wave Animation

```cpp
void animateWave(Cube& cube, double time) {
    int n = cube.getSubdivisionLevels();

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                SubCell& cell = cube.getSubCellMutable(i, j, k);

                // Apply wave displacement to all vertices
                double wave = 0.1 * sin(time + i * 0.5 + j * 0.3);
                for (int v = 0; v < 8; ++v) {
                    cell.vertices[v] += Vector3D(0, 0, wave);
                }
            }
        }
    }
    cube.refreshTriangulation();
}
```

#### 2. Selective Layer Display

```cpp
void showLayer(Cube& cube, int current_layer) {
    int n = cube.getSubdivisionLevels();

    // Hide all subcells
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                cube.setSubCellActive(i, j, k, false);
            }
        }
    }

    // Show only current layer
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cube.setSubCellActive(i, j, current_layer, true);
        }
    }
    cube.refreshTriangulation();
}
```

### Performance Notes

- **Memory Usage**: `n³` subcells with 8 vertices each
- **Triangulation**: Call `refreshTriangulation()` after vertex modifications
- **Selective Rendering**: Use `getSubCellFacets()` and `getPlaneFacets()` for efficiency
- **Real-Time Updates**: Vertex modifications are O(1), retriangulation is O(n³)

### Integration with Rendering Pipeline

The cube integrates seamlessly with the existing rendering system:
- Uses same `Vector3D`, `Quaternion`, and `FacetBox` classes
- Compatible with `drawFacet()` function for OpenGL rendering
- Supports same transformation and manipulation methods as `Dodecahedron`

---

## 🎯 NEW: String-Based Coordinate System

### Enhanced `getPlaneFacets()` Method

**Update**: The `getPlaneFacets` method now supports an intuitive string-based coordinate system alongside the traditional numeric approach.

#### Method Signatures

```cpp
// Traditional method (still supported)
FacetBox getPlaneFacets(int axis, int layer);

// NEW: String-based method
FacetBox getPlaneFacets(const std::string& coord1,
                       const std::string& coord2,
                       const std::string& coord3);
```

#### How It Works

The new method accepts three string parameters representing coordinates in **X, Y, Z order**:
- **Two parameters** must be axis names: `"x"`, `"y"`, or `"z"`
- **One parameter** must be a numeric string: `"0"`, `"-1"`, `"2"`, etc.

The numeric coordinate determines:
1. **Which axis is fixed** (perpendicular to the plane)
2. **Which layer** to extract from that axis

#### Usage Examples

```cpp
// Center planes (most common usage)
FacetBox xy_plane = cube.getPlaneFacets("x", "y", "0");    // XY plane through center (Z=0)
FacetBox xz_plane = cube.getPlaneFacets("x", "0", "z");    // XZ plane through center (Y=0)
FacetBox yz_plane = cube.getPlaneFacets("0", "y", "z");    // YZ plane through center (X=0)

// Offset planes - positive values
FacetBox xy_top   = cube.getPlaneFacets("x", "y", "2");    // XY plane above center (Z=+2)
FacetBox xz_right = cube.getPlaneFacets("x", "2", "z");    // XZ plane right of center (Y=+2)
FacetBox yz_front = cube.getPlaneFacets("2", "y", "z");    // YZ plane in front of center (X=+2)

// Offset planes - negative values
FacetBox xy_bottom = cube.getPlaneFacets("x", "y", "-2");   // XY plane below center (Z=-2)
FacetBox xz_left   = cube.getPlaneFacets("x", "-2", "z");   // XZ plane left of center (Y=-2)
FacetBox yz_back   = cube.getPlaneFacets("-2", "y", "z");   // YZ plane behind center (X=-2)
```

#### Coordinate Range

For an n×n×n subdivision, numeric coordinates must be in range **[-n/2, n/2]**:

| Subdivision | Valid Range | Center | Examples |
|-------------|-------------|---------|----------|
| 4×4×4       | [-2, 1]     | 0       | "-2", "-1", "0", "1" |
| 8×8×8       | [-4, 3]     | 0       | "-4", "-3", ..., "2", "3" |
| 16×16×16    | [-8, 7]     | 0       | "-8", "-7", ..., "6", "7" |

#### Error Handling

The method provides comprehensive validation with descriptive error messages:

```cpp
try {
    // ❌ Invalid: two numeric coordinates
    cube.getPlaneFacets("0", "1", "z");
} catch (const std::invalid_argument& e) {
    // Error: "Exactly one coordinate must be numeric, got 2"
}

try {
    // ❌ Invalid: coordinate out of range
    cube.getPlaneFacets("x", "y", "10");  // Assuming 8×8×8 subdivision
} catch (const std::out_of_range& e) {
    // Error: "Coordinate 10 out of range [-4, 4]"
}

try {
    // ❌ Invalid: wrong axis names
    cube.getPlaneFacets("x", "0", "y");  // Should be "z" for Z coordinate
} catch (const std::invalid_argument& e) {
    // Error: "For fixed Y coordinate, other coordinates must be 'x' and 'z'"
}
```

#### Advantages of String-Based System

1. **More Intuitive**: `"x", "y", "0"` clearly means "XY plane at Z=0"
2. **Self-Documenting**: Code is more readable than `getPlaneFacets(2, 4)`
3. **Flexible Ordering**: Coordinates can be specified in logical order
4. **Better Error Messages**: Validation provides specific guidance
5. **Logical Center**: `"0"` always means center, regardless of subdivision level

#### Migration Guide

```cpp
// OLD: Traditional method
FacetBox plane1 = cube.getPlaneFacets(0, 4);  // YZ plane, layer 4
FacetBox plane2 = cube.getPlaneFacets(1, 4);  // XZ plane, layer 4
FacetBox plane3 = cube.getPlaneFacets(2, 4);  // XY plane, layer 4

// NEW: String-based method (assuming 8×8×8 subdivision where layer 4 = center)
FacetBox plane1 = cube.getPlaneFacets("0", "y", "z");    // YZ plane, center
FacetBox plane2 = cube.getPlaneFacets("x", "0", "z");    // XZ plane, center
FacetBox plane3 = cube.getPlaneFacets("x", "y", "0");    // XY plane, center
```

## ⚠️ Critical Coordinate System Guidelines

### Understanding Coordinate Systems

The Cube class uses **three different coordinate systems** that must not be confused:

#### 1. **Physical Coordinates** (Internal Storage)
- **Range**: `[0, n)` where n is the subdivision level
- **Usage**: Internal grid storage `subcells_[i][j][k]`
- **Direct Access**: Use `getSubCellPhysical(i, j, k)` (internal method)

#### 2. **Logical Coordinates** (Public API)
- **Range**: `[-center, center)` where `center = n/2`
- **Center**: `(0, 0, 0)` represents the center of the cube
- **Usage**: Most public methods expect logical coordinates
- **Example**: For n=8, logical coordinates range from -4 to +3

#### 3. **String Coordinates** (NEW - Intuitive API)
- **Range**: Numeric strings in `[-n/2, n/2]` plus axis names `"x"`, `"y"`, `"z"`
- **Center**: `"0"` always represents center, regardless of subdivision level
- **Usage**: New `getPlaneFacets(string, string, string)` method
- **Example**: `"x", "y", "0"` means XY plane through center

### Coordinate Conversion

The class automatically converts between systems:

```cpp
// Logical to Physical conversion (internal)
int center = subdivision_levels_ / 2;
int physical_i = center + logical_x;  // e.g., center=4, logical_x=1 → physical_i=5

// Physical to Logical conversion (internal)
int logical_x = physical_i - center;  // e.g., physical_i=5, center=4 → logical_x=1
```

### ❌ Common Out-of-Range Errors

#### **Error 1: Using Physical Coordinates in Public Methods**

```cpp
// ❌ WRONG: Using physical coordinates (0 to n-1)
int n = cube.getSubdivisionLevels();  // e.g., n = 8
for (int i = 0; i < n; i++) {         // i goes 0,1,2,3,4,5,6,7
    for (int j = 0; j < n; j++) {
        cube.getSubCell(i, j, k);     // ERROR: expects logical coordinates!
    }
}
```

**Error Message**: `std::out_of_range: Cube::getSubCell: logical coordinates out of range`

#### **Error 2: Wrong Center Calculation**

```cpp
// ❌ WRONG: Assuming (1,1,1) is always the center
cube.getSubCell(1, 1, 1);  // Only works for n=2 or n=3!
```

### ✅ Correct Usage Patterns

#### **Pattern 1: Proper Loop Iteration**

```cpp
// ✅ CORRECT: Use logical coordinate ranges
int n = cube.getSubdivisionLevels();
int center = n / 2;  // Calculate center offset

for (int i = -center; i < center; i++) {      // Logical coordinates
    for (int j = -center; j < center; j++) {
        for (int k = -center; k < center; k++) {
            cube.getSubCell(i, j, k);         // Now within valid range
            cube.setSubCellActive(i, j, k, true);
        }
    }
}
```

#### **Pattern 2: Center Access**

```cpp
// ✅ CORRECT: Always use (0,0,0) for center in logical coordinates
const SubCell& center_cell = cube.getSubCell(0, 0, 0);  // Always valid
Vector3D center_pos = cube.getSubcellCenter(0, 0, 0);   // Always valid
```

#### **Pattern 3: Coordinate Range Validation**

```cpp
bool isValidLogicalCoordinate(const Cube& cube, int x, int y, int z) {
    int n = cube.getSubdivisionLevels();
    int center = n / 2;
    return (x >= -center && x < center &&
            y >= -center && y < center &&
            z >= -center && z < center);
}

// Use before accessing
if (isValidLogicalCoordinate(cube, x, y, z)) {
    cube.getSubCell(x, y, z);
}
```

### 📊 Coordinate Reference Tables

#### For n=8 Subdivision (8³=512 subcells):
| Coordinate Type | X Range | Y Range | Z Range | Center Location |
|----------------|---------|---------|---------|-----------------|
| Physical       | [0, 7]  | [0, 7]  | [0, 7]  | (4, 4, 4)       |
| Logical        | [-4, 3] | [-4, 3] | [-4, 3] | (0, 0, 0)       |

#### For n=4 Subdivision (4³=64 subcells):
| Coordinate Type | X Range | Y Range | Z Range | Center Location |
|----------------|---------|---------|---------|-----------------|
| Physical       | [0, 3]  | [0, 3]  | [0, 3]  | (2, 2, 2)       |
| Logical        | [-2, 1] | [-2, 1] | [-2, 1] | (0, 0, 0)       |

### 🛠️ Debugging Tips

1. **Enable Debug Output**: The `logicalToPhysical()` method prints conversion details:
```cpp
// Debug output shows coordinate conversion
cube.getSubCell(x, y, z);  // Will print: i = 5, j = 5, k = 5
```

2. **Check Subdivision Level**:
```cpp
if (!cube.hasSubdivision()) {
    cerr << "Error: No subdivision data available!" << endl;
    return;
}
int n = cube.getSubdivisionLevels();
cout << "Subdivision level: " << n << "³ subcells" << endl;
```

3. **Validate Before Access**:
```cpp
try {
    cube.getSubCell(x, y, z);
} catch (const std::out_of_range& e) {
    cerr << "Coordinate error: " << e.what() << endl;
    cerr << "Attempted coordinates: (" << x << "," << y << "," << z << ")" << endl;
}
```

### 🚀 Compilation Instructions

```bash
# Compile with OpenGL libraries (essential for linking)
g++ -o main11 main11.cpp -lGL -lGLU -lglut

# Run the program
./main11
```

**Note**: The `-lGL -lGLU -lglut` flags are **required** for all OpenGL programs to avoid linking errors.

