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
// Get all triangles in a plane
FacetBox plane_faces = cube.getPlaneFacets(2, 0);  // XY plane at Z=0

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

