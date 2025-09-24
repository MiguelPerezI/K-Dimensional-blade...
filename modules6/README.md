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

