/*
┌─────────────────────────────────────────────────────────────────────────────────┐
│                              CUBE OBJECT LIFECYCLE                              │
└─────────────────────────────────────────────────────────────────────────────────┘

CREATION PHASE:
┌─────────────────────┐    ┌─────────────────────┐    ┌─────────────────────┐
│   Cube() Default    │    │  Cube(r, center)    │    │ Cube(r,center,n)    │
│                     │    │                     │    │                     │
│ • Empty cube        │    │ • Basic 8 vertices  │    │ • Basic 8 vertices  │
│ • No vertices       │    │ • 12 triangles      │    │ • Then subdivide    │
│ • subdivision=0     │    │ • subdivision=0     │    │ • n³ subcells       │
│ • built=false       │    │ • built=false       │    │ • n³×12 triangles   │
└─────────────────────┘    └─────────────────────┘    └─────────────────────┘
         │                           │                           │
         │                           │                           │
         ▼                           ▼                           ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│                              CUBE OBJECT STATE                                  │
│                                                                                 │
│  Basic Cube (n=0):          Subdivided Cube (n>0):                              │
│  ┌─────────────────┐        ┌─────────────────────────────────────────────┐     │
│  │ verts_[8]       │        │ verts_[8]           subcells_[n][n][n]      │     │
│  │┌─┬─┬─┬─┬─┬─┬─┬─┐│        │ ┌─┬─┬─┬─┬─┬─┬─┬─┐    ┌─────┬─────┬─────┐    │     │
│  ││0│1│2│3│4│5│6│7││        │ │0│1│2│3│4│5│6│7│    │[0]  │[1]  │...  │    │     │
│  │└─┴─┴─┴─┴─┴─┴─┴─┘│        │ └─┴─┴─┴─┴─┴─┴─┴─┘    │     │     │     │    │     │
│  │                 │        │                      │ ┌─┬─┬─┐ ┌─┬─┬─┐ │    │     │
│  │ facets_[12]     │        │ facets_[n³×12]       │ │ │ │ │ │ │ │ │ │    │     │
│  │ 12 triangles    │        │ Many triangles       │ └─┴─┴─┘ └─┴─┴─┘ │    │     │
│  └─────────────────┘        └─────────────────────────────────────────────┘     │
└─────────────────────────────────────────────────────────────────────────────────┘

AVAILABLE METHODS AFTER CREATION:

┌─────────────────────────────────────────────────────────────────────────────────┐
│                              METHOD CATEGORIES                                  │
└─────────────────────────────────────────────────────────────────────────────────┘

1. BASIC CUBE OPERATIONS (Always Available):
   ┌─────────────────────────────────────────────────────────────────────────┐
   │ ACCESS & QUERY:                    GEOMETRIC OPERATIONS:                │
   │ • faceCount()       → size_t       • center()          → Vector3D       │
   │ • operator[](k)     → Facet&       • translate(offset) → void           │
   │ • getFacets()       → FacetBox&    • scale(s, pivot)   → void           │
   │                                                                         │
   │ EXAMPLE USAGE:                                                          │
   │ cube.faceCount()           // Returns 12 for basic, n³×12 for subdivided│
   │ cube[5]                    // Get 6th triangle                          │
   │ cube.translate(Vector3D(1,0,0))  // Move cube 1 unit in X               │
   └─────────────────────────────────────────────────────────────────────────┘

2. SUBDIVISION OPERATIONS:
   ┌─────────────────────────────────────────────────────────────────────────┐
   │ CREATE SUBDIVISION:                CHECK SUBDIVISION:                   │
   │ • subdivide(n)      → void         • hasSubdivision()    → bool         │
   │                                    • getSubdivisionLevels() → int       │
   │                                                                         │
   │ FLOW DIAGRAM:                                                           │
   │ Basic Cube ──subdivide(3)──→ 3×3×3 = 27 subcells                        │
   │                               Each subcell: 8 vertices + 12 triangles   │
   │                               Total: 27×12 = 324 triangles              │
   └─────────────────────────────────────────────────────────────────────────┘

3. SUBDIVISION-ONLY METHODS (Require hasSubdivision() == true):
   ┌─────────────────────────────────────────────────────────────────────────┐
   │ COORDINATE CONVERSION:                                                  │
   │ • logicalToPhysical(x,y,z)  → tuple<int,int,int>                        │
   │ • physicalToLogical(i,j,k)  → tuple<int,int,int>                        │
   │                                                                         │
   │ LOGICAL COORDS (User-friendly):    PHYSICAL COORDS (Internal):          │
   │ ┌─────────────────────────────┐    ┌─────────────────────────────┐      │
   │ │(-1,1,0)  │(0,1,0) │(1,1,0)  │    │(0,2,1)  │(1,2,1) │(2,2,1)   │      │
   │ │─────────┼────────┼───────── │    │─────────┼────────┼───────── │      │
   │ │(-1,0,0) │(0,0,0) │(1,0,0)   │ ←→ │(0,1,1)  │(1,1,1) │(2,1,1)   │      │
   │ │─────────┼────────┼───────── │    │─────────┼────────┼───────── │      │
   │ │(-1,-1,0)│(0,-1,0)│(1,-1,0)  │    │(0,0,1)  │(1,0,1) │(2,0,1)   │      │
   │ └─────────────────────────────┘    └─────────────────────────────┘      │
   │              n=5, Z=0 plane                    Array indices            │
   └─────────────────────────────────────────────────────────────────────────┘

4. SUBCELL ACCESS (Subdivision Required):
   ┌─────────────────────────────────────────────────────────────────────────┐
   │ READ-ONLY ACCESS:                  MODIFIABLE ACCESS:                   │
   │ • getSubCell(x,y,z)         → const SubCell&                           │
   │ • getSubCellPhysical(i,j,k) → const SubCell&                           │
   │ • getSubCellMutable(x,y,z)  → SubCell&                                 │
   │ • getSubCellMutablePhysical(i,j,k) → SubCell&                          │
   │                                                                         │
   │ SUBCELL PROPERTIES ACCESS:                                              │
   │ • getSubcellCenter(x,y,z)   → Vector3D                                 │
   │ • getSubcellRadius(x,y,z)   → double                                   │
   │ • setSubCellActive(x,y,z,active) → void                                │
   └─────────────────────────────────────────────────────────────────────────┘

5. VERTEX MANIPULATION (Subdivision Required):
   ┌─────────────────────────────────────────────────────────────────────────┐
   │ MODIFY VERTICES:                   REBUILD TRIANGULATION:               │
   │ • updateSubCellVertex(x,y,z,       • refreshTriangulation() → void     │
   │   vertex_idx, new_pos) → void                                           │
   │                                                                         │
   │ WORKFLOW:                                                               │
   │ 1. updateSubCellVertex(0,0,0, 3, newPos)  // Modify vertex 3 of center │
   │ 2. updateSubCellVertex(1,0,0, 1, newPos2) // Modify another vertex     │
   │ 3. refreshTriangulation()                  // Rebuild all triangles    │
   └─────────────────────────────────────────────────────────────────────────┘

6. PATTERN & PLANE EXTRACTION (Subdivision Required):
   ┌─────────────────────────────────────────────────────────────────────────┐
   │ PLANE EXTRACTION:                                                       │
   │ • getPlane(axis, layer)     → vector<SubCell&>                         │
   │ • getPlaneFacets(axis, layer) → FacetBox                               │
   │ • getPlaneFacets("x","y","0") → FacetBox  // String interface          │
   │                                                                         │
   │ PATTERN SELECTION:                                                      │
   │ • getCheckerboardSubcells() → vector<SubCell&>                         │
   │ • getCheckerboardFacets()   → FacetBox                                 │
   │                                                                         │
   │ INDIVIDUAL SUBCELL FACETS:                                              │
   │ • getSubCellFacets(x,y,z)   → FacetBox                                 │
   └─────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────────┐
│                            TYPICAL USAGE PATTERNS                              │
└─────────────────────────────────────────────────────────────────────────────────┘

PATTERN 1: Basic Cube Operations
┌───────────────────────────────────────────┐
│ Cube cube(1.0, Vector3D(0,0,0));         │ // Create basic cube
│ cout << cube.faceCount();  // 12          │ // Check face count
│ cube.translate(Vector3D(5,0,0));          │ // Move cube
│ auto center = cube.center();              │ // Get new center
│ cube.scale(2.0, center);                  │ // Double size
└───────────────────────────────────────────┘

PATTERN 2: Subdivision and Access
┌───────────────────────────────────────────┐
│ Cube cube(1.0, Vector3D(0,0,0), 5);      │ // Create 5×5×5 subdivision
│ if (cube.hasSubdivision()) {              │ // Always check first
│   auto center_cell = cube.getSubCell(0,0,0); // Get center subcell
│   auto corner_cell = cube.getSubCell(2,2,2); // Get corner subcell
│   cube.setSubCellActive(1,1,1, false);    │ // Disable one subcell
│ }                                         │
└───────────────────────────────────────────┘

PATTERN 3: Vertex Modification
┌───────────────────────────────────────────┐
│ cube.updateSubCellVertex(0,0,0, 0,        │ // Modify center subcell
│                         Vector3D(1,1,1)); │ // vertex 0 position
│ cube.updateSubCellVertex(0,0,0, 7,        │ // Modify center subcell  
│                         Vector3D(2,2,2)); │ // vertex 7 position
│ cube.refreshTriangulation();              │ // REQUIRED after changes
└───────────────────────────────────────────┘

PATTERN 4: Plane Extraction
┌───────────────────────────────────────────┐
│ // Get center XY plane (Z=0)             │
│ auto xy_plane = cube.getPlaneFacets("x","y","0");                    │
│                                           │
│ // Get all subcells in YZ plane (X=1)    │
│ auto yz_subcells = cube.getPlane(0, 3);  │
│                                           │
│ // Get checkerboard pattern              │
│ auto checker_facets = cube.getCheckerboardFacets();                 │
└───────────────────────────────────────────┘

ERROR CONDITIONS TO HANDLE:
┌─────────────────────────────────────────────────────────────────────────────────┐
│ • std::out_of_range: Invalid coordinates or indices                            │
│ • std::runtime_error: Calling subdivision methods on non-subdivided cube       │
│ • std::invalid_argument: Invalid string coordinates in getPlaneFacets()        │
│                                                                                 │
│ ALWAYS CHECK: hasSubdivision() before using subdivision-specific methods       │
└─────────────────────────────────────────────────────────────────────────────────┘
*/

//This diagram shows:
//
//1. **Creation Phase**: Three different constructors and their resulting states
//2. **Object Structure**: Internal data layout for basic vs subdivided cubes
//3. **Method Categories**: Organized by functionality and requirements
//4. **Coordinate Systems**: Visual comparison of logical vs physical coordinates
//5. **Usage Patterns**: Common workflows with example code
//6. **Error Handling**: What to watch out for


#ifndef CUBE_H  // Header guard to prevent multiple inclusions
#define CUBE_H

using namespace std;  // Use standard namespace for convenience

// Include standard library headers for various functionality
#include <stdexcept>   // For exception handling (std::out_of_range, std::runtime_error)
#include <vector>      // For dynamic arrays (std::vector)
#include <array>       // For fixed-size arrays (std::array)
#include <functional>  // For reference_wrapper in plane operations
#include <tuple>       // For coordinate transformation return values
#include <algorithm>   // For std::all_of in string validation
#include <string>      // For string coordinate parsing
#include <cctype>      // For character type checking (isdigit)
#include <math.h>      // For mathematical functions (sqrt)
#include "Vector3D.hpp"    // Custom 3D vector class for spatial coordinates
#include "Quaternion.hpp"  // Custom quaternion class for vertex storage
#include "FacetBox.hpp"    // Custom container for triangular faces

/**
 * @brief Represents a cube in 3D space with advanced subdivision and triangulation capabilities.
 *
 * This class stores a cube as 8 vertices (using Quaternion points) and a collection of 
 * triangular faces (FacetBox). The cube can be subdivided into n×n×n smaller subcubes,
 * each generating 12 triangular faces using standard cube triangulation. Supports both
 * basic cube operations (translate, scale) and advanced subdivision features including
 * coordinate transformation, plane extraction, and selective subcell manipulation.
 * 
 * Features:
 * - Basic cube geometry with 8 vertices and 12 triangular faces
 * - Subdivision into n×n×n smaller cubes for detailed modeling
 * - Coordinate system conversion (logical ↔ physical)
 * - Plane extraction and pattern-based subcell selection
 * - Individual vertex manipulation with retriangulation
 * - Selective subcell activation/deactivation
 * 
 * Uses Rule-of-Zero for memory management with std::vector internally.
 */
class Cube {
public:
    // Friend function for stream output - allows direct access to private members
    friend std::ostream& operator<<(std::ostream& os, const Cube& cube);

    // === SUBDIVISION DATA STRUCTURES ===
    /**
     * @brief Represents a single subcell within a subdivided cube.
     * 
     * Each subcell is essentially a smaller cube with its own complete geometry,
     * spatial properties, and state. Subcells can be individually manipulated,
     * activated/deactivated, and have their vertices modified independently.
     * The subcell knows its position in both 3D space and within the subdivision grid.
     */
    struct SubCell {
        std::array<Vector3D, 8> vertices;  // 8 cube vertices in standard order (0-7) for consistent triangulation
        Vector3D center;                   // Geometric center point of subcell in 3D world space
        double radius;                     // Half side length (distance from center to face, not corner)
        int i, j, k;                      // Physical grid coordinates [0, n) within the subdivision matrix
        bool active = true;               // Visibility/processing flag - can be disabled for selective rendering

        // Default constructor - creates uninitialized subcell (used for vector resizing)
        SubCell() = default;

        // Parameterized constructor - initializes all subcell properties from provided data
        SubCell(const std::array<Vector3D, 8>& verts, const Vector3D& c, double r, int x, int y, int z)
            : vertices(verts), center(c), radius(r), i(x), j(y), k(z) {}  // Member initializer list for efficiency
    };

public:
    /* === CONSTRUCTORS AND RULE OF ZERO === */
    
    // Default constructor - creates empty cube with no geometry or subdivision data
    Cube() : subdivision_levels_(0), subdivision_built_(false) {}

    // Destructor - compiler-generated cleanup handles all std::vector cleanup automatically
    ~Cube() = default;

    // Copy constructor - performs deep copy of all cube data including subdivision grid
    Cube(const Cube&) = default;

    // Move constructor - efficient transfer of resources without copying (C++11 optimization)
    Cube(Cube&&) noexcept = default;

    // Copy assignment operator - assigns complete cube state from another cube instance
    Cube& operator=(const Cube&) = default;

    // Move assignment operator - efficient assignment via resource transfer (C++11 optimization)
    Cube& operator=(Cube&&) noexcept = default;

    /**
     * @brief Construct a basic cube with specified radius and center position.
     * 
     * Creates a simple cube with 8 vertices and 12 triangular faces. The radius
     * parameter represents half the side length (distance from center to face).
     * This constructor generates the standard cube triangulation immediately.
     * 
     * @param radius Distance from center to face (half of side length).
     * @param center Geometric center position of the cube in 3D space.
     */
    Cube(double radius, const Vector3D& center) : subdivision_levels_(0), subdivision_built_(false) {
        initVertices(radius, center);  // Generate the 8 corner vertices using standard cube layout
        buildFacets();                 // Create 12 triangular faces from vertices using lookup table
    }

    /**
     * @brief Construct a subdivided cube with specified subdivision levels.
     * 
     * Creates a cube and immediately subdivides it into n×n×n smaller subcubes.
     * Each subcube maintains its own 8 vertices and generates 12 triangular faces,
     * resulting in n³×12 total triangular faces for the complete structure.
     * 
     * @param radius Distance from center to face of the main cube.
     * @param center Geometric center of the main cube.
     * @param subdivisions Number of subdivision levels (0 = basic cube, n = n×n×n subcubes).
     */
    Cube(double radius, const Vector3D& center, int subdivisions) : subdivision_levels_(0), subdivision_built_(false) {
        initVertices(radius, center);  // Generate the 8 corner vertices of the main cube boundary
        buildFacets();                 // Create basic 12 triangular faces for the main cube
        if (subdivisions > 0) {        // Only perform subdivision for positive subdivision counts
            subdivide(subdivisions);   // Replace basic triangulation with n×n×n subcube triangulation
        }
    }

    /* === FACE COUNT AND ACCESS === */
    
    /**
     * @brief Get the total number of triangular faces in the cube.
     * 
     * For basic cubes: returns 12 (standard cube triangulation).
     * For subdivided cubes: returns n³×12 where n is subdivision levels.
     * 
     * @return Total count of triangular faces available for rendering.
     */
    size_t faceCount() const noexcept { return facets_.size(); }

    /**
     * @brief Access a specific triangular face by index.
     * 
     * Provides array-style access to individual triangular faces with automatic
     * bounds checking. Each face is a Facet object containing 3 vertices.
     * 
     * @param k Face index in range [0, faceCount()).
     * @return Const reference to the k-th triangular face.
     * @throws std::out_of_range if k is invalid.
     */
    const Facet& operator[](size_t k) const {
        return facets_[k];  // Delegate to FacetBox's bounds-checked access
    }

    /* === GEOMETRIC OPERATIONS === */
    
    /**
     * @brief Calculate the geometric center of all vertices.
     * 
     * Computes the centroid by averaging the positions of all cube vertices.
     * For basic cubes, this should match the original center parameter.
     * For subdivided cubes, this represents the center of the overall structure.
     * 
     * @return Vector3D representing the average position of all vertices.
     * @throws std::runtime_error if no vertices have been initialized.
     */
    Vector3D center() const {
        if (verts_.empty()) {  // Validate that cube has been properly initialized
            throw std::runtime_error("Cube::center: no vertices initialized");
        }
        Vector3D sum{0,0,0};  // Initialize accumulator vector with zero components
        for (auto const& q : verts_) {  // Iterate through all stored vertices
            sum += q.V();  // Extract Vector3D from Quaternion and add to running sum
        }
        return sum / static_cast<double>(verts_.size());  // Return arithmetic mean of all positions
    }

    /**
     * @brief Translate the entire cube by a specified offset vector.
     * 
     * Moves all vertices by adding the offset to their current positions.
     * This affects both the main cube vertices and any subdivision data.
     * Automatically rebuilds triangulation to reflect new vertex positions.
     * 
     * @param offset Vector3D displacement to apply to all vertices.
     */
    void translate(const Vector3D& offset) {
        for (auto& q : verts_) {  // Process each vertex stored as Quaternion
            Vector3D p = q.V() + offset;  // Calculate new position by adding offset
            q = Quaternion(0.0, p);       // Update quaternion with new position (scalar=0)
        }
        buildFacets();  // Regenerate triangular faces to reflect new vertex positions
    }

    /**
     * @brief Scale all vertices relative to a specified pivot point.
     * 
     * Applies uniform scaling by multiplying distances from the pivot point
     * by the scaling factor. Vertices closer to pivot move less, creating
     * natural scaling behavior. Automatically rebuilds triangulation.
     * 
     * @param s Scaling factor (1.0=no change, >1.0=larger, <1.0=smaller).
     * @param pivot The fixed point to scale around.
     */
    void scale(double s, const Vector3D& pivot) {
        for (auto& q : verts_) {  // Process each vertex individually
            Vector3D p = pivot + s * (q.V() - pivot);  // Scale distance from pivot by factor s
            q = Quaternion(0.0, p);  // Update quaternion with scaled position
        }
        buildFacets();  // Regenerate triangular faces with new scaled geometry
    }

    /**
     * @brief Get read-only access to the underlying triangular face collection.
     * 
     * Provides direct access to the FacetBox containing all triangular faces.
     * Useful for rendering systems that need to iterate through all faces
     * or perform batch operations on the complete face set.
     * 
     * @return Const reference to the internal FacetBox container.
     */
    const FacetBox& getFacets() const noexcept {
        return facets_;  // Direct access to internal triangular face storage
    }

    /* === SUBDIVISION OPERATIONS === */
    
    /**
     * @brief Subdivide the cube into n×n×n smaller subcubes with full triangulation.
     * 
     * Replaces the basic cube triangulation with a complex subdivision where each
     * subcube maintains its own 8 vertices and generates 12 triangular faces.
     * The subdivision creates a 3D grid of smaller cubes that can be individually
     * manipulated, activated/deactivated, or modified. This enables detailed
     * geometric modeling and selective rendering operations.
     * 
     * @param n Number of subdivision levels (n×n×n total subcubes created).
     */
    void subdivide(int n) {
        if (n <= 0) return;  // Early exit for invalid subdivision counts

        // Store current cube properties before rebuilding
        Vector3D cubeCenter = center();    // Get current geometric center
        double cubeRadius = getRadius();   // Get current cube radius (half side length)

        // Update subdivision state tracking
        subdivision_levels_ = n;

        // Replace simple triangulation with complex subdivided structure
        facets_.clear();  // Remove all existing triangular faces
        buildSubdividedCube(cubeCenter, cubeRadius, n);  // Generate n³ subcubes with triangulation
        subdivision_built_ = true;  // Mark subdivision data as available and valid
    }

    /* === SUBDIVISION ACCESS AND QUERIES === */
    
    /**
     * @brief Check if subdivision data is available and valid.
     * 
     * Verifies both that subdivision has been performed and that the subdivision
     * data structures are properly populated. Required before accessing subcells.
     * 
     * @return True if subdivision data exists and is accessible.
     */
    bool hasSubdivision() const noexcept {
        return subdivision_built_ && !subcells_.empty();  // Both flag set and data structures populated
    }

    /**
     * @brief Get the current subdivision level.
     * 
     * Returns the subdivision parameter used to create the current subcube grid.
     * For n subdivision levels, there are n×n×n total subcubes in the structure.
     * 
     * @return Number of subdivision levels (n means n×n×n subcubes total).
     */
    int getSubdivisionLevels() const noexcept {
        return subdivision_levels_;  // Direct access to stored subdivision count
    }

    /* === COORDINATE TRANSFORMATION === */
    
    /**
     * @brief Convert user-friendly logical coordinates to internal grid indices.
     * 
     * Transforms logical coordinates (centered at origin) to physical array indices
     * for internal subcell access. Logical coordinates allow users to think of
     * the cube center as (0,0,0) with positive/negative indices extending outward.
     * 
     * Examples for n=5: logical (0,0,0) → physical (2,2,2), logical (2,2,2) → physical (4,4,4)
     * Examples for n=4: logical (0,0,0) → physical (2,2,2), logical (1,1,1) → physical (3,3,3)
     * 
     * @param x Logical X coordinate (can be negative, 0=center).
     * @param y Logical Y coordinate (can be negative, 0=center).
     * @param z Logical Z coordinate (can be negative, 0=center).
     * @return Tuple of (i,j,k) physical grid indices for array access.
     * @throws std::runtime_error if no subdivision available.
     */
    std::tuple<int,int,int> logicalToPhysical(int x, int y, int z) const {
        if (!hasSubdivision()) {  // Ensure subdivision data exists and is valid
            throw std::runtime_error("Cube::logicalToPhysical: no subdivision available");
        }

        int center = subdivision_levels_ / 2;  // Integer division gives center offset for grid
        int i = center + x;  // Convert logical X to physical array index
        int j = center + y;  // Convert logical Y to physical array index
        int k = center + z;  // Convert logical Z to physical array index

        return std::make_tuple(i, j, k);  // Return as tuple for structured binding support
    }

    /**
     * @brief Convert internal grid indices to user-friendly logical coordinates.
     * 
     * Transforms physical array indices back to logical coordinates centered at origin.
     * This is the inverse operation of logicalToPhysical() and is useful for
     * converting internal processing results back to user coordinate space.
     * 
     * @param i Physical X grid index in range [0, n).
     * @param j Physical Y grid index in range [0, n).
     * @param k Physical Z grid index in range [0, n).
     * @return Tuple of (x,y,z) logical coordinates centered at origin.
     * @throws std::runtime_error if no subdivision available.
     */
    std::tuple<int,int,int> physicalToLogical(int i, int j, int k) const {
        if (!hasSubdivision()) {  // Ensure subdivision data exists and is valid
            throw std::runtime_error("Cube::physicalToLogical: no subdivision available");
        }

        int center = subdivision_levels_ / 2;  // Calculate center offset for conversion
        int x = i - center;  // Convert physical array index to logical X
        int y = j - center;  // Convert physical array index to logical Y
        int z = k - center;  // Convert physical array index to logical Z

        return std::make_tuple(x, y, z);  // Return as tuple for structured binding support
    }

    /* === SUBCELL ACCESS === */
    
    /**
     * @brief Access a specific subcell using logical coordinates.
     * 
     * Provides read-only access to a subcell using user-friendly coordinates
     * where (0,0,0) represents the center subcell. Automatically handles
     * coordinate transformation and bounds checking.
     * 
     * @param x Logical X-coordinate (can be negative, 0=center).
     * @param y Logical Y-coordinate (can be negative, 0=center).
     * @param z Logical Z-coordinate (can be negative, 0=center).
     * @return Const reference to the specified subcell.
     * @throws std::out_of_range if coordinates invalid or no subdivision.
     */
    const SubCell& getSubCell(int x, int y, int z) const {
        auto [i, j, k] = logicalToPhysical(x, y, z);  // Convert to internal grid indices

        // Validate that converted coordinates are within valid grid bounds
        if (i < 0 || i >= subdivision_levels_ ||
            j < 0 || j >= subdivision_levels_ ||
            k < 0 || k >= subdivision_levels_) {
            throw std::out_of_range("Cube::getSubCell: logical coordinates out of range");
        }
        return subcells_[i][j][k];  // Return reference to subcell in 3D grid
    }

    /**
     * @brief Access a specific subcell using physical grid indices.
     * 
     * Provides read-only access to a subcell using direct array indices.
     * This is primarily for internal use or when working directly with
     * the physical grid layout. Includes full bounds checking.
     * 
     * @param i Physical X-coordinate in grid range [0, n).
     * @param j Physical Y-coordinate in grid range [0, n).
     * @param k Physical Z-coordinate in grid range [0, n).
     * @return Const reference to the specified subcell.
     * @throws std::out_of_range if coordinates invalid or no subdivision.
     */
    const SubCell& getSubCellPhysical(int i, int j, int k) const {
        if (!hasSubdivision()) {  // Ensure subdivision data exists and is valid
            throw std::runtime_error("Cube::getSubCellPhysical: no subdivision available");
        }
        // Validate that all coordinates are within valid grid bounds [0, n)
        if (i < 0 || i >= subdivision_levels_ ||
            j < 0 || j >= subdivision_levels_ ||
            k < 0 || k >= subdivision_levels_) {
            throw std::out_of_range("Cube::getSubCellPhysical: coordinates out of range");
        }
        return subcells_[i][j][k];  // Direct array access to subcell in 3D grid
    }

    /**
     * @brief Get modifiable access to a subcell using logical coordinates.
     * 
     * Provides read-write access to a subcell for modifications such as
     * vertex updates, activation state changes, or other property modifications.
     * Changes may require calling refreshTriangulation() afterward.
     * 
     * @param x Logical X-coordinate (can be negative, 0=center).
     * @param y Logical Y-coordinate (can be negative, 0=center).
     * @param z Logical Z-coordinate (can be negative, 0=center).
     * @return Mutable reference to the specified subcell.
     * @throws std::out_of_range if coordinates invalid or no subdivision.
     */
    SubCell& getSubCellMutable(int x, int y, int z) {
        auto [i, j, k] = logicalToPhysical(x, y, z);  // Convert to internal grid indices

        // Validate that converted coordinates are within valid grid bounds
        if (i < 0 || i >= subdivision_levels_ ||
            j < 0 || j >= subdivision_levels_ ||
            k < 0 || k >= subdivision_levels_) {
            throw std::out_of_range("Cube::getSubCellMutable: logical coordinates out of range");
        }
        return subcells_[i][j][k];  // Return modifiable reference to subcell
    }

    /**
     * @brief Get modifiable access to a subcell using physical grid indices.
     * 
     * Provides read-write access to a subcell using direct array indices.
     * This is primarily for internal use or when working directly with
     * the physical grid layout during batch operations.
     * 
     * @param i Physical X-coordinate in grid range [0, n).
     * @param j Physical Y-coordinate in grid range [0, n).
     * @param k Physical Z-coordinate in grid range [0, n).
     * @return Mutable reference to the specified subcell.
     * @throws std::out_of_range if coordinates invalid or no subdivision.
     */
    SubCell& getSubCellMutablePhysical(int i, int j, int k) {
        if (!hasSubdivision()) {  // Ensure subdivision data exists and is valid
            throw std::runtime_error("Cube::getSubCellMutablePhysical: no subdivision available");
        }
        // Validate that all coordinates are within valid grid bounds [0, n)
        if (i < 0 || i >= subdivision_levels_ ||
            j < 0 || j >= subdivision_levels_ ||
            k < 0 || k >= subdivision_levels_) {
            throw std::out_of_range("Cube::getSubCellMutablePhysical: coordinates out of range");
        }
        return subcells_[i][j][k];  // Direct modifiable array access to subcell
    }

    /* === PLANE EXTRACTION === */
    
    /**
     * @brief Extract all subcells from a specific 2D plane within the 3D subdivision.
     * 
     * Returns a collection of subcell references that form a 2D slice through the
     * 3D subdivision grid. Useful for layer-by-layer processing, visualization,
     * or analysis of specific cross-sections through the subdivided cube.
     * 
     * @param axis Plane orientation: 0=YZ plane (varying X), 1=XZ plane (varying Y), 2=XY plane (varying Z).
     * @param layer Layer index within the specified axis [0, n).
     * @return Vector of const references to all subcells in the specified plane.
     * @throws std::out_of_range if axis or layer parameters are invalid.
     */
    std::vector<std::reference_wrapper<const SubCell>> getPlane(int axis, int layer) const {
        if (!hasSubdivision()) {  // Ensure subdivision data exists and is valid
            throw std::runtime_error("Cube::getPlane: no subdivision available");
        }
        // Validate axis (must be 0-2) and layer (must be within grid bounds)
        if (axis < 0 || axis > 2 || layer < 0 || layer >= subdivision_levels_) {
            throw std::out_of_range("Cube::getPlane: invalid axis or layer");
        }

        std::vector<std::reference_wrapper<const SubCell>> plane;  // Container for plane subcells
        int n = subdivision_levels_;  // Cache subdivision count for loop efficiency

        switch (axis) {
            case 0: // YZ plane (fixed X = layer) - slice perpendicular to X-axis
                for (int j = 0; j < n; ++j) {      // Iterate through all Y coordinates
                    for (int k = 0; k < n; ++k) {  // Iterate through all Z coordinates
                        plane.emplace_back(subcells_[layer][j][k]);  // Add subcell reference to plane
                    }
                }
                break;
            case 1: // XZ plane (fixed Y = layer) - slice perpendicular to Y-axis
                for (int i = 0; i < n; ++i) {      // Iterate through all X coordinates
                    for (int k = 0; k < n; ++k) {  // Iterate through all Z coordinates
                        plane.emplace_back(subcells_[i][layer][k]);  // Add subcell reference to plane
                    }
                }
                break;
            case 2: // XY plane (fixed Z = layer) - slice perpendicular to Z-axis
                for (int i = 0; i < n; ++i) {      // Iterate through all X coordinates
                    for (int j = 0; j < n; ++j) {  // Iterate through all Y coordinates
                        plane.emplace_back(subcells_[i][j][layer]);  // Add subcell reference to plane
                    }
                }
                break;
        }
        return plane;  // Return complete collection of plane subcell references
    }

    /* === VERTEX MANIPULATION === */
    
    /**
     * @brief Update a specific vertex within a subcell and mark for retriangulation.
     * 
     * Modifies the position of one vertex within a specified subcell. Since vertices
     * define the geometry of triangular faces, this operation typically requires
     * calling refreshTriangulation() afterward to update the visual representation.
     * 
     * @param x Logical X-coordinate of target subcell (0=center).
     * @param y Logical Y-coordinate of target subcell (0=center).
     * @param z Logical Z-coordinate of target subcell (0=center).
     * @param vertex_idx Vertex index within subcell [0,7] (standard cube vertex ordering).
     * @param new_pos New 3D position for the specified vertex.
     * @throws std::out_of_range if coordinates invalid, vertex_idx invalid, or no subdivision.
     */
    void updateSubCellVertex(int x, int y, int z, int vertex_idx, const Vector3D& new_pos) {
        if (vertex_idx < 0 || vertex_idx > 7) {  // Validate vertex index (cube has 8 vertices: 0-7)
            throw std::out_of_range("Cube::updateSubCellVertex: vertex_idx must be [0,7]");
        }
        SubCell& cell = getSubCellMutable(x, y, z);  // Get modifiable reference to target subcell
        cell.vertices[vertex_idx] = new_pos;         // Update the specified vertex position

        // Note: Triangulation refresh will be handled by explicit refreshTriangulation() call
        // This allows batch vertex updates before expensive triangulation regeneration
    }

    /**
     * @brief Regenerate all triangular faces after vertex modifications.
     * 
     * Completely rebuilds the triangulation from current subcell vertex positions.
     * This is necessary after any vertex position changes to ensure the visual
     * representation matches the current geometry. Can be expensive for large
     * subdivisions, so batch vertex updates before calling this method.
     */
    void refreshTriangulation() {
        if (!hasSubdivision()) return;  // Exit early if no subdivision exists

        facets_.clear();  // Remove all existing triangular faces from collection

        // Rebuild triangles from current subcell vertices using physical grid traversal
        for (int i = 0; i < subdivision_levels_; ++i) {      // Iterate through X dimension
            for (int j = 0; j < subdivision_levels_; ++j) {  // Iterate through Y dimension
                for (int k = 0; k < subdivision_levels_; ++k) {  // Iterate through Z dimension
                    const SubCell& cell = subcells_[i][j][k];  // Get reference to current subcell
                    if (!cell.active) continue;  // Skip disabled subcells for selective rendering

                    // Generate 12 triangular faces from current vertex positions
                    for (int t = 0; t < 12; ++t) {  // Each cube generates 12 triangular faces
                        auto const& tri = cube_triangles_[t];  // Get vertex indices for triangle t
                        // Create triangle from current vertices and add to face collection
                        facets_.push(cell.vertices[tri[0]], cell.vertices[tri[1]], cell.vertices[tri[2]]);
                    }
                }
            }
        }
    }

    /* === FACET GENERATION === */
    
    /**
     * @brief Generate triangular faces for a specific subcell.
     * 
     * Creates a FacetBox containing the 12 triangular faces that represent
     * the visual geometry of a single subcell. Useful for rendering individual
     * subcells or creating partial visualizations of the subdivision.
     * 
     * @param x Logical X-coordinate of target subcell (0=center).
     * @param y Logical Y-coordinate of target subcell (0=center).
     * @param z Logical Z-coordinate of target subcell (0=center).
     * @return FacetBox containing 12 triangular faces (or empty if subcell inactive).
     * @throws std::out_of_range if coordinates invalid or no subdivision.
     */
    FacetBox getSubCellFacets(int x, int y, int z) const {
        const SubCell& cell = getSubCell(x, y, z);  // Get reference to target subcell
        FacetBox subcell_facets;  // Container for this subcell's triangular faces

        if (!cell.active) return subcell_facets; // Return empty collection if subcell disabled

        // Generate 12 triangular faces from the subcell's current 8 vertices
        for (int t = 0; t < 12; ++t) {  // Each cube generates 12 triangular faces
            auto const& tri = cube_triangles_[t];  // Get vertex indices from lookup table
            // Create triangle from 3 vertices and add to subcell face collection
            subcell_facets.push(cell.vertices[tri[0]], cell.vertices[tri[1]], cell.vertices[tri[2]]);
        }
        return subcell_facets;  // Return complete set of subcell triangular faces
    }

    /**
     * @brief Generate triangular faces for all subcells within a 2D plane.
     * 
     * Creates a FacetBox containing all triangular faces from subcells that
     * lie within the specified plane. Useful for layer-by-layer rendering
     * or creating cross-sectional views of the subdivided cube.
     * 
     * @param axis Plane orientation: 0=YZ plane, 1=XZ plane, 2=XY plane.
     * @param layer Layer index within the specified axis [0, n).
     * @return FacetBox containing triangular faces from all active subcells in the plane.
     * @throws std::out_of_range if axis or layer parameters are invalid.
     */
    FacetBox getPlaneFacets(int axis, int layer) const {
        auto plane_subcells = getPlane(axis, layer);  // Get all subcells within the specified plane
        FacetBox plane_facets;  // Container for all triangular faces in the plane

        // Process each subcell within the plane
        for (const SubCell& cell : plane_subcells) {
            if (!cell.active) continue;  // Skip disabled subcells for selective rendering

            // Generate 12 triangular faces for each active subcell in the plane
            for (int t = 0; t < 12; ++t) {  // Each cube generates 12 triangular faces
                auto const& tri = cube_triangles_[t];  // Get vertex indices from lookup table
                // Create triangle from 3 vertices and add to plane face collection
                plane_facets.push(cell.vertices[tri[0]], cell.vertices[tri[1]], cell.vertices[tri[2]]);
            }
        }
        return plane_facets;  // Return complete set of plane triangular faces
    }

    /* === PATTERN-BASED SUBCELL SELECTION === */
    
    /**
     * @brief Extract subcells following a 3D checkerboard pattern.
     * 
     * Returns subcells where the sum of physical coordinates (i+j+k) is divisible by 4,
     * creating a 3D checkerboard pattern similar to alternating squares on a 2D board
     * but extended to three dimensions. This pattern is useful for creating visual
     * effects, selective processing, or algorithmic operations on subcell subsets.
     * 
     * The checkerboard pattern ensures spatial distribution while reducing density,
     * making it useful for sampling, performance optimization, or artistic effects.
     * 
     * @return Vector of const references to subcells matching the checkerboard pattern.
     * @throws std::runtime_error if no subdivision available.
     */
    std::vector<std::reference_wrapper<const SubCell>> getCheckerboardSubcells() const {
        if (!hasSubdivision()) {  // Ensure subdivision data exists and is valid
            throw std::runtime_error("Cube::getCheckerboardSubcells: no subdivision available");
        }

        std::vector<std::reference_wrapper<const SubCell>> checkerboard;  // Container for pattern subcells
        int n = subdivision_levels_;  // Cache subdivision count for loop efficiency

        // Iterate through all subcells using physical grid coordinates
        for (int i = 0; i < n; ++i) {      // X dimension physical coordinates
            for (int j = 0; j < n; ++j) {  // Y dimension physical coordinates
                for (int k = 0; k < n; ++k) {  // Z dimension physical coordinates
                    // Select subcells where coordinate sum is divisible by 4 (checkerboard pattern)
                    if ((i + j + k) % 4 == 0) {
                        checkerboard.emplace_back(subcells_[i][j][k]);  // Add subcell reference to pattern
                    }
                }
            }
        }
        return checkerboard;  // Return complete collection of pattern-matching subcells
    }

    /**
     * @brief Generate triangular faces for all subcells in the checkerboard pattern.
     * 
     * Creates a FacetBox containing triangular faces from all subcells that match
     * the 3D checkerboard pattern. This creates a visually distinctive rendering
     * with spatially distributed geometry, useful for performance testing, visual
     * effects, or demonstrating pattern-based operations.
     * 
     * @return FacetBox containing triangular faces from all active checkerboard subcells.
     * @throws std::runtime_error if no subdivision available.
     */
    FacetBox getCheckerboardFacets() const {
        auto checkerboard_subcells = getCheckerboardSubcells();  // Get all pattern-matching subcells
        FacetBox checkerboard_facets;  // Container for pattern triangular faces

        // Process each subcell in the checkerboard pattern
        for (const SubCell& cell : checkerboard_subcells) {
            if (!cell.active) continue;  // Skip disabled subcells for selective rendering

            // Generate 12 triangular faces for each active pattern subcell
            for (int t = 0; t < 12; ++t) {  // Each cube generates 12 triangular faces
                auto const& tri = cube_triangles_[t];  // Get vertex indices from lookup table
                // Create triangle from 3 vertices and add to pattern face collection
                checkerboard_facets.push(cell.vertices[tri[0]], cell.vertices[tri[1]], cell.vertices[tri[2]]);
            }
        }
        return checkerboard_facets;  // Return complete set of pattern triangular faces
    }

    /* === STRING-BASED PLANE ACCESS === */
    
    /**
     * @brief Generate plane facets using string-based coordinate specification.
     * 
     * Provides a user-friendly interface for specifying planes using string coordinates.
     * Two coordinates should be axis names ("x", "y", "z") and one should be a numeric
     * value specifying the plane position in logical coordinates.
     * 
     * This method enables intuitive plane specification:
     * - getPlaneFacets("x", "y", "0") → center XY plane (perpendicular to Z axis)
     * - getPlaneFacets("x", "0", "z") → center XZ plane (perpendicular to Y axis)  
     * - getPlaneFacets("0", "y", "z") → center YZ plane (perpendicular to X axis)
     * 
     * @param coord1 First coordinate: axis name ("x", "y", "z") or numeric string.
     * @param coord2 Second coordinate: axis name ("x", "y", "z") or numeric string.
     * @param coord3 Third coordinate: axis name ("x", "y", "z") or numeric string.
     * @return FacetBox containing triangular faces from the specified plane.
     * @throws std::invalid_argument if coordinate specification is invalid.
     * @throws std::out_of_range if numeric coordinate is outside valid range.
     */
    FacetBox getPlaneFacets(const std::string& coord1, const std::string& coord2, const std::string& coord3) const {
        if (!hasSubdivision()) {  // Ensure subdivision data exists and is valid
            throw std::runtime_error("Cube::getPlaneFacets: no subdivision available");
        }

        // Helper lambda to validate if string contains only numeric characters
        auto isNumeric = [](const std::string& str) -> bool {
            if (str.empty()) return false;  // Empty strings are not numeric
            size_t start = (str[0] == '-') ? 1 : 0;  // Skip negative sign if present
            return start < str.size() && std::all_of(str.begin() + start, str.end(), ::isdigit);
        };

        // Helper lambda to convert string to integer with range validation
        auto parseCoord = [this](const std::string& str) -> int {
            int val = std::stoi(str);  // Convert string to integer (may throw)
            int max_coord = subdivision_levels_ / 2;  // Maximum valid logical coordinate
            int min_coord = -max_coord;               // Minimum valid logical coordinate
            if (val < min_coord || val > max_coord) {  // Validate coordinate is within bounds
                throw std::out_of_range("Coordinate " + str + " out of range [" +
                                      std::to_string(min_coord) + ", " + std::to_string(max_coord) + "]");
            }
            return val;  // Return validated coordinate value
        };

        std::vector<std::string> coords = {coord1, coord2, coord3};  // Collect all coordinates for analysis
        int numeric_count = 0;   // Count of numeric coordinates (should be exactly 1)
        int numeric_idx = -1;    // Index of the numeric coordinate
        int numeric_value = 0;   // Parsed value of the numeric coordinate

        // Analyze coordinate specification to find the single numeric coordinate
        for (int i = 0; i < 3; ++i) {
            if (isNumeric(coords[i])) {  // Check if current coordinate is numeric
                numeric_count++;         // Increment numeric coordinate counter
                numeric_idx = i;         // Store index of numeric coordinate
                numeric_value = parseCoord(coords[i]);  // Parse and validate numeric value
            }
        }

        if (numeric_count != 1) {  // Validate exactly one numeric coordinate exists
            throw std::invalid_argument("Exactly one coordinate must be numeric, got " + std::to_string(numeric_count));
        }

        // Determine axis orientation and validate coordinate specification
        int axis;  // Will store axis perpendicular to the desired plane
        if (numeric_idx == 0) {
            // X is fixed (numeric), varying Y and Z → YZ plane (axis = 0, perpendicular to X)
            if (coords[1] != "y" || coords[2] != "z") {
                throw std::invalid_argument("For fixed X coordinate, other coordinates must be 'y' and 'z'");
            }
            axis = 0;  // YZ plane (perpendicular to X-axis)
        } else if (numeric_idx == 1) {
            // Y is fixed (numeric), varying X and Z → XZ plane (axis = 1, perpendicular to Y)
            if (coords[0] != "x" || coords[2] != "z") {
                throw std::invalid_argument("For fixed Y coordinate, other coordinates must be 'x' and 'z'");
            }
            axis = 1;  // XZ plane (perpendicular to Y-axis)
        } else {
            // Z is fixed (numeric), varying X and Y → XY plane (axis = 2, perpendicular to Z)
            if (coords[0] != "x" || coords[1] != "y") {
                throw std::invalid_argument("For fixed Z coordinate, other coordinates must be 'x' and 'y'");
            }
            axis = 2;  // XY plane (perpendicular to Z-axis)
        }

        // Convert logical coordinate to physical layer index for internal access
        int center = subdivision_levels_ / 2;    // Calculate center offset for coordinate system
        int layer = center + numeric_value;      // Convert logical to physical coordinate

        if (layer < 0 || layer >= subdivision_levels_) {  // Validate layer is within grid bounds
            throw std::out_of_range("Calculated layer " + std::to_string(layer) + " out of range [0, " +
                                  std::to_string(subdivision_levels_) + ")");
        }

        return getPlaneFacets(axis, layer);  // Delegate to numeric plane access method
    }

    /* === SUBCELL PROPERTY ACCESS === */
    
    /**
     * @brief Get the geometric center of a specific subcell.
     * 
     * Returns the 3D position of the subcell's center point, which is useful
     * for positioning, distance calculations, or spatial analysis operations.
     * 
     * @param x Logical X-coordinate of target subcell (0=center).
     * @param y Logical Y-coordinate of target subcell (0=center).
     * @param z Logical Z-coordinate of target subcell (0=center).
     * @return Vector3D representing the exact center position of the subcell.
     * @throws std::out_of_range if coordinates invalid or no subdivision.
     */
    Vector3D getSubcellCenter(int x, int y, int z) const {
        return getSubCell(x, y, z).center;  // Return center point directly from subcell data
    }

    /**
     * @brief Get the corner-to-center distance of a specific subcell.
     * 
     * Calculates the distance from the subcell center to any of its corners.
     * This is useful for collision detection, bounding sphere calculations,
     * or determining the maximum extent of a subcell from its center.
     * 
     * Note: This returns the 3D diagonal distance, not the face-to-center distance.
     * 
     * @param x Logical X-coordinate of target subcell (0=center).
     * @param y Logical Y-coordinate of target subcell (0=center).
     * @param z Logical Z-coordinate of target subcell (0=center).
     * @return double representing distance from center to corner.
     * @throws std::out_of_range if coordinates invalid or no subdivision.
     */
    double getSubcellRadius(int x, int y, int z) const {
        const SubCell& cell = getSubCell(x, y, z);  // Get reference to target subcell
        // The radius stored is half the side length, but corner distance requires 3D diagonal
        // For a cube: corner_distance = sqrt(3) * half_side_length
        return cell.radius * sqrt(3.0);  // Convert from half-side-length to corner distance
    }

    /* === SUBCELL STATE MANAGEMENT === */
    
    /**
     * @brief Enable or disable a specific subcell for rendering and processing.
     * 
     * Controls whether a subcell participates in triangulation, rendering, and
     * other operations. Disabled subcells are skipped during facet generation,
     * enabling selective visualization and processing of subdivision regions.
     * 
     * @param x Logical X-coordinate of target subcell (0=center).
     * @param y Logical Y-coordinate of target subcell (0=center).
     * @param z Logical Z-coordinate of target subcell (0=center).
     * @param active True to enable subcell, false to disable.
     * @throws std::out_of_range if coordinates invalid or no subdivision.
     */
    void setSubCellActive(int x, int y, int z, bool active) {
        getSubCellMutable(x, y, z).active = active;  // Set the active/inactive state flag
    }

private:
    // === CORE DATA MEMBERS ===
    std::vector<Quaternion> verts_;  // 8 corner vertices of the main cube (using quaternions for storage)
    FacetBox facets_;               // Collection of all triangular faces for rendering/processing

    // === SUBDIVISION DATA MEMBERS ===
    // 3D grid storage: subcells_[i][j][k] - indexed by physical coordinates [0,n) × [0,n) × [0,n)
    std::vector<std::vector<std::vector<SubCell>>> subcells_;
    int subdivision_levels_;   // Current subdivision level (n means n×n×n total subcells)
    bool subdivision_built_;   // Flag indicating if subdivision data has been generated and is valid

    // === TRIANGULATION LOOKUP TABLE ===
    // Static lookup table defining the 12 triangular faces of a standard cube
    // Each entry contains 3 vertex indices that form one triangle
    // Based on standard cube vertex ordering: vertices 0-3 front face, 4-7 back face
    static constexpr int cube_triangles_[12][3] = {
        // Front face (z = +radius) - divided into 2 triangles from vertices 0,1,2,3
        {0, 1, 2}, {2, 3, 0},  // Triangle 1: vertices (0,1,2), Triangle 2: vertices (2,3,0)
        // Back face (z = -radius) - divided into 2 triangles from vertices 4,5,6,7  
        {4, 6, 5}, {6, 4, 7},  // Triangle 3: vertices (4,6,5), Triangle 4: vertices (6,4,7)
        // Left face (x = -radius) - divided into 2 triangles from vertices 0,3,4,7
        {4, 0, 3}, {3, 7, 4},  // Triangle 5: vertices (4,0,3), Triangle 6: vertices (3,7,4)
        // Right face (x = +radius) - divided into 2 triangles from vertices 1,2,5,6
        {1, 5, 6}, {6, 2, 1},  // Triangle 7: vertices (1,5,6), Triangle 8: vertices (6,2,1)
        // Bottom face (y = -radius) - divided into 2 triangles from vertices 0,1,4,5
        {4, 5, 1}, {1, 0, 4},  // Triangle 9: vertices (4,5,1), Triangle 10: vertices (1,0,4)
        // Top face (y = +radius) - divided into 2 triangles from vertices 2,3,6,7
        {3, 2, 6}, {6, 7, 3}   // Triangle 11: vertices (3,2,6), Triangle 12: vertices (6,7,3)
    };

    /* === PRIVATE HELPER METHODS === */

    /**
     * @brief Initialize the 8 fundamental vertices of a cube.
     * 
     * Creates the 8 corner points that define basic cube geometry using standard
     * vertex ordering. The ordering is critical for consistent triangulation
     * across all cubes and subcubes in the system. Vertices are stored as
     * quaternions with zero scalar part for compatibility with the rotation system.
     * 
     * @param r Half the side length (distance from center to face).
     * @param center The center position of the cube in 3D space.
     */
    void initVertices(double r, const Vector3D& center) {
        verts_.clear();     // Remove any existing vertex data
        verts_.reserve(8);  // Pre-allocate space for exactly 8 vertices for efficiency

        // Generate all 8 vertices using standard cube vertex ordering
        // This ordering ensures consistent triangulation and normal direction
        verts_.emplace_back(0.0, Vector3D(center.x() - r, center.y() - r, center.z() + r)); // 0: front-bottom-left
        verts_.emplace_back(0.0, Vector3D(center.x() + r, center.y() - r, center.z() + r)); // 1: front-bottom-right  
        verts_.emplace_back(0.0, Vector3D(center.x() + r, center.y() + r, center.z() + r)); // 2: front-top-right
        verts_.emplace_back(0.0, Vector3D(center.x() - r, center.y() + r, center.z() + r)); // 3: front-top-left
        verts_.emplace_back(0.0, Vector3D(center.x() - r, center.y() - r, center.z() - r)); // 4: back-bottom-left
        verts_.emplace_back(0.0, Vector3D(center.x() + r, center.y() - r, center.z() - r)); // 5: back-bottom-right
        verts_.emplace_back(0.0, Vector3D(center.x() + r, center.y() + r, center.z() - r)); // 6: back-top-right
        verts_.emplace_back(0.0, Vector3D(center.x() - r, center.y() + r, center.z() - r)); // 7: back-top-left
    }

    /**
     * @brief Generate the 12 triangular faces for basic cube rendering.
     * 
     * Converts the 8 stored vertices into 12 triangular faces using the static
     * triangulation lookup table. This creates the standard cube triangulation
     * where each face is divided into 2 triangles, providing consistent geometry
     * for rendering systems that work with triangular primitives.
     */
    void buildFacets() {
        facets_.clear();  // Remove any existing triangular faces

        // Generate all 12 triangular faces using the static lookup table
        for (int i = 0; i < 12; ++i) {  // Process each triangular face
            auto const& t = cube_triangles_[i];  // Get vertex indices for triangle i
            // Create triangular face from 3 vertices and add to collection
            Facet face(verts_[t[0]], verts_[t[1]], verts_[t[2]]);
            facets_.push(face);  // Add completed triangle to rendering collection
        }
    }

    /**
     * @brief Core subdivision algorithm that generates n×n×n subcubes.
     * 
     * This is the heart of the subdivision system. It replaces the simple cube
     * triangulation with a complex grid of smaller cubes, each maintaining its
     * own geometry and properties. Each subcube generates 12 triangular faces,
     * and the complete subdivision creates a detailed mesh suitable for complex
     * geometric operations, deformation, and selective rendering.
     * 
     * @param center Center position of the overall cube structure.
     * @param radius Half side length of the overall cube structure.
     * @param n Number of subdivisions per axis (creates n³ total subcubes).
     */
    void buildSubdividedCube(const Vector3D& center, double radius, int n) {
        double step = (2.0 * radius) / n; // Calculate size of each subcube (total cube divided by n)
        double subRadius = step / 2.0;     // Half side length of each subcube

        // Initialize the 3D subcell grid with proper dimensions
        subcells_.clear();  // Remove any existing subdivision data
        subcells_.resize(n, std::vector<std::vector<SubCell>>(n, std::vector<SubCell>(n)));

        // Generate n³ subcubes in a systematic grid pattern
        for (int i = 0; i < n; ++i) {      // X dimension: left to right
            for (int j = 0; j < n; ++j) {  // Y dimension: bottom to top
                for (int k = 0; k < n; ++k) {  // Z dimension: back to front
                    // Calculate center position for this specific subcube
                    Vector3D subCenter(
                        center.x() - radius + step * (i + 0.5),  // X: start at left edge + half-step offset
                        center.y() - radius + step * (j + 0.5),  // Y: start at bottom edge + half-step offset  
                        center.z() - radius + step * (k + 0.5)   // Z: start at back edge + half-step offset
                    );

                    // Generate 8 vertices for this subcube using standard ordering
                    std::array<Vector3D, 8> subVerts = {
                        Vector3D(subCenter.x() - subRadius, subCenter.y() - subRadius, subCenter.z() + subRadius), // 0
                        Vector3D(subCenter.x() + subRadius, subCenter.y() - subRadius, subCenter.z() + subRadius), // 1
                        Vector3D(subCenter.x() + subRadius, subCenter.y() + subRadius, subCenter.z() + subRadius), // 2
                        Vector3D(subCenter.x() - subRadius, subCenter.y() + subRadius, subCenter.z() + subRadius), // 3
                        Vector3D(subCenter.x() - subRadius, subCenter.y() - subRadius, subCenter.z() - subRadius), // 4
                        Vector3D(subCenter.x() + subRadius, subCenter.y() - subRadius, subCenter.z() - subRadius), // 5
                        Vector3D(subCenter.x() + subRadius, subCenter.y() + subRadius, subCenter.z() - subRadius), // 6
                        Vector3D(subCenter.x() - subRadius, subCenter.y() + subRadius, subCenter.z() - subRadius)  // 7
                    };

                    // Store complete subcell data in the 3D grid for later access
                    subcells_[i][j][k] = SubCell(subVerts, subCenter, subRadius, i, j, k);

                    // Generate immediate triangulation for rendering - 12 triangular faces per subcube
                    for (int t = 0; t < 12; ++t) {  // Process each triangular face
                        auto const& tri = cube_triangles_[t];  // Get vertex indices from lookup table
                        // Create triangle from subcube vertices and add to main face collection
                        facets_.push(subVerts[tri[0]], subVerts[tri[1]], subVerts[tri[2]]);
                    }
                }
            }
        }
    }

    /**
     * @brief Calculate the radius (half side length) of the current cube.
     * 
     * Determines the cube's size by measuring the distance between adjacent vertices.
     * This is used internally for subdivision calculations and geometric operations.
     * 
     * @return Half the side length of the cube (distance from center to face).
     */
    double getRadius() const {
        if (verts_.size() < 2) return 0.0;  // Return zero for uninitialized cube
        Vector3D v0 = verts_[0].V();        // Get first vertex position
        Vector3D v1 = verts_[1].V();        // Get second vertex position  
        return abs(v1.x() - v0.x()) / 2.0;  // Calculate half the distance between adjacent vertices
    }
};

//-----------------------------------------------------------------------------
// Stream output operator for Cube debugging and inspection
//-----------------------------------------------------------------------------
/**
 * @brief Stream output operator for comprehensive cube information display.
 * 
 * Provides detailed debugging output showing all vertices and triangular faces
 * in the cube. This is invaluable for development, testing, and troubleshooting
 * geometric operations. The output includes vertex coordinates and face definitions.
 */
inline std::ostream& operator<<(std::ostream& os, const Cube& cube) {
    os << "--- Cube Vertices ---\n";  // Header for vertex information section
    for (size_t i = 0; i < cube.verts_.size(); ++i) {  // Display each vertex with index
        os << "  [" << i << "] " << cube.verts_[i].V() << "\n";  // Show index and 3D position
    }
    os << "--- Cube Facets (" << cube.facets_.size() << " triangles) ---\n";  // Header with count
    for (size_t i = 0; i < cube.facets_.size(); ++i) {  // Display each triangular face
        os << "  Face " << i << ": " << cube.facets_[i] << "\n";  // Show face index and geometry
    }
    return os;  // Return stream for chaining
}

#endif // CUBE_H
