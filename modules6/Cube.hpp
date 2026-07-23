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

    /**
     * @brief Write cube geometry to STL file format for 3D printing or CAD applications.
     *
     * Exports all triangular faces of the cube (basic or subdivided) to an STL file.
     * The STL format stores triangular mesh data with surface normals, making it
     * suitable for 3D printing, CAD software, and mesh processing applications.
     * Only active subcells are included in subdivided cubes.
     *
     * @param filename Path to output STL file (will be created or overwritten).
     * @throws std::runtime_error if file cannot be created or written.
     */
    void writeSTL(const std::string& filename) const {
        std::ofstream stl(filename);  // Create output file stream
        if (!stl.is_open()) {  // Validate file was successfully opened
            throw std::runtime_error("Cube::writeSTL: Cannot create file " + filename);
        }
    
        stl << "solid Cube\n";  // STL header with object name
    
        // Write all triangular faces from the cube's facet collection
        for (size_t i = 0; i < facets_.size(); ++i) {
            const Facet& face = facets_[i];  // Get current triangular face
    
            // Extract the three vertices from the facet using operator[]
            Vector3D v1 = face[0];  // Vertex A (index 0)
            Vector3D v2 = face[1];  // Vertex B (index 1)
            Vector3D v3 = face[2];  // Vertex C (index 2)
    
            // Get the precomputed normal or calculate it
            Vector3D normal = face.getNormal();  // Use the stored normal
    
            // Write STL facet in standard format
            stl << "facet normal ";
            stl << std::scientific << normal.x() << " " << normal.y() << " " << normal.z() << "\n";
            stl << "\touter loop\n";
            stl << "\t\tvertex " << v1.x() << " " << v1.y() << " " << v1.z() << "\n";
            stl << "\t\tvertex " << v2.x() << " " << v2.y() << " " << v2.z() << "\n";
            stl << "\t\tvertex " << v3.x() << " " << v3.y() << " " << v3.z() << "\n";
            stl << "\tendloop\n";
            stl << "endfacet\n";
        }
    
        stl << "endsolid Cube\n";  // STL footer
        stl.close();  // Ensure file is properly closed
    }

    /**
     * @brief Write multiple cube geometries to a single STL file.
     * 
     * Exports all triangular faces from a collection of cubes into one STL file.
     * This is useful for creating complex multi-cube models, architectural structures,
     * or batch processing of cube collections for 3D printing applications.
     * 
     * @param filename Path to output STL file (will be created or overwritten).
     * @param cubes Vector of cube objects to include in the STL file.
     * @param objectName Name to use in STL header (default: "MultiCube").
     * @throws std::runtime_error if file cannot be created or written.
     */
    static void writeMultiSTL(const std::string& filename, 
                             const std::vector<Cube>& cubes, 
                             const std::string& objectName = "MultiCube") {
        std::ofstream stl(filename);  // Create output file stream
        if (!stl.is_open()) {  // Validate file was successfully opened
            throw std::runtime_error("Cube::writeMultiSTL: Cannot create file " + filename);
        }
        
        stl << "solid " << objectName << "\n";  // STL header with custom object name
        
        // Write all triangular faces from each cube in the collection
        for (size_t cubeIndex = 0; cubeIndex < cubes.size(); ++cubeIndex) {
            const Cube& cube = cubes[cubeIndex];  // Get current cube
            
            // Process all faces in this cube
            for (size_t i = 0; i < cube.facets_.size(); ++i) {
                const Facet& face = cube.facets_[i];  // Get current triangular face
                
                // Extract the three vertices from the facet using operator[]
                Vector3D v1 = face[0];  // Vertex A (index 0)
                Vector3D v2 = face[1];  // Vertex B (index 1)
                Vector3D v3 = face[2];  // Vertex C (index 2)
                
                // Get the precomputed normal
                Vector3D normal = face.getNormal();  // Use the stored normal
                
                // Write STL facet in standard format
                stl << "facet normal ";
                stl << std::scientific << normal.x() << " " << normal.y() << " " << normal.z() << "\n";
                stl << "\touter loop\n";
                stl << "\t\tvertex " << v1.x() << " " << v1.y() << " " << v1.z() << "\n";
                stl << "\t\tvertex " << v2.x() << " " << v2.y() << " " << v2.z() << "\n";
                stl << "\t\tvertex " << v3.x() << " " << v3.y() << " " << v3.z() << "\n";
                stl << "\tendloop\n";
                stl << "endfacet\n";
            }
        }
        
        stl << "endsolid " << objectName << "\n";  // STL footer
        stl.close();  // Ensure file is properly closed
    }

    /**
     * @brief Write cube geometry to STL file with pattern support.
     *
     * Exports triangular faces from the current cube into an STL file.
     * Supports different extraction modes including full geometry, checkerboard patterns,
     * and plane slices for creating complex geometric models.
     *
     * @param filename Path to output STL file (will be created or overwritten).
     * @param objectName Name to use in STL header (default: "Cube").
     * @param mode Extraction mode: "full", "checkerboard", "plane_xy", "plane_xz", "plane_yz".
     * @param layer Layer index for plane modes (ignored for other modes, default: 0).
     * @throws std::runtime_error if file cannot be created or written.
     * @throws std::invalid_argument if mode is invalid or plane layer is out of range.
     */
    void writeSTL_s(const std::string& filename,
                    const std::string& objectName = "Cube",
                    const std::string& mode = "full",
                    int layer = 0) const {
        std::ofstream stl(filename);
        if (!stl.is_open()) {
            throw std::runtime_error("Cube::writeSTL_s: Cannot create file " + filename);
        }
    
        stl << "solid " << objectName << "\n";
    
        FacetBox facetsToWrite;  // Container for facets to export
    
        // Select facets based on extraction mode
        if (mode == "full") {
            // Use all facets from the cube
            facetsToWrite = getFacets();
    
        } else if (mode == "checkerboard") {
            // Use checkerboard pattern (only works with subdivided cubes)
            if (!hasSubdivision()) {
                // For non-subdivided cubes, fall back to full geometry
                facetsToWrite = getFacets();
            } else {
                facetsToWrite = getCheckerboardFacets();
            }
    
        } else if (mode == "plane_xy") {
            // Extract XY plane (perpendicular to Z-axis)
            if (!hasSubdivision()) {
                throw std::invalid_argument("Plane extraction requires subdivided cube");
            }
            facetsToWrite = getPlaneFacets(2, layer);  // axis=2 (Z), layer
    
        } else if (mode == "plane_xz") {
            // Extract XZ plane (perpendicular to Y-axis)
            if (!hasSubdivision()) {
                throw std::invalid_argument("Plane extraction requires subdivided cube");
            }
            facetsToWrite = getPlaneFacets(1, layer);  // axis=1 (Y), layer
    
        } else if (mode == "plane_yz") {
            // Extract YZ plane (perpendicular to X-axis)
            if (!hasSubdivision()) {
                throw std::invalid_argument("Plane extraction requires subdivided cube");
            }
            facetsToWrite = getPlaneFacets(0, layer);  // axis=0 (X), layer
    
        } else {
            throw std::invalid_argument("Invalid mode: " + mode +
                ". Valid modes: full, checkerboard, plane_xy, plane_xz, plane_yz");
        }
    
        // Write all selected facets to STL
        for (size_t i = 0; i < facetsToWrite.size(); ++i) {
            const Facet& face = facetsToWrite[i];
    
            // Extract vertices
            Vector3D v1 = face[0];  // Vertex A
            Vector3D v2 = face[1];  // Vertex B
            Vector3D v3 = face[2];  // Vertex C
    
            // Get normal
            Vector3D normal = face.getNormal();
    
            // Write STL facet
            stl << "facet normal ";
            stl << std::scientific << normal.x() << " " << normal.y() << " " << normal.z() << "\n";
            stl << "\touter loop\n";
            stl << "\t\tvertex " << v1.x() << " " << v1.y() << " " << v1.z() << "\n";
            stl << "\t\tvertex " << v2.x() << " " << v2.y() << " " << v2.z() << "\n";
            stl << "\t\tvertex " << v3.x() << " " << v3.y() << " " << v3.z() << "\n";
            stl << "\tendloop\n";
            stl << "endfacet\n";
        }
    
        stl << "endsolid " << objectName << "\n";
        stl.close();
    }

    /**
     * @brief Write cube geometry to STL file with pattern support including parameterized checkerboard.
     *
     * Exports triangular faces from the current cube into an STL file.
     * Supports different extraction modes including full geometry, parameterized checkerboard patterns,
     * and plane slices for creating complex geometric models.
     *
     * @param filename Path to output STL file (will be created or overwritten).
     * @param objectName Name to use in STL header (default: "Cube").
     * @param mode Extraction mode: "full", "checkerboard", "plane_xy", "plane_xz", "plane_yz".
     * @param x X-axis modulo parameter for checkerboard (default: 2), or layer for plane modes.
     * @param y Y-axis modulo parameter for checkerboard (default: 9).
     * @param z Z-axis modulo parameter for checkerboard (default: 2).
     * @throws std::runtime_error if file cannot be created or written.
     * @throws std::invalid_argument if mode is invalid or parameters are out of range.
     */
    void writeSTL_s(const std::string& filename,
                    const std::string& objectName,
                    const std::string& mode,
                    int x,
                    int y = 9,
                    int z = 2) const {
        std::ofstream stl(filename);
        if (!stl.is_open()) {
            throw std::runtime_error("Cube::writeSTL_s: Cannot create file " + filename);
        }
    
        stl << "solid " << objectName << "\n";
    
        FacetBox facetsToWrite;
    
        if (mode == "full") {
            facetsToWrite = getFacets();
    
        } else if (mode == "checkerboard") {
            if (!hasSubdivision()) {
                facetsToWrite = getFacets();
            } else {
                // Use parameterized checkerboard with x, y, z as modulo parameters
                facetsToWrite = getCheckerboardFacets(x, y, z);
            }
    
        } else if (mode == "plane_xy") {
            if (!hasSubdivision()) {
                throw std::invalid_argument("Plane extraction requires subdivided cube");
            }
            facetsToWrite = getPlaneFacets(2, x);  // x is layer parameter for plane modes
    
        } else if (mode == "plane_xz") {
            if (!hasSubdivision()) {
                throw std::invalid_argument("Plane extraction requires subdivided cube");
            }
            facetsToWrite = getPlaneFacets(1, x);  // x is layer parameter for plane modes
    
        } else if (mode == "plane_yz") {
            if (!hasSubdivision()) {
                throw std::invalid_argument("Plane extraction requires subdivided cube");
            }
            facetsToWrite = getPlaneFacets(0, x);  // x is layer parameter for plane modes
    
        } else {
            throw std::invalid_argument("Invalid mode: " + mode +
                ". Valid modes: full, checkerboard, plane_xy, plane_xz, plane_yz");
        }
    
        // Write all selected facets to STL
        for (size_t i = 0; i < facetsToWrite.size(); ++i) {
            const Facet& face = facetsToWrite[i];
    
            Vector3D v1 = face[0];
            Vector3D v2 = face[1];
            Vector3D v3 = face[2];
    
            Vector3D normal = face.getNormal();
    
            stl << "facet normal ";
            stl << std::scientific << normal.x() << " " << normal.y() << " " << normal.z() << "\n";
            stl << "\touter loop\n";
            stl << "\t\tvertex " << v1.x() << " " << v1.y() << " " << v1.z() << "\n";
            stl << "\t\tvertex " << v2.x() << " " << v2.y() << " " << v2.z() << "\n";
            stl << "\t\tvertex " << v3.x() << " " << v3.y() << " " << v3.z() << "\n";
            stl << "\tendloop\n";
            stl << "endfacet\n";
        }
    
        stl << "endsolid " << objectName << "\n";
        stl.close();
    }

    /**
     * @brief Write multiple cube geometries to a single STL file with pattern support.
     *
     * Exports triangular faces from a collection of cubes into one STL file.
     * Supports different extraction modes including full geometry, checkerboard patterns,
     * and plane slices for creating complex multi-cube models.
     *
     * @param filename Path to output STL file (will be created or overwritten).
     * @param cubes Vector of cube objects to include in the STL file.
     * @param objectName Name to use in STL header (default: "MultiCube").
     * @param mode Extraction mode: "full", "checkerboard", "plane_xy", "plane_xz", "plane_yz".
     * @param layer Layer index for plane modes (ignored for other modes, default: 0).
     * @throws std::runtime_error if file cannot be created or written.
     * @throws std::invalid_argument if mode is invalid or plane layer is out of range.
     */
    static void writeMultiSTL_m(const std::string& filename,
                             const std::vector<Cube>& cubes,
                             const std::string& objectName = "MultiCube",
                             const std::string& mode = "full",
                             int layer = 0) {
        std::ofstream stl(filename);
        if (!stl.is_open()) {
            throw std::runtime_error("Cube::writeMultiSTL: Cannot create file " + filename);
        }
    
        stl << "solid " << objectName << "\n";
    
        // Process each cube in the collection
        for (size_t cubeIndex = 0; cubeIndex < cubes.size(); ++cubeIndex) {
            const Cube& cube = cubes[cubeIndex];
            FacetBox facetsToWrite;  // Container for facets from current cube
    
            // Select facets based on extraction mode
            if (mode == "full") {
                // Use all facets from the cube
                facetsToWrite = cube.getFacets();
    
            } else if (mode == "checkerboard") {
                // Use checkerboard pattern (only works with subdivided cubes)
                if (!cube.hasSubdivision()) {
                    // For non-subdivided cubes, fall back to full geometry
                    facetsToWrite = cube.getFacets();
                } else {
                    facetsToWrite = cube.getCheckerboardFacets();
                }
    
            } else if (mode == "plane_xy") {
                // Extract XY plane (perpendicular to Z-axis)
                if (!cube.hasSubdivision()) {
                    throw std::invalid_argument("Plane extraction requires subdivided cube");
                }
                facetsToWrite = cube.getPlaneFacets(2, layer);  // axis=2 (Z), layer
    
            } else if (mode == "plane_xz") {
                // Extract XZ plane (perpendicular to Y-axis)
                if (!cube.hasSubdivision()) {
                    throw std::invalid_argument("Plane extraction requires subdivided cube");
                }
                facetsToWrite = cube.getPlaneFacets(1, layer);  // axis=1 (Y), layer
    
            } else if (mode == "plane_yz") {
                // Extract YZ plane (perpendicular to X-axis)
                if (!cube.hasSubdivision()) {
                    throw std::invalid_argument("Plane extraction requires subdivided cube");
                }
                facetsToWrite = cube.getPlaneFacets(0, layer);  // axis=0 (X), layer
    
            } else {
                throw std::invalid_argument("Invalid mode: " + mode +
                    ". Valid modes: full, checkerboard, plane_xy, plane_xz, plane_yz");
            }
    
            // Write all selected facets to STL
            for (size_t i = 0; i < facetsToWrite.size(); ++i) {
                const Facet& face = facetsToWrite[i];
    
                // Extract vertices
                Vector3D v1 = face[0];  // Vertex A
                Vector3D v2 = face[1];  // Vertex B
                Vector3D v3 = face[2];  // Vertex C
    
                // Get normal
                Vector3D normal = face.getNormal();
    
                // Write STL facet
                stl << "facet normal ";
                stl << std::scientific << normal.x() << " " << normal.y() << " " << normal.z() << "\n";
                stl << "\touter loop\n";
                stl << "\t\tvertex " << v1.x() << " " << v1.y() << " " << v1.z() << "\n";
                stl << "\t\tvertex " << v2.x() << " " << v2.y() << " " << v2.z() << "\n";
                stl << "\t\tvertex " << v3.x() << " " << v3.y() << " " << v3.z() << "\n";
                stl << "\tendloop\n";
                stl << "endfacet\n";
            }
        }
    
        stl << "endsolid " << objectName << "\n";
        stl.close();
    }


    /**
     * @brief Write multiple cube geometries to a single STL file with parameterized pattern support.
     *
     * Exports triangular faces from a collection of cubes into one STL file.
     * Supports different extraction modes including full geometry, parameterized checkerboard patterns,
     * and plane slices for creating complex multi-cube models.
     *
     * @param filename Path to output STL file (will be created or overwritten).
     * @param cubes Vector of cube objects to include in the STL file.
     * @param objectName Name to use in STL header (default: "MultiCube").
     * @param mode Extraction mode: "full", "checkerboard", "plane_xy", "plane_xz", "plane_yz".
     * @param x X-axis modulo parameter for checkerboard (default: 2), or layer for plane modes.
     * @param y Y-axis modulo parameter for checkerboard (default: 9).
     * @param z Z-axis modulo parameter for checkerboard (default: 2).
     * @throws std::runtime_error if file cannot be created or written.
     * @throws std::invalid_argument if mode is invalid or parameters are out of range.
     */
    static void writeMultiSTL_m(const std::string& filename,
                             const std::vector<Cube>& cubes,
                             const std::string& objectName,
                             const std::string& mode,
                             int x,
                             int y = 9,
                             int z = 2) {
        std::ofstream stl(filename);
        if (!stl.is_open()) {
            throw std::runtime_error("Cube::writeMultiSTL_m: Cannot create file " + filename);
        }
    
        stl << "solid " << objectName << "\n";
    
        for (size_t cubeIndex = 0; cubeIndex < cubes.size(); ++cubeIndex) {
            const Cube& cube = cubes[cubeIndex];
            FacetBox facetsToWrite;
    
            if (mode == "full") {
                facetsToWrite = cube.getFacets();
    
            } else if (mode == "checkerboard") {
                if (!cube.hasSubdivision()) {
                    facetsToWrite = cube.getFacets();
                } else {
                    // Use parameterized checkerboard with x, y, z as modulo parameters
                    facetsToWrite = cube.getCheckerboardFacets(x, y, z);
                }
    
            } else if (mode == "plane_xy") {
                if (!cube.hasSubdivision()) {
                    throw std::invalid_argument("Plane extraction requires subdivided cube");
                }
                facetsToWrite = cube.getPlaneFacets(2, x);  // x is layer parameter
    
            } else if (mode == "plane_xz") {
                if (!cube.hasSubdivision()) {
                    throw std::invalid_argument("Plane extraction requires subdivided cube");
                }
                facetsToWrite = cube.getPlaneFacets(1, x);  // x is layer parameter
    
            } else if (mode == "plane_yz") {
                if (!cube.hasSubdivision()) {
                    throw std::invalid_argument("Plane extraction requires subdivided cube");
                }
                facetsToWrite = cube.getPlaneFacets(0, x);  // x is layer parameter
    
            } else {
                throw std::invalid_argument("Invalid mode: " + mode +
                    ". Valid modes: full, checkerboard, plane_xy, plane_xz, plane_yz");
            }
    
            // Write all selected facets to STL
            for (size_t i = 0; i < facetsToWrite.size(); ++i) {
                const Facet& face = facetsToWrite[i];
    
                Vector3D v1 = face[0];
                Vector3D v2 = face[1];
                Vector3D v3 = face[2];
    
                Vector3D normal = face.getNormal();
    
                stl << "facet normal ";
                stl << std::scientific << normal.x() << " " << normal.y() << " " << normal.z() << "\n";
                stl << "\touter loop\n";
                stl << "\t\tvertex " << v1.x() << " " << v1.y() << " " << v1.z() << "\n";
                stl << "\t\tvertex " << v2.x() << " " << v2.y() << " " << v2.z() << "\n";
                stl << "\t\tvertex " << v3.x() << " " << v3.y() << " " << v3.z() << "\n";
                stl << "\tendloop\n";
                stl << "endfacet\n";
            }
        }
    
        stl << "endsolid " << objectName << "\n";
        stl.close();
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
                    if (/*k%2 == 0 && i%2 == 0 && j%2 == 0*/ (i%2 == 0 && k%2 == 0) || j%9 == 0) {
                        checkerboard.emplace_back(subcells_[i][j][k]);  // Add subcell reference to pattern
                    }
                }
            }
        }
        return checkerboard;  // Return complete collection of pattern-matching subcells
    }

    /**
     * @brief Extract subcells following a parameterized 3D checkerboard pattern.
     * 
     * Returns subcells where the pattern (i%x == 0 && k%z == 0) || j%y == 0 is satisfied.
     * This allows customizable checkerboard patterns with different modulo values for
     * each axis, enabling flexible geometric sampling and visual effects.
     * 
     * @param x Modulo divisor for i-coordinate (X-axis pattern frequency).
     * @param y Modulo divisor for j-coordinate (Y-axis pattern frequency).
     * @param z Modulo divisor for k-coordinate (Z-axis pattern frequency).
     * @return Vector of const references to subcells matching the custom pattern.
     * @throws std::runtime_error if no subdivision available.
     * @throws std::invalid_argument if x, y, or z are less than 1.
     */
    std::vector<std::reference_wrapper<const SubCell>> getCheckerboardSubcells(int x, int y, int z) const {
        if (!hasSubdivision()) {
            throw std::runtime_error("Cube::getCheckerboardSubcells: no subdivision available");
        }
        
        if (x < 1 || y < 1 || z < 1) {
            throw std::invalid_argument("Cube::getCheckerboardSubcells: x, y, z must be >= 1");
        }
    
        std::vector<std::reference_wrapper<const SubCell>> checkerboard;
        int n = subdivision_levels_;
    
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    // Apply parameterized checkerboard pattern
                    if ((i%x == 0 && k%z == 0) || j%y == 0) {
                        checkerboard.emplace_back(subcells_[i][j][k]);
                    }
                }
            }
        }
        return checkerboard;
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

    /**
     * @brief Generate triangular faces for parameterized checkerboard pattern.
     * 
     * Creates a FacetBox containing triangular faces from all subcells that match
     * the custom checkerboard pattern defined by (i%x == 0 && k%z == 0) || j%y == 0.
     * 
     * @param x Modulo divisor for i-coordinate (X-axis pattern frequency).
     * @param y Modulo divisor for j-coordinate (Y-axis pattern frequency).
     * @param z Modulo divisor for k-coordinate (Z-axis pattern frequency).
     * @return FacetBox containing triangular faces from matching active subcells.
     * @throws std::runtime_error if no subdivision available.
     * @throws std::invalid_argument if x, y, or z are less than 1.
     */
    FacetBox getCheckerboardFacets(int x, int y, int z) const {
        auto checkerboard_subcells = getCheckerboardSubcells(x, y, z);
        FacetBox checkerboard_facets;
    
        for (const SubCell& cell : checkerboard_subcells) {
            if (!cell.active) continue;
    
            for (int t = 0; t < 12; ++t) {
                auto const& tri = cube_triangles_[t];
                checkerboard_facets.push(cell.vertices[tri[0]], cell.vertices[tri[1]], cell.vertices[tri[2]]);
            }
        }
        return checkerboard_facets;
    }

    /**
     * @brief Flatten the entire subdivided lattice into a contiguous position array.
     *
     * Fills `positions` with n^3 * 8 * 3 floats, iterating subcells in physical
     * grid order [i][j][k] and local cube vertex 0..7. Each subcell contributes its
     * own 8 vertices (vertices are duplicated across neighbours, matching the
     * existing storage). This is the static vertex source for a GPU vertex buffer
     * uploaded once after subdivision is finalized.
     *
     * @param positions Output buffer (cleared and filled). Length = n*n*n*8*3.
     * @throws std::runtime_error if no subdivision available.
     */
    void fillVertexLattice(std::vector<float>& positions) const {
        if (!hasSubdivision()) {
            throw std::runtime_error("Cube::fillVertexLattice: no subdivision available");
        }
        const int n = subdivision_levels_;
        positions.clear();
        positions.reserve((size_t)n * n * n * 8 * 3);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    const SubCell& cell = subcells_[i][j][k];
                    for (int l = 0; l < 8; ++l) {
                        const Vector3D& v = cell.vertices[l];
                        positions.push_back((float)v.x());
                        positions.push_back((float)v.y());
                        positions.push_back((float)v.z());
                    }
                }
            }
        }
    }

    /**
     * @brief Build a triangle index buffer for the parameterized checkerboard pattern.
     *
     * Mirrors getCheckerboardFacets(x,y,z) but emits raw vertex indices into the
     * lattice produced by fillVertexLattice() instead of constructing Facet objects.
     * The predicate (i%x==0 && k%z==0) || j%y== 0, the active-cell filter, and the
     * cube_triangles_ lookup table are identical to the Facet-generating path, so the
     * rendered geometry is byte-for-byte the same — only the per-frame CPU cost
     * (cross product + sqrt per Facet) is eliminated. Uses unsigned int (no GL types)
     * so Cube stays GL-agnostic. Divisors < 1 are clamped to 1 to avoid modulo-by-zero.
     *
     * @param x Modulo divisor for i (X-axis), clamped to >= 1.
     * @param y Modulo divisor for j (Y-axis), clamped to >= 1.
     * @param z Modulo divisor for k (Z-axis), clamped to >= 1.
     * @param indices Output index buffer (cleared and filled). 12 triangles per active cell.
     * @throws std::runtime_error if no subdivision available.
     */
    void fillCheckerboardIndices(int x, int y, int z, std::vector<unsigned int>& indices) const {
        if (!hasSubdivision()) {
            throw std::runtime_error("Cube::fillCheckerboardIndices: no subdivision available");
        }
        if (x < 1) x = 1;
        if (y < 1) y = 1;
        if (z < 1) z = 1;

        const int n = subdivision_levels_;
        indices.clear();
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    if (!((i % x == 0 && k % z == 0) || j % y == 0)) continue;
                    const SubCell& cell = subcells_[i][j][k];
                    if (!cell.active) continue;
                    const unsigned int base = (unsigned int)(((i * n) + j) * n + k) * 8u;
                    for (int t = 0; t < 12; ++t) {
                        const auto& tri = cube_triangles_[t];
                        indices.push_back(base + (unsigned int)tri[0]);
                        indices.push_back(base + (unsigned int)tri[1]);
                        indices.push_back(base + (unsigned int)tri[2]);
                    }
                }
            }
        }
    }

    /**
     * @brief Build a flat octagonal (truncated-cube) vertex lattice for the checkerboard
     *        selection — the VBO counterpart of getCheckerboardFacetsOctagonal().
     *
     * Mirrors getCheckerboardFacetsOctagonal(x,y,z,hollow) but, instead of pushing
     * 152 Facet objects, expands each selected subcell into a fixed block of 126
     * generated vertices (the same outer/inner/centroid/corner points that
     * pushOctagonalCubeFacets computes) and appends them to a flat float buffer.
     * Use fillCheckerboardIndicesOctagonal() to get the matching triangle index
     * buffer; the two share the identical iteration order, predicate, active-cell
     * filter, and a per-cell 126-vertex layout so vertex-block base = selOrd*126
     * aligns with index base = selOrd*126.
     *
     * Per-cell vertex layout (126 verts, 3 floats each):
     *   faces f in 0..5, block at offset f*17 (17 verts/face):
     *     outer[0..7] at f*17+0..7, inner[0..7] at f*17+8..15, centroid at f*17+16
     *   corners c in 0..7, block at offset 102 + c*3: cut-points 0,1,2
     * (hollow does NOT change the vertex set — it only drops the center-fan
     * triangles in the index buffer — so 126 verts/cell regardless of hollow.)
     *
     * @param x Modulo divisor for i (clamped to >= 1).
     * @param y Modulo divisor for j (clamped to >= 1).
     * @param z Modulo divisor for k (clamped to >= 1).
     * @param positions Output flat vertex buffer (cleared and filled), 126*3 floats
     *                  per selected active cell. Pass through to glBufferSubData().
     * @param hollow Currently unused for the vertex set (kept for API symmetry with
     *               fillCheckerboardIndicesOctagonal and getCheckerboardFacetsOctagonal).
     * @throws std::runtime_error if no subdivision available.
     */
    void fillOctagonalVertexLattice(int x, int y, int z,
                                    std::vector<float>& positions, bool hollow) const {
        if (!hasSubdivision()) {
            throw std::runtime_error("Cube::fillOctagonalVertexLattice: no subdivision available");
        }
        if (x < 1) x = 1;
        if (y < 1) y = 1;
        if (z < 1) z = 1;

        const int n = subdivision_levels_;
        const double trunc      = default_trunc_;
        const double innerScale = default_inner_scale_;
        positions.clear();
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    if (!((i % x == 0 && k % z == 0) || j % y == 0)) continue;
                    const SubCell& cell = subcells_[i][j][k];
                    if (!cell.active) continue;
                    const std::array<Vector3D, 8>& v = cell.vertices;

                    // 6 faces: outer[8], inner[8], centroid[1] = 17 verts/face.
                    for (int f = 0; f < 6; ++f) {
                        const int* F = octagonal_faces_[f];
                        Vector3D P[4] = { v[F[0]], v[F[1]], v[F[2]], v[F[3]] };
                        Vector3D C = (P[0] + P[1] + P[2] + P[3]) / 4.0;   // face centroid
                        Vector3D outer[8];
                        for (int a = 0; a < 4; ++a) {
                            int nxt = (a + 1) & 3;
                            outer[2*a]     = P[a]   + trunc * (P[nxt] - P[a]);    // F[a]
                            outer[2*a + 1] = P[nxt] + trunc * (P[a]   - P[nxt]); // B[nxt]
                        }
                        Vector3D inner[8];
                        for (int a = 0; a < 8; ++a)
                            inner[a] = C + innerScale * (outer[a] - C);

                        for (int a = 0; a < 8; ++a) {
                            positions.push_back((float)outer[a].x());
                            positions.push_back((float)outer[a].y());
                            positions.push_back((float)outer[a].z());
                        }
                        for (int a = 0; a < 8; ++a) {
                            positions.push_back((float)inner[a].x());
                            positions.push_back((float)inner[a].y());
                            positions.push_back((float)inner[a].z());
                        }
                        positions.push_back((float)C.x());
                        positions.push_back((float)C.y());
                        positions.push_back((float)C.z());
                    }
                    // 8 corners: 3 cut-points each (order = corner_neighbors_[c]).
                    for (int c = 0; c < 8; ++c) {
                        const int* nn = corner_neighbors_[c];
                        Vector3D Vc = v[c];
                        Vector3D cp[3] = {
                            Vc + trunc * (v[nn[0]] - Vc),
                            Vc + trunc * (v[nn[1]] - Vc),
                            Vc + trunc * (v[nn[2]] - Vc)
                        };
                        for (int p = 0; p < 3; ++p) {
                            positions.push_back((float)cp[p].x());
                            positions.push_back((float)cp[p].y());
                            positions.push_back((float)cp[p].z());
                        }
                    }
                }
            }
        }
        (void)hollow;  // vertex set is independent of hollow (only triangle set changes)
    }

    /**
     * @brief Build the triangle index buffer for the octagonal checkerboard lattice.
     *
     * Emits 152 (solid) or 104 (hollow) triangle indices per selected active cell,
     * referencing the 126-vertex-per-cell layout produced by
     * fillOctagonalVertexLattice(). Iteration order, predicate
     * (i%x==0 && k%z==0) || j%y==0, and the active-cell filter are identical to
     * fillOctagonalVertexLattice, and a running selected-cell ordinal sets
     * base = selOrd*126 so the indices line up with the vertex buffer. Winding
     * matches pushOctagonalCubeFacets exactly (same fb.push() order), so the
     * rendered geometry is byte-for-byte the same as getCheckerboardFacetsOctagonal
     * — only the per-frame Facet construction is eliminated.
     *
     * Independent of deformation (depends only on selection + hollow), so rebuild
     * only when ii/jj/kk or hollow changes.
     *
     * @param indices Output index buffer (cleared and filled). 456 (solid) or
     *                312 (hollow) unsigned ints per selected active cell.
     * @throws std::runtime_error if no subdivision available.
     */
    void fillCheckerboardIndicesOctagonal(int x, int y, int z,
                                           std::vector<unsigned int>& indices,
                                           bool hollow) const {
        if (!hasSubdivision()) {
            throw std::runtime_error("Cube::fillCheckerboardIndicesOctagonal: no subdivision available");
        }
        if (x < 1) x = 1;
        if (y < 1) y = 1;
        if (z < 1) z = 1;

        const int n = subdivision_levels_;
        const unsigned int VPV  = 126;   // verts per cell
        const unsigned int FACE = 17;    // verts per face block (8 outer + 8 inner + 1 centroid)
        indices.clear();
        unsigned int sel = 0;            // running selected-cell ordinal (matches fillOctagonalVertexLattice)
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    if (!((i % x == 0 && k % z == 0) || j % y == 0)) continue;
                    const SubCell& cell = subcells_[i][j][k];
                    if (!cell.active) continue;
                    const unsigned int base = sel * VPV;

                    // Per face: 16 annular + (8 center-fan unless hollow).
                    for (int f = 0; f < 6; ++f) {
                        unsigned int fb = base + (unsigned int)f * FACE;  // outer @ +0..7, inner @ +8..15, C @ +16
                        for (int a = 0; a < 8; ++a) {
                            int b = (a + 1) & 7;
                            unsigned int oa = fb + (unsigned int)a;
                            unsigned int ob = fb + (unsigned int)b;
                            unsigned int ia = fb + 8u + (unsigned int)a;
                            unsigned int ib = fb + 8u + (unsigned int)b;
                            // push(outer[a], outer[b], inner[b])
                            indices.push_back(oa); indices.push_back(ob); indices.push_back(ib);
                            // push(outer[a], inner[b], inner[a])
                            indices.push_back(oa); indices.push_back(ib); indices.push_back(ia);
                        }
                        if (!hollow) {
                            unsigned int cc = fb + 16u;
                            for (int a = 0; a < 8; ++a) {
                                int b = (a + 1) & 7;
                                unsigned int ia = fb + 8u + (unsigned int)a;
                                unsigned int ib = fb + 8u + (unsigned int)b;
                                // push(C, inner[a], inner[b])
                                indices.push_back(cc); indices.push_back(ia); indices.push_back(ib);
                            }
                        }
                    }
                    // 8 corner triangles: cut-points at base+102 + c*3 + 0..2.
                    for (int c = 0; c < 8; ++c) {
                        unsigned int cb = base + 102u + (unsigned int)c * 3u;
                        indices.push_back(cb + 0u);
                        indices.push_back(cb + 1u);
                        indices.push_back(cb + 2u);
                    }
                    ++sel;
                }
            }
        }
    }

    /**
     * @brief Build a flat truncated-octahedron vertex lattice for the checkerboard
     *        selection — the VBO counterpart of getCheckerboardFacetsTruncatedOctahedron().
     *
     * Mirrors getCheckerboardFacetsTruncatedOctahedron(x,y,z,s) but, instead of
     * pushing 44 Facet objects, expands each selected active subcell into a fixed
     * block of 24 generated vertices (the faceCenter + s*(edgeMid-faceCenter)
     * points) and appends them to a flat float buffer. Use
     * fillCheckerboardIndicesTruncatedOctahedron() to get the matching triangle
     * index buffer; the two share the identical iteration order, predicate,
     * active-cell filter, and a per-cell 24-vertex layout so vertex-block base =
     * selOrd*24 aligns with index base = selOrd*24.
     *
     * Per-cell vertex layout (24 verts, 3 floats each): face f in 0..5, edge e in
     * 0..3, at flat index f*4+e: faceCenter(f)+s*(edgeMid(f,e)-faceCenter(f)).
     * Positions depend on s (re-upload the VBO when s changes); the index buffer
     * is independent of s.
     *
     * @param x Modulo divisor for i (clamped to >= 1).
     * @param y Modulo divisor for j (clamped to >= 1).
     * @param z Modulo divisor for k (clamped to >= 1).
     * @param positions Output flat vertex buffer (cleared and filled), 24*3 floats
     *                  per selected active cell. Pass through to glBufferSubData().
     * @param s Morph parameter in (0,1] (clamped); default regular TO.
     * @param hollow If true, emit 96 verts/cell (24 outer + 24 inner-square + 48
     *               inner-hex) for the hollow frame; false emits 24 outer verts/cell.
     * @param inset Hollow border inset ratio in (0,1) (clamped); only used when hollow.
     * @throws std::runtime_error if no subdivision available.
     */
    void fillTruncatedOctahedronVertexLattice(int x, int y, int z,
                                              std::vector<float>& positions,
                                              double s = default_to_scale_,
                                              bool hollow = false,
                                              double inset = default_to_inset_) const {
        if (!hasSubdivision()) {
            throw std::runtime_error("Cube::fillTruncatedOctahedronVertexLattice: no subdivision available");
        }
        if (x < 1) x = 1;
        if (y < 1) y = 1;
        if (z < 1) z = 1;
        if (!(s > 0.0)) s = 0.0001;
        if (s > 1.0)    s = 1.0;
        if (hollow) {
            if (!(inset > 0.0)) inset = 0.0001;
            if (inset >= 1.0)  inset = 0.9999;
        }

        const int n = subdivision_levels_;
        positions.clear();
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    if (!((i % x == 0 && k % z == 0) || j % y == 0)) continue;
                    const SubCell& cell = subcells_[i][j][k];
                    if (!cell.active) continue;
                    const std::array<Vector3D, 8>& v = cell.vertices;

                    // 24 surface verts, flat index = face*4 + edge.
                    Vector3D vert[24];
                    for (int f = 0; f < 6; ++f) {
                        const int* F = octagonal_faces_[f];
                        Vector3D P[4] = { v[F[0]], v[F[1]], v[F[2]], v[F[3]] };
                        Vector3D C = (P[0] + P[1] + P[2] + P[3]) / 4.0;   // face center
                        for (int e = 0; e < 4; ++e) {
                            int nxt = (e + 1) & 3;
                            Vector3D em = (P[e] + P[nxt]) / 2.0;        // edge midpoint
                            vert[f * 4 + e] = C + s * (em - C);
                        }
                    }
                    // emit the 24 outer verts (solid and hollow both).
                    for (int idx = 0; idx < 24; ++idx) {
                        positions.push_back((float)vert[idx].x());
                        positions.push_back((float)vert[idx].y());
                        positions.push_back((float)vert[idx].z());
                    }
                    if (!hollow) continue;

                    // HOLLOW: append inner verts. Per-cell layout after the 24 outer:
                    //   24 inner-square verts at 24 + f*4 + e   (f in 0..5, e in 0..3)
                    //   48 inner-hex    verts at 48 + c*6 + i   (c in 0..7, i in 0..5)
                    // inner = faceCentroid + inset*(outer - faceCentroid), faceCentroid =
                    // centroid of the face's own outer verts. 96 verts/cell total.
                    for (int f = 0; f < 6; ++f) {
                        int b = f * 4;
                        Vector3D o[4] = { vert[b], vert[b + 1], vert[b + 2], vert[b + 3] };
                        Vector3D Cf = (o[0] + o[1] + o[2] + o[3]) / 4.0;
                        for (int e = 0; e < 4; ++e) {
                            Vector3D p = Cf + inset * (o[e] - Cf);
                            positions.push_back((float)p.x());
                            positions.push_back((float)p.y());
                            positions.push_back((float)p.z());
                        }
                    }
                    for (int c = 0; c < 8; ++c) {
                        const int* h = to_hex_faces_[c];
                        Vector3D o[6] = { vert[h[0]], vert[h[1]], vert[h[2]], vert[h[3]], vert[h[4]], vert[h[5]] };
                        Vector3D Ch = (o[0] + o[1] + o[2] + o[3] + o[4] + o[5]) / 6.0;
                        for (int i = 0; i < 6; ++i) {
                            Vector3D p = Ch + inset * (o[i] - Ch);
                            positions.push_back((float)p.x());
                            positions.push_back((float)p.y());
                            positions.push_back((float)p.z());
                        }
                    }
                }
            }
        }
    }

    /**
     * @brief Build the triangle index buffer for the truncated-octahedron checkerboard lattice.
     *
     * Emits 132 (solid) or 432 (hollow) triangle indices per selected active cell,
     * referencing the 24- (solid) or 96- (hollow) vertex-per-cell layout produced by
     * fillTruncatedOctahedronVertexLattice(). Solid: 44 triangles (6 squares*2 +
     * 8 hexagons*4). Hollow: 144 frame triangles (6 squares*8 + 8 hexagons*12), inner
     * face skipped. Iteration order, predicate (i%x==0 && k%z==0) || j%y==0, and the
     * active-cell filter are identical to fillTruncatedOctahedronVertexLattice, and a
     * running selected-cell ordinal sets base = selOrd*VPV (VPV=24 solid, 96 hollow)
     * so the indices line up with the vertex buffer. Winding matches
     * pushTruncatedOctahedronFacets exactly (same fb.push() order), so the rendered
     * geometry is byte-for-byte the same as the CPU path — only the per-frame Facet
     * construction (cross product + sqrt per Facet) is eliminated.
     *
     * Independent of s and of the inset ratio (depends only on selection and hollow),
     * so rebuild only when ii/jj/kk or hollow change.
     *
     * @param indices Output index buffer (cleared and filled). 132 (solid) or 432
     *                (hollow) unsigned ints per selected active cell.
     * @param hollow If true, emit the 432 hollow-frame indices (VPV=96); else 132 solid.
     * @param inset Unused (indices are inset-independent); kept for API symmetry.
     * @throws std::runtime_error if no subdivision available.
     */
    void fillCheckerboardIndicesTruncatedOctahedron(int x, int y, int z,
                                                    std::vector<unsigned int>& indices,
                                                    bool hollow = false,
                                                    double inset = default_to_inset_) const {
        if (!hasSubdivision()) {
            throw std::runtime_error("Cube::fillCheckerboardIndicesTruncatedOctahedron: no subdivision available");
        }
        if (x < 1) x = 1;
        if (y < 1) y = 1;
        if (z < 1) z = 1;
        (void)inset;   // indices depend only on hollow (topology), not the inset ratio

        const int n = subdivision_levels_;
        const unsigned int VPV = hollow ? 96u : 24u;   // verts per cell (matches fillTruncatedOctahedronVertexLattice)
        indices.clear();
        unsigned int sel = 0;          // running selected-cell ordinal (matches fillTruncatedOctahedronVertexLattice)
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    if (!((i % x == 0 && k % z == 0) || j % y == 0)) continue;
                    const SubCell& cell = subcells_[i][j][k];
                    if (!cell.active) continue;
                    const unsigned int base = sel * VPV;

                    if (!hollow) {
                        // SOLID: 6 squares (2 tris) + 8 hexagons (4 fan tris) = 132 indices.
                        for (int f = 0; f < 6; ++f) {
                            unsigned int b = base + (unsigned int)f * 4u;
                            indices.push_back(b + 0u); indices.push_back(b + 1u); indices.push_back(b + 2u);
                            indices.push_back(b + 0u); indices.push_back(b + 2u); indices.push_back(b + 3u);
                        }
                        for (int c = 0; c < 8; ++c) {
                            const int* h = to_hex_faces_[c];
                            indices.push_back(base + (unsigned int)h[0]);
                            indices.push_back(base + (unsigned int)h[1]);
                            indices.push_back(base + (unsigned int)h[2]);
                            indices.push_back(base + (unsigned int)h[0]);
                            indices.push_back(base + (unsigned int)h[2]);
                            indices.push_back(base + (unsigned int)h[3]);
                            indices.push_back(base + (unsigned int)h[0]);
                            indices.push_back(base + (unsigned int)h[3]);
                            indices.push_back(base + (unsigned int)h[4]);
                            indices.push_back(base + (unsigned int)h[0]);
                            indices.push_back(base + (unsigned int)h[4]);
                            indices.push_back(base + (unsigned int)h[5]);
                        }
                    } else {
                        // HOLLOW: 6 squares (8 tris) + 8 hexagons (12 tris) = 432 indices.
                        // Layout: outer at base+f*4+e (squares) / base+to_hex_faces_[c][i] (hex);
                        // inner-square at base+24+f*4+e; inner-hex at base+48+c*6+i.
                        // Winding (o0,o1,i1),(o0,i1,i0) matches pushTruncatedOctahedronFacets.
                        for (int f = 0; f < 6; ++f) {
                            unsigned int of = base + (unsigned int)f * 4u;
                            unsigned int iff = base + 24u + (unsigned int)f * 4u;
                            for (int e = 0; e < 4; ++e) {
                                unsigned int o0 = of + (unsigned int)e;
                                unsigned int o1 = of + (unsigned int)((e + 1) & 3);
                                unsigned int i0 = iff + (unsigned int)e;
                                unsigned int i1 = iff + (unsigned int)((e + 1) & 3);
                                indices.push_back(o0); indices.push_back(o1); indices.push_back(i1);
                                indices.push_back(o0); indices.push_back(i1); indices.push_back(i0);
                            }
                        }
                        for (int c = 0; c < 8; ++c) {
                            const int* h = to_hex_faces_[c];
                            unsigned int ih = base + 48u + (unsigned int)c * 6u;
                            for (int i = 0; i < 6; ++i) {
                                unsigned int o0 = base + (unsigned int)h[i];
                                unsigned int o1 = base + (unsigned int)h[(i + 1) % 6];
                                unsigned int i0 = ih + (unsigned int)i;
                                unsigned int i1 = ih + (unsigned int)((i + 1) % 6);
                                indices.push_back(o0); indices.push_back(o1); indices.push_back(i1);
                                indices.push_back(o0); indices.push_back(i1); indices.push_back(i0);
                            }
                        }
                    }
                    ++sel;
                }
            }
        }
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

    /* === OCTAGONAL (TRUNCATED-CUBE) MESH — a second method for building a cube ===
     * Each face becomes an octagon (outer ring + inner ring + center fan = 24
     * triangles) and 8 corner triangles close the solid => 152 triangles/cube.
     * Derived on the fly from the same 8 corner vertices the plain 12-triangle
     * path uses, so the slow-inversion animation is identical. The truncation
     * (default_trunc_) and inner scale (default_inner_scale_) are fixed; edit
     * those constexpr to tune the look. FacetBox results are returned by value.
     */

    /**
     * @brief Rebuild this cube's stored facets as an octagonal (truncated-cube) mesh.
     *
     * Second method for building a cube's surface, parallel to buildFacets():
     * populates facets_ from the 8 main-cube vertices (verts_) using the
     * truncated-cube octagonal mesh (152 triangles). For a basic (non-subdivided)
     * cube only; for subdivided cubes prefer the on-the-fly getters below, since
     * storing 152*n^3 facets is expensive.
     */
    void buildFacetsOctagonal(bool hollow = false) {
        facets_.clear();
        std::array<Vector3D, 8> v;
        for (int i = 0; i < 8; ++i) v[i] = verts_[i].V();   // Quaternion -> Vector3D
        pushOctagonalCubeFacets(facets_, v, default_trunc_, default_inner_scale_, hollow);
    }

    /**
     * @brief Build a fresh FacetBox of this cube as an octagonal (truncated-cube) mesh.
     *
     * For a non-subdivided cube, builds from the 8 main vertices (verts_). For a
     * subdivided cube, iterates all active subcells. Returns by value (freshly
     * built), matching getCheckerboardFacets' contract (not the const-ref of getFacets).
     */
    FacetBox getFacetsOctagonal(bool hollow = false) const {
        FacetBox fb;
        if (!hasSubdivision()) {
            std::array<Vector3D, 8> v;
            for (int i = 0; i < 8; ++i) v[i] = verts_[i].V();
            pushOctagonalCubeFacets(fb, v, default_trunc_, default_inner_scale_, hollow);
        } else {
            const int n = subdivision_levels_;
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    for (int k = 0; k < n; ++k) {
                        const SubCell& cell = subcells_[i][j][k];
                        if (!cell.active) continue;
                        pushOctagonalCubeFacets(fb, cell.vertices, default_trunc_, default_inner_scale_, hollow);
                    }
        }
        return fb;
    }

    /**
     * @brief Checkerboard subcells as an octagonal (truncated-cube) mesh.
     *
     * Same subcell selection as getCheckerboardFacets(x,y,z) (unchanged predicate),
     * but each selected active subcell is triangulated as a truncated cube via
     * pushOctagonalCubeFacets instead of the 12-triangle cube_triangles_ table.
     * This is the render path used by main15_slow_inversion_truncated_cube.
     */
    FacetBox getCheckerboardFacetsOctagonal(int x, int y, int z, bool hollow = false) const {
        auto cells = getCheckerboardSubcells(x, y, z);
        FacetBox fb;
        for (const SubCell& cell : cells) {
            if (!cell.active) continue;
            pushOctagonalCubeFacets(fb, cell.vertices, default_trunc_, default_inner_scale_, hollow);
        }
        return fb;
    }

    /**
     * @brief A 2D plane of subcells as an octagonal (truncated-cube) mesh.
     *
     * Same selection as getPlane(axis, layer); each active subcell is triangulated
     * as a truncated cube.
     */
    FacetBox getPlaneFacetsOctagonal(int axis, int layer, bool hollow = false) const {
        auto cells = getPlane(axis, layer);
        FacetBox fb;
        for (const SubCell& cell : cells) {
            if (!cell.active) continue;
            pushOctagonalCubeFacets(fb, cell.vertices, default_trunc_, default_inner_scale_, hollow);
        }
        return fb;
    }

    /**
     * @brief One subcell as an octagonal (truncated-cube) mesh (152 triangles).
     */
    FacetBox getSubCellFacetsOctagonal(int x, int y, int z, bool hollow = false) const {
        const SubCell& cell = getSubCell(x, y, z);
        FacetBox fb;
        if (!cell.active) return fb;
        pushOctagonalCubeFacets(fb, cell.vertices, default_trunc_, default_inner_scale_, hollow);
        return fb;
    }

    /**
     * @brief Rebuild stored facets from subcells as an octagonal mesh.
     *
     * WARNING: stores 152*n^3 Facets (~1 GB at n=37). Not used by the slow-inversion
     * live path (which reads subcells on the fly via getCheckerboardFacetsOctagonal);
     * provided only for API parity with refreshTriangulation.
     */
    void refreshTriangulationOctagonal(bool hollow = false) {
        if (!hasSubdivision()) return;
        facets_.clear();
        const int n = subdivision_levels_;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                for (int k = 0; k < n; ++k) {
                    const SubCell& cell = subcells_[i][j][k];
                    if (!cell.active) continue;
                    pushOctagonalCubeFacets(facets_, cell.vertices, default_trunc_, default_inner_scale_, hollow);
                }
    }

    /**
     * @brief Write this cube's octagonal (truncated-cube) mesh to an STL file.
     *
     * Parallel to writeSTL_s (6-arg): same modes and ASCII format. "checkerboard"
     * -> getCheckerboardFacetsOctagonal(x,y,z) (or getFacetsOctagonal if not
     * subdivided); "full" -> getFacetsOctagonal; "plane_xy/xz/yz" ->
     * getPlaneFacetsOctagonal. Only active subcells are included.
     *
     * @param filename Path to output STL file (created or overwritten).
     * @param objectName Name to use in the STL header.
     * @param mode Extraction mode: "full", "checkerboard", "plane_xy", "plane_xz", "plane_yz".
     * @param x Checkerboard X modulo, or plane layer for plane modes.
     * @param y Checkerboard Y modulo (default 9).
     * @param z Checkerboard Z modulo (default 2).
     * @param hollow If true, omit the 8 center-fan triangles per face so each face
     *               has an octagonal hole (hollow truncated cubes); default false (solid).
     * @throws std::runtime_error if the file cannot be created.
     * @throws std::invalid_argument if mode is invalid or a plane mode is used without subdivision.
     */
    void writeSTL_s_octagonal(const std::string& filename,
                              const std::string& objectName,
                              const std::string& mode,
                              int x,
                              int y = 9,
                              int z = 2,
                              bool hollow = false) const {
        std::ofstream stl(filename);
        if (!stl.is_open())
            throw std::runtime_error("Cube::writeSTL_s_octagonal: Cannot create file " + filename);

        stl << "solid " << objectName << "\n";

        FacetBox facetsToWrite;
        if (mode == "full") {
            facetsToWrite = getFacetsOctagonal(hollow);
        } else if (mode == "checkerboard") {
            if (!hasSubdivision()) facetsToWrite = getFacetsOctagonal(hollow);
            else                  facetsToWrite = getCheckerboardFacetsOctagonal(x, y, z, hollow);
        } else if (mode == "plane_xy") {
            if (!hasSubdivision()) throw std::invalid_argument("Plane extraction requires subdivided cube");
            facetsToWrite = getPlaneFacetsOctagonal(2, x, hollow);
        } else if (mode == "plane_xz") {
            if (!hasSubdivision()) throw std::invalid_argument("Plane extraction requires subdivided cube");
            facetsToWrite = getPlaneFacetsOctagonal(1, x, hollow);
        } else if (mode == "plane_yz") {
            if (!hasSubdivision()) throw std::invalid_argument("Plane extraction requires subdivided cube");
            facetsToWrite = getPlaneFacetsOctagonal(0, x, hollow);
        } else {
            throw std::invalid_argument("Invalid mode: " + mode +
                ". Valid modes: full, checkerboard, plane_xy, plane_xz, plane_yz");
        }

        for (size_t i = 0; i < facetsToWrite.size(); ++i) {
            const Facet& face = facetsToWrite[i];
            Vector3D v1 = face[0], v2 = face[1], v3 = face[2];
            Vector3D normal = face.getNormal();
            stl << "facet normal " << std::scientific
                << normal.x() << " " << normal.y() << " " << normal.z() << "\n";
            stl << "\touter loop\n";
            stl << "\t\tvertex " << v1.x() << " " << v1.y() << " " << v1.z() << "\n";
            stl << "\t\tvertex " << v2.x() << " " << v2.y() << " " << v2.z() << "\n";
            stl << "\t\tvertex " << v3.x() << " " << v3.y() << " " << v3.z() << "\n";
            stl << "\tendloop\n";
            stl << "endfacet\n";
        }

        stl << "endsolid " << objectName << "\n";
        stl.close();
    }

    /* === TRUNCATED-OCTAHEDRON (8 HEXAGONS + 6 SQUARES) MESH — public API ===
     * Parallel to the octagonal (truncated-cube) family above, but each cube
     * corner becomes a hexagon and each cube face a square (44 tris/cube). Every
     * getter takes the morph parameter s (default default_to_scale_ = 0.5, the
     * regular truncated octahedron); s=1.0 is the cuboctahedron (cube-like end,
     * corner faces = triangles), s<0.5 stretches the corner hexagons larger.
     */

    /**
     * @brief Rebuild this cube's stored facets as a truncated-octahedron mesh.
     *
     * Third surface method, parallel to buildFacets()/buildFacetsOctagonal():
     * populates facets_ from the 8 main-cube vertices (verts_) using the
     * truncated-octahedron mesh (44 triangles). For a basic (non-subdivided)
     * cube only; for subdivided cubes prefer the on-the-fly getters below.
     *
     * @param s Morph parameter in (0,1] (default regular TO).
     */
    void buildFacetsTruncatedOctahedron(double s = default_to_scale_,
                                        bool hollow = false,
                                        double inset = default_to_inset_) {
        facets_.clear();
        std::array<Vector3D, 8> v;
        for (int i = 0; i < 8; ++i) v[i] = verts_[i].V();   // Quaternion -> Vector3D
        pushTruncatedOctahedronFacets(facets_, v, s, hollow, inset);
    }

    /**
     * @brief Build a fresh FacetBox of this cube as a truncated-octahedron mesh.
     *
     * For a non-subdivided cube, builds from the 8 main vertices (verts_). For a
     * subdivided cube, iterates all active subcells. Returns by value (freshly
     * built), matching getFacetsOctagonal's contract.
     *
     * @param s Morph parameter in (0,1] (default regular TO).
     */
    FacetBox getFacetsTruncatedOctahedron(double s = default_to_scale_,
                                          bool hollow = false,
                                          double inset = default_to_inset_) const {
        FacetBox fb;
        if (!hasSubdivision()) {
            std::array<Vector3D, 8> v;
            for (int i = 0; i < 8; ++i) v[i] = verts_[i].V();
            pushTruncatedOctahedronFacets(fb, v, s, hollow, inset);
        } else {
            const int n = subdivision_levels_;
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    for (int k = 0; k < n; ++k) {
                        const SubCell& cell = subcells_[i][j][k];
                        if (!cell.active) continue;
                        pushTruncatedOctahedronFacets(fb, cell.vertices, s, hollow, inset);
                    }
        }
        return fb;
    }

    /**
     * @brief Checkerboard subcells as a truncated-octahedron mesh.
     *
     * Same subcell selection as getCheckerboardFacets(x,y,z) (unchanged
     * predicate), but each selected active subcell is triangulated as a
     * truncated octahedron via pushTruncatedOctahedronFacets.
     *
     * @param s Morph parameter in (0,1] (default regular TO).
     */
    FacetBox getCheckerboardFacetsTruncatedOctahedron(int x, int y, int z,
                                                      double s = default_to_scale_,
                                                      bool hollow = false,
                                                      double inset = default_to_inset_) const {
        auto cells = getCheckerboardSubcells(x, y, z);
        FacetBox fb;
        for (const SubCell& cell : cells) {
            if (!cell.active) continue;
            pushTruncatedOctahedronFacets(fb, cell.vertices, s, hollow, inset);
        }
        return fb;
    }

    /**
     * @brief A 2D plane of subcells as a truncated-octahedron mesh.
     *
     * Same selection as getPlane(axis, layer); each active subcell is triangulated
     * as a truncated octahedron.
     *
     * @param s Morph parameter in (0,1] (default regular TO).
     */
    FacetBox getPlaneFacetsTruncatedOctahedron(int axis, int layer,
                                               double s = default_to_scale_,
                                               bool hollow = false,
                                               double inset = default_to_inset_) const {
        auto cells = getPlane(axis, layer);
        FacetBox fb;
        for (const SubCell& cell : cells) {
            if (!cell.active) continue;
            pushTruncatedOctahedronFacets(fb, cell.vertices, s, hollow, inset);
        }
        return fb;
    }

    /**
     * @brief One subcell as a truncated-octahedron mesh (44 triangles).
     *
     * @param s Morph parameter in (0,1] (default regular TO).
     */
    FacetBox getSubCellFacetsTruncatedOctahedron(int x, int y, int z,
                                                 double s = default_to_scale_,
                                                 bool hollow = false,
                                                 double inset = default_to_inset_) const {
        const SubCell& cell = getSubCell(x, y, z);
        FacetBox fb;
        if (!cell.active) return fb;
        pushTruncatedOctahedronFacets(fb, cell.vertices, s, hollow, inset);
        return fb;
    }

    /**
     * @brief Rebuild stored facets from subcells as a truncated-octahedron mesh.
     *
     * WARNING: stores 44*n^3 Facets (large at high n). Not used by single-shape
     * demos or the on-the-fly render path (which reads subcells via
     * getCheckerboardFacetsTruncatedOctahedron); provided only for API parity with
     * refreshTriangulation / refreshTriangulationOctagonal.
     *
     * @param s Morph parameter in (0,1] (default regular TO).
     */
    void refreshTriangulationTruncatedOctahedron(double s = default_to_scale_,
                                                 bool hollow = false,
                                                 double inset = default_to_inset_) {
        if (!hasSubdivision()) return;
        facets_.clear();
        const int n = subdivision_levels_;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                for (int k = 0; k < n; ++k) {
                    const SubCell& cell = subcells_[i][j][k];
                    if (!cell.active) continue;
                    pushTruncatedOctahedronFacets(facets_, cell.vertices, s, hollow, inset);
                }
    }

    /**
     * @brief Write this cube's truncated-octahedron mesh to an STL file.
     *
     * Parallel to writeSTL_s_octagonal: same modes and ASCII format. "checkerboard"
     * -> getCheckerboardFacetsTruncatedOctahedron(x,y,z,s) (or getFacetsTruncatedOctahedron
     * if not subdivided); "full" -> getFacetsTruncatedOctahedron; "plane_xy/xz/yz"
     * -> getPlaneFacetsTruncatedOctahedron. Only active subcells are included.
     *
     * @param filename Path to output STL file (created or overwritten).
     * @param objectName Name to use in the STL header.
     * @param mode Extraction mode: "full", "checkerboard", "plane_xy", "plane_xz", "plane_yz".
     * @param x Checkerboard X modulo, or plane layer for plane modes.
     * @param y Checkerboard Y modulo (default 9).
     * @param z Checkerboard Z modulo (default 2).
     * @param s Morph parameter in (0,1] (default regular TO).
     * @param hollow If true, write the hollow frame mesh (144 tris/cube); false=solid (44).
     * @param inset Hollow border inset ratio in (0,1) (default 0.5); only used when hollow=true.
     * @throws std::runtime_error if the file cannot be created.
     * @throws std::invalid_argument if mode is invalid or a plane mode is used without subdivision.
     */
    void writeSTL_s_truncated_octahedron(const std::string& filename,
                                         const std::string& objectName,
                                         const std::string& mode,
                                         int x,
                                         int y = 9,
                                         int z = 2,
                                         double s = default_to_scale_,
                                         bool hollow = false,
                                         double inset = default_to_inset_) const {
        std::ofstream stl(filename);
        if (!stl.is_open())
            throw std::runtime_error("Cube::writeSTL_s_truncated_octahedron: Cannot create file " + filename);

        stl << "solid " << objectName << "\n";

        FacetBox facetsToWrite;
        if (mode == "full") {
            facetsToWrite = getFacetsTruncatedOctahedron(s, hollow, inset);
        } else if (mode == "checkerboard") {
            if (!hasSubdivision()) facetsToWrite = getFacetsTruncatedOctahedron(s, hollow, inset);
            else                   facetsToWrite = getCheckerboardFacetsTruncatedOctahedron(x, y, z, s, hollow, inset);
        } else if (mode == "plane_xy") {
            if (!hasSubdivision()) throw std::invalid_argument("Plane extraction requires subdivided cube");
            facetsToWrite = getPlaneFacetsTruncatedOctahedron(2, x, s, hollow, inset);
        } else if (mode == "plane_xz") {
            if (!hasSubdivision()) throw std::invalid_argument("Plane extraction requires subdivided cube");
            facetsToWrite = getPlaneFacetsTruncatedOctahedron(1, x, s, hollow, inset);
        } else if (mode == "plane_yz") {
            if (!hasSubdivision()) throw std::invalid_argument("Plane extraction requires subdivided cube");
            facetsToWrite = getPlaneFacetsTruncatedOctahedron(0, x, s, hollow, inset);
        } else {
            throw std::invalid_argument("Invalid mode: " + mode +
                ". Valid modes: full, checkerboard, plane_xy, plane_xz, plane_yz");
        }

        for (size_t i = 0; i < facetsToWrite.size(); ++i) {
            const Facet& face = facetsToWrite[i];
            Vector3D v1 = face[0], v2 = face[1], v3 = face[2];
            Vector3D normal = face.getNormal();
            stl << "facet normal " << std::scientific
                << normal.x() << " " << normal.y() << " " << normal.z() << "\n";
            stl << "\touter loop\n";
            stl << "\t\tvertex " << v1.x() << " " << v1.y() << " " << v1.z() << "\n";
            stl << "\t\tvertex " << v2.x() << " " << v2.y() << " " << v2.z() << "\n";
            stl << "\t\tvertex " << v3.x() << " " << v3.y() << " " << v3.z() << "\n";
            stl << "\tendloop\n";
            stl << "endfacet\n";
        }

        stl << "endsolid " << objectName << "\n";
        stl.close();
    }


    /* === TRUNCATED-CUBOCTAHEDRON (6 OCTAGONS + 8 HEXAGONS + 12 SQUARES) MESH — public API ===
     * Parallel to the truncated-octahedron family above, but each cube face becomes
     * an octagon, each cube corner a hexagon and each cube edge a square (92 tris/cube
     * solid, 288 hollow). No morph parameter: the truncated cuboctahedron is regular
     * only at the fixed a/b encoded in tc_vert_norm_ (see the private tables below);
     * the interactive knobs are hollow on/off and the inset ratio.
     */

    /**
     * @brief Rebuild this cube's stored facets as a truncated-cuboctahedron mesh.
     *
     * Fourth surface method, parallel to buildFacetsTruncatedOctahedron():
     * populates facets_ from the 8 main-cube vertices (verts_) using the
     * truncated-cuboctahedron mesh (92 triangles). For a basic (non-subdivided)
     * cube only; for subdivided cubes prefer the on-the-fly getters below.
     */
    void buildFacetsTruncatedCuboctahedron(bool hollow = false,
                                            double inset = default_tc_inset_) {
        facets_.clear();
        std::array<Vector3D, 8> v;
        for (int i = 0; i < 8; ++i) v[i] = verts_[i].V();   // Quaternion -> Vector3D
        pushTruncatedCuboctahedronFacets(facets_, v, hollow, inset);
    }

    /**
     * @brief Build a fresh FacetBox of this cube as a truncated-cuboctahedron mesh.
     *
     * For a non-subdivided cube, builds from the 8 main vertices (verts_). For a
     * subdivided cube, iterates all active subcells. Returns by value (freshly
     * built), matching getFacetsTruncatedOctahedron's contract.
     */
    FacetBox getFacetsTruncatedCuboctahedron(bool hollow = false,
                                             double inset = default_tc_inset_) const {
        FacetBox fb;
        if (!hasSubdivision()) {
            std::array<Vector3D, 8> v;
            for (int i = 0; i < 8; ++i) v[i] = verts_[i].V();
            pushTruncatedCuboctahedronFacets(fb, v, hollow, inset);
        } else {
            const int n = subdivision_levels_;
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    for (int k = 0; k < n; ++k) {
                        const SubCell& cell = subcells_[i][j][k];
                        if (!cell.active) continue;
                        pushTruncatedCuboctahedronFacets(fb, cell.vertices, hollow, inset);
                    }
        }
        return fb;
    }

    /**
     * @brief Checkerboard subcells as a truncated-cuboctahedron mesh.
     *
     * Same subcell selection as getCheckerboardFacets(x,y,z) (unchanged
     * predicate), but each selected active subcell is triangulated as a
     * truncated cuboctahedron via pushTruncatedCuboctahedronFacets.
     */
    FacetBox getCheckerboardFacetsTruncatedCuboctahedron(int x, int y, int z,
                                                         bool hollow = false,
                                                         double inset = default_tc_inset_) const {
        auto cells = getCheckerboardSubcells(x, y, z);
        FacetBox fb;
        for (const SubCell& cell : cells) {
            if (!cell.active) continue;
            pushTruncatedCuboctahedronFacets(fb, cell.vertices, hollow, inset);
        }
        return fb;
    }

    /**
     * @brief A 2D plane of subcells as a truncated-cuboctahedron mesh.
     *
     * Same selection as getPlane(axis, layer); each active subcell is triangulated
     * as a truncated cuboctahedron.
     */
    FacetBox getPlaneFacetsTruncatedCuboctahedron(int axis, int layer,
                                                  bool hollow = false,
                                                  double inset = default_tc_inset_) const {
        auto cells = getPlane(axis, layer);
        FacetBox fb;
        for (const SubCell& cell : cells) {
            if (!cell.active) continue;
            pushTruncatedCuboctahedronFacets(fb, cell.vertices, hollow, inset);
        }
        return fb;
    }

    /**
     * @brief One subcell as a truncated-cuboctahedron mesh (92 triangles).
     */
    FacetBox getSubCellFacetsTruncatedCuboctahedron(int x, int y, int z,
                                                    bool hollow = false,
                                                    double inset = default_tc_inset_) const {
        const SubCell& cell = getSubCell(x, y, z);
        FacetBox fb;
        if (!cell.active) return fb;
        pushTruncatedCuboctahedronFacets(fb, cell.vertices, hollow, inset);
        return fb;
    }

    /**
     * @brief Rebuild stored facets from subcells as a truncated-cuboctahedron mesh.
     *
     * WARNING: stores 92*n^3 Facets (large at high n). Provided only for API parity
     * with refreshTriangulation / refreshTriangulationOctagonal /
     * refreshTriangulationTruncatedOctahedron; not used by the on-the-fly demos.
     */
    void refreshTriangulationTruncatedCuboctahedron(bool hollow = false,
                                                     double inset = default_tc_inset_) {
        if (!hasSubdivision()) return;
        facets_.clear();
        const int n = subdivision_levels_;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                for (int k = 0; k < n; ++k) {
                    const SubCell& cell = subcells_[i][j][k];
                    if (!cell.active) continue;
                    pushTruncatedCuboctahedronFacets(facets_, cell.vertices, hollow, inset);
                }
    }

    /**
     * @brief Write this cube's truncated-cuboctahedron mesh to an STL file.
     *
     * Parallel to writeSTL_s_truncated_octahedron: same modes and ASCII format.
     * "checkerboard" -> getCheckerboardFacetsTruncatedCuboctahedron(x,y,z,hollow,inset)
     * (or getFacetsTruncatedCuboctahedron if not subdivided); "full" ->
     * getFacetsTruncatedCuboctahedron; "plane_xy/xz/yz" ->
     * getPlaneFacetsTruncatedCuboctahedron. Only active subcells are included.
     *
     * @param filename Path to output STL file (created or overwritten).
     * @param objectName Name to use in the STL header.
     * @param mode Extraction mode: "full", "checkerboard", "plane_xy", "plane_xz", "plane_yz".
     * @param x Checkerboard X modulo, or plane layer for plane modes.
     * @param y Checkerboard Y modulo (default 9).
     * @param z Checkerboard Z modulo (default 2).
     * @param hollow If true, write the hollow frame mesh (288 tris/cube); false=solid (92).
     * @param inset Hollow border inset ratio in (0,1) (default 0.5); only used when hollow=true.
     * @throws std::runtime_error if the file cannot be created.
     * @throws std::invalid_argument if mode is invalid or a plane mode is used without subdivision.
     */
    void writeSTL_s_truncated_cuboctahedron(const std::string& filename,
                                            const std::string& objectName,
                                            const std::string& mode,
                                            int x,
                                            int y = 9,
                                            int z = 2,
                                            bool hollow = false,
                                            double inset = default_tc_inset_) const {
        std::ofstream stl(filename);
        if (!stl.is_open())
            throw std::runtime_error("Cube::writeSTL_s_truncated_cuboctahedron: Cannot create file " + filename);

        stl << "solid " << objectName << "\n";

        FacetBox facetsToWrite;
        if (mode == "full") {
            facetsToWrite = getFacetsTruncatedCuboctahedron(hollow, inset);
        } else if (mode == "checkerboard") {
            if (!hasSubdivision()) facetsToWrite = getFacetsTruncatedCuboctahedron(hollow, inset);
            else                   facetsToWrite = getCheckerboardFacetsTruncatedCuboctahedron(x, y, z, hollow, inset);
        } else if (mode == "plane_xy") {
            if (!hasSubdivision()) throw std::invalid_argument("Plane extraction requires subdivided cube");
            facetsToWrite = getPlaneFacetsTruncatedCuboctahedron(2, x, hollow, inset);
        } else if (mode == "plane_xz") {
            if (!hasSubdivision()) throw std::invalid_argument("Plane extraction requires subdivided cube");
            facetsToWrite = getPlaneFacetsTruncatedCuboctahedron(1, x, hollow, inset);
        } else if (mode == "plane_yz") {
            if (!hasSubdivision()) throw std::invalid_argument("Plane extraction requires subdivided cube");
            facetsToWrite = getPlaneFacetsTruncatedCuboctahedron(0, x, hollow, inset);
        } else {
            throw std::invalid_argument("Invalid mode: " + mode +
                ". Valid modes: full, checkerboard, plane_xy, plane_xz, plane_yz");
        }

        for (size_t i = 0; i < facetsToWrite.size(); ++i) {
            const Facet& face = facetsToWrite[i];
            Vector3D v1 = face[0], v2 = face[1], v3 = face[2];
            Vector3D normal = face.getNormal();
            stl << "facet normal " << std::scientific
                << normal.x() << " " << normal.y() << " " << normal.z() << "\n";
            stl << "\touter loop\n";
            stl << "\t\tvertex " << v1.x() << " " << v1.y() << " " << v1.z() << "\n";
            stl << "\t\tvertex " << v2.x() << " " << v2.y() << " " << v2.z() << "\n";
            stl << "\t\tvertex " << v3.x() << " " << v3.y() << " " << v3.z() << "\n";
            stl << "\tendloop\n";
            stl << "endfacet\n";
        }

        stl << "endsolid " << objectName << "\n";
        stl.close();
    }


    /* === CANTITRUNCATED CUBIC HONEYCOMB — mixed-cell builder ===
     * The cantitruncated cubic honeycomb (uniform space-filling tessellation,
     * vertex figure mirrored sphenoid) is made of truncated cuboctahedra, truncated
     * octahedra and cubes in the ratio 1 : 1 : 3.
     *
     * Arrangement (from the omnitruncation t_{0,1,2}{4,3,4}): in the half-integer
     * grid, a point is a cell center iff it has 0, 1 or 3 half-integer coordinates
     * (2 = an edge center, no cell). Mapping that onto this Cube's simple-cubic
     * subcell grid by index parity, in each 2x2x2 block `odd = #(i,j,k odd)` gives
     *   odd == 0 -> 1 truncated cuboctahedron  (orig. cube centers)
     *   odd == 1 -> 3 cubes                     (orig. face centers, 3 orientations)
     *   odd == 3 -> 1 truncated octahedron      (orig. vertices)
     *   odd == 2 -> 3 gaps (skip)               (orig. edge centers, no cell)
     * which is exactly 1 : 1 : 3 per fundamental block, tiling periodically.
     *
     * This is the framework's stylized per-subcell version (each cell is built at
     * subcell scale from that subcell's 8 corners, like the existing truncated-cube
     * and truncated-octahedron honeycombs); the odd==2 gaps are the analog of the
     * checkerboard gaps those demos leave. Because every cell is rebuilt from its
     * (possibly deformed) subcell corners, the existing sphere-inversion machinery
     * (sigmaCenter / composeSigmaBlended deforming subcell vertices before this
     * call) deforms all three cell types together with no special handling.
     */

    /**
     * @brief Build the cantitruncated cubic honeycomb (1:1:3) over all active
     *        subcells and return it as one merged FacetBox.
     *
     * @param hollow If true, the truncated-cuboctahedron and truncated-octahedron
     *               cells are emitted as hollow frames (288 / 144 tris); cubes stay
     *               solid 12 tris (a cube has no face to inset meaningfully here).
     * @param inset  Hollow border inset ratio for the two truncated cell types.
     * @param toS    Morph parameter for the truncated-octahedron cells (default
     *               regular truncated octahedron, 0.5).
     * @throws std::runtime_error if this cube has no subdivision.
     */
    FacetBox getCantitruncatedHoneycombFacets(bool hollow = false,
                                              double inset = default_tc_inset_,
                                              double toS = default_to_scale_) const {
        FacetBox fb;
        if (!hasSubdivision())
            throw std::runtime_error("Cube::getCantitruncatedHoneycombFacets: no subdivision available");
        const int n = subdivision_levels_;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                for (int k = 0; k < n; ++k) {
                    const SubCell& cell = subcells_[i][j][k];
                    if (!cell.active) continue;
                    const int odd = (i & 1) + (j & 1) + (k & 1);
                    if      (odd == 0) pushTruncatedCuboctahedronFacets(fb, cell.vertices, hollow, inset);
                    else if (odd == 3) pushTruncatedOctahedronFacets(fb, cell.vertices, toS, hollow, inset);
                    else if (odd == 1) pushCubeFacets(fb, cell.vertices);
                    // odd == 2: edge-center gap -> skip (analog of checkerboard gaps)
                }
        return fb;
    }

    /**
     * @brief Write the cantitruncated cubic honeycomb to an ASCII STL file.
     *
     * @param filename Path to output STL file (created or overwritten).
     * @param objectName Name to use in the STL header.
     * @param hollow If true, the truncated cells are hollow frames; cubes stay solid.
     * @param inset  Hollow border inset ratio (only used when hollow=true).
     * @param toS    Morph parameter for the truncated-octahedron cells (default 0.5).
     * @throws std::runtime_error if the file cannot be created or there is no subdivision.
     */
    void writeSTL_s_cantitruncated_honeycomb(const std::string& filename,
                                             const std::string& objectName,
                                             bool hollow = false,
                                             double inset = default_tc_inset_,
                                             double toS = default_to_scale_) const {
        std::ofstream stl(filename);
        if (!stl.is_open())
            throw std::runtime_error("Cube::writeSTL_s_cantitruncated_honeycomb: Cannot create file " + filename);

        FacetBox facetsToWrite = getCantitruncatedHoneycombFacets(hollow, inset, toS);

        stl << "solid " << objectName << "\n";
        for (size_t i = 0; i < facetsToWrite.size(); ++i) {
            const Facet& face = facetsToWrite[i];
            Vector3D v1 = face[0], v2 = face[1], v3 = face[2];
            Vector3D normal = face.getNormal();
            stl << "facet normal " << std::scientific
                << normal.x() << " " << normal.y() << " " << normal.z() << "\n";
            stl << "\touter loop\n";
            stl << "\t\tvertex " << v1.x() << " " << v1.y() << " " << v1.z() << "\n";
            stl << "\t\tvertex " << v2.x() << " " << v2.y() << " " << v2.z() << "\n";
            stl << "\t\tvertex " << v3.x() << " " << v3.y() << " " << v3.z() << "\n";
            stl << "\tendloop\n";
            stl << "endfacet\n";
        }
        stl << "endsolid " << objectName << "\n";
        stl.close();
    }


    /* === OCTAGON-CONNECTED CUBOCTAHEDRON LATTICE ===
     * A truncated cuboctahedron at EVERY active subcell. Each cell's 6 octagons
     * sit on the 6 subcell faces, so two face-adjacent cuboctahedra share a
     * coincident octagon on the shared face — the octagons become the "windows"
     * that connect the lattice. The hexagon + square frames bound the octahedral
     * voids (at lattice vertices) and cubic voids (at lattice edges) left open
     * between the cells. This is the cantitruncated honeycomb skeleton with the
     * filler cells (truncated octahedra + cubes) removed.
     *
     * To avoid z-fighting and wasted triangles on the coincident shared octagons,
     * each shared octagon is drawn ONCE: a cell emits its +x,+y,+z octagons
     * always and skips its -x,-y,-z octagon when the neighbour on that side
     * exists (that neighbour draws it via its + face). Hexagons and squares are
     * never coincident (distinct void-boundary planes), so they are always
     * emitted. Inversion works unchanged (deforms subcell corners before rebuild).
     */

    /**
     * @brief Build the octagon-connected cuboctahedron lattice (one hollow
     *        cuboctahedron per active subcell) and return it as one FacetBox.
     *
     * @param hollow If true (default) each cell is a hollow frame -> the octagons
     *               read as windows connecting neighbours and you see through the
     *               lattice; false -> solid cuboctahedra stuck face-to-face.
     * @param inset  Hollow border inset ratio (only used when hollow=true).
     * @throws std::runtime_error if this cube has no subdivision.
     */
    FacetBox getCuboctahedronLatticeFacets(bool hollow = true,
                                           double inset = default_tc_inset_) const {
        FacetBox fb;
        if (!hasSubdivision())
            throw std::runtime_error("Cube::getCuboctahedronLatticeFacets: no subdivision available");
        const int n = subdivision_levels_;
        // Octagon face indices: 0=Front z+, 1=Back z-, 2=Right x+, 3=Left x-,
        // 4=Top y+, 5=Bottom y-. Skip the negative faces when the neighbour on
        // that side exists (i>0/j>0/k>0) -> that neighbour's + face draws it.
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                for (int k = 0; k < n; ++k) {
                    const SubCell& cell = subcells_[i][j][k];
                    if (!cell.active) continue;
                    const int mask = (i > 0 ? (1 << 3) : 0)   // Left  x-
                                   | (j > 0 ? (1 << 5) : 0)   // Bottom y-
                                   | (k > 0 ? (1 << 1) : 0);  // Back  z-
                    pushTruncatedCuboctahedronFacets(fb, cell.vertices, hollow, inset, mask);
                }
        return fb;
    }

    /**
     * @brief Octagon-connected cuboctahedron lattice on a modular sublattice.
     *
     * Same selection predicate as getCheckerboardSubcells(x,y,z)
     * ((i%x==0 && k%z==0) || j%y==0): only the matching active subcells are
     * triangulated as (hollow) truncated cuboctahedra, and the lattice is
     * "stuck together" through the octagonal faces. To avoid drawing each shared
     * octagon twice, a cell skips its -x/-y/-z octagon ONLY when the neighbour on
     * that side is ALSO selected+active (that neighbour draws the shared octagon
     * via its + face); octagons facing a non-selected neighbour are kept as the
     * cell's outer shell. Hexagons/squares are never coincident and are always
     * emitted. Selection of a sparse modulus therefore yields a set of cuboctahedra
     * (possibly disconnected) whose shared octagons are still drawn once.
     *
     * @param x,y,z Modular periods (>=1) for the i, j, k axes respectively.
     * @param hollow If true, hollow frames (octagons read as windows); false solid.
     * @param inset  Hollow border inset ratio (only used when hollow=true).
     * @throws std::runtime_error if no subdivision; std::invalid_argument if x/y/z < 1.
     */
    FacetBox getCheckerboardFacetsCuboctahedronLattice(int x, int y, int z,
                                                       bool hollow = true,
                                                       double inset = default_tc_inset_) const {
        FacetBox fb;
        if (!hasSubdivision())
            throw std::runtime_error("Cube::getCheckerboardFacetsCuboctahedronLattice: no subdivision available");
        if (x < 1 || y < 1 || z < 1)
            throw std::invalid_argument("Cube::getCheckerboardFacetsCuboctahedronLattice: x, y, z must be >= 1");
        const int n = subdivision_levels_;
        // Same predicate as getCheckerboardSubcells.
        auto selected = [=](int i, int j, int k) {
            return ((i % x == 0 && k % z == 0) || (j % y == 0));
        };
        // Octagon face indices: 0=Front z+, 1=Back z-, 2=Right x+, 3=Left x-,
        // 4=Top y+, 5=Bottom y-. Cull a negative octagon only when the neighbour
        // on that side is also selected+active (it draws the shared octagon).
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                for (int k = 0; k < n; ++k) {
                    const SubCell& cell = subcells_[i][j][k];
                    if (!cell.active || !selected(i, j, k)) continue;
                    const int mask = (i > 0 && subcells_[i-1][j][k].active   && selected(i-1, j, k) ? (1 << 3) : 0)  // Left  x-
                                   | (j > 0 && subcells_[i][j-1][k].active   && selected(i, j-1, k) ? (1 << 5) : 0)  // Bottom y-
                                   | (k > 0 && subcells_[i][j][k-1].active   && selected(i, j, k-1) ? (1 << 1) : 0); // Back  z-
                    pushTruncatedCuboctahedronFacets(fb, cell.vertices, hollow, inset, mask);
                }
        return fb;
    }

    /**
     * @brief Write the octagon-connected cuboctahedron lattice to an ASCII STL.
     *
     * @param filename Path to output STL file (created or overwritten).
     * @param objectName Name to use in the STL header.
     * @param hollow If true (default) hollow frames; false solid cuboctahedra.
     * @param inset  Hollow border inset ratio (only used when hollow=true).
     * @throws std::runtime_error if the file cannot be created or there is no subdivision.
     */
    void writeSTL_s_cuboctahedron_lattice(const std::string& filename,
                                          const std::string& objectName,
                                          bool hollow = true,
                                          double inset = default_tc_inset_) const {
        std::ofstream stl(filename);
        if (!stl.is_open())
            throw std::runtime_error("Cube::writeSTL_s_cuboctahedron_lattice: Cannot create file " + filename);

        FacetBox facetsToWrite = getCuboctahedronLatticeFacets(hollow, inset);

        stl << "solid " << objectName << "\n";
        for (size_t i = 0; i < facetsToWrite.size(); ++i) {
            const Facet& face = facetsToWrite[i];
            Vector3D v1 = face[0], v2 = face[1], v3 = face[2];
            Vector3D normal = face.getNormal();
            stl << "facet normal " << std::scientific
                << normal.x() << " " << normal.y() << " " << normal.z() << "\n";
            stl << "\touter loop\n";
            stl << "\t\tvertex " << v1.x() << " " << v1.y() << " " << v1.z() << "\n";
            stl << "\t\tvertex " << v2.x() << " " << v2.y() << " " << v2.z() << "\n";
            stl << "\t\tvertex " << v3.x() << " " << v3.y() << " " << v3.z() << "\n";
            stl << "\tendloop\n";
            stl << "endfacet\n";
        }
        stl << "endsolid " << objectName << "\n";
        stl.close();
    }

    /**
     * @brief Write the modular (checkerboard) cuboctahedron lattice to ASCII STL.
     *
     * Same selection as getCheckerboardFacetsCuboctahedronLattice(x,y,z,hollow,inset);
     * the STL therefore matches what the render shows for the current <ii,jj,kk>.
     *
     * @throws std::runtime_error if the file cannot be created or there is no subdivision.
     */
    void writeSTL_s_cuboctahedron_lattice_checkerboard(const std::string& filename,
                                                       const std::string& objectName,
                                                       int x, int y, int z,
                                                       bool hollow = true,
                                                       double inset = default_tc_inset_) const {
        std::ofstream stl(filename);
        if (!stl.is_open())
            throw std::runtime_error("Cube::writeSTL_s_cuboctahedron_lattice_checkerboard: Cannot create file " + filename);

        FacetBox facetsToWrite = getCheckerboardFacetsCuboctahedronLattice(x, y, z, hollow, inset);

        stl << "solid " << objectName << "\n";
        for (size_t i = 0; i < facetsToWrite.size(); ++i) {
            const Facet& face = facetsToWrite[i];
            Vector3D v1 = face[0], v2 = face[1], v3 = face[2];
            Vector3D normal = face.getNormal();
            stl << "facet normal " << std::scientific
                << normal.x() << " " << normal.y() << " " << normal.z() << "\n";
            stl << "\touter loop\n";
            stl << "\t\tvertex " << v1.x() << " " << v1.y() << " " << v1.z() << "\n";
            stl << "\t\tvertex " << v2.x() << " " << v2.y() << " " << v2.z() << "\n";
            stl << "\t\tvertex " << v3.x() << " " << v3.y() << " " << v3.z() << "\n";
            stl << "\tendloop\n";
            stl << "endfacet\n";
        }
        stl << "endsolid " << objectName << "\n";
        stl.close();
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

    // === OCTAGONAL (TRUNCATED-CUBE) MESH TABLES ===
    // A second way to build a cube's surface: each of the 6 square faces becomes
    // an octagon (a "truncated cube"), triangulated as an outer octagon ring + an
    // inner octagon ring + a center fan (24 triangles per face) plus 8 corner
    // triangles that close the solid => 152 triangles per cube. It is built on
    // the fly from the SAME 8 corner vertices the plain cube_triangles_ path uses,
    // so it deforms with the lattice under inversion exactly like the 12-triangle
    // mesh and needs no extra storage.
    //
    // octagonal_faces_: the 4 corner indices of each face in CCW-OUTWARD order
    // (the same winding convention as cube_triangles_ above).
    static constexpr int octagonal_faces_[6][4] = {
        {0, 1, 2, 3},  // Front  (z = +radius, outward +z)
        {4, 7, 6, 5},  // Back   (z = -radius, outward -z)
        {1, 5, 6, 2},  // Right  (x = +radius, outward +x)
        {0, 3, 7, 4},  // Left   (x = -radius, outward -x)
        {3, 2, 6, 7},  // Top    (y = +radius, outward +y)
        {0, 4, 5, 1}   // Bottom (y = -radius, outward -y)
    };
    // corner_neighbors_: for each cube corner, its 3 edge-neighbours ordered so the
    // corner triangle (Pa, Pb, Pc) with Pa = Vc + t*(n0 - Vc), Pb = Vc + t*(n1 - Vc),
    // Pc = Vc + t*(n2 - Vc) has an OUTWARD normal (along Vc - cubeCenter). Winding is
    // fixed by this table, so no runtime dot-check is needed and the result does not
    // depend on cell.center (which goes stale under inversion).
    static constexpr int corner_neighbors_[8][3] = {
        {1, 3, 4},  // 0: (-,-,+)
        {0, 5, 2},  // 1: (+,-,+)
        {3, 1, 6},  // 2: (+,+,+)   <-- NOT {1,3,6}; that winding is inward
        {2, 7, 0},  // 3: (-,+,+)
        {5, 0, 7},  // 4: (-,-,-)
        {4, 6, 1},  // 5: (+,-,-)
        {7, 2, 5},  // 6: (+,+,-)
        {6, 4, 3}   // 7: (-,+,-)
    };
    // Regular octagon: truncation t = 1/(2+sqrt(2)) makes the 4 long outer edges
    // (2r(1-2t)) equal to the 4 short cut-corner edges (2rt*sqrt(2)). Hardcoded as
    // a literal because std::sqrt is not constexpr until C++26 and M_SQRT2 is not
    // guaranteed under -std=c++17 (strict). innerScale shrinks the inner ring
    // toward the face centroid (a homothety, so the inner octagon stays regular).
    static constexpr double default_trunc_       = 0.2928932188134524;  // = 1/(2+sqrt(2))
    static constexpr double default_inner_scale_  = 0.5;

    /**
     * @brief Push the 152 triangles of one truncated cube (octagonal mesh) into a FacetBox.
     *
     * Given the 8 corner vertices of a (possibly deformed) cube, emits:
     *   - 6 octagonal faces, 24 triangles each: an outer octagon (truncated square)
     *     and an inner octagon scaled toward the face centroid by innerScale,
     *     joined by 16 annular triangles, plus an 8-triangle center fan.
     *   - 8 corner triangles (one per truncated corner) to close the solid.
     * All winding is fixed by octagonal_faces_ and corner_neighbors_ (outward),
     * matching the cube_triangles_ convention; the Facet constructor computes the
     * deformed normal eagerly from the current vertex positions.
     *
     * @param fb Destination FacetBox (appended to; FacetBox has no reserve, so
     *            reallocations occur — fine for the sizes used here).
     * @param v  The 8 cube corners in standard order (same as initVertices/SubCell).
     * @param trunc Corner-cut fraction in (0, 0.5); pass default_trunc_ for a regular octagon.
     * @param innerScale Inner ring scale toward the face centroid in (0, 1).
     */
    static void pushOctagonalCubeFacets(FacetBox& fb,
                                        const std::array<Vector3D, 8>& v,
                                        double trunc,
                                        double innerScale,
                                        bool hollow) {
        for (int f = 0; f < 6; ++f) {
            const int* F = octagonal_faces_[f];
            Vector3D P[4] = { v[F[0]], v[F[1]], v[F[2]], v[F[3]] };
            Vector3D C = (P[0] + P[1] + P[2] + P[3]) / 4.0;   // face centroid

            // Outer octagon ring in CCW-outward order: [F0, B1, F1, B2, F2, B3, F3, B0]
            // where F[i] = P[i] + t*(P[i+1] - P[i]) (forward along edge i -> i+1)
            // and   B[i] = P[i] + t*(P[i-1] - P[i]) (backward along edge i -> i-1).
            Vector3D outer[8];
            for (int i = 0; i < 4; ++i) {
                int nxt = (i + 1) & 3;
                outer[2 * i]     = P[i]   + trunc * (P[nxt] - P[i]);    // F[i]
                outer[2 * i + 1] = P[nxt] + trunc * (P[i]   - P[nxt]); // B[nxt]
            }
            // Inner ring = homothety of the outer ring toward the face centroid C.
            Vector3D inner[8];
            for (int i = 0; i < 8; ++i)
                inner[i] = C + innerScale * (outer[i] - C);

            // 16 annular triangles between the outer and inner rings (outward).
            for (int i = 0; i < 8; ++i) {
                int j = (i + 1) & 7;
                fb.push(outer[i], outer[j], inner[j]);
                fb.push(outer[i], inner[j], inner[i]);
            }
            // 8 center-fan triangles filling the inner octagon (outward). Skipped
            // when hollow=true, leaving an octagonal hole in each face so the
            // truncated cube is hollow (you can see through the face centers).
            if (!hollow) {
                for (int i = 0; i < 8; ++i) {
                    int j = (i + 1) & 7;
                    fb.push(C, inner[i], inner[j]);
                }
            }
        }
        // 8 corner triangles closing the truncated cube (outward, fixed winding).
        for (int c = 0; c < 8; ++c) {
            const int* n = corner_neighbors_[c];
            Vector3D Vc = v[c];
            fb.push(Vc + trunc * (v[n[0]] - Vc),
                    Vc + trunc * (v[n[1]] - Vc),
                    Vc + trunc * (v[n[2]] - Vc));
        }
    }

    /* === TRUNCATED-OCTAHEDRON MESH — a third method for building a cube ===
     * Each cube FACE becomes a square and each cube CORNER becomes a hexagon
     * (8 hexagons + 6 squares = 14 faces, 24 verts, 36 edges, vertex config
     * 4.6.6) => 44 triangles/cube. Built on the fly from the SAME 8 corner
     * vertices the plain and octagonal paths use, so it deforms with the lattice
     * under inversion exactly like them and needs no extra storage.
     *
     * Construction: the 24 surface vertices lie on the segments from each
     * cube-face center to each cube-face-edge midpoint:
     *   vert(f,e) = faceCenter(f) + s*(edgeMid(f,e) - faceCenter(f))
     * with morph parameter s in (0,1] (default_to_scale_):
     *   s=1.0  -> cuboctahedron (rectified cube; corner faces are small
     *             triangles, the most cube-like end)
     *   s=0.5  -> regular truncated octahedron (8 regular hexagons + 6 squares)
     *   s<0.5  -> corner hexagons stretch larger than regular (squares shrink)
     *   s>1   -> self-intersects (clamped to 1.0)
     * Faces stay planar + convex and the polyhedron stays globally convex for
     * every s in (0,1]. Winding is fixed by octagonal_faces_ (CCW-outward
     * squares) and to_hex_faces_ (CCW-outward hexagons); no runtime dot-check,
     * and the result does not depend on cell.center (which goes stale under
     * inversion). FacetBox results are returned by value.
     */

    // to_hex_faces_: the 6 vertex indices (flat idx = face*4+edge, face order =
    // octagonal_faces_ [Front,Back,Right,Left,Top,Bottom], edge = consecutive
    // corner edge of that face) of each corner hexagon, in CCW-OUTWARD order.
    // Six of the eight required an explicit winding flip vs the natural edge
    // walk (e.g. H2 must be {11,10,17,16,2,1}, not {1,2,16,17,10,11} which is
    // inward) — the same "corner-2 gotcha" category as corner_neighbors_
    // (Cube.hpp:2063). The table stores the verified outward order, so no
    // runtime dot-check is needed and the result does not depend on cell.center.
    static constexpr int to_hex_faces_[8][6] = {
        { 0,  3, 12, 15, 20, 23},  // H0 corner 0 (-,-,+)
        {23, 22,  8, 11,  1,  0},  // H1 corner 1 (+,-,+)
        {11, 10, 17, 16,  2,  1},  // H2 corner 2 (+,+,+)
        {16, 19, 13, 12,  3,  2},  // H3 corner 3 (-,+,+)
        { 4,  7, 21, 20, 15, 14},  // H4 corner 4 (-,-,-)
        { 9,  8, 22, 21,  7,  6},  // H5 corner 5 (+,-,-)
        {18, 17, 10,  9,  6,  5},  // H6 corner 6 (+,+,-)
        {14, 13, 19, 18,  5,  4}   // H7 corner 7 (-,+,-)
    };

    // Morph parameter s (see block comment above). 0.5 = regular truncated
    // octahedron; 1.0 = cuboctahedron (cube-like); (0,0.5) stretches hexagons.
    static constexpr double default_to_scale_ = 0.5;

    // Hollow border inset ratio in (0,1) (used only when hollow=true): each face is
    // inset toward its own centroid by this ratio to form a border frame, and the
    // inner face (the hole) is skipped. inner vertex = faceCentroid + inset*(outer -
    // faceCentroid). Larger inset = bigger hole / thinner border. 0.5 = medium.
    static constexpr double default_to_inset_ = 0.5;

    /**
     * @brief Push the triangles of one truncated octahedron into a FacetBox.
     *
     * Given the 8 corner vertices of a (possibly deformed) cube, emits 6 squares
     * (2 tris each) + 8 hexagons (4 fan tris each) = 44 solid triangles, or — when
     * hollow=true — 144 frame triangles (each face inset toward its own centroid to
     * form a border, the inner face skipped as a hole). The 24 surface vertices are
     * derived as faceCenter + s*(edgeMid - faceCenter) for each of the 6 faces x 4
     * edges, using octagonal_faces_ for the face corner order. Winding is
     * CCW-outward (fixed by octagonal_faces_ and to_hex_faces_; the hollow homothety
     * toward each face centroid preserves that winding, so no frame winding table is
     * needed).
     *
     * @param fb Destination FacetBox (appended to).
     * @param v  The 8 cube corners in standard order (same as initVertices/SubCell).
     * @param s  Morph parameter in (0,1]; clamped. 0.5=regular TO, 1.0=cuboctahedron.
     * @param hollow If true, hollow out each face (inset frame, inner face skipped ->
     *               144 tris); if false, emit the solid 44-tri mesh.
     * @param inset Hollow border inset ratio in (0,1) (only used when hollow=true):
     *              inner vertex = faceCentroid + inset*(outer - faceCentroid).
     *              Larger inset = bigger hole / thinner border. default_to_inset_=0.5.
     */
    static void pushTruncatedOctahedronFacets(FacetBox& fb,
                                              const std::array<Vector3D, 8>& v,
                                              double s,
                                              bool hollow,
                                              double inset) {
        if (!(s > 0.0)) s = 0.0001;   // guard: squares collapse to points at s=0
        if (s > 1.0)    s = 1.0;      // clamp: self-intersects for s>1
        if (hollow) {
            if (!(inset > 0.0)) inset = 0.0001;   // guard: border vanishes at inset=0
            if (inset >= 1.0)  inset = 0.9999;   // clamp: hole vanishes at inset>=1
        }

        // 24 surface vertices, flat index = face*4 + edge.
        Vector3D vert[24];
        for (int f = 0; f < 6; ++f) {
            const int* F = octagonal_faces_[f];
            Vector3D P[4] = { v[F[0]], v[F[1]], v[F[2]], v[F[3]] };
            Vector3D C = (P[0] + P[1] + P[2] + P[3]) / 4.0;        // face center
            for (int e = 0; e < 4; ++e) {
                int nxt = (e + 1) & 3;
                Vector3D em = (P[e] + P[nxt]) / 2.0;              // edge midpoint
                vert[f * 4 + e] = C + s * (em - C);
            }
        }

        if (!hollow) {
            // SOLID: 6 squares (2 tris each) + 8 corner hexagons (4 fan tris each) = 44 tris.
            for (int f = 0; f < 6; ++f) {
                int b = f * 4;
                fb.push(vert[b + 0], vert[b + 1], vert[b + 2]);
                fb.push(vert[b + 0], vert[b + 2], vert[b + 3]);
            }
            for (int c = 0; c < 8; ++c) {
                const int* h = to_hex_faces_[c];
                fb.push(vert[h[0]], vert[h[1]], vert[h[2]]);
                fb.push(vert[h[0]], vert[h[2]], vert[h[3]]);
                fb.push(vert[h[0]], vert[h[3]], vert[h[4]]);
                fb.push(vert[h[0]], vert[h[4]], vert[h[5]]);
            }
            return;
        }

        // HOLLOW: for each face, inset its outer polygon toward that face's own centroid
        // to form a border frame, and skip the inner face (the hole). 144 tris total
        // (6 squares*8 + 8 hexagons*12), all CCW-outward — the homothety toward the
        // centroid preserves each face's outward winding (verified), so no frame
        // winding table is needed.
        // 6 squares: inner[e] = Cf + inset*(outer[e] - Cf), Cf = centroid of the 4 verts.
        for (int f = 0; f < 6; ++f) {
            int b = f * 4;
            Vector3D o[4] = { vert[b + 0], vert[b + 1], vert[b + 2], vert[b + 3] };
            Vector3D Cf = (o[0] + o[1] + o[2] + o[3]) / 4.0;        // face center
            Vector3D in[4];
            for (int e = 0; e < 4; ++e) in[e] = Cf + inset * (o[e] - Cf);
            for (int e = 0; e < 4; ++e) {
                int nxt = (e + 1) & 3;
                fb.push(o[e], o[nxt], in[nxt]);
                fb.push(o[e], in[nxt], in[e]);
            }
        }
        // 8 corner hexagons: inner[i] = Ch + inset*(outer[i] - Ch), Ch = centroid of the 6.
        for (int c = 0; c < 8; ++c) {
            const int* h = to_hex_faces_[c];
            Vector3D o[6] = { vert[h[0]], vert[h[1]], vert[h[2]], vert[h[3]], vert[h[4]], vert[h[5]] };
            Vector3D Ch = (o[0] + o[1] + o[2] + o[3] + o[4] + o[5]) / 6.0;   // face center
            Vector3D in[6];
            for (int i = 0; i < 6; ++i) in[i] = Ch + inset * (o[i] - Ch);
            for (int i = 0; i < 6; ++i) {
                int nxt = (i + 1) % 6;
                fb.push(o[i], o[nxt], in[nxt]);
                fb.push(o[i], in[nxt], in[i]);
            }
        }
    }

    /* === TRUNCATED-CUBOCTAHEDRON MESH — a fourth method for building a cube ===
     * The truncated cuboctahedron (great rhombicuboctahedron) is the Archimedean
     * solid with vertex configuration 4.6.8: 6 octagons (one per cube FACE),
     * 8 hexagons (one per cube CORNER) and 12 squares (one per cube EDGE) = 26
     * faces, 48 vertices, 72 edges => 92 triangles/cube (solid). It is the
     * cantitruncated (omnitruncated) cube, t_{0,1,2}{4,3}.
     *
     * Construction (inversion-compatible, like the other three methods): the 48
     * surface vertices are the canonical positions (permutations of (±a, ±b, ±1)
     * in normalized cube space, a = 1/(1+2√2) the corner-cut coord, b = (1+√2)/(1+2√2)
     * the edge-cut coord, ±1 selects the cube face). Each vertex is computed ONCE
     * by trilinear interpolation of the SAME 8 cube corners the other methods use,
     * so it is a fixed affine combination of the 8 corners -> it deforms with the
     * lattice under inversion exactly like them and needs no extra storage.
     * Because every vertex is computed once and referenced from all three of its
     * faces (octagon+hexagon+square), the mesh is watertight under deformation.
     *
     * The four tables below were generated AND verified by
     * truncated_cuboctahedron_verify.py (48 verts, 92 tris, V−E+F=2, every mesh
     * edge shared by 2, all 72 polyhedron edges length 2, every face regular, all
     * normals outward). Winding is CCW-outward, fixed by the tables, so no runtime
     * dot-check is needed and the result does not depend on cell.center (which
     * goes stale under inversion). FacetBox results are returned by value.
     *
     * Hollow (mirrors the truncated-octahedron hollow convention): each face is
     * inset toward its own centroid by `inset` (a homothety, which preserves
     * winding and planarity) and the inner face is skipped as a hole -> 288 tris
     * (6 octagons*16 + 8 hexagons*12 + 12 squares*8).
     */

    // tc_corner_sign_: sign pattern (sx,sy,sz) of each cube corner in the
    // standard initVertices/SubCell order 0..7. Used by the trilinear interp.
    static constexpr int tc_corner_sign_[8][3] = {
        {-1, -1, +1},  // 0: (-,-,+)
        {+1, -1, +1},  // 1: (+,-,+)
        {+1, +1, +1},  // 2: (+,+,+)
        {-1, +1, +1},  // 3: (-,+,+)
        {-1, -1, -1},  // 4: (-,-,-)
        {+1, -1, -1},  // 5: (+,-,-)
        {+1, +1, -1},  // 6: (+,+,-)
        {-1, +1, -1}   // 7: (-,+,-)
    };

    // tc_vert_norm_: the 48 surface-vertex normalized positions (perms of
    // (±a, ±b, ±1) with a≈0.2612038750 = 1/(1+2√2), b≈0.6306019375 = (1+√2)/(1+2√2)).
    // Vertex g = Σ_i w_i * corner_i, w_i = (1/8)·∏_axis (1 + s_i·coord) — trilinear
    // interpolation of the 8 corners at tc_vert_norm_[g] (cube half-side = 1 here).
    static constexpr double tc_vert_norm_[48][3] = {
        { 0.2612038750, 0.6306019375, 1.0000000000 },  // 0
        { 0.2612038750, 0.6306019375, -1.0000000000 },  // 1
        { 0.2612038750, -0.6306019375, 1.0000000000 },  // 2
        { 0.2612038750, -0.6306019375, -1.0000000000 },  // 3
        { -0.2612038750, 0.6306019375, 1.0000000000 },  // 4
        { -0.2612038750, 0.6306019375, -1.0000000000 },  // 5
        { -0.2612038750, -0.6306019375, 1.0000000000 },  // 6
        { -0.2612038750, -0.6306019375, -1.0000000000 },  // 7
        { 0.2612038750, 1.0000000000, 0.6306019375 },  // 8
        { 0.2612038750, 1.0000000000, -0.6306019375 },  // 9
        { 0.2612038750, -1.0000000000, 0.6306019375 },  // 10
        { 0.2612038750, -1.0000000000, -0.6306019375 },  // 11
        { -0.2612038750, 1.0000000000, 0.6306019375 },  // 12
        { -0.2612038750, 1.0000000000, -0.6306019375 },  // 13
        { -0.2612038750, -1.0000000000, 0.6306019375 },  // 14
        { -0.2612038750, -1.0000000000, -0.6306019375 },  // 15
        { 0.6306019375, 0.2612038750, 1.0000000000 },  // 16
        { 0.6306019375, 0.2612038750, -1.0000000000 },  // 17
        { 0.6306019375, -0.2612038750, 1.0000000000 },  // 18
        { 0.6306019375, -0.2612038750, -1.0000000000 },  // 19
        { -0.6306019375, 0.2612038750, 1.0000000000 },  // 20
        { -0.6306019375, 0.2612038750, -1.0000000000 },  // 21
        { -0.6306019375, -0.2612038750, 1.0000000000 },  // 22
        { -0.6306019375, -0.2612038750, -1.0000000000 },  // 23
        { 0.6306019375, 1.0000000000, 0.2612038750 },  // 24
        { 0.6306019375, 1.0000000000, -0.2612038750 },  // 25
        { 0.6306019375, -1.0000000000, 0.2612038750 },  // 26
        { 0.6306019375, -1.0000000000, -0.2612038750 },  // 27
        { -0.6306019375, 1.0000000000, 0.2612038750 },  // 28
        { -0.6306019375, 1.0000000000, -0.2612038750 },  // 29
        { -0.6306019375, -1.0000000000, 0.2612038750 },  // 30
        { -0.6306019375, -1.0000000000, -0.2612038750 },  // 31
        { 1.0000000000, 0.2612038750, 0.6306019375 },  // 32
        { 1.0000000000, 0.2612038750, -0.6306019375 },  // 33
        { 1.0000000000, -0.2612038750, 0.6306019375 },  // 34
        { 1.0000000000, -0.2612038750, -0.6306019375 },  // 35
        { -1.0000000000, 0.2612038750, 0.6306019375 },  // 36
        { -1.0000000000, 0.2612038750, -0.6306019375 },  // 37
        { -1.0000000000, -0.2612038750, 0.6306019375 },  // 38
        { -1.0000000000, -0.2612038750, -0.6306019375 },  // 39
        { 1.0000000000, 0.6306019375, 0.2612038750 },  // 40
        { 1.0000000000, 0.6306019375, -0.2612038750 },  // 41
        { 1.0000000000, -0.6306019375, 0.2612038750 },  // 42
        { 1.0000000000, -0.6306019375, -0.2612038750 },  // 43
        { -1.0000000000, 0.6306019375, 0.2612038750 },  // 44
        { -1.0000000000, 0.6306019375, -0.2612038750 },  // 45
        { -1.0000000000, -0.6306019375, 0.2612038750 },  // 46
        { -1.0000000000, -0.6306019375, -0.2612038750 }  // 47
    };

    // tc_oct_faces_ / tc_hex_faces_ / tc_sq_faces_: global vertex indices (into
    // the 48 above) of each face, in CCW-OUTWARD order (verified by the script).
    // Octagon order matches octagonal_faces_ [Front z+, Back z-, Right x+, Left x-, Top y+, Bottom y-].
    static constexpr int tc_oct_faces_[6][8] = {
        { 22,  6,  2, 18, 16,  0,  4, 20 },  // octagon 0  Front  (z = +radius)
        { 21,  5,  1, 17, 19,  3,  7, 23 },  // octagon 1  Back   (z = -radius)
        { 43, 35, 33, 41, 40, 32, 34, 42 },  // octagon 2  Right  (x = +radius)
        { 46, 38, 36, 44, 45, 37, 39, 47 },  // octagon 3  Left   (x = -radius)
        { 28, 12,  8, 24, 25,  9, 13, 29 },  // octagon 4  Top    (y = +radius)
        { 31, 15, 11, 27, 26, 10, 14, 30 }   // octagon 5  Bottom (y = -radius)
    };
    static constexpr int tc_hex_faces_[8][6] = {
        {  0, 16, 32, 40, 24,  8 },  // hexagon 0  corner 0 (-,-,+)
        {  9, 25, 41, 33, 17,  1 },  // hexagon 1  corner 1 (+,-,+)
        { 10, 26, 42, 34, 18,  2 },  // hexagon 2  corner 2 (+,+,+)
        {  3, 19, 35, 43, 27, 11 },  // hexagon 3  corner 3 (-,+,+)
        { 36, 20,  4, 12, 28, 44 },  // hexagon 4  corner 4 (-,-,-)
        { 45, 29, 13,  5, 21, 37 },  // hexagon 5  corner 5 (+,-,-)
        { 46, 30, 14,  6, 22, 38 },  // hexagon 6  corner 6 (+,+,-)
        { 39, 23,  7, 15, 31, 47 }   // hexagon 7  corner 7 (-,+,-)
    };
    static constexpr int tc_sq_faces_[12][4] = {
        { 24, 40, 41, 25 },  // square 0   edge +x,+y
        { 27, 43, 42, 26 },  // square 1   edge +x,-y
        { 44, 28, 29, 45 },  // square 2   edge -x,+y
        { 47, 31, 30, 46 },  // square 3   edge -x,-y
        { 18, 34, 32, 16 },  // square 4   edge +x,+z
        { 17, 33, 35, 19 },  // square 5   edge +x,-z
        { 38, 22, 20, 36 },  // square 6   edge -x,+z
        { 37, 21, 23, 39 },  // square 7   edge -x,-z
        {  4,  0,  8, 12 },  // square 8   edge +y,+z
        { 13,  9,  1,  5 },  // square 9   edge +y,-z
        { 14, 10,  2,  6 },  // square 10  edge -y,+z
        {  7,  3, 11, 15 }   // square 11  edge -y,-z
    };

    // Hollow border inset ratio in (0,1) (used only when hollow=true): each face
    // is inset toward its own centroid by this ratio and the inner face is skipped.
    // inner = centroid + inset*(outer - centroid). Larger inset = bigger hole.
    static constexpr double default_tc_inset_ = 0.5;

    /**
     * @brief Push the triangles of one truncated cuboctahedron into a FacetBox.
     *
     * Given the 8 corner vertices of a (possibly deformed) cube, emits 6 octagons
     * (6 fan tris each) + 8 hexagons (4 each) + 12 squares (2 each) = 92 solid
     * triangles, or — when hollow=true — 288 frame triangles (each face inset
     * toward its own centroid, inner face skipped). The 48 surface vertices are
     * derived by trilinear interpolation of the 8 corners at tc_vert_norm_[g].
     * Winding is CCW-outward, fixed by tc_oct_faces_/tc_hex_faces_/tc_sq_faces_;
     * the hollow homothety toward each face centroid preserves that winding.
     *
     * @param fb Destination FacetBox (appended to).
     * @param v  The 8 cube corners in standard order (same as initVertices/SubCell).
     * @param hollow If true, emit the 288-tri hollow frame; false, the 92-tri solid mesh.
     * @param inset Hollow border inset ratio in (0,1) (only used when hollow=true).
     */
    static void pushTruncatedCuboctahedronFacets(FacetBox& fb,
                                                  const std::array<Vector3D, 8>& v,
                                                  bool hollow,
                                                  double inset,
                                                  int cullOctMask = 0) {
        if (hollow) {
            if (!(inset > 0.0)) inset = 0.0001;   // guard: border vanishes at inset=0
            if (inset >= 1.0)  inset = 0.9999;   // clamp: hole vanishes at inset>=1
        }

        // 48 surface vertices via trilinear interpolation of the 8 cube corners
        // at the normalized positions tc_vert_norm_[g] (perms of (±a, ±b, ±1)).
        Vector3D vert[48];
        for (int g = 0; g < 48; ++g) {
            const double x = tc_vert_norm_[g][0];
            const double y = tc_vert_norm_[g][1];
            const double z = tc_vert_norm_[g][2];
            Vector3D p;                                   // zero
            for (int i = 0; i < 8; ++i) {
                const double w = (1.0 / 8.0)
                    * (1.0 + tc_corner_sign_[i][0] * x)
                    * (1.0 + tc_corner_sign_[i][1] * y)
                    * (1.0 + tc_corner_sign_[i][2] * z);
                p += w * v[i];
            }
            vert[g] = p;
        }

        if (!hollow) {
            // SOLID: fan each face CCW-outward. 6*6 + 8*4 + 12*2 = 92 tris.
            // cullOctMask bit i skips octagon face i (used by the cuboctahedron
            // lattice to draw each shared octagon once; 0 = all 6, the default).
            for (int i = 0; i < 6; ++i) {
                if (cullOctMask & (1 << i)) continue;
                const int* f = tc_oct_faces_[i];
                for (int t = 1; t < 7; ++t)
                    fb.push(vert[f[0]], vert[f[t]], vert[f[t + 1]]);
            }
            for (int i = 0; i < 8; ++i) {
                const int* f = tc_hex_faces_[i];
                for (int t = 1; t < 5; ++t)
                    fb.push(vert[f[0]], vert[f[t]], vert[f[t + 1]]);
            }
            for (int i = 0; i < 12; ++i) {
                const int* f = tc_sq_faces_[i];
                fb.push(vert[f[0]], vert[f[1]], vert[f[2]]);
                fb.push(vert[f[0]], vert[f[2]], vert[f[3]]);
            }
            return;
        }

        // HOLLOW: each face inset toward its own centroid by `inset`; inner face
        // skipped. Homothety preserves CCW-outward winding. 288 tris total.
        auto hollowFace = [&](const int* f, int n) {
            Vector3D o[8];
            Vector3D C;                                   // zero -> centroid
            for (int i = 0; i < n; ++i) { o[i] = vert[f[i]]; C += o[i]; }
            C = C / static_cast<double>(n);
            Vector3D in[8];
            for (int i = 0; i < n; ++i) in[i] = C + inset * (o[i] - C);
            for (int i = 0; i < n; ++i) {
                int j = (i + 1) % n;
                fb.push(o[i], o[j], in[j]);
                fb.push(o[i], in[j], in[i]);
            }
        };
        for (int i = 0; i < 6; ++i)  { if (cullOctMask & (1 << i)) continue; hollowFace(tc_oct_faces_[i], 8); }
        for (int i = 0; i < 8; ++i)  hollowFace(tc_hex_faces_[i], 6);
        for (int i = 0; i < 12; ++i) hollowFace(tc_sq_faces_[i], 4);
    }

    /**
     * @brief Push the 12 plain-cube triangles into a FacetBox (cube_triangles_ walk).
     *
     * Used by the cantitruncated honeycomb builder for the cube cells (3 per 2x2x2
     * block), so the per-cell loop allocates no temporary FacetBox.
     */
    static void pushCubeFacets(FacetBox& fb, const std::array<Vector3D, 8>& v) {
        for (int i = 0; i < 12; ++i) {
            const int* t = cube_triangles_[i];
            fb.push(v[t[0]], v[t[1]], v[t[2]]);
        }
    }

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
