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
 * @brief Represents a cube in 3D space with triangulation capabilities.
 *
 * Stores 8 vertices (as Quaternion points) and a FacetBox of triangular faces.
 * Can be subdivided using the standard cube triangulation with 12 triangles per subcube.
 * Uses Rule-of-Zero for special members and std::vector internally.
 */

class Cube {
public:
    friend std::ostream& operator<<(std::ostream& os, const Cube& cube);

    // === SUBDIVISION DATA STRUCTURES ===
    /**
     * @brief Represents a single subcell within a subdivided cube
     * Each subcell is a smaller cube with its own geometry and properties
     */
    struct SubCell {
        std::array<Vector3D, 8> vertices;  // 8 cube vertices in standard order (0-7)
        Vector3D center;                   // Center point of subcell in 3D space
        double radius;                     // Half side length (distance from center to face)
        int i, j, k;                      // Physical grid coordinates [0, n) in subdivision
        bool active = true;               // Can be disabled for selective rendering/processing

        // Default constructor - creates empty subcell
        SubCell() = default;

        // Parameterized constructor - initializes all subcell properties
        SubCell(const std::array<Vector3D, 8>& verts, const Vector3D& c, double r, int x, int y, int z)
            : vertices(verts), center(c), radius(r), i(x), j(y), k(z) {}  // Initialize all members
    };

public:
    /* - Constructors */
    /* - Rule of Zero: default special members */
    // Default constructor - creates empty cube with no subdivision
    Cube() : subdivision_levels_(0), subdivision_built_(false) {}

    // Destructor - compiler-generated cleanup is sufficient
    ~Cube() = default;

    // Copy constructor - deep copy of all cube data
    Cube(const Cube&) = default;

    // Move constructor - efficient transfer of resources without copying
    Cube(Cube&&) noexcept = default;

    // Copy assignment operator - assigns from another cube
    Cube& operator=(const Cube&) = default;

    // Move assignment operator - efficient assignment via resource transfer
    Cube& operator=(Cube&&) noexcept = default;

    /**
     * @brief Construct a cube with given radius (half side length) centered at 'center'.
     * @param radius Distance from center to face (half of side length).
     * @param center Geometric center of the cube.
     */
    // Parameterized constructor - creates basic cube at specified location and size
    Cube(double radius, const Vector3D& center) : subdivision_levels_(0), subdivision_built_(false) {
        initVertices(radius, center);  // Generate the 8 corner vertices of the cube
        buildFacets();                 // Create 12 triangular faces from vertices
    }

    /**
     * @brief Construct a subdivided cube with n subdivision levels.
     * @param radius Distance from center to face.
     * @param center Geometric center of the cube.
     * @param subdivisions Number of subdivision levels (0 = basic cube).
     */
    // Constructor with subdivision - creates cube with n levels of internal division
    Cube(double radius, const Vector3D& center, int subdivisions) : subdivision_levels_(0), subdivision_built_(false) {
        initVertices(radius, center);  // Generate the 8 corner vertices of the main cube
        buildFacets();                 // Create basic 12 triangular faces
        if (subdivisions > 0) {        // Only subdivide if positive subdivision count
            subdivide(subdivisions);   // Break cube into n×n×n smaller subcubes
        }
    }

    /* - Face Count */
    /**
     * @brief Number of triangular faces.
     */
    // Returns total number of triangular faces in the cube (12 for basic, more for subdivided)
    size_t faceCount() const noexcept { return facets_.size(); }

    /* - Facet Access */
    /**
     * @brief Access the k-th triangular face as a Facet.
     * @param k Index in range [0, faceCount()).
     * @throws std::out_of_range if k is invalid.
     */
    // Array-style access to individual triangular faces with bounds checking
    const Facet& operator[](size_t k) const {
        return facets_[k];  // Return reference to the k-th triangular face
    }

    /* - Geometric Center */
    /**
     * @brief Compute and return the geometric center of all vertices.
     */
    // Calculate geometric center by averaging all vertex positions
    Vector3D center() const {
        if (verts_.empty()) {  // Check if cube has been initialized
            throw std::runtime_error("Cube::center: no vertices initialized");
        }
        Vector3D sum{0,0,0};  // Initialize accumulator vector to zero
        for (auto const& q : verts_) {  // Iterate through all vertices
            sum += q.V();  // Add each vertex position to the sum
        }
        return sum / static_cast<double>(verts_.size());  // Return average position
    }

    /* - Translate by offset */
    /**
     * @brief Translate the entire cube by 'offset'.
     */
    // Move the entire cube by adding offset to each vertex position
    void translate(const Vector3D& offset) {
        for (auto& q : verts_) {  // Iterate through all vertices
            Vector3D p = q.V() + offset;  // Calculate new position with offset
            q = Quaternion(0.0, p);       // Update quaternion with new position
        }
        buildFacets();  // Rebuild triangular faces with new vertex positions
    }

    /* - Scale */
    /**
     * @brief Scale all vertices relative to a pivot point.
     * @param s Scaling factor.
     * @param pivot The point to scale around.
     */
    // Scale cube relative to a pivot point by factor s
    void scale(double s, const Vector3D& pivot) {
        for (auto& q : verts_) {  // Process each vertex
            Vector3D p = pivot + s * (q.V() - pivot);  // Scale distance from pivot
            q = Quaternion(0.0, p);  // Update quaternion with scaled position
        }
        buildFacets();  // Rebuild triangular faces with scaled vertices
    }

    /**
     * @brief Retrieve the underlying FacetBox.
     */
    // Return reference to the container holding all triangular faces
    const FacetBox& getFacets() const noexcept {
        return facets_;  // Direct access to internal facet storage
    }

    /**
     * @brief Subdivide the cube into n levels of smaller cubes, each triangulated.
     * @param n Number of subdivision levels.
     */
    // Break cube into n×n×n smaller subcubes, each with 12 triangular faces
    void subdivide(int n) {
        if (n <= 0) return;  // Exit early for invalid subdivision count

        // Clear current facets and rebuild with subdivision
        Vector3D cubeCenter = center();    // Get current cube center
        double cubeRadius = getRadius();   // Get current cube radius

        // Store subdivision info for later access
        subdivision_levels_ = n;

        // Rebuild with subdivision - replace simple cube with complex subdivided one
        facets_.clear();  // Remove existing triangular faces
        buildSubdividedCube(cubeCenter, cubeRadius, n);  // Generate n³ subcubes
        subdivision_built_ = true;  // Mark subdivision as completed
    }

    // === SUBDIVISION ACCESS METHODS ===

    /**
     * @brief Check if subdivision data is available.
     */
    // Check if the cube has been subdivided and subdivision data is available
    bool hasSubdivision() const noexcept {
        return subdivision_built_ && !subcells_.empty();  // Both flag set and data present
    }

    // === COORDINATE TRANSFORMATION ===

    /**
     * @brief Convert logical coordinates (centered at origin) to physical grid indices.
     * @param x Logical X coordinate (can be negative)
     * @param y Logical Y coordinate (can be negative)
     * @param z Logical Z coordinate (can be negative)
     * @returns Tuple of (i,j,k) physical grid indices
     *
     * For n=5: logical (0,0,0) -> physical (2,2,2), logical (2,2,2) -> physical (4,4,4)
     * For n=4: logical (0,0,0) -> physical (2,2,2), logical (1,1,1) -> physical (3,3,3)
     */
    // Convert user-friendly centered coordinates to internal grid indices
    std::tuple<int,int,int> logicalToPhysical(int x, int y, int z) const {
        if (!hasSubdivision()) {  // Ensure subdivision data exists
            throw std::runtime_error("Cube::logicalToPhysical: no subdivision available");
        }

        int center = subdivision_levels_ / 2;  // Integer division gives us the center offset
        int i = center + x;  // Convert logical X to physical grid index
        int j = center + y;  // Convert logical Y to physical grid index
        int k = center + z;  // Convert logical Z to physical grid index

        return std::make_tuple(i, j, k);  // Return as tuple for structured binding
    }

    /**
     * @brief Convert physical grid indices to logical coordinates (centered at origin).
     * @param i Physical X grid index [0, n)
     * @param j Physical Y grid index [0, n)
     * @param k Physical Z grid index [0, n)
     * @returns Tuple of (x,y,z) logical coordinates
     */
    // Convert internal grid indices to user-friendly centered coordinates
    std::tuple<int,int,int> physicalToLogical(int i, int j, int k) const {
        if (!hasSubdivision()) {  // Ensure subdivision data exists
            throw std::runtime_error("Cube::physicalToLogical: no subdivision available");
        }

        int center = subdivision_levels_ / 2;  // Calculate center offset
        int x = i - center;  // Convert physical grid index to logical X
        int y = j - center;  // Convert physical grid index to logical Y
        int z = k - center;  // Convert physical grid index to logical Z

        return std::make_tuple(x, y, z);  // Return as tuple for structured binding
    }

    /**
     * @brief Get subdivision levels (n means n³ subcells).
     */
    // Return the current subdivision level (n means n×n×n subcubes total)
    int getSubdivisionLevels() const noexcept {
        return subdivision_levels_;  // Direct access to subdivision count
    }

    /**
     * @brief Access a specific subcell by 3D logical coordinates (centered at origin).
     * @param x Logical X-coordinate (can be negative, 0=center)
     * @param y Logical Y-coordinate (can be negative, 0=center)
     * @param z Logical Z-coordinate (can be negative, 0=center)
     * @throws std::out_of_range if coordinates invalid or no subdivision
     */
    // Get read-only access to a specific subcell using logical coordinates
    const SubCell& getSubCell(int x, int y, int z) const {
        auto [i, j, k] = logicalToPhysical(x, y, z);  // Convert to physical grid indices

        // Validate that physical coordinates are within grid bounds
        if (i < 0 || i >= subdivision_levels_ ||
            j < 0 || j >= subdivision_levels_ ||
            k < 0 || k >= subdivision_levels_) {
            throw std::out_of_range("Cube::getSubCell: logical coordinates out of range");
        }
        return subcells_[i][j][k];  // Return reference to specific subcell
    }

    /**
     * @brief Access a specific subcell by physical grid indices (for internal use).
     * @param i Physical X-coordinate in grid [0, n)
     * @param j Physical Y-coordinate in grid [0, n)
     * @param k Physical Z-coordinate in grid [0, n)
     * @throws std::out_of_range if coordinates invalid or no subdivision
     */
    // Get read-only access to subcell using direct physical grid indices
    const SubCell& getSubCellPhysical(int i, int j, int k) const {
        if (!hasSubdivision()) {  // Ensure subdivision data exists
            throw std::runtime_error("Cube::getSubCellPhysical: no subdivision available");
        }
        // Validate that all coordinates are within grid bounds [0, n)
        if (i < 0 || i >= subdivision_levels_ ||
            j < 0 || j >= subdivision_levels_ ||
            k < 0 || k >= subdivision_levels_) {
            throw std::out_of_range("Cube::getSubCellPhysical: coordinates out of range");
        }
        return subcells_[i][j][k];  // Direct array access to subcell
    }

    /**
     * @brief Get mutable access to a subcell for modifications (using logical coordinates).
     */
    // Get modifiable access to a subcell using logical coordinates
    SubCell& getSubCellMutable(int x, int y, int z) {
        auto [i, j, k] = logicalToPhysical(x, y, z);  // Convert to physical indices

        // Validate that physical coordinates are within grid bounds
        if (i < 0 || i >= subdivision_levels_ ||
            j < 0 || j >= subdivision_levels_ ||
            k < 0 || k >= subdivision_levels_) {
            throw std::out_of_range("Cube::getSubCellMutable: logical coordinates out of range");
        }
        return subcells_[i][j][k];  // Return modifiable reference to subcell
    }

    /**
     * @brief Get mutable access to a subcell by physical indices (for internal use).
     */
    // Get modifiable access to subcell using direct physical grid indices
    SubCell& getSubCellMutablePhysical(int i, int j, int k) {
        if (!hasSubdivision()) {  // Ensure subdivision data exists
            throw std::runtime_error("Cube::getSubCellMutablePhysical: no subdivision available");
        }
        // Validate that all coordinates are within grid bounds [0, n)
        if (i < 0 || i >= subdivision_levels_ ||
            j < 0 || j >= subdivision_levels_ ||
            k < 0 || k >= subdivision_levels_) {
            throw std::out_of_range("Cube::getSubCellMutablePhysical: coordinates out of range");
        }
        return subcells_[i][j][k];  // Direct modifiable array access to subcell
    }

    /**
     * @brief Get all subcells in a specific plane.
     * @param axis 0=YZ plane (varying X), 1=XZ plane (varying Y), 2=XY plane (varying Z)
     * @param layer Layer index [0, n)
     * @returns Vector of subcell references
     */
    // Extract all subcells from a 2D slice of the 3D subdivision grid
    std::vector<std::reference_wrapper<const SubCell>> getPlane(int axis, int layer) const {
        if (!hasSubdivision()) {  // Ensure subdivision data exists
            throw std::runtime_error("Cube::getPlane: no subdivision available");
        }
        // Validate axis (0-2) and layer (0 to n-1)
        if (axis < 0 || axis > 2 || layer < 0 || layer >= subdivision_levels_) {
            throw std::out_of_range("Cube::getPlane: invalid axis or layer");
        }

        std::vector<std::reference_wrapper<const SubCell>> plane;  // Container for plane subcells
        int n = subdivision_levels_;  // Cache subdivision count for loops

        switch (axis) {
            case 0: // YZ plane (fixed X = layer) - slice perpendicular to X-axis
                for (int j = 0; j < n; ++j) {      // Iterate through Y coordinates
                    for (int k = 0; k < n; ++k) {  // Iterate through Z coordinates
                        plane.emplace_back(subcells_[layer][j][k]);  // Add subcell reference
                    }
                }
                break;
            case 1: // XZ plane (fixed Y = layer) - slice perpendicular to Y-axis
                for (int i = 0; i < n; ++i) {      // Iterate through X coordinates
                    for (int k = 0; k < n; ++k) {  // Iterate through Z coordinates
                        plane.emplace_back(subcells_[i][layer][k]);  // Add subcell reference
                    }
                }
                break;
            case 2: // XY plane (fixed Z = layer) - slice perpendicular to Z-axis
                for (int i = 0; i < n; ++i) {      // Iterate through X coordinates
                    for (int j = 0; j < n; ++j) {  // Iterate through Y coordinates
                        plane.emplace_back(subcells_[i][j][layer]);  // Add subcell reference
                    }
                }
                break;
        }
        return plane;  // Return vector of subcell references
    }

    /**
     * @brief Update a specific vertex of a subcell and mark for retriangulation.
     * @param x,y,z Logical subcell coordinates (centered at origin)
     * @param vertex_idx Vertex index [0,7] within the subcell
     * @param new_pos New position for the vertex
     */
    // Modify a single vertex position within a specific subcell
    void updateSubCellVertex(int x, int y, int z, int vertex_idx, const Vector3D& new_pos) {
        if (vertex_idx < 0 || vertex_idx > 7) {  // Validate vertex index (cube has 8 vertices)
            throw std::out_of_range("Cube::updateSubCellVertex: vertex_idx must be [0,7]");
        }
        SubCell& cell = getSubCellMutable(x, y, z);  // Get modifiable reference to subcell
        cell.vertices[vertex_idx] = new_pos;         // Update the specified vertex position

        // Mark as needing retriangulation
        // (will be handled in refreshTriangulation())
    }

    /**
     * @brief Rebuild all triangulation after vertex modifications.
     * Call this after updating vertices to refresh the visual mesh.
     */
    // Regenerate all triangular faces after subcell vertex modifications
    void refreshTriangulation() {
        if (!hasSubdivision()) return;  // Exit if no subdivision exists

        facets_.clear();  // Remove all existing triangular faces

        // Rebuild triangles from current subcell vertices (using physical indices internally)
        for (int i = 0; i < subdivision_levels_; ++i) {      // Iterate through X dimension
            for (int j = 0; j < subdivision_levels_; ++j) {  // Iterate through Y dimension
                for (int k = 0; k < subdivision_levels_; ++k) {  // Iterate through Z dimension
                    const SubCell& cell = subcells_[i][j][k];  // Get reference to current subcell
                    if (!cell.active) continue;  // Skip disabled subcells for selective rendering

                    // Generate 12 triangles from current vertices
                    for (int t = 0; t < 12; ++t) {  // Each cube has 12 triangular faces
                        auto const& tri = cube_triangles_[t];  // Get vertex indices for triangle t
                        // Create triangle from 3 vertices and add to facet collection
                        facets_.push(cell.vertices[tri[0]], cell.vertices[tri[1]], cell.vertices[tri[2]]);
                    }
                }
            }
        }
    }

    /**
     * @brief Get facets for a specific subcell.
     * @param x,y,z Logical coordinates (centered at origin)
     * @returns FacetBox containing the 12 triangles of the specified subcell
     */
    // Generate triangular faces for a single subcell at logical coordinates
    FacetBox getSubCellFacets(int x, int y, int z) const {
        const SubCell& cell = getSubCell(x, y, z);  // Get reference to target subcell
        FacetBox subcell_facets;  // Container for this subcell's triangular faces

        if (!cell.active) return subcell_facets; // Return empty if disabled

        // Generate 12 triangular faces from the subcell's 8 vertices
        for (int t = 0; t < 12; ++t) {  // Each cube has 12 triangular faces
            auto const& tri = cube_triangles_[t];  // Get vertex indices for triangle t
            // Create triangle from 3 vertices and add to face collection
            subcell_facets.push(cell.vertices[tri[0]], cell.vertices[tri[1]], cell.vertices[tri[2]]);
        }
        return subcell_facets;  // Return complete set of subcell faces
    }

    /**
     * @brief Get facets for all subcells in a plane.
     */
    // Generate triangular faces for all subcells within a specified 2D plane
    FacetBox getPlaneFacets(int axis, int layer) const {
        auto plane_subcells = getPlane(axis, layer);  // Get all subcells in the plane
        FacetBox plane_facets;  // Container for all plane triangular faces

        // Process each subcell in the plane
        for (const SubCell& cell : plane_subcells) {
            if (!cell.active) continue;  // Skip disabled subcells

            // Generate 12 triangular faces for each active subcell
            for (int t = 0; t < 12; ++t) {  // Each cube has 12 triangular faces
                auto const& tri = cube_triangles_[t];  // Get vertex indices for triangle t
                // Create triangle from 3 vertices and add to plane collection
                plane_facets.push(cell.vertices[tri[0]], cell.vertices[tri[1]], cell.vertices[tri[2]]);
            }
        }
        return plane_facets;  // Return complete set of plane faces
    }

    /**
     * @brief Get subcells that form a checkerboard pattern where (i+j+k)%2 == 0.
     * @returns Vector of subcell references following the 3D checkerboard pattern
     *
     * This creates a 3D checkerboard pattern similar to alternating black/white squares
     * on a 2D checkerboard, but extended to 3 dimensions. Useful for creating patterns,
     * selective rendering, or algorithmic processing of subcells.
     */
    // Extract subcells following a 3D checkerboard pattern for selective operations
    std::vector<std::reference_wrapper<const SubCell>> getCheckerboardSubcells() const {
        if (!hasSubdivision()) {  // Ensure subdivision data exists
            throw std::runtime_error("Cube::getCheckerboardSubcells: no subdivision available");
        }

        std::vector<std::reference_wrapper<const SubCell>> checkerboard;  // Container for pattern subcells
        int n = subdivision_levels_;  // Cache subdivision count for loops

        // Iterate through all subcells using physical coordinates
        for (int i = 0; i < n; ++i) {      // X dimension
            for (int j = 0; j < n; ++j) {  // Y dimension
                for (int k = 0; k < n; ++k) {  // Z dimension
                    // Select subcells where sum of coordinates is divisible by 4 (checkerboard pattern)
                    if ((i + j + k) % 4 == 0) {
                        checkerboard.emplace_back(subcells_[i][j][k]);  // Add subcell reference
                    }
                }
            }
        }
        return checkerboard;  // Return vector of pattern-matching subcells
    }

    /**
     * @brief Get facets for all subcells that form a checkerboard pattern where (i+j+k)%2 == 0.
     * @returns FacetBox containing triangles for all checkerboard-pattern subcells
     */
    // Generate triangular faces for all subcells in the checkerboard pattern
    FacetBox getCheckerboardFacets() const {
        auto checkerboard_subcells = getCheckerboardSubcells();  // Get pattern subcells
        FacetBox checkerboard_facets;  // Container for pattern triangular faces

        // Process each subcell in the checkerboard pattern
        for (const SubCell& cell : checkerboard_subcells) {
            if (!cell.active) continue;  // Skip disabled subcells

            // Generate 12 triangular faces for each active pattern subcell
            for (int t = 0; t < 12; ++t) {  // Each cube has 12 triangular faces
                auto const& tri = cube_triangles_[t];  // Get vertex indices for triangle t
                // Create triangle from 3 vertices and add to pattern collection
                checkerboard_facets.push(cell.vertices[tri[0]], cell.vertices[tri[1]], cell.vertices[tri[2]]);
            }
        }
        return checkerboard_facets;  // Return complete set of pattern faces
    }

    /**
     * @brief Get facets for all subcells in a plane using string coordinates.
     * @param coord1 First coordinate - either "x", "y", "z" or a numeric string like "0", "-1", "2"
     * @param coord2 Second coordinate - either "x", "y", "z" or a numeric string like "0", "-1", "2"
     * @param coord3 Third coordinate - either "x", "y", "z" or a numeric string like "0", "-1", "2"
     *
     * Examples:
     * - getPlaneFacets("x", "y", "0") gets center XY plane (perpendicular to Z axis)
     * - getPlaneFacets("x", "0", "z") gets center XZ plane (perpendicular to Y axis)
     * - getPlaneFacets("0", "y", "z") gets center YZ plane (perpendicular to X axis)
     *
     * The numeric coordinate must be in range [-n/2, n/2] for n×n×n subdivision.
     */
    FacetBox getPlaneFacets(const std::string& coord1, const std::string& coord2, const std::string& coord3) const {
        if (!hasSubdivision()) {
            throw std::runtime_error("Cube::getPlaneFacets: no subdivision available");
        }

        // Helper lambda to check if string is numeric
        auto isNumeric = [](const std::string& str) -> bool {
            if (str.empty()) return false;
            size_t start = (str[0] == '-') ? 1 : 0;
            return start < str.size() && std::all_of(str.begin() + start, str.end(), ::isdigit);
        };

        // Helper lambda to convert string to int with validation
        auto parseCoord = [this](const std::string& str) -> int {
            int val = std::stoi(str);
            int max_coord = subdivision_levels_ / 2;
            int min_coord = -max_coord;
            if (val < min_coord || val > max_coord) {
                throw std::out_of_range("Coordinate " + str + " out of range [" +
                                      std::to_string(min_coord) + ", " + std::to_string(max_coord) + "]");
            }
            return val;
        };

        std::vector<std::string> coords = {coord1, coord2, coord3};
        int numeric_count = 0;
        int numeric_idx = -1;
        int numeric_value = 0;

        // Find which coordinate is numeric
        for (int i = 0; i < 3; ++i) {
            if (isNumeric(coords[i])) {
                numeric_count++;
                numeric_idx = i;
                numeric_value = parseCoord(coords[i]);
            }
        }

        if (numeric_count != 1) {
            throw std::invalid_argument("Exactly one coordinate must be numeric, got " + std::to_string(numeric_count));
        }

        // Determine axis and layer
        int axis;
        if (numeric_idx == 0) {
            // X is fixed, varying Y and Z -> YZ plane (axis = 0)
            if (coords[1] != "y" || coords[2] != "z") {
                throw std::invalid_argument("For fixed X coordinate, other coordinates must be 'y' and 'z'");
            }
            axis = 0;
        } else if (numeric_idx == 1) {
            // Y is fixed, varying X and Z -> XZ plane (axis = 1)
            if (coords[0] != "x" || coords[2] != "z") {
                throw std::invalid_argument("For fixed Y coordinate, other coordinates must be 'x' and 'z'");
            }
            axis = 1;
        } else {
            // Z is fixed, varying X and Y -> XY plane (axis = 2)
            if (coords[0] != "x" || coords[1] != "y") {
                throw std::invalid_argument("For fixed Z coordinate, other coordinates must be 'x' and 'y'");
            }
            axis = 2;
        }

        // Convert logical coordinate to physical layer
        int center = subdivision_levels_ / 2;
        int layer = center + numeric_value;

        if (layer < 0 || layer >= subdivision_levels_) {
            throw std::out_of_range("Calculated layer " + std::to_string(layer) + " out of range [0, " +
                                  std::to_string(subdivision_levels_) + ")");
        }

        return getPlaneFacets(axis, layer);
    }

    /**
     * @brief Get the center of a specific subcell.
     * @param x Logical X-coordinate (0=center)
     * @param y Logical Y-coordinate (0=center)
     * @param z Logical Z-coordinate (0=center)
     * @returns Vector3D representing the exact center of the subcell
     * @throws std::out_of_range if coordinates invalid or no subdivision
     */
    // Get the geometric center point of a specific subcell
    Vector3D getSubcellCenter(int x, int y, int z) const {
        return getSubCell(x, y, z).center;  // Return center from subcell data structure
    }

    /**
     * @brief Get the radius of a specific subcell.
     * The radius is the distance from the center to any corner of the subcell.
     * @param x Logical X-coordinate (0=center)
     * @param y Logical Y-coordinate (0=center)
     * @param z Logical Z-coordinate (0=center)
     * @returns double representing the radius (distance from center to corner)
     * @throws std::out_of_range if coordinates invalid or no subdivision
     */
    // Get the radius (distance from center to corner) of a specific subcell
    double getSubcellRadius(int x, int y, int z) const {
        const SubCell& cell = getSubCell(x, y, z);  // Get reference to target subcell
        // The radius stored is half the side length, but for distance to corner
        // we need sqrt(3) * radius for a cube (3D diagonal calculation)
        return cell.radius * sqrt(3.0);  // Convert from half-side-length to corner distance
    }

    /**
     * @brief Enable/disable a specific subcell for rendering.
     * @param x,y,z Logical coordinates (centered at origin)
     */
    // Enable or disable rendering/processing for a specific subcell
    void setSubCellActive(int x, int y, int z, bool active) {
        getSubCellMutable(x, y, z).active = active;  // Set the active flag in subcell
    }

private:
    // === CORE DATA MEMBERS ===
    std::vector<Quaternion> verts_;  // 8 corner vertices of the main cube (stored as quaternions)
    FacetBox facets_;               // Collection of all triangular faces for rendering

    // === SUBDIVISION DATA MEMBERS ===
    // 3D Grid of subcells: subcells_[i][j][k] - indexed by physical coordinates
    std::vector<std::vector<std::vector<SubCell>>> subcells_;
    int subdivision_levels_;   // Current subdivision level (n means n×n×n subcells total)
    bool subdivision_built_;   // Flag tracking if subdivision data has been generated

    // === TRIANGULATION LOOKUP TABLE ===
    // Static lookup table defining vertex indices for the 12 triangular faces of a cube
    // Each cube face is divided into 2 triangles, resulting in 12 total triangles
    static constexpr int cube_triangles_[12][3] = {
        // Front face (z = +radius) - 2 triangles from vertices 0,1,2,3
        {0, 1, 2}, {2, 3, 0},  // Triangle 1: (0,1,2), Triangle 2: (2,3,0)
        // Back face (z = -radius) - 2 triangles from vertices 4,5,6,7
        {4, 6, 5}, {6, 4, 7},  // Triangle 3: (4,6,5), Triangle 4: (6,4,7)
        // Left face (x = -radius) - 2 triangles from vertices 0,3,4,7
        {4, 0, 3}, {3, 7, 4},  // Triangle 5: (4,0,3), Triangle 6: (3,7,4)
        // Right face (x = +radius) - 2 triangles from vertices 1,2,5,6
        {1, 5, 6}, {6, 2, 1},  // Triangle 7: (1,5,6), Triangle 8: (6,2,1)
        // Bottom face (y = -radius) - 2 triangles from vertices 0,1,4,5
        {4, 5, 1}, {1, 0, 4},  // Triangle 9: (4,5,1), Triangle 10: (1,0,4)
        // Top face (y = +radius) - 2 triangles from vertices 2,3,6,7
        {3, 2, 6}, {6, 7, 3}   // Triangle 11: (3,2,6), Triangle 12: (6,7,3)
    };


    /**
     * @brief Initialize the 8 vertices of a cube.
     * Creates the fundamental 8 corner points that define the cube geometry
     */
    void initVertices(double r, const Vector3D& center) {
        verts_.clear();     // Remove any existing vertices
        verts_.reserve(8);  // Optimize memory allocation for 8 vertices

        // Generate 8 vertices of cube: all combinations of ±r around center
        // Vertex ordering follows standard cube convention for consistent triangulation
        verts_.emplace_back(0.0, Vector3D(center.x() - r, center.y() - r, center.z() + r)); // 0: (-x, -y, +z)
        verts_.emplace_back(0.0, Vector3D(center.x() + r, center.y() - r, center.z() + r)); // 1: (+x, -y, +z)
        verts_.emplace_back(0.0, Vector3D(center.x() + r, center.y() + r, center.z() + r)); // 2: (+x, +y, +z)
        verts_.emplace_back(0.0, Vector3D(center.x() - r, center.y() + r, center.z() + r)); // 3: (-x, +y, +z)
        verts_.emplace_back(0.0, Vector3D(center.x() - r, center.y() - r, center.z() - r)); // 4: (-x, -y, -z)
        verts_.emplace_back(0.0, Vector3D(center.x() + r, center.y() - r, center.z() - r)); // 5: (+x, -y, -z)
        verts_.emplace_back(0.0, Vector3D(center.x() + r, center.y() + r, center.z() - r)); // 6: (+x, +y, -z)
        verts_.emplace_back(0.0, Vector3D(center.x() - r, center.y() + r, center.z() - r)); // 7: (-x, +y, -z)
    }

    /**
     * @brief Build basic cube facets (12 triangles).
     * Converts the 8 vertices into 12 triangular faces for rendering
     */
    void buildFacets() {
        facets_.clear();  // Remove any existing triangular faces

        // Generate 12 triangles using the static triangulation lookup table
        for (int i = 0; i < 12; ++i) {  // Process each of the 12 triangular faces
            auto const& t = cube_triangles_[i];  // Get vertex indices for triangle i
            // Create triangle from 3 vertices referenced by lookup table
            Facet face(verts_[t[0]], verts_[t[1]], verts_[t[2]]);
            facets_.push(face);  // Add completed triangle to facet collection
        }
    }

    /**
     * @brief Build subdivided cube with n subdivision levels and populate subcells data.
     * Core subdivision algorithm that creates n×n×n smaller cubes from the original
     */
    void buildSubdividedCube(const Vector3D& center, double radius, int n) {
        double step = (2.0 * radius) / n; // Size of each subcube (full cube / n subdivisions)
        double subRadius = step / 2.0;     // Radius of each subcube (half its side length)

        // Initialize 3D subcells grid - creates n×n×n empty subcells
        subcells_.clear();  // Remove any existing subdivision data
        subcells_.resize(n, std::vector<std::vector<SubCell>>(n, std::vector<SubCell>(n)));

        // Generate n³ subcubes and store them in the grid
        for (int i = 0; i < n; ++i) {      // X dimension subdivision
            for (int j = 0; j < n; ++j) {  // Y dimension subdivision
                for (int k = 0; k < n; ++k) {  // Z dimension subdivision
                    // Calculate subcube center using grid coordinates
                    Vector3D subCenter(
                        center.x() - radius + step * (i + 0.5),  // X: start at left edge, offset by grid position
                        center.y() - radius + step * (j + 0.5),  // Y: start at bottom edge, offset by grid position
                        center.z() - radius + step * (k + 0.5)   // Z: start at back edge, offset by grid position
                    );

                    // Create 8 vertices for this subcube (following standard cube vertex ordering)
                    std::array<Vector3D, 8> subVerts = {
                        Vector3D(subCenter.x() - subRadius, subCenter.y() - subRadius, subCenter.z() + subRadius), // 0: (-x,-y,+z)
                        Vector3D(subCenter.x() + subRadius, subCenter.y() - subRadius, subCenter.z() + subRadius), // 1: (+x,-y,+z)
                        Vector3D(subCenter.x() + subRadius, subCenter.y() + subRadius, subCenter.z() + subRadius), // 2: (+x,+y,+z)
                        Vector3D(subCenter.x() - subRadius, subCenter.y() + subRadius, subCenter.z() + subRadius), // 3: (-x,+y,+z)
                        Vector3D(subCenter.x() - subRadius, subCenter.y() - subRadius, subCenter.z() - subRadius), // 4: (-x,-y,-z)
                        Vector3D(subCenter.x() + subRadius, subCenter.y() - subRadius, subCenter.z() - subRadius), // 5: (+x,-y,-z)
                        Vector3D(subCenter.x() + subRadius, subCenter.y() + subRadius, subCenter.z() - subRadius), // 6: (+x,+y,-z)
                        Vector3D(subCenter.x() - subRadius, subCenter.y() + subRadius, subCenter.z() - subRadius)  // 7: (-x,+y,-z)
                    };

                    // Store complete subcell data in the 3D grid
                    subcells_[i][j][k] = SubCell(subVerts, subCenter, subRadius, i, j, k);

                    // Generate 12 triangles for immediate rendering using the standard triangulation
                    for (int t = 0; t < 12; ++t) {  // Each subcube generates 12 triangular faces
                        auto const& tri = cube_triangles_[t];  // Get vertex indices for triangle t
                        // Create triangle from subcube vertices and add to main facet collection
                        facets_.push(subVerts[tri[0]], subVerts[tri[1]], subVerts[tri[2]]);
                    }
                }
            }
        }
    }

    /**
     * @brief Get the radius (half side length) of the cube.
     */
    double getRadius() const {
        if (verts_.size() < 2) return 0.0;
        Vector3D v0 = verts_[0].V();
        Vector3D v1 = verts_[1].V();
        return abs(v1.x() - v0.x()) / 2.0;
    }
};

//-----------------------------------------------------------------------------
// Stream output for Cube
//-----------------------------------------------------------------------------
inline std::ostream& operator<<(std::ostream& os, const Cube& cube) {
    os << "--- Cube Vertices ---\n";
    for (size_t i = 0; i < cube.verts_.size(); ++i) {
        os << "  [" << i << "] " << cube.verts_[i].V() << "\n";
    }
    os << "--- Cube Facets (" << cube.facets_.size() << " triangles) ---\n";
    for (size_t i = 0; i < cube.facets_.size(); ++i) {
        os << "  Face " << i << ": " << cube.facets_[i] << "\n";
    }
    return os;
}

#endif // CUBE_H
