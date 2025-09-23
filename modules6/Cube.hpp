#ifndef CUBE_H
#define CUBE_H

using namespace std;

#include <stdexcept>
#include <vector>
#include <array>
#include <functional>
#include <math.h>
#include "Vector3D.hpp"
#include "Quaternion.hpp"
#include "FacetBox.hpp"

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
    struct SubCell {
        std::array<Vector3D, 8> vertices;  // 8 cube vertices in standard order
        Vector3D center;                   // Center point of subcell
        double radius;                     // Half side length
        int i, j, k;                      // Grid coordinates
        bool active = true;               // Can be disabled for selective rendering

        SubCell() = default;
        SubCell(const std::array<Vector3D, 8>& verts, const Vector3D& c, double r, int x, int y, int z)
            : vertices(verts), center(c), radius(r), i(x), j(y), k(z) {}
    };

public:
    /* - Constructors */
    /* - Rule of Zero: default special members */
    Cube() : subdivision_levels_(0), subdivision_built_(false) {}
    ~Cube() = default;
    Cube(const Cube&) = default;
    Cube(Cube&&) noexcept = default;
    Cube& operator=(const Cube&) = default;
    Cube& operator=(Cube&&) noexcept = default;

    /**
     * @brief Construct a cube with given radius (half side length) centered at 'center'.
     * @param radius Distance from center to face (half of side length).
     * @param center Geometric center of the cube.
     */
    Cube(double radius, const Vector3D& center) : subdivision_levels_(0), subdivision_built_(false) {
        initVertices(radius, center);
        buildFacets();
    }

    /**
     * @brief Construct a subdivided cube with n subdivision levels.
     * @param radius Distance from center to face.
     * @param center Geometric center of the cube.
     * @param subdivisions Number of subdivision levels (0 = basic cube).
     */
    Cube(double radius, const Vector3D& center, int subdivisions) : subdivision_levels_(0), subdivision_built_(false) {
        initVertices(radius, center);
        buildFacets();
        if (subdivisions > 0) {
            subdivide(subdivisions);
        }
    }

    /* - Face Count */
    /**
     * @brief Number of triangular faces.
     */
    size_t faceCount() const noexcept { return facets_.size(); }

    /* - Facet Access */
    /**
     * @brief Access the k-th triangular face as a Facet.
     * @param k Index in range [0, faceCount()).
     * @throws std::out_of_range if k is invalid.
     */
    const Facet& operator[](size_t k) const {
        return facets_[k];
    }

    /* - Geometric Center */
    /**
     * @brief Compute and return the geometric center of all vertices.
     */
    Vector3D center() const {
        if (verts_.empty()) {
            throw std::runtime_error("Cube::center: no vertices initialized");
        }
        Vector3D sum{0,0,0};
        for (auto const& q : verts_) {
            sum += q.V();
        }
        return sum / static_cast<double>(verts_.size());
    }

    /* - Translate by offset */
    /**
     * @brief Translate the entire cube by 'offset'.
     */
    void translate(const Vector3D& offset) {
        for (auto& q : verts_) {
            Vector3D p = q.V() + offset;
            q = Quaternion(0.0, p);
        }
        buildFacets();
    }

    /* - Scale */
    /**
     * @brief Scale all vertices relative to a pivot point.
     * @param s Scaling factor.
     * @param pivot The point to scale around.
     */
    void scale(double s, const Vector3D& pivot) {
        for (auto& q : verts_) {
            Vector3D p = pivot + s * (q.V() - pivot);
            q = Quaternion(0.0, p);
        }
        buildFacets();
    }

    /**
     * @brief Retrieve the underlying FacetBox.
     */
    const FacetBox& getFacets() const noexcept {
        return facets_;
    }

    /**
     * @brief Subdivide the cube into n levels of smaller cubes, each triangulated.
     * @param n Number of subdivision levels.
     */
    void subdivide(int n) {
        if (n <= 0) return;

        // Clear current facets and rebuild with subdivision
        Vector3D cubeCenter = center();
        double cubeRadius = getRadius();

        // Store subdivision info
        subdivision_levels_ = n;

        // Rebuild with subdivision
        facets_.clear();
        buildSubdividedCube(cubeCenter, cubeRadius, n);
        subdivision_built_ = true;
    }

    // === SUBDIVISION ACCESS METHODS ===

    /**
     * @brief Check if subdivision data is available.
     */
    bool hasSubdivision() const noexcept {
        return subdivision_built_ && !subcells_.empty();
    }

    /**
     * @brief Get subdivision levels (n means n³ subcells).
     */
    int getSubdivisionLevels() const noexcept {
        return subdivision_levels_;
    }

    /**
     * @brief Access a specific subcell by 3D grid coordinates.
     * @param i X-coordinate in grid [0, n)
     * @param j Y-coordinate in grid [0, n)
     * @param k Z-coordinate in grid [0, n)
     * @throws std::out_of_range if coordinates invalid or no subdivision
     */
    const SubCell& getSubCell(int i, int j, int k) const {
        if (!hasSubdivision()) {
            throw std::runtime_error("Cube::getSubCell: no subdivision available");
        }
        if (i < 0 || i >= subdivision_levels_ ||
            j < 0 || j >= subdivision_levels_ ||
            k < 0 || k >= subdivision_levels_) {
            throw std::out_of_range("Cube::getSubCell: coordinates out of range");
        }
        return subcells_[i][j][k];
    }

    /**
     * @brief Get mutable access to a subcell for modifications.
     */
    SubCell& getSubCellMutable(int i, int j, int k) {
        if (!hasSubdivision()) {
            throw std::runtime_error("Cube::getSubCellMutable: no subdivision available");
        }
        if (i < 0 || i >= subdivision_levels_ ||
            j < 0 || j >= subdivision_levels_ ||
            k < 0 || k >= subdivision_levels_) {
            throw std::out_of_range("Cube::getSubCellMutable: coordinates out of range");
        }
        return subcells_[i][j][k];
    }

    /**
     * @brief Get all subcells in a specific plane.
     * @param axis 0=YZ plane (varying X), 1=XZ plane (varying Y), 2=XY plane (varying Z)
     * @param layer Layer index [0, n)
     * @returns Vector of subcell references
     */
    std::vector<std::reference_wrapper<const SubCell>> getPlane(int axis, int layer) const {
        if (!hasSubdivision()) {
            throw std::runtime_error("Cube::getPlane: no subdivision available");
        }
        if (axis < 0 || axis > 2 || layer < 0 || layer >= subdivision_levels_) {
            throw std::out_of_range("Cube::getPlane: invalid axis or layer");
        }

        std::vector<std::reference_wrapper<const SubCell>> plane;
        int n = subdivision_levels_;

        switch (axis) {
            case 0: // YZ plane (fixed X = layer)
                for (int j = 0; j < n; ++j) {
                    for (int k = 0; k < n; ++k) {
                        plane.emplace_back(subcells_[layer][j][k]);
                    }
                }
                break;
            case 1: // XZ plane (fixed Y = layer)
                for (int i = 0; i < n; ++i) {
                    for (int k = 0; k < n; ++k) {
                        plane.emplace_back(subcells_[i][layer][k]);
                    }
                }
                break;
            case 2: // XY plane (fixed Z = layer)
                for (int i = 0; i < n; ++i) {
                    for (int j = 0; j < n; ++j) {
                        plane.emplace_back(subcells_[i][j][layer]);
                    }
                }
                break;
        }
        return plane;
    }

    /**
     * @brief Update a specific vertex of a subcell and mark for retriangulation.
     * @param i,j,k Subcell coordinates
     * @param vertex_idx Vertex index [0,7] within the subcell
     * @param new_pos New position for the vertex
     */
    void updateSubCellVertex(int i, int j, int k, int vertex_idx, const Vector3D& new_pos) {
        if (vertex_idx < 0 || vertex_idx > 7) {
            throw std::out_of_range("Cube::updateSubCellVertex: vertex_idx must be [0,7]");
        }
        SubCell& cell = getSubCellMutable(i, j, k);
        cell.vertices[vertex_idx] = new_pos;

        // Mark as needing retriangulation
        // (will be handled in refreshTriangulation())
    }

    /**
     * @brief Rebuild all triangulation after vertex modifications.
     * Call this after updating vertices to refresh the visual mesh.
     */
    void refreshTriangulation() {
        if (!hasSubdivision()) return;

        facets_.clear();

        // Rebuild triangles from current subcell vertices
        for (int i = 0; i < subdivision_levels_; ++i) {
            for (int j = 0; j < subdivision_levels_; ++j) {
                for (int k = 0; k < subdivision_levels_; ++k) {
                    const SubCell& cell = subcells_[i][j][k];
                    if (!cell.active) continue;  // Skip disabled subcells

                    // Generate 12 triangles from current vertices
                    for (int t = 0; t < 12; ++t) {
                        auto const& tri = cube_triangles_[t];
                        facets_.push(cell.vertices[tri[0]], cell.vertices[tri[1]], cell.vertices[tri[2]]);
                    }
                }
            }
        }
    }

    /**
     * @brief Get facets for a specific subcell.
     * @returns FacetBox containing the 12 triangles of the specified subcell
     */
    FacetBox getSubCellFacets(int i, int j, int k) const {
        const SubCell& cell = getSubCell(i, j, k);
        FacetBox subcell_facets;

        if (!cell.active) return subcell_facets; // Return empty if disabled

        for (int t = 0; t < 12; ++t) {
            auto const& tri = cube_triangles_[t];
            subcell_facets.push(cell.vertices[tri[0]], cell.vertices[tri[1]], cell.vertices[tri[2]]);
        }
        return subcell_facets;
    }

    /**
     * @brief Get facets for all subcells in a plane.
     */
    FacetBox getPlaneFacets(int axis, int layer) const {
        auto plane_subcells = getPlane(axis, layer);
        FacetBox plane_facets;

        for (const SubCell& cell : plane_subcells) {
            if (!cell.active) continue;

            for (int t = 0; t < 12; ++t) {
                auto const& tri = cube_triangles_[t];
                plane_facets.push(cell.vertices[tri[0]], cell.vertices[tri[1]], cell.vertices[tri[2]]);
            }
        }
        return plane_facets;
    }

    /**
     * @brief Enable/disable a specific subcell for rendering.
     */
    void setSubCellActive(int i, int j, int k, bool active) {
        getSubCellMutable(i, j, k).active = active;
    }

private:
    std::vector<Quaternion> verts_;
    FacetBox facets_;

    // 3D Grid of subcells: subcells_[i][j][k]
    std::vector<std::vector<std::vector<SubCell>>> subcells_;
    int subdivision_levels_;  // Current subdivision level (n means n³ subcells)
    bool subdivision_built_;  // Track if subdivision data is available

    // Basic cube triangulation: 12 triangles (2 per face)
    static constexpr int cube_triangles_[12][3] = {
        // Front face (z = +radius)
        {0, 1, 2}, {2, 3, 0},
        // Back face (z = -radius)  
        {4, 6, 5}, {6, 4, 7},
        // Left face (x = -radius)
        {4, 0, 3}, {3, 7, 4},
        // Right face (x = +radius)
        {1, 5, 6}, {6, 2, 1},
        // Bottom face (y = -radius)
        {4, 5, 1}, {1, 0, 4},
        // Top face (y = +radius)
        {3, 2, 6}, {6, 7, 3}
    };


    /**
     * @brief Initialize the 8 vertices of a cube.
     */
    void initVertices(double r, const Vector3D& center) {
        verts_.clear();
        verts_.reserve(8);
        
        // 8 vertices of cube: all combinations of ±r around center
        verts_.emplace_back(0.0, Vector3D(center.x() - r, center.y() - r, center.z() + r)); // 0
        verts_.emplace_back(0.0, Vector3D(center.x() + r, center.y() - r, center.z() + r)); // 1
        verts_.emplace_back(0.0, Vector3D(center.x() + r, center.y() + r, center.z() + r)); // 2
        verts_.emplace_back(0.0, Vector3D(center.x() - r, center.y() + r, center.z() + r)); // 3
        verts_.emplace_back(0.0, Vector3D(center.x() - r, center.y() - r, center.z() - r)); // 4
        verts_.emplace_back(0.0, Vector3D(center.x() + r, center.y() - r, center.z() - r)); // 5
        verts_.emplace_back(0.0, Vector3D(center.x() + r, center.y() + r, center.z() - r)); // 6
        verts_.emplace_back(0.0, Vector3D(center.x() - r, center.y() + r, center.z() - r)); // 7
    }

    /**
     * @brief Build basic cube facets (12 triangles).
     */
    void buildFacets() {
        facets_.clear();
        for (int i = 0; i < 12; ++i) {
            auto const& t = cube_triangles_[i];
            Facet face(verts_[t[0]], verts_[t[1]], verts_[t[2]]);
            facets_.push(face);
        }
    }

    /**
     * @brief Build subdivided cube with n subdivision levels and populate subcells data.
     */
    void buildSubdividedCube(const Vector3D& center, double radius, int n) {
        double step = (2.0 * radius) / n; // Size of each subcube
        double subRadius = step / 2.0;     // Radius of each subcube

        // Initialize 3D subcells grid
        subcells_.clear();
        subcells_.resize(n, std::vector<std::vector<SubCell>>(n, std::vector<SubCell>(n)));

        // Generate n³ subcubes and store them in the grid
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    // Calculate subcube center
                    Vector3D subCenter(
                        center.x() - radius + step * (i + 0.5),
                        center.y() - radius + step * (j + 0.5),
                        center.z() - radius + step * (k + 0.5)
                    );

                    // Create vertices for this subcube (standard cube vertex order)
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

                    // Store subcell in the grid
                    subcells_[i][j][k] = SubCell(subVerts, subCenter, subRadius, i, j, k);

                    // Generate triangles for immediate rendering
                    for (int t = 0; t < 12; ++t) {
                        auto const& tri = cube_triangles_[t];
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