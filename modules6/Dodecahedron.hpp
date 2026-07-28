#ifndef DODECAHEDRON_H
#define DODECAHEDRON_H

using namespace std;

#include <stdexcept>
#include <vector>
#include <array>    // for our temporary pts list
#include <math.h>
#include "Vector3D.hpp"
#include "Quaternion.hpp"
#include "FacetBox.hpp"

/**
 *  —————————————————————————————————————————————————————————————
 * @brief Represents a regular dodecahedron in 3D space.
 *
 * Stores 20 vertices (as Quaternion points) and a FacetBox of triangular faces.
 * If built with FaceMode::Pentagons, also stores the 12 pentagon face centers.
 * Uses Rule-of-Zero for special members and std::array/std::vector internally. 
 *  —————————————————————————————————————————————————————————————
 */

class Dodecahedron {
public:
    // allow operator<< to access private members
    friend std::ostream& operator<<(std::ostream& os, const Dodecahedron& dd);

public:
    /* - Face Mode Enum ————————————————————————————————————————————————————*/
    enum class FaceMode { Triangles, Pentagons };
    
    /* - Constructors ———————————————————————————————————————————————————————*/
    /* - Rule of Zero: default special members ------------------------------*/
    Dodecahedron() = default;
    ~Dodecahedron() = default;
    Dodecahedron(const Dodecahedron&) = default;
    Dodecahedron(Dodecahedron&&) noexcept = default;
    Dodecahedron& operator=(const Dodecahedron&) = default;
    Dodecahedron& operator=(Dodecahedron&&) noexcept = default;

    /**
     * @brief Construct a dodecahedron of given radius centered at 'center'.
     * @param radius Distance from center to each vertex.
     * @param center Geometric center of the dodecahedron.
     * Triangular version (36 facets)
     */
    Dodecahedron(double radius, const Vector3D& center) : faceMode_(FaceMode::Triangles) {
        initVertices(radius, center);
        buildFacets();
    }

    /**
     * @brief Construct a dodecahedron of given radius centered at 'center'.
     * @param radius Distance from center to each vertex.
     * @param center Geometric center of the dodecahedron.
     * @param mode FaceMode::Triangles (36 faces) or FaceMode::Pentagons (60 faces)
     */
    Dodecahedron(double radius, const Vector3D& center, FaceMode mode) : faceMode_(mode) {
        initVertices(radius, center);
        if (mode == FaceMode::Pentagons) buildFacetsPenta();
        else                              buildFacets();
    }

    /* - Face Count —————————————————————————————————————————————————————————*/
    /**
     * @brief Number of triangular faces (36 for Triangles mode, 60 for Pentagons mode).
     */
    size_t faceCount() const noexcept { return facets_.size(); }

    /**
     * @brief Number of pentagon face centers (12 if FaceMode::Pentagons, 0 otherwise).
     */
    size_t pentagonCenterCount() const noexcept { return pentagonCenters_.size(); }

    /* - Facet Access ———————————————————————————————————————————————————————*/
    /**
     * @brief Access the k-th triangular face as a Facet.
     * @param k Index in range [0, faceCount()).
     * @throws std::out_of_range if k is invalid.
     */
    const Facet& operator[](size_t k) const {
        return facets_[k];
    }

    /* - Pentagon Center Access —————————————————————————————————————————————*/
    /**
     * @brief Access the k-th pentagon face center.
     * @param k Index in range [0, pentagonCenterCount()).
     * @throws std::out_of_range if k is invalid or not in Pentagon mode.
     */
    const Vector3D& getPentagonCenter(size_t k) const {
        if (faceMode_ != FaceMode::Pentagons) {
            throw std::runtime_error("Pentagon centers only available in FaceMode::Pentagons");
        }
        if (k >= pentagonCenters_.size()) {
            throw std::out_of_range("Pentagon center index out of range");
        }
        return pentagonCenters_[k];
    }

    /**
     * @brief Get all pentagon centers (empty if not in Pentagon mode).
     */
    const std::vector<Vector3D>& getPentagonCenters() const noexcept {
        return pentagonCenters_;
    }

    /**
     * @brief Get current face mode.
     */
    FaceMode getFaceMode() const noexcept { return faceMode_; }

    /* - Geometric Center ———————————————————————————————————————————————————*/
    /**
     * @brief Compute and return the geometric center of all vertices.
     */
    Vector3D center() const {
        if (verts_.empty()) {
            throw std::runtime_error("Dodecahedron::center: no vertices initialized");
        }
        Vector3D sum{0,0,0};
        for (auto const& q : verts_) {
            sum += q.V();
        }
        return sum / 20.0;
    }

    /* - Translate by offset ————————————————————————————————————————————————*/
    /**
     * @brief Translate the entire dodecahedron by 'offset'.
     */
    void translate(const Vector3D& offset) {
        // Translate vertices
        for (auto& q : verts_) {
            Vector3D p = q.V() + offset;
            q = Quaternion(0.0, p);
        }
        
        // Translate pentagon centers if they exist
        for (auto& center : pentagonCenters_) {
            center += offset;
        }
        
        // Rebuild facets
        if (faceMode_ == FaceMode::Pentagons) buildFacetsPenta();
        else                                  buildFacets();
    }

    /* - Scale ——————————————————————————————————————————————————————————————*/
    /**
     * @brief Scale all vertices relative to a pivot point.
     * @param s Scaling factor.
     * @param pivot The point to scale around.
     */
    void scale(double s, const Vector3D& pivot) {
        // Scale vertices
        for (auto& q : verts_) {
            Vector3D p = pivot + s * (q.V() - pivot);
            q = Quaternion(0.0, p);
        }
        
        // Scale pentagon centers if they exist
        for (auto& center : pentagonCenters_) {
            center = pivot + s * (center - pivot);
        }
        
        // Rebuild facets
        if (faceMode_ == FaceMode::Pentagons) buildFacetsPenta();
        else                                  buildFacets();
    }

    /**
     * @brief Rotate the dodecahedron around an axis through a point.
     * @param axis The rotation axis (will be normalized).
     * @param angle Rotation angle in radians.
     * @param pivot Point on the rotation axis.
     */
    void rotate(const Vector3D& axis, double angle, const Vector3D& pivot) {
        // Create rotation quaternion
        Vector3D normalizedAxis = unit(axis);  // Instead of axis.normalized()
        Quaternion rotQuat(std::cos(angle/2.0), std::sin(angle/2.0) * normalizedAxis);
        
        // Rotate vertices
        for (auto& q : verts_) {
            Vector3D p = q.V() - pivot;
            Quaternion pQuat(0.0, p);
            Quaternion rotated = rotQuat * pQuat * rotQuat.conjugate();
            q = Quaternion(0.0, rotated.V() + pivot);
        }
        
        // Rotate pentagon centers if they exist
        for (auto& center : pentagonCenters_) {
            Vector3D p = center - pivot;
            Quaternion pQuat(0.0, p);
            Quaternion rotated = rotQuat * pQuat * rotQuat.conjugate();
            center = rotated.V() + pivot;
        }
        
        // Rebuild facets
        if (faceMode_ == FaceMode::Pentagons) buildFacetsPenta();
        else                                  buildFacets();
    }

    /**
     * @brief Retrieve the underlying FacetBox of triangular facets.
     */
    const FacetBox& getFacets() const noexcept {
        return facets_;
    }

    /**
     * @brief Build a fresh FacetBox as the solid dodecahedron OR a hollow frame.
     *
     * When hollow=false this returns the existing facets_ (the 36-triangle solid
     * mesh). When hollow=true each of the 12 pentagonal faces (penta_) becomes an
     * inset border ring (120 triangles): inner[e] = C + inset*(o[e]-C) about the
     * face centroid C, 2 triangles per edge, inner face skipped. `inset` is the
     * border-thickness ratio in (0,1); only used when hollow. Same homothety +
     * 2-tris/edge pattern as Cube's frame pushers, so the winding is outward under
     * GL_CCW (penta_ is the same CCW-outward order buildFacetsPenta renders).
     */
    FacetBox getFacetsFrame(bool hollow = false, double inset = 0.5) const {
        if (!hollow) return facets_;                 // solid: existing 36 tris
        if (!(inset > 0.0)) inset = 0.0001;
        if (inset >= 1.0)   inset = 0.9999;
        FacetBox fb;
        for (auto const& idxs : penta_) {
            auto [i0, i1, i2, i3, i4] = idxs;
            Vector3D o[5] = { verts_[i0].V(), verts_[i1].V(), verts_[i2].V(),
                             verts_[i3].V(), verts_[i4].V() };
            Vector3D C = (o[0] + o[1] + o[2] + o[3] + o[4]) / 5.0;   // face centroid
            Vector3D in[5];
            for (int e = 0; e < 5; ++e) in[e] = C + inset * (o[e] - C);
            for (int e = 0; e < 5; ++e) {
                int j = (e + 1) % 5;
                fb.push(o[e], o[j], in[j]);
                fb.push(o[e], in[j], in[e]);
            }
        }   // 12 pentagons * 5 edges * 2 = 120 tris
        return fb;
    }

/**
 * @brief Apply sigma (sphere inversion) transformation to all components of the dodecahedron.
 * Transforms all vertices and pentagon centers, then rebuilds the entire structure.
 * @param center Center of the inversion sphere.
 * @param radius Radius of the inversion sphere.
 */
void sigma(const Vector3D& center, double radius) {
    // Apply sigma to all vertices
    for (auto& vertex : verts_) {
        Vector3D transformedPoint = ::sigma(vertex.V(), center, radius);
        vertex = Quaternion(0.0, transformedPoint);
    }
    
    // Apply sigma to pentagon centers if they exist
    for (auto& pentCenter : pentagonCenters_) {
        pentCenter = ::sigma(pentCenter, center, radius);
    }
    
    // Completely rebuild the dodecahedron structure
    //if (faceMode_ == FaceMode::Pentagons) {
    //    buildFacetsPenta();
    //} else {
    //    buildFacets();
    //}
}

private:
    std::vector<Quaternion> verts_;      ///< 20 dodecahedron vertices
    FacetBox facets_;                    ///< Triangular facets (36 or 60 depending on mode)
    std::vector<Vector3D> pentagonCenters_; ///< 12 pentagon face centers (only for Pentagon mode)
    FaceMode faceMode_ = FaceMode::Triangles; ///< Current face construction mode

    /* - Mesh data ——————————————————————————————————————————————————————————*/
    // Indices of vertices forming each of the 36 triangles
    static constexpr int tri_[36][3] = {
        {9, 15, 14},    {8, 9, 14},     {13, 8, 14},
        {1, 7, 5},      {0, 1, 5},      {6, 0, 5},
        {12, 8, 13},    {16, 12, 13},   {19, 16, 13},
        {17, 18, 7},    {2, 17, 7},     {1, 2, 7},
        {18, 17, 15},   {10, 18, 15},   {9, 10, 15},
        {4, 0, 6},      {19, 4, 6},     {16, 19, 6},
        {11, 10, 9},    {12, 11, 9},    {8, 12, 9},
        {2, 1, 0},      {3, 2, 0},      {4, 3, 0},
        {11, 12, 16},   {5, 11, 16},    {6, 5, 16},
        {14, 15, 17},   {3, 14, 17},    {2, 3, 17},
        {7, 18, 10},    {5, 7, 10},     {11, 5, 10},
        {14, 3, 4},     {13, 14, 4},    {19, 13, 4}
    };

    // Indices of vertices forming each of the 12 pentagons                                                                                                                             
    static constexpr int penta_[12][5] = {
        { 7, 18, 10, 11,  5}, {14, 15, 17,  2,  3},
        {14,  3,  4, 19, 13}, {11, 12, 16,  6,  5},
        { 9, 15, 14, 13,  8}, { 1,  7,  5,  6,  0},
        {12,  8, 13, 19, 16}, {17, 18,  7,  1,  2},
        {17, 15,  9, 10, 18}, { 0,  6, 16, 19,  4},
        {11, 10,  9,  8, 12}, { 1,  0,  4,  3,  2}
    };

    /* - Initialize vertices ————————————————————————————————————————————————*/
    /**
     * @brief Compute the 20 vertex positions for a regular dodecahedron.
     */
    void initVertices(double r, const Vector3D& center) {
        verts_.clear();
        verts_.reserve(20);
        double gold = 0.5 * (1.0 + std::sqrt(5.0));
        double g1 = 1.0 / gold;
        double g2 = 1.0 / (gold * gold);
        
        // explicit 20 vertices:
        verts_.emplace_back(0.0, Vector3D(center.x() +  g2*r, center.y() + 0.0*r, center.z() + 1.0*r));
        verts_.emplace_back(0.0, Vector3D(center.x() -  g2*r, center.y() + 0.0*r, center.z() + 1.0*r));
        verts_.emplace_back(0.0, Vector3D(center.x() -  g1*r, center.y() +  g1*r, center.z() +  g1*r));
        verts_.emplace_back(0.0, Vector3D(center.x() +  0.0*r, center.y() + 1.0*r, center.z() +  g2*r));
        verts_.emplace_back(0.0, Vector3D(center.x() +  g1*r, center.y() +  g1*r, center.z() +  g1*r));
        verts_.emplace_back(0.0, Vector3D(center.x() +  0.0*r, center.y() - 1.0*r, center.z() +  g2*r));
        verts_.emplace_back(0.0, Vector3D(center.x() +  g1*r, center.y() -  g1*r, center.z() +  g1*r));
        verts_.emplace_back(0.0, Vector3D(center.x() -  g1*r, center.y() -  g1*r, center.z() +  g1*r));
        verts_.emplace_back(0.0, Vector3D(center.x() +  g2*r, center.y() + 0.0*r, center.z() - 1.0*r));
        verts_.emplace_back(0.0, Vector3D(center.x() -  g2*r, center.y() + 0.0*r, center.z() - 1.0*r));
        verts_.emplace_back(0.0, Vector3D(center.x() -  g1*r, center.y() -  g1*r, center.z() -  g1*r));
        verts_.emplace_back(0.0, Vector3D(center.x() +  0.0*r, center.y() - 1.0*r, center.z() -  g2*r));
        verts_.emplace_back(0.0, Vector3D(center.x() +  g1*r, center.y() -  g1*r, center.z() -  g1*r));
        verts_.emplace_back(0.0, Vector3D(center.x() +  g1*r, center.y() +  g1*r, center.z() -  g1*r));
        verts_.emplace_back(0.0, Vector3D(center.x() +  0.0*r, center.y() + 1.0*r, center.z() -  g2*r));
        verts_.emplace_back(0.0, Vector3D(center.x() -  g1*r, center.y() +  g1*r, center.z() -  g1*r));
        verts_.emplace_back(0.0, Vector3D(center.x() + 1.0*r, center.y() -  g2*r, center.z() +  0.0*r));
        verts_.emplace_back(0.0, Vector3D(center.x() - 1.0*r, center.y() +  g2*r, center.z() +  0.0*r));
        verts_.emplace_back(0.0, Vector3D(center.x() - 1.0*r, center.y() -  g2*r, center.z() +  0.0*r));
        verts_.emplace_back(0.0, Vector3D(center.x() + 1.0*r, center.y() +  g2*r, center.z() +  0.0*r));
    }

    /* - Build Mesh —————————————————————————————————————————————————————————*/
    /**
     * @brief Rebuilds the FacetBox from current vertices (Triangle mode).
     */
    void buildFacets() {
        facets_.clear();
        pentagonCenters_.clear(); // No pentagon centers in triangle mode
        
        for (int i = 0; i < 36; ++i) {
            auto const& t = tri_[i];
            Facet face( verts_[t[0]], verts_[t[1]], verts_[t[2]] );
            facets_.push(face);
        }
    }

    /* - Build Mesh from Pentagons ——————————————————————————————————————————*/
    /**
     * @brief Rebuilds the FacetBox from current vertices using pentagon decomposition.
     */
    void buildFacetsPenta() {
        facets_.clear();
        pentagonCenters_.clear();
        pentagonCenters_.reserve(12);
    
        // For each of the 12 pentagon index‐tuples in penta_…
        for (auto const& idxs : penta_) {
            // Unpack the five vertex‐indices:
            auto [i0,i1,i2,i3,i4] = idxs;
    
            // Grab each point:
            const Vector3D P0 = verts_[i0].V();
            const Vector3D P1 = verts_[i1].V();
            const Vector3D P2 = verts_[i2].V();
            const Vector3D P3 = verts_[i3].V();
            const Vector3D P4 = verts_[i4].V();
    
            // Calculate the true geometric center of the pentagon
            Vector3D pentagonCenter = (P0 + P1 + P2 + P3 + P4) / 5.0;
            
            // Store the pentagon center
            pentagonCenters_.push_back(pentagonCenter);
    
            // Helper to push the triangle (centroid, A, B)
            auto pushEdge = [&](const Vector3D& A, const Vector3D& B) {
                facets_.push(pentagonCenter, A, B);
            };
    
            // Build the 5 triangular facets around the pentagon
            pushEdge(P0, P1);
            pushEdge(P1, P2);
            pushEdge(P2, P3);
            pushEdge(P3, P4);
            pushEdge(P4, P0);
        }
    }
};

//-----------------------------------------------------------------------------
// Stream output for Dodecahedron
//-----------------------------------------------------------------------------
/**
 * @brief Print all vertices, facets, and pentagon centers of a dodecahedron.
 */
inline std::ostream& operator<<(std::ostream& os, const Dodecahedron& dd) {
    // Output vertices
    os << "--- Dodecahedron Vertices ---\n";
    for (size_t i = 0; i < dd.verts_.size(); ++i) {
        os << "  [" << i << "] " << dd.verts_[i].V() << "\n";
    }
    
    // Output pentagon centers if available
    if (!dd.pentagonCenters_.empty()) {
        os << "--- Pentagon Face Centers ---\n";
        for (size_t i = 0; i < dd.pentagonCenters_.size(); ++i) {
            os << "  Pentagon " << i << " center: " << dd.pentagonCenters_[i] << "\n";
        }
    }
    
    // Output facets
    os << "--- Dodecahedron Facets (" << dd.facets_.size() << " triangles) ---\n";
    for (size_t i = 0; i < dd.facets_.size(); ++i) {
        os << "  Face " << i << ": " << dd.facets_[i] << "\n";
    }
    
    return os;
}

using FaceMode = Dodecahedron::FaceMode;

#endif // DODECAHEDRON_H
