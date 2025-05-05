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
 * Stores 20 vertices (as Quaternion points) and a FacetBox of 36 triangular faces.
 * Uses Rule-of-Zero for special members and std::array/std::vector internally. 
 *  —————————————————————————————————————————————————————————————
 */

/**
      _----------_,
    ,"__         _-:,
   /    ""--_--""...:\
  /         |.........\
 /          |..........\
/,         _'_........./.
! -,    _-"   "-_... ,;;:
\   -_-"         "-_/;;;.
 \   \             /;;;.
  \   \           /;;;.                                                                                                                                                                       
   '.  \         /;;;'
     "-_\_______/;;'        ~[Dodecahedron Class]

**/

class Dodecahedron {
public:
    // allow operator<< to access private members
    friend std::ostream& operator<<(std::ostream& os, const Dodecahedron& dd);

public:
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
     */
    Dodecahedron(double radius, const Vector3D& center) {
        initVertices(radius, center);
        buildFacets();
    }

    /* - Face Count —————————————————————————————————————————————————————————*/
    /**
     * @brief Number of triangular faces (12 faces × 3 triangles per pentagon = 36).
     */
    size_t faceCount() const noexcept { return facets_.size(); }

    /* - Facet Access ———————————————————————————————————————————————————————*/
    /**
     * @brief Access the k-th triangular face as a Facet.
     * @param k Index in range [0, faceCount()).
     * @throws std::out_of_range if k is invalid.
     */
    const Facet& operator[](size_t k) const {
        return facets_[k];
    }

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
            cout << "--> " << q.V();
            sum += q.V();
            cout << "   sum = " << sum << "\n";
        }
        return sum / 20.0;
    }

    /* - Translate by offset ————————————————————————————————————————————————*/
    /**
     * @brief Translate the entire dodecahedron by 'offset'.
     */
    void translate(const Vector3D& offset) {
        for (auto& q : verts_) {
            Vector3D p = q.V() + offset;
            q = Quaternion(0.0, p);
        }
        buildFacets();
    }

    /* - Scale ——————————————————————————————————————————————————————————————*/
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


private:
    std::vector<Quaternion> verts_; ///< 20 dodecahedron vertices
    FacetBox facets_;                ///< 36 triangular facets

    /* - Mesh data ——————————————————————————————————————————————————————————*/
    // Indices of vertices forming each of the 36 triangles
    static constexpr int tri_[36][3] = {
        {19,14, 3}, {19,13,14}, {19, 4, 3},
        {19,13, 8}, {19, 8,12}, {19,12,16},
        { 4,19,16}, {16, 0, 4}, { 0,16, 6},
        { 0, 2, 3}, { 0, 1, 2}, { 0, 3, 4},
        { 2,17,15}, { 2,15,14}, { 2,14, 3},
        { 8,13,14}, { 8,14,15}, { 8,15, 9},
        { 5, 6,12}, {11, 5,12}, { 6,16,12},
        { 8, 9,12}, { 9,10,12}, {10,11,12},
        { 9,15,10}, {15,17,10}, {17,18,10},
        {18, 7,10}, { 7, 5,10}, { 5,11,10},
        {18,17, 2}, { 7,18, 2}, { 1, 7, 2},
        { 0, 7, 1}, { 0, 5, 7}, { 0, 6, 5}
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
     * @brief Rebuilds the FacetBox from current vertices.
     */
    void buildFacets() {
        facets_.clear();
        for (int i = 0; i < 36; ++i) {
            auto const& t = tri_[i];
                    // Construct each triangular face explicitly
        Facet face( verts_[t[0]], verts_[t[1]], verts_[t[2]] );
        facets_.push(face);
        }
    }
};



//-----------------------------------------------------------------------------
// Stream output for Dodecahedron
//-----------------------------------------------------------------------------
/**
 * @brief Print all vertices and facets of a dodecahedron.
 */
inline std::ostream& operator<<(std::ostream& os, const Dodecahedron& dd) {
    // Output vertices
    os << "--- Dodecahedron Vertices ---\n";
    for (size_t i = 0; i < dd.verts_.size(); ++i) {
        os << "  [" << i << "] " << dd.verts_[i].V() << "\n";
    }
    // Output facets
    os << "--- Dodecahedron Facets (" << dd.facets_.size() << " triangles) ---\n";
    for (size_t i = 0; i < dd.facets_.size(); ++i) {
        os << "  Face " << i << ": " << dd.facets_[i] << "\n";
    }
    return os;
}

//class Torus {
//
//	private:
//		QuaternionBox q;
//		int n;
//		FacetBox f;
//	public:
//		Torus(double R, double r, const Vector3D& c, int N);
//		int getN() const;
//		Facet  operator [] (int k) const;
//		double G1(double u, double v, double R, double r);
//        	double G2(double u, double v, double R, double r);
//        	double G3(double u, double v, double R, double r);
//		int getBoxSize() const;	
//};
//
//istream& operator >> (istream& is, Facet& a);
//ostream& operator << (ostream& os, const Facet& a);

#endif // DODECAHEDRON_H
