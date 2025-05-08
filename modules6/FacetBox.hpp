#ifndef FACETBOX_H
#define FACETBOX_H

using namespace std;

#include <vector>
#include <iostream>
#include "Facet.hpp"

/**
 *—————————————————————————————————————————————————————————————
 * @brief Container for managing a dynamic collection of Facet objects.
 * 
 * FacetBox leverages std::vector to store facets and provides
 * methods for adding, accessing, replacing, and transforming them.
 * It also computes the geometric center of all vertices.
 *—————————————————————————————————————————————————————————————
 *
 **/


 /*—————————————————————————————————————————————————————|*
 *| [A] _________[B]  [A] _________[B]   [A] _________[B]|        
 *|    \^========/       \^========/        \^========/  |    
 *|     \^=[N]==/         \^=[N]==/          \^=[N]==/   |    
 *|      \^====/           \^====/            \^====/    |        
 *|       \^==/             \^==/              \^==/     | ~ [FacetBox Class]
 *|        \^/               \^/                \^/      |        
 *|        [C]               [C]                [C]      |            
 **—————————————————————————————————————————————————————|*              
 */

class FacetBox {
public:
    /* - From Facet method ——————————————————————————————————————————————————*/
    /**
     * @brief Subdivide a single Facet about its centroid into three smaller facets.
     *
     * Given Facet f with vertices A,B,C and centroid D=f.getCenter(),
     * returns a FacetBox containing triangles (D,A,B), (D,B,C), and (D,C,A).
     */
    static FacetBox fromFacet(const Facet& f) {
        FacetBox box;
        Vector3D D = f.getCenter();
        Vector3D A = f[0];
        Vector3D B = f[1];
        Vector3D C = f[2];
        box.push(D, A, B);
        box.push(D, B, C);
        box.push(D, C, A);
        return box;
    }

    /* — From Facet — midpoint subdivision into 4 facets ————————————————————— */
    /**
     * FacetBox method that subdivides a given Facet into 4 parts by midpoints.
     */
    static FacetBox subdivide4(const Facet& f) {
        FacetBox box;
        Vector3D A =  f[0];
        Vector3D B =  f[1];
        Vector3D C =  f[2];
        Vector3D A0 = line(0.5, A, B);
        Vector3D B0 = line(0.5, B, C);
        Vector3D C0 = line(0.5, C, A);
        box.push(A, A0, C0);
        box.push(B, B0, A0);
        box.push(C, C0, B0);
        box.push(A0, B0, C0);
        return box;
    }

    /* — From Facet — midpoint subdivision into 6 facets ————————————————————— */
    /**
     * FacetBox method that subdivides a given Facet into 4 parts by midpoints.
     */
    static FacetBox subdivide6(const Facet& f) {
        FacetBox box;
        Vector3D D =  f.getCenter();
        Vector3D A =  f[0];
        Vector3D B =  f[1];
        Vector3D C =  f[2];
        Vector3D A0 = line(0.5, A, B);
        Vector3D B0 = line(0.5, B, C);
        Vector3D C0 = line(0.5, C, A);
        box.push(D,  A, A0);
        box.push(D, A0,  B);
        box.push(D,  B, B0);
        box.push(D, B0,  C);
        box.push(D,  C, C0);
        box.push(D, C0,  A);
        return box;
    }

    /* — From Facet — Sierpinski subdivision————————————————————————————————— */
    /**
     * FacetBox method that subdivides a given Facet into 4 parts by midpoints.
     */
    static FacetBox sierpinski(const Facet& f) {
        FacetBox box;
        Vector3D A = f[0];
        Vector3D B = f[1];
        Vector3D C = f[2];
        Vector3D A0 = line(0.5, A, B);
        Vector3D B0 = line(0.5, B, C);
        Vector3D C0 = line(0.5, C, A);
        box.push(A, A0, C0);
        box.push(B, B0, A0);
        box.push(C, C0, B0);
        //box.push(A0, B0, C0);
        return box;
    }

public:
    /* - Constructors ———————————————————————————————————————————————————————*/
    /* - Rule of Zero: let the compiler generate default special members-----*/
    FacetBox() = default;
    ~FacetBox() = default;
    FacetBox(const FacetBox&) = default;
    FacetBox(FacetBox&&) noexcept = default;
    FacetBox& operator=(const FacetBox&) = default;
    FacetBox& operator=(FacetBox&&) noexcept = default;

    /* - Basic constructor---------------------------------------------------*/
    /**
     * @brief Construct a FacetBox from a single facet
     * @param facet The facet to add
     */
    explicit FacetBox(const Facet& facet) {
        push(facet);
    }

    /* - Element Access —————————————————————————————————————————————————————*/
    /**
     * @brief Access a facet by index (bounds-checked)
     * @throws std::out_of_range
     */
    const Facet& operator[](size_t idx) const {
        return facets_.at(idx);
    }

    Facet& operator[](size_t idx) {
        return facets_.at(idx);
    }

    /* - Element Push ———————————————————————————————————————————————————————*/
    /**
     * @brief Add a facet to the collection
     * @param facet The facet to append
     */
    void push(const Facet& facet) {
        facets_.push_back(facet);
    }

    
    /**
     * @brief Construct and append a facet from three points
     */
    void push(const Vector3D& A, const Vector3D& B, const Vector3D& C) {
        facets_.emplace_back(A, B, C);
    }

    /* - FacetBox algebra ———————————————————————————————————————————————————*/
    /**
     * @brief Merge another FacetBox into this one (union of facets).
     */
    void merge(const FacetBox& other) {
        facets_.reserve(facets_.size() + other.facets_.size());
        for (auto const& f : other.facets_) {
            facets_.push_back(f);
        }
    }

    /**
     * @brief Combine operator: this = this ∪ other.
     */
    FacetBox& operator+=(const FacetBox& other) {
        merge(other);
        return *this;
    }


    /* - Element Replace ————————————————————————————————————————————————————*/
    /**
     * @brief Replace the facet at index idx
     * @throws std::out_of_range
     */
    void replace(size_t idx, const Facet& facet) {
        facets_.at(idx) = facet;
    }
    void replace(size_t idx,
                 const Vector3D& A,
                 const Vector3D& B,
                 const Vector3D& C)
    {
        facets_.at(idx) = Facet(A, B, C);
    }


    /* - Element Count ——————————————————————————————————————————————————————*/
    /**
     * @brief Number of facets stored
     */
    size_t size() const noexcept {
        return facets_.size();
    }


    /* - Element Delete —————————————————————————————————————————————————————*/
    /**
     * @brief Remove all facets
     */
    void clear() noexcept {
        facets_.clear();
    }

    /* - Recompute center of mass ———————————————————————————————————————————*/
    /**
     * @brief Compute the centroid of all vertices across all facets
     * @return The average position of every corner point
     */
    Vector3D center() const {
        Vector3D sum(0,0,0);
        size_t count = 0;
        for (auto const& f : facets_) {
            sum += f[0]; sum += f[1]; sum += f[2];
            count += 3;
        }
        return count ? sum / static_cast<double>(count)
                     : Vector3D(0,0,0);
    }

    /* - Translate Facets ———————————————————————————————————————————————————*/
    /**
     * @brief Translate every facet by an offset vector
     */
    void translate(const Vector3D& offset) {
        for (auto& f : facets_)
            f.translate(offset);
    }

    /* - Crunch Facets ——————————————————————————————————————————————————————*/
    /**
     * @brief Scale (crunch) every facet relative to pivot
     * @param t Scaling factor
     * @param pivot The point to scale around
     */
    void crunch(double t, const Vector3D& pivot) {
        for (auto& f : facets_)
            f.crunch(t, pivot);
    }

    /* Hyperbolic on entire mesh ———————————————————————————————————————————*/
    void applyHyperboloid() {
        for(auto& f: facets_) f.applyHyperboloid();
    }
    FacetBox hyperboloid() const {
        FacetBox out = *this;
        out.applyHyperboloid();
        return out;
    }

    /* — Refinement — subdivide all triangles n times around centroids ————————— */
    /**
     * @brief Choose a subdivision strategy for refining.
     */
    enum class SubdivisionMode { Centroid3, Midpoint4, Midpoint6, Sierpinski};

    /**
     * @brief Return a new FacetBox with each facet subdivided n times.
     * @param n Number of refinement iterations
     * @param mode Which subdivision pattern to use: Centroid3 or Midpoint4
     */

    FacetBox refine(int n, SubdivisionMode mode = SubdivisionMode::Centroid3) const {
        FacetBox curr = *this;
        FacetBox next;
        for (int pass = 0; pass < n; ++pass) {
            next.clear();
            for (auto const& f : curr.facets_) {
                switch(mode) {
                    case SubdivisionMode::Centroid3:
                        next += fromFacet(f);
                        break;
                    case SubdivisionMode::Midpoint4:
                        next += subdivide4(f);
                        break;
                    case SubdivisionMode::Midpoint6:
                        next += subdivide6(f);
                        break;
                    case SubdivisionMode::Sierpinski:
                        next += sierpinski(f);
                        break;
                }
            }
            std::swap(curr, next);
        }
        return curr;
    }

    /* - Example usage */
    // mesh3 = box.refine(2, FacetBox::SubdivisionMode::Centroid3);
    // mesh4 = box.refine(2, FacetBox::SubdivisionMode::Midpoint4);


private:
    std::vector<Facet> facets_;  ///< underlying storage

};


/**
 * @brief Stream output for debugging / serialization
 */
inline std::ostream& operator<<(std::ostream& os, const FacetBox& fb) {
    os << "FacetBox(size=" << fb.size() << ")\n";
    for (size_t i = 0; i < fb.size(); ++i) {
        os << "  [" << i << "] " << fb[i] << "\n";
    }
    return os;
}

/**
 * @brief Free-function operator to union two FacetBoxes.
 */
inline FacetBox operator+(FacetBox lhs, const FacetBox& rhs) {
    lhs += rhs;
    return lhs;
}

#endif // FACETBOX_H
