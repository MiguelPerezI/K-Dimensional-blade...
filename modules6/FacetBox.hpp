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

#endif // FACETBOX_H
