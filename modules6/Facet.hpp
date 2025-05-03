#ifndef FACET_H
#define FACET_H

using namespace std;

#include <iostream>
#include "Quaternion.hpp"
#include <math.h>

/*—————————————————————————————————————————————————————————————
 * Facet represents a triangle defined by three vertices (as Quaternions)
 * and its associated normal vector.
 *—————————————————————————————————————————————————————————————*/
class Facet {

/*
 * [A] _________[B]
 *    \^========/
 *     \^=[N]==/
 *      \^====/
 *       \^==/ [Facet Class]
 *        \^/
 *        [C]
 */

    private:
        Quaternion A{}; // Vertex A
        Quaternion B{}; // Vertex B
        Quaternion C{}; // Vertex C
        Quaternion N{}; // Normal (computed as cross(B - A, C - A))

    public:
    /*————— element access ——————————————————————————————————————————————————*/
    /*-Element access: 0 → A, 1 → B, 2 → C ----------------------------------*/
    Vector3D  operator [] (int k) const;
    
    /*Accessor for the normal -----------------------------------------------*/
    Vector3D getNormal() const { return N.V(); }

    /*———— constructors —————————————————————————————————————————————————————*/
     // Rule of Zero: default operations are sufficient
    Facet() = default;
    ~Facet() = default;
    Facet(const Facet&) = default;
    Facet(Facet&&) noexcept = default;
    Facet& operator=(const Facet&) = default;
    Facet& operator=(Facet&&) noexcept = default;
    
    /* Constructor from three Quaternions -----------------------------------*/
    Facet(const Quaternion& a, const Quaternion& b, const Quaternion& c);

    /* Constructor from three Vector3D points -------------------------------*/
    Facet(const Vector3D& a, const Vector3D& b, const Vector3D& c);

    /*———— update method ————————————————————————————————————————————————————*/
    void updateFacet(const Vector3D& a, const Vector3D& b, const Vector3D& c);
    
    /*———— geometry helpers ——————————————————————————————————————————————————*/
    /* Get the center (centroid) of the triangle -----------------------------*/
    Vector3D getCenter() const;
    
    /* Translate the facet by a vector ---------------------------------------*/
    void translate(const Vector3D& offset);

    /* Scale (crunch) the triangle with respect to a point -------------------*/
    void crunch(double t, const Vector3D& a);

};

istream& operator >> (istream& is, Facet& a);
ostream& operator << (ostream& os, const Facet& a);

#endif
