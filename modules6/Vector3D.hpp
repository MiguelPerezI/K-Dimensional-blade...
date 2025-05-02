#ifndef VECTOR3D_H
#define VECTOR3D_H

using namespace std;

#include <iostream>
#include <cmath>
#include <stdexcept>

/*—————————————————————————————————————————————————————————————
 * Simple 3-component vector-v⃗ for geometry / linear-algebra work
 *—————————————————————————————————————————————————————————————*/
class Vector3D {
    protected:
        double u[4]; // x, y, z   (size 3 is enough)

    public:
        /*————— element access ——————————————————————————————————————————————————*/
        double& operator [] (int k);        // write-access
        double  operator [] (int k) const;  // read-access

        double x() const noexcept { return u[0]; };  // Access 'x' component
        double y() const noexcept { return u[1]; };  // Access 'y' component
        double z() const noexcept { return u[2]; };  // Access 'z' component


        /*———— constructors —————————————————————————————————————————————————————*/
        /* Construct a Vector3D by providing three real numbers -----------------*/        
        Vector3D(double xx = 0, double yy = 0, double zz = 0);
        
        /* Contruct a Vector3D by copying component data from another vector ----*/
        Vector3D(const Vector3D& a) = default;

        /* Construct a Vector3D by copying another vector using the '=' operator */
        Vector3D& operator =  (const Vector3D& a) = default;

        /*----- compound operator  ----------------------------------------------*/
        Vector3D& operator+=(const Vector3D& a) noexcept;
        Vector3D& operator-=(const Vector3D& rhs) noexcept;
        Vector3D& operator/=(double scalar);
        /*----- free-function operators (declared below) ------------------------*/
};

/* —————————————————————————————— free functions ————————————————————————————————*/
Vector3D  operator+(const Vector3D& a, const Vector3D& b) noexcept; // Addition
Vector3D  operator-(const Vector3D& a, const Vector3D& b) noexcept; // Subtraction
Vector3D  operator-(const Vector3D& a) noexcept;                    // unary minus
Vector3D  operator*(double s, const Vector3D& v) noexcept;          // scalar × vector
Vector3D  operator/(const Vector3D& v, double s);                   // scalar / vector 
double    operator*(const Vector3D& a, const Vector3D& b) noexcept; // dot product
Vector3D  operator%(const Vector3D& a, const Vector3D& b) noexcept; // cross product
bool      operator==(const Vector3D&, const Vector3D& ) noexcept; // Compare



double   abs   (const Vector3D& v) noexcept;                       // Euclidean length
double   infty (const Vector3D& v) noexcept;                       // ∞-norm
Vector3D unit  (const Vector3D& v);                       // normalized copy



/* helper for linear interpolation on a segment */
Vector3D line(double t, const Vector3D& p, const Vector3D& q) noexcept;

/* stream operators */
std::istream& operator>>(std::istream& is, Vector3D& v);
std::ostream& operator<<(std::ostream& os, const Vector3D& v);

/* center of 4 vectors*/
Vector3D centerM4(
                    const Vector3D& a0, const Vector3D& a1, 
                    const Vector3D& a2, const Vector3D& a3) noexcept;

/* center of 5 vectors*/
Vector3D cPenta(
                    Vector3D a0, Vector3D a1, 
                    Vector3D a2, Vector3D a3, Vector3D a4) noexcept;

/* center for 20 vectors*/
Vector3D cDodeca(Vector3D a[]) noexcept;

/* center for 8 vectors*/
Vector3D centerM8(
                  const Vector3D& a0, const Vector3D& a1, 
                  const Vector3D& a2, const Vector3D& a3,
                  const Vector3D& a4, const Vector3D& a5, 
                  const Vector3D& a6, const Vector3D& a7) noexcept;

/* cross product as function*/
Vector3D cruz(const Vector3D& a, const Vector3D& b) noexcept;

#endif
