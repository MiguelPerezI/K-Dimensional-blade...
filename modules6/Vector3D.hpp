#ifndef VECTOR3D_H
#define VECTOR3D_H

using namespace std;

#include <iostream>
#include <cmath>//<math.h>

/*--------------------------------------------------------------
 * Simple 3-component vector for geometry / linear-algebra work
 *-------------------------------------------------------------*/
class Vector3D {
    protected:
        double u[4]; // x, y, z   (size 3 is enough)

    public:
        /*----- element access --------------------------------------------------*/
        double& operator [] (int k);        // write-access
        double  operator [] (int k) const;  // read-access

        double x() const { return u[0]; };
        double y() const { return u[1]; };
        double z() const { return u[2]; };


        /*----- constructors ----------------------------------------------------*/
        Vector3D(double xx = 0, double yy = 0, double zz = 0);
        Vector3D(const Vector3D& a);
        Vector3D& operator =  (const Vector3D& a);

        /*----- compound operator  ------------------------------------------------------*/
        Vector3D& operator += (const Vector3D& a);

        /*----- free-function operators (declared below) ------------------------*/
};

                               /*--- free functions ---*/
Vector3D operator+(const Vector3D& a, const Vector3D& b);
Vector3D operator-(const Vector3D& a, const Vector3D& b);
Vector3D operator-(const Vector3D& a);                    // unary minus
Vector3D operator*(double s, const Vector3D& v);          // scalar × vector
Vector3D operator/(const Vector3D& v, double s);
double   operator*(const Vector3D& a, const Vector3D& b); // dot product
Vector3D operator%(const Vector3D& a, const Vector3D& b); // cross product


bool     operator == (const Vector3D&, const Vector3D& );

double   abs   (const Vector3D& v);                       // Euclidean length
double   infty (const Vector3D& v);                       // ∞-norm
Vector3D unit  (const Vector3D& v);                       // normalized copy


/* helper for linear interpolation on a segment */
Vector3D line(double t, const Vector3D& p, const Vector3D& q);

/* stream operators */
std::istream& operator>>(std::istream& is, Vector3D& v);
std::ostream& operator<<(std::ostream& os, const Vector3D& v);


Vector3D centerM4(const Vector3D& a0, const Vector3D& a1, const Vector3D& a2, const Vector3D& a3);
Vector3D cPenta(Vector3D a0, Vector3D a1, Vector3D a2, Vector3D a3, Vector3D a4);
Vector3D cDodeca(Vector3D a[]);
Vector3D centerM8(const Vector3D& a0, const Vector3D& a1, const Vector3D& a2, const Vector3D& a3,
                  const Vector3D& a4, const Vector3D& a5, const Vector3D& a6, const Vector3D& a7);
Vector3D cruz(const Vector3D& a, const Vector3D& b);
istream& operator >> (istream& is, Vector3D& a);
ostream& operator << (ostream& os, const Vector3D& a);

Vector3D line(double t, Vector3D& b, Vector3D& e);

#endif
