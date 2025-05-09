#include <cmath>
#include <fstream>
#include "Vector3D.hpp"

using namespace std;


/*————— element access ——————————————————————————————————————————————————*/
/* Write access ---------------------------------------------------------*/
double& Vector3D::operator[](int k)
{
    if (k < 0 || k > 2)
        throw std::out_of_range("Vector3D::operator[] index must be 0,1,2");
    return u[k];
}

/* Read access ----------------------------------------------------------*/
double Vector3D::operator[](int k) const
{
    if (k < 0 || k > 2)
        throw std::out_of_range("Vector3D::operator[] index must be 0,1,2");
    return u[k];
}
/*———————————————————————————————————————————————————————————————————————*/



/* Constructors —————————————————————————————————————————————————————————*/
/* Construct a Vector3D by providing three real numbers -----------------*/
/* Explicitly initialize all components via member-initializer list */
Vector3D::Vector3D(double xx, double yy, double zz)
    : u{ xx, yy, zz }
{ }
/*———————————————————————————————————————————————————————————————————————*/



/*----- compound operator  ----------------------------------------------*/
Vector3D& Vector3D::operator+= (const Vector3D& a) noexcept
{
    return *this= *this+a;
}

Vector3D& Vector3D::operator-= (const Vector3D& a) noexcept
{
    return *this= *this-a;
}

/* Unitary scalar / Vector ----------------------------------------------*/
Vector3D& Vector3D::operator/= (double scalar)
{
    if (std::abs(scalar) < 1e-12)
        throw std::runtime_error("Vector3D::operator/= division by zero");

    u[0] /= scalar;
    u[1] /= scalar;
    u[2] /= scalar;

    return *this;
}

/*———————————————————————————————————————————————————————————————————————*/

/* Free functions ———————————————————————————————————————————————————————*/
/* Adition ---------------------------------------------------------------*/
Vector3D operator + (const Vector3D& a, const Vector3D& b) noexcept
{
    return Vector3D(a.x() + b.x(),
                    a.y() + b.y(),
                    a.z() + b.z());
}

/* Subtraction -----------------------------------------------------------*/
Vector3D operator - (const Vector3D& a, const Vector3D& b) noexcept
{
    return Vector3D(a.x() - b.x(),
                    a.y() - b.y(),
                    a.z() - b.z());
}

/* Unitary minus ---------------------------------------------------------*/
Vector3D operator - (const Vector3D& a) noexcept
{
   return Vector3D(-a.x(), -a.y(), -a.z());
}

/* Scalar x Vector -------------------------------------------------------*/
Vector3D operator * (const double s, const Vector3D& b) noexcept
{
   return Vector3D( s * b.x(), s * b.y(), s * b.z() );
}

/* Scalar / Vector -------------------------------------------------------*/
Vector3D operator / (const Vector3D& a, const double s)
{
    if (abs(s) < 1e-12)
        throw std::runtime_error("Vector3D::operator/ division by zero");
    return Vector3D(
        a.x() / s, 
        a.y() / s, 
        a.z() / s
    );
}


/* Dot Product -----------------------------------------------------------*/
double operator * (const Vector3D& a, const Vector3D& b) noexcept
{
   return a.x() * b.x() + a.y() * b.y() + a.z() * b.z();
}

/* Cross Product ---------------------------------------------------------*/
Vector3D operator % (const Vector3D& a, const Vector3D& b) noexcept
{
   return Vector3D(a.y() * b.z() - a.z() * b.y(),
                  -a.x() * b.z() + a.z() * b.x(),
                   a.x() * b.y() - a.y() * b.x() );
}

/* Compare Vector closeness upto radius of 1E-8 --------------------------*/
bool     operator == (const Vector3D& a, const Vector3D& b) noexcept
{
    return abs(a-b)<1E-12;
}

/* Euclidean length ------------------------------------------------------*/
double abs(const Vector3D& a) noexcept
{
   double x= fabs(a.x());
   double y= fabs(a.y());
   double z= fabs(a.z());
   double t= x>y ? x>z ? x : y>z ? y : z : y>z ? y : z;
   if (t == 0)
      return 0.0;
   x/= t;
   y/= t;
   z/= t;
   return t*sqrt(x*x + y*y + z*z);
}

/* ∞-norm ---------------------------------------------------------------*/                                                                                          
double infty(const Vector3D& a) noexcept
{
   return fabs(a.x()) + fabs(a.y()) + fabs(a.z());
}

/* normalized copy ------------------------------------------------------*/
Vector3D unit(const Vector3D& v)
{
    double n= abs(v);
    return (n < 1e-12) ? Vector3D(0,0,0)
                       : v / n;
}

/* helper for linear interpolation on a segment -------------------------*/
Vector3D line(double t, const Vector3D& b, const Vector3D& e) noexcept
{
   return b + t * (e-b);
}

/*—————————————————————————————————————————————————————————-------------
 * @brief Check if three 3D points are colinear:
 *  returns true if (q-p)×(r-p) ≈ (0,0,0).
 *—————————————————————————————————————————————————————————-------------*/
bool areColinear(
    const Vector3D& p,
    const Vector3D& q,
    const Vector3D& r
) noexcept {
    // Form the two direction vectors
    Vector3D v1 = q - p;
    Vector3D v2 = r - p;
    // their cross-product is zero ⇔ they are colinear
    return (v1 % v2) == Vector3D(0,0,0);
}

/* stream operators -----------------------------------------------------*/
istream& operator >> (istream& is, Vector3D& a)
{
   while (!is.eof() && is.get() != '(')
       ;

   if (is.eof())
      {
       a= Vector3D(0,0,0);
       return is;
      }

   double x;
   is >> x;

   while (!is.eof() && is.get() != ',')
       ;

   if (is.eof())
      {
       a= Vector3D(x,0,0);
       return is;
      }

   double y;
   is >> y;

   while (!is.eof() && is.get() != ',')
       ;
   if (is.eof())
      {
       a= Vector3D(x,y,0);
       return is;
      }

   double z;
   is >> z;

   while (!is.eof() && is.get() != ')')
       ;
   a= Vector3D(x,y,z);

   return is;
}

ostream& operator << (ostream& os, const Vector3D& a)
{
   return os << "(" << a.x() << ", " << a.y() << ", " << a.z() << ") ";
}

/*———— Geometric helpers ————————————————————————————————————————————————*/
/* center for 4 vectors -------------------------------------------------*/
Vector3D centerM4(
                    const Vector3D& a0, const Vector3D& a1,
                    const Vector3D& a2, const Vector3D& a3) noexcept
{
	return (a0 + a1 + a2 + a3) * 0.25;
}

/*center for 5 vectors --------------------------------------------------*/
Vector3D cPenta(
                    Vector3D a0, Vector3D a1, 
                    Vector3D a2, Vector3D a3, Vector3D a4) noexcept
{
	return (a0 + a1 + a2 + a3 + a4) * 0.20;
}

/* center for 20 vectors ------------------------------------------------*/
Vector3D cDodeca(Vector3D a[]) noexcept
{
    Vector3D sum{0,0,0};
    for (int i = 0; i < 20; i++)
        sum += a[i];
    return sum * (1.0/20.0);
}

/* center for 8 vectors -------------------------------------------------*/
Vector3D centerM8(
                    const Vector3D& a0, const Vector3D& a1, 
                    const Vector3D& a2, const Vector3D& a3,
                    const Vector3D& a4, const Vector3D& a5, 
                    const Vector3D& a6, const Vector3D& a7) noexcept
{
    Vector3D sum = a0 + a1 + a2 + a3 + a4 + a5 + a6 + a7;
    return sum * (1.0/8.0);
}

/* cross product as function --------------------------------------------*/
Vector3D cruz(const Vector3D& a, const Vector3D& b) noexcept
{
  return Vector3D(a.y() * b.z() - a.z() * b.y(),
                 -a.x() * b.z() + a.z() * b.x(),
                  a.x() * b.y() - a.y() * b.x() );
}
