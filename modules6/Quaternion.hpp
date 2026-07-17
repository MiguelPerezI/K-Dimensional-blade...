#ifndef Quaternion_H
#define Quaternion_H

using namespace std;

#include <iostream>
#include "Vector3D.hpp"
#include "Vector4D.hpp"
#include <iosfwd>          // forward-declare std::istream / ostream
#include <cmath>

/*—————————————————————————————————————————————————————————————
 * Simple Quaternion module - uses Vector3D-v⃗/Vector4D-v⃗ class
 *—————————————————————————————————————————————————————————————*/
class Quaternion {
    private:
        /*————————————————————————————————————————————————————————————————————————————
         * Quaternion q = u + v⃗
         *   where u is the scalar part, v⃗ is the 3D vector part
         *——————————————————————————————————————————————————————————————————————————*/
        double  u{};            // scalar part (initialized to 0)
        Vector3D v{};           // vector part

    public:
        /* ---------- Rule of Zero compliant -------------------------------------- */
        Quaternion() = default;
        ~Quaternion() = default;
        Quaternion(Quaternion&&) noexcept = default;
        Quaternion& operator=(const Quaternion&) = default;
        Quaternion& operator=(Quaternion&&) noexcept = default;
        Quaternion(const Quaternion&) = default;
         
        /* ---------- main constructors  ------------------------------------------ */
        explicit Quaternion(double s, const Vector3D& vv = Vector3D{});// body in .cpp
        /* pure-vector constructor: q = (0 , a⃗) */
        Quaternion(const Vector3D& a) : u{0.0}, v{a} {}
        
        /*--- element-style access ------------------------------------------------ */
        /*  k == 0  → vector part   (reference / copy)
            k == 1  → scalar broadcast (u,u,u)
            else    → (0,0,0)*/
        Vector3D& operator[](int k);            // non-const → reference
        Vector3D  operator[](int k) const;      // const     → value

        /*------------- inspectors -------------------------------------------------*/
        /*   (cheap, noexcept, never modify state)                                  */
        double      r() const noexcept { return u; }
        Vector3D    V() const noexcept { return v; }
        double      i() const noexcept { return v.x(); };
	    double      j() const noexcept { return v.y(); };
	    double      k() const noexcept { return v.z(); };
    
        /*---- Conversions ---------------------------------------------------------*/
        Vector4D v4() const noexcept {                        // (u,x,y,z)
            return { u, v.x(), v.y(), v.z() };
        }

        /*----- Mutating operators -------------------------------------------------*/
         Quaternion& operator+=(const Quaternion&);
         Quaternion& operator-=(const Quaternion&);
         Quaternion& operator/=(double);


        /*----- Utilities ----------------------------------------------------------*/
        Quaternion conjugate() const noexcept                // u  −v⃗
        { return Quaternion{ u, -v }; }

        Quaternion toHyperboloid() const;
/*----- free-function operators (declared below) ------------------------*/
};

/* —————————————————————————————— free functions ————————————————————————————————*/
Quaternion operator + (const Quaternion& a, const Quaternion& b);
Quaternion operator - (const Quaternion& a, const Quaternion& b);
Quaternion operator * (const Quaternion& a, const Quaternion& b);
Quaternion operator * (const double a, const Quaternion& b);
bool operator==(const Quaternion& q1, const Quaternion& q2) noexcept;

Quaternion sigma(const Quaternion& x, const Quaternion& a, double r);
Quaternion Qan(double theta, const Vector3D& n);
Quaternion cross(const Quaternion& n, const Quaternion& z);
Quaternion rotate(const Quaternion& p, const Quaternion& normal, const Quaternion& J);
Quaternion lerp(double t, const Quaternion& p, Quaternion& q) noexcept;

/* --- rotation / parallel-transport helpers --------------------------------
 * q-prefixed names avoid clashing with Vector3D::unit / Vector3D::abs.
 * Used by the surface-following car camera in main17_gpu.cpp.
 */
double    qabs(const Quaternion& q) noexcept;          // norm sqrt(r² + ‖v‖²)
Quaternion qunit(const Quaternion& q);                  // normalized copy; identity if ~0
Quaternion qFromToRotation(const Vector3D& from,       // shortest-arc quaternion that
                           const Vector3D& to);         // rotates unit(from) onto unit(to)
Vector3D   qRotateVec(const Quaternion& q,             // q v q⁻¹ sandwich (q auto-normalized),
                      const Vector3D& v);               // returning a Vector3D directly
Quaternion qslerp(const Quaternion& a,                 // shortest-path spherical lerp
                  const Quaternion& b, double t);
Quaternion qFromBasis(const Vector3D& forward,         // unit quaternion R with R*(0,0,1)=forward
                      const Vector3D& up);              // and R*(0,1,0)=up (up re-orthogonalized)

istream& operator >> (istream& is, Quaternion& a);
ostream& operator << (ostream& os, const Quaternion& a);

#endif
