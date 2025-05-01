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
        /*  q = u + v⃗   where  u is scalar part,  v⃗ = (x,y,z) the vector part  */
        double  u {};          // scalar part
        Vector3D v;            // vector part

    public:
        /*--- element-style access ------------------------------------------------*/
        /*  k == 0  → vector part   (reference / copy)  
            k == 1  → scalar broadcast (u,u,u)  
            else    → (0,0,0)                                               */
        Vector3D& operator[](int k);            // non-const → reference
        Vector3D  operator[](int k) const;      // const     → value 


        /* ---------- constructors  ------------------------------------------ */
        Quaternion(double s = 0.0,
                const Vector3D& vv = Vector3D{});          // body in .cpp

        Quaternion(const Vector4D& a);
        Quaternion(const Quaternion&);


        /*------------- simple accessors --------------------------------------*/
        /*— 1. Inspectors ————————————————————————————————————————————*/
        /*   (cheap, noexcept, never modify state)                           */
        double   r() const noexcept { return u; }
        Vector3D V() const noexcept { return v; }
        double i() const noexcept { return v.x(); };
	    double j() const noexcept { return v.y(); };
	    double k() const noexcept { return v.z(); };
        //Vector4D v4() const noexcept {return Vector4D(u, v.x(), v.y(), v.z());};
    
        /*— 2. Conversions ————————————————————————————————————————————*/
        Vector4D v4() const noexcept {                        // (u,x,y,z)
            return { u, v.x(), v.y(), v.z() };
        }


        /*— 3. Mutating operators ————————————————————————————————————*/
         Quaternion& operator=(const Quaternion&);
         Quaternion& operator+=(const Quaternion&);
         Quaternion& operator/=(double);


        /*— 4. Utility functions ————————————————————————————————————*/
        Quaternion conjugate() const noexcept                // u  −v⃗
        { return { u, -v }; }

        /* pure-vector constructor: q = (0 , a⃗) */
        Quaternion(const Vector3D& a) : u{0.0}, v{a} {}

/*----- free-function operators (declared below) ------------------------*/
};

/* —————————————————————————————— free functions ————————————————————————————————*/
Quaternion operator + (const Quaternion& a, const Quaternion& b);
Quaternion operator - (const Quaternion& a, const Quaternion& b);
Quaternion operator * (const Quaternion& a, const Quaternion& b);
Quaternion Qan(double theta, const Vector3D& n);
Quaternion cross(const Quaternion& n, const Quaternion& z);
Quaternion rotate(const Quaternion& p, const Quaternion& normal, const Quaternion& J);

Quaternion operator * (const double a, const Quaternion& b);
istream& operator >> (istream& is, Quaternion& a);
ostream& operator << (ostream& os, const Quaternion& a);

#endif
