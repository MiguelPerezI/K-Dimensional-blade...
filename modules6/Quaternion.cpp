
using namespace std;

#include "Quaternion.hpp"
#include <cmath>
#include <iomanip>
#include <fstream>
#include <stdexcept>

/* ---------- Constructors -------------------------------------------------- */

// Explicit constructor with scalar and vector parts (using initializer list)
Quaternion::Quaternion(double scalar, const Vector3D& vector) 
    : u(scalar), v(vector) {}

/* ---------- Element Access ------------------------------------------------ */
// Non-const operator[]
/*---------------------------------------------------------------------------
 * Non-const overload (returns a reference so the caller can
 * modify v’s components directly when k == 0).
 *--------------------------------------------------------------------------- */
Vector3D& Quaternion::operator[](int k)
{
    if (k == 0)
        return v;  // vector part
    throw std::out_of_range("Quaternion non-const operator[]: valid index is 0 only");
}

/*---------------------------------------------------------------------------
 * Const overload (returns a value).
 *  k == 0 → vector part        (x, y, z)
 *  k == 1 → scalar broadcast   (u, u, u)
 *  else   → zero vector        (0, 0, 0)
 *--------------------------------------------------------------------------- */
Vector3D Quaternion::operator[](int k) const
{
    switch (k) {
        case 0:  return v;                         // (x, y, z)
        case 1:  return Vector3D(u, u, u);         // (u, u, u)
        default: return Vector3D(0.0, 0.0, 0.0);   // out-of-range → zero
    }
}

/* ---------- Mutating Operators -------------------------------------------- */

// Compound addition (modifies object)
Quaternion& Quaternion::operator += (const Quaternion& other)
{
   return *this= *this + other;
}

// Compound addition (modifies object)
Quaternion& Quaternion::operator -= (const Quaternion& other)
{
   return *this= *this - other;
}

// Compound division by scalar (with exception safety)
Quaternion& Quaternion::operator/=(double scalar) {
    if (std::abs(scalar) < 1e-12)
        throw std::invalid_argument("Quaternion division by zero");
    u /= scalar;
    v /= scalar;
    return *this;
}

/*----- Hyperbolic projection ———————————————————————————————————————————————*/
/*----- Hyperbolic projection — subgroup toHyperboloid ------------------
 * Project a pure-vector quaternion (r()==0) from the unit ball into
 * the hyperboloid model via µ(x) = (1/√(1−‖x‖²), x/√(1−‖x‖²)),
 * then re-project to the gnomic disk: (0, x') where x' = x/(1+u').
 * @throws std::domain_error if ‖x‖ ≥ 1
 */
Quaternion Quaternion::toHyperboloid() const {
    double sq = v * v;  // squared norm of vector part
    if (sq >= 1.0) throw std::domain_error("toHyperboloid: vector norm must be < 1");
    double s = 1.0 / std::sqrt(1.0 - sq);
    Quaternion ret{ s, s * v };
    // final gnomic re-projection
    return Quaternion{ 0.0, (1.0 / (1.0 + ret.r())) * ret.V() };
}

/* Quaternion spherical inversion (sigma method) -------------------------*/
Quaternion sigma(const Quaternion& x, const Quaternion& a, double r)
{
    // For quaternions, we need to define "distance"
    // Option 1: Treat as 4D vectors (u, v.x, v.y, v.z)
    // Option 2: Use quaternion-specific distance measure
    
    // Using 4D approach:
    Quaternion diff = x - a;
    
    // Distance squared in 4D space: ||(u, x, y, z)||²
    double dist_squared = diff.r() * diff.r() + (diff.V() * diff.V());
    
    if (std::abs(dist_squared) < 1e-12)
        throw std::runtime_error("sigma: quaternion coincides with sphere center");
    
    double r_squared = r * r;
    double factor = r_squared / dist_squared;
    
    return a + (factor * diff);
}

/* ---------- Free-function Operators --------------------------------------- */
// Quaternion addition
Quaternion operator + (const Quaternion& a, const Quaternion& b) {
    return Quaternion(a.r() + b.r(), a.V() + b.V());
}

// Quaternion subtraction
Quaternion operator - (const Quaternion& a, const Quaternion& b) {
    return Quaternion(a.r() - b.r(), a.V() - b.V());
}

// Quaternion multiplication (Hamilton product)
Quaternion operator * (const Quaternion& a, const Quaternion& b) {
    return Quaternion(
        (a.r() * b.r()) - (a.V() * b.V()),
        (a.r() * b.V()) + (b.r() * a.V()) + cruz(a.V(), b.V()));		
}

// Scalar multiplication
Quaternion operator * (const double a, const Quaternion& b) {
   return Quaternion( b.r(), a * Vector3D(b.V()) );
}

// Quaternion comparison
bool operator==(const Quaternion& q1, const Quaternion& q2) noexcept {
    return std::abs(q1.r() - q2.r()) < 1e-12 && q1.V() == q2.V();
}


/* ---------- Utility and Helper Functions ----------------------------------*/
// Quaternion from angle-axis representation
Quaternion Qan(double theta, const Vector3D& axis) {
    constexpr double epsilon = 1e-10;
    if (std::abs(theta) < epsilon)
        return Quaternion(0.0, axis);
    else
        return Quaternion(std::cos(0.5 * theta), std::sin(0.5 * theta) * unit(axis));
};

// Cross product as quaternion (pure vector quaternion)
Quaternion cross(const Quaternion& a, const Quaternion& b) {
    Vector3D result = (1.0/(a.V() * b.V())) * (a.V() % b.V());
    return Quaternion(0.0, result);
}

/**
 * @brief Rotate a point-quaternion p from position a to position b around a given axis.
 *
 * We treat p as a “point” in 3D space encoded as a pure-vector quaternion,
 * then:
 *   1. Determine the rotation axis (τ) and angle (φ) needed to turn the ray from b→a 
 *      onto the global “up” direction (0,0,1).
 *   2. Build the corresponding rotation quaternion Q = Qan(φ, τ).
 *   3. Apply that rotation to p:  p1 = Q * p * Q⁻¹.
 *   4. Finally translate p1 so that it’s moved from origin-centered coordinates back to
 *      “around point b” in world space.
 *
 * @param p      A pure-vector Quaternion representing the point to rotate.
 * @param a      The original 3D position (Vector3D) of the point.
 * @param b      The target 3D position (Vector3D) to rotate around.
 * @param normal A unit-quaternion whose vector part is the desired rotation axis.
 * @return       A new Quaternion encoding the rotated (and translated) point.
 */
Quaternion rotate(
    const Quaternion& p, 
    const Vector3D& a, 
    const Vector3D& b, 
    const Quaternion& normal
) {

	// 1) Compute unit direction from b to a as a pure-vector quaternion:
    //    difference = (0, normalize(a - b))
	Quaternion difference{ 0.0, unit(a-b) };

    // 2) Quick colinearity check: does (a - b) already align with the z-axis?
    //    If line(...) returns 1, the vector (a-b) lies in the z-direction.
	int isColinearWithZ = areColinear(unit(a-b), {0, 0, 1}, {0, 0, 0});
	if (isColinearWithZ == true) {
        // Already aligned—no rotation needed.
        return p;
	}

    // 3) Compute the rotation axis τ = normalize(normal × difference):
    //    - normal.V(): axis we want to rotate around (as Vector3D)
    //    - difference.V(): direction to align with z
    Vector3D tau = unit(unit(normal.V()) % unit(difference.V()));

    // 4) Compute the rotation angle φ between normal and difference:
    //    φ = arccos( normal∙difference )
    double phi =   acos(unit(normal.V()) * unit(difference.V()));

    // 5) Build the unit quaternion representing rotation of φ about τ:
    //    τQ = [ cos(φ/2),  sin(φ/2) * τ ]
    //    τ = eigenvector for space
    Quaternion tauQ = Qan(phi, tau);
   
    // 6) Rotate the original point-quaternion p:  p1 = τQ * p * τQ.conjugate()
    // Rotate p with respect to the eigenvector τ
    Quaternion p1 = tauQ * p * tauQ.conjugate();

    // 7) Shift p1 so that the center moves from origin to point b:
    //    Add pure-vector quaternion (0, b) 
    //    i.e. translate p1
    return  (p1 + Quaternion(0, b));
		


}

/*———————————————————————————————————————————————————————————————————————————
 * @brief Linearly interpolate between two quaternions:
 *  component-wise lerp: (1–t)*p + t*q.
 *———————————————————————————————————————————————————————————————————————————*/
Quaternion lerp(
    double t,
    const Quaternion& p,
    const Quaternion& q
) noexcept {
    // lerp the real (scalar) parts
    double u = p.r() + t * (q.r() - p.r());
    // lerp the vector parts via the Vector3D::line helper
    Vector3D v = line(t, p.V(), q.V());
    return Quaternion(u, v);
}


ostream& operator << (ostream& os, const Quaternion& a) {

   int w= os.width();
   int p= os.precision();
   os << setw(0) << "(" 
      << setw(w) << setprecision(p) << a.r() << setw(0) << ",  (" 
      << setw(w) << setprecision(p) << a.i() << setw(0) << ", " 
      << setw(w) << setprecision(p) << a.j() << setw(0) << ", " 
      << setw(w) << setprecision(p) << a.k() << setw(0) << ") )";
   os.width(w);
   os.precision(p);
   return os;
}



