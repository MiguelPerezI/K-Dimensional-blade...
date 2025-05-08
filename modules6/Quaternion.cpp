
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

/*----- Hyperbolic projection ——————————————————————————————————————————*/
    /*----- Hyperbolic projection — subgroup toHyperboloid ----------------------
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


/* ---------- Utility and Helper Functions ---------------------------------- */
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

int  line(const Vector3D& r, const Vector3D& a, const Vector3D& b) {

                        Vector3D cro = (a-r) % (b-r); 
                        if (cro == Vector3D(0, 0, 0)) return 1;
                        else return 0;
}

Quaternion rotate(const Quaternion& p, const Vector3D& a, const Vector3D& b, const Quaternion& normal) {

	//The diference Vector = a - b;
	Quaternion difference = Quaternion( 0.0, unit(a-b));

	int ll = line(unit(a-b), Vector3D(0, 0, 1), Vector3D(0, 0, 0));

	if (ll != 1) {
		

		Vector3D tau = unit(unit(normal.V()) % unit(difference.V()));
		double phi =   acos(unit(normal.V()) * unit(difference.V()));
	
		//tau = eigenvector for space
		Quaternion tauQ = Qan(phi, tau);

		//Rotate p with respect to the eigenvector
		Quaternion p1 = tauQ * p * tauQ.conjugate();
		
		//translate p1
		return  (p1 + Quaternion(0, b));
		
	} else { return p;}


}





