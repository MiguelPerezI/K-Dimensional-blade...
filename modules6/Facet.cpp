
using namespace std;

#include "Facet.hpp"
#include <cmath>
#include <iomanip>
#include <fstream>
#include <stdexcept>

/*——————————————————————————————————————————————————————————————————————————————*/
// Constructor: three Quaternion vertices
/*——————————————————————————————————————————————————————————————————————————————*/
Facet::Facet(const Quaternion& a, const Quaternion& b, const Quaternion& c)
    : A(a), B(b), C(c)
{
    // Compute normalized facet normal from cross product of edges
    Vector3D edge1 = B.V() - A.V();
    Vector3D edge2 = C.V() - A.V();
    Vector3D crossVec = edge1 % edge2;
    N = Quaternion(0.0, unit(crossVec));
}


/*——————————————————————————————————————————————————————————————————————————————*/
// Constructor: three Vector3D points
/*——————————————————————————————————————————————————————————————————————————————*/
Facet::Facet(const Vector3D& a, const Vector3D& b, const Vector3D& c)
    : A(0.0, a), B(0.0, b), C(0.0, c)
{
    Vector3D edge1 = b - a;
    Vector3D edge2 = c - a;
    Vector3D crossVec = edge1 % edge2;
    N = Quaternion(0.0, unit(crossVec));
}


/*——————————————————————————————————————————————————————————————————————————————*/
// Element access: 0→A, 1→B, 2→C, 3→normal
/*——————————————————————————————————————————————————————————————————————————————*/
Vector3D Facet::operator[](int k) const
{
    switch (k) {
        case 0: return A.V();
        case 1: return B.V();
        case 2: return C.V();
        case 3: return N.V();
        default: throw std::out_of_range("Facet index must be 0..3");
    }
}


/*——————————————————————————————————————————————————————————————————————————————*/
// Update facet vertices and recompute normal
/*——————————————————————————————————————————————————————————————————————————————*/
void Facet::updateFacet(const Vector3D& a, const Vector3D& b, const Vector3D& c)
{
    A = Quaternion(0.0, a);
    B = Quaternion(0.0, b);
    C = Quaternion(0.0, c);
    Vector3D crossVec = (b - a) % (c - a);
    N = Quaternion(0.0, unit(crossVec));
}


/*——————————————————————————————————————————————————————————————————————————————*/
// Get triangle centroid
/*——————————————————————————————————————————————————————————————————————————————*/
Vector3D Facet::getCenter() const
{
    return (A.V() + B.V() + C.V()) / 3.0;
}


/*——————————————————————————————————————————————————————————————————————————————*/
// Output stream
/*——————————————————————————————————————————————————————————————————————————————*/
std::ostream& operator<<(std::ostream& os, const Facet& f)
{
    int w = os.width();
    int p = os.precision();
    os << std::fixed << std::setprecision(p)
       << "Facet< A=" << f[0]
       << ", B=" << f[1]
       << ", C=" << f[2]
       << ", N=" << f[3]
       << ">";
    os.width(w);
    os.precision(p);
    return os;
}

/*——————————————————————————————————————————————————————————————————————————————*/
// Translate facet by offset vector
/*——————————————————————————————————————————————————————————————————————————————*/
void Facet::translate(const Vector3D& offset)
{
    A = A + Quaternion(0.0, offset);
    B = B + Quaternion(0.0, offset);
    C = C + Quaternion(0.0, offset);
    // Recompute normal after translation
    Vector3D edge1 = B.V() - A.V();
    Vector3D edge2 = C.V() - A.V();
    N = Quaternion(0.0, unit(edge1 % edge2));
}


/*——————————————————————————————————————————————————————————————————————————————*/
// Scale (crunch) the triangle towards point a by factor t
/*——————————————————————————————————————————————————————————————————————————————*/
void Facet::crunch(double t, const Vector3D& pivot)
{
    auto scalePoint = [&](const Quaternion& Q) {
        Vector3D P = Q.V();
        return pivot + t * (P - pivot);
    };
    Vector3D a = scalePoint(A);
    Vector3D b = scalePoint(B);
    Vector3D c = scalePoint(C);
    updateFacet(a, b, c);
}
