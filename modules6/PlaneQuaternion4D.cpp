
using namespace std;

#include "PlaneQuaternion4D.hpp"
#include <cmath>
#include <iomanip>
#include <fstream>


PlaneQuaternion4D::PlaneQuaternion4D() {}

PlaneQuaternion4D::PlaneQuaternion4D(const Vector4D& a, const Vector4D& b) {
	
	//DEFINE THE NORMAL
	this->normal = Vector4D(0, 0, 0, 1);
        this->point[0] = Vector4D( 1, 0, 0, 0);
	this->point[1] = Vector4D( 0, 1, 0, 0);
	this->point[2] = Vector4D( 0, 0, 1, 0);
	this->point[3] = Vector4D(-1, 0, 0, 0);
	this->point[4] = Vector4D( 0,-1, 0, 0);
	this->point[5] = Vector4D( 0, 0,-1, 0);

	this->base = Vector4D(b);
	//DEFINE THE DIFFERENCE VECTOR (v = a-b)
	Vector4D J = unit(a-b);

	//CHECK OF v coincides with (0, 0, 1)
	int ll = linePointIntersection4D(J, Vector4D(0, 0, 0, 1), Vector4D(0, 0, 0, 0)); 
	
	//Rotate
	//if (ll != 1) {
	//	
	//	base = Vector3D(b);
	//	Quaternion N = Quaternion(0, normal.V());
	//	normal = rotate(normal, a, b, N);
	//	Quaternion p0 = rotate(Quaternion(point[0]), a, b, N);
	//	Quaternion p1 = rotate(Quaternion(point[1]), a, b, N);
	//	Quaternion p2 = rotate(Quaternion(point[2]), a, b, N);
	//	Quaternion p3 = rotate(Quaternion(point[3]), a, b, N);

	//	point[0] = Vector3D(p0.V());
        //        point[1] = Vector3D(p1.V());
        //        point[2] = Vector3D(p2.V());
        //        point[3] = Vector3D(p3.V());

	//} else {
	//	normal = Quaternion(0.0, Vector3D(1e-20, 0, 1));
	//}

		

}	

//Fetch PlaneQuaternion Data
//Vector3D PlaneQuaternion::operator [] (int k) const {
//   if (k > 5)
//      return Vector3D(0, 0, 0);
//   else {
//   
//   	if (k == 0) return normal.V();
//	if (k == 1) return base;
//	if (k == 2) return point[0];
//	if (k == 3) return point[1];
//	if (k == 4) return point[2];
//	if (k == 5) return point[3];
//   }
//}


int  PlaneQuaternion4D::linePointIntersection4D(const Vector4D& r, const Vector4D& a, const Vector4D& b) {

	Vector4D cro = Cross(a, b, r);
	if (cro == Vector4D(0, 0, 0, 0)) return 1;
	else return 0;
}

int PlaneQuaternion4D::equal4R(double a, double b) {

	int epsilon = 1e-20;
	double cc = (a - b) * (a - b);

	if (sqrt(cc) < epsilon)
		return 1;
	else
		return 0;
}

//double PlaneQuaternion4D::intersectionLine4D(const Vector4D& a, const Vector4D& b) {
//
//
//	double D = normal * base;
//	double Dpr = normal * b;
//	double Dps = normal *(a-b);
//
//	double t = (D - Dpr) / Dps;
//
//	if (0.0 < t && t < 1.0 ) return t;
//        else {
//                if (equal4R(t, 0.0) == 1) return t;
//                if (equal4R(t, 1.0) == 1) return t;
//                if (t < 0.0 || 1.0 < t) return -10000.0;
//        }
//}

double PlaneQuaternion4D::intersectionLine4D(const Vector4D& a,
                                             const Vector4D& b)
{
    const double D   = normal * base;      // n·P0   (plane constant)
    const double Dpr = normal * b;         // n·P1
    const double Dps = normal * (a - b);   // n·(P0 − P1)

    if (std::abs(Dps) < 1e-12)               // line is parallel to plane
        return std::numeric_limits<double>::infinity();  // or throw

    const double t = (D - Dpr) / Dps;        // parametric distance

    if (0.0 < t && t < 1.0)          return t;      // interior hit
    if (equal4R(t, 0.0) == 1)        return 0.0;    // hits at a
    if (equal4R(t, 1.0) == 1)        return 1.0;    // hits at b

    return -10000.0;                              // no intersection
}
