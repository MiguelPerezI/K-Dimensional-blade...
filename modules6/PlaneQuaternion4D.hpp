#ifndef PLANE4D_H
#define PLANE4D_H

using namespace std;

#include <iostream>
#include "Quaternion.hpp"
#include <math.h>

class PlaneQuaternion4D {
       	private:
	 	Vector4D normal;
		Vector4D point[6], base;
		int m;
	public:
        	Vector4D  operator [] (int k) const;
        	PlaneQuaternion4D();
		PlaneQuaternion4D(const Vector4D& a, const Vector4D& b);
		double intersectionLine4D(const Vector4D& a, const Vector4D& b);
		int  linePointIntersection4D(const Vector4D& r, const Vector4D& a, const Vector4D& b);
		int equal4R(double a, double b);
};


#endif
