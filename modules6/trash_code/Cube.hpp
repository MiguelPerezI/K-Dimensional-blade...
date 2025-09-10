#ifndef CUBE_H
#define CUBE_H

using namespace std;

#include <iostream>
#include "FacetBox.hpp"
#include <math.h>

class Cube {
       	private:
		FacetBox f;
	public:
        	FacetBox  operator [] (int k) const;
        	Cube(double radius, const Vector3D& center);
		Cube();
};

ostream& operator << (ostream& os, const Facet& a);

#endif
