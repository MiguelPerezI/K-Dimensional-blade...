#ifndef DODECAHEDRON_H
#define DODECAHEDRON_H

using namespace std;

#include <iostream>
#include "Facet.hpp"
#include <math.h>

class Dodecahedron {
       	private:
	 	Quaternion v[20];
		Facet f[8];
	public:
        	Facet  operator [] (int k) const;
        	Dodecahedron(double radius, const Vector3D& center);
		Dodecahedron();
};

istream& operator >> (istream& is, Facet& a);
ostream& operator << (ostream& os, const Facet& a);

#endif
