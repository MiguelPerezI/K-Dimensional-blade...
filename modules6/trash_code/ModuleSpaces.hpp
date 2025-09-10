#ifndef MODULESPACES_H
#define MODULESPACES_H

using namespace std;

#include <iostream>
#include "FacetBox.hpp"
#include "FacetGash.hpp"
#include <math.h>

class ModuleSpaces {
       	private:
	 	FacetBox box0, box1;
		FacetGash In;
		int cut;
	public:
        	FacetBox  operator [] (int k) const;
        	ModuleSpaces(const Vector3D& a, const Vector3D& b, const FacetBox& D);
		ModuleSpaces() {};
		Vector3D piecewise(double t, const Vector3D& a, const Vector3D& b);
		void restart();
		int getCut() const {return cut;}
};

#endif
