using namespace std;

#include "ModuleSpaces.hpp"
#include <cmath>
#include <iomanip>
#include <fstream>



FacetBox  ModuleSpaces::operator [] (int k) const {
	if (k == 0) return box0;
	else return box1;
}

ModuleSpaces::ModuleSpaces(const Vector3D& a, const Vector3D& b, const FacetBox& D) {

	//Prepare Blade
	In = FacetGash(a, b);

	//First half
	cut = 0;
	for (int i = 0; i < D.getN(); i++) {
		cut = In.cutFacet(D[i]);
		if (In.checkFacet(D[i]) == 1)
			box0.push(D[i]);
	}

	if (cut != 0) {
		In.readList(&box0);

		In.restart();
		In.updateOrientation();
		//Second half
        	
		for (int i = 0; i < D.getN(); i++) {
        	        In.cutFacet(D[i]);
        	        if (In.checkFacet(D[i]) == 1)
        	                box1.push(D[i]);
        	}

        		In.readList(&box1);
			In.restart();

	}
}

Vector3D ModuleSpaces::piecewise(double t, const Vector3D& a, const Vector3D& b) {return (t*(a-b)) + b;}
void ModuleSpaces::restart() {
	
	box0.empty();
	box1.empty();
	cut = 0;
}

