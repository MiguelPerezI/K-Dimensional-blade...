using namespace std;

#include "FacetBox.hpp"
#include <cmath>
#include <iomanip>
#include <fstream>


FacetBox::FacetBox() {};
FacetBox::FacetBox(const Facet& facet) {
	
	F = (Facet *) malloc (1 * sizeof(Facet));
	n = 1;
	F[0] = Facet(facet);

	Vector3D tmp = Vector3D(facet[0] + facet[1] + facet[2]);
	center = Vector3D(0.333333333333333 * tmp);

}

Facet FacetBox::operator [] (int k) const {
   if (k > n-1)
      	return F[0];
   else {
	return F[k];
   }
}

Vector3D FacetBox::getCenter() const {return center;}
int FacetBox::getN() const {return n;}
void FacetBox::updateCenter() {
	
	//center = 0.0 * center;
	//for (int i = 0; i < getN(); i++) 
	//	center = center + F[i].getNormal().V();
	//
	//center = (1.0/((double) getN())) * center;

}


