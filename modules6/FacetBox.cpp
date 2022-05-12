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
	
	center = 0.0 * center;
	for (int i = 0; i < getN(); i++) 
		center = center + F[i].getCenter();
	
	center = (1.0/((double) getN())) * center;

}

void FacetBox::pushFacet(const Facet& facet) {

	if (n == 0) {
	
		F[0] = Facet(facet[0], facet[1], facet[2]);
		n = 1;
	} else {
	
		F = (Facet *) realloc(F, sizeof(Facet) * (n+1));
		n = n + 1;
		F[n-1] = Facet(facet[0], facet[1], facet[2]);
		updateCenter();
	}
}

void FacetBox::pushFacet(const Vector3D& A, const Vector3D& B, const Vector3D& C) {

	if (n == 0) {
		
		 F[0] = Facet(A, B, C);
		 n = 1;
	} else {
	
		F = (Facet *) realloc(F, sizeof(Facet) * (n+1));
		n = n + 1;
		F[n-1] = Facet(A, B, C);
	}

	updateCenter();

}
