#ifndef FACETBOX_H
#define FACETBOX_H

using namespace std;

#include <iostream>
#include "Facet.hpp"
#include <math.h>

class FacetBox {
       	private:
	 	Facet * F;
		int n;
		Vector3D center;
	public:
        	Facet  operator [] (int k) const;
        	FacetBox();
		FacetBox(const Facet& facet);
		Vector3D getCenter() const;
		void updateCenter();
		void pushFacet(const Facet& facet);
		void pushFacet(const Vector3D& A, const Vector3D& B, const Vector3D& C);
		int getN() const;
		void translate(const Vector3D& t);

};


#endif
