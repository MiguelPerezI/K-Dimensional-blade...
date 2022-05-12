#ifndef FACETGASH_H
#define FACETGASH_H

using namespace std;

#include <iostream>
#include "PlaneQuaternion.hpp"
#include "Facet.hpp"
#include "FacetBox.hpp"
#include <math.h>

class FacetGash {
       	private:
	 	Quaternion * * M;
		Quaternion head, base;
		PlaneQuaternion blade;
		int n, m;
		Vector3D inter0, inter1, inter2;
		Facet facet;
		Vector3D * V;
		Facet * facets;

	public:
        	Vector3D  operator [] (int k) const;
		FacetGash() {};
		FacetGash(const Vector3D& planeH, const Vector3D& planeB);
		void updateFacetGash(const Vector3D& newHead, const Vector3D& newBase);
		PlaneQuaternion getSaw() const;
                Vector3D getNormal () const;
                Quaternion getM(int i, int j) const;
                int getN() const;
                int getM() const;
                double  getT(int nn) const;
                double getT1(int nn) const;
                double getT2(int nn) const;
                Facet getFacet(int n0) const;
                Vector3D getDir(double t) const;
                int checkPoint(const Quaternion& p, const Quaternion& J);
                void updateOrientation();
                void intersectFacet(const Facet& facet0);
		void restart();
		void readListC(const Facet& facet, FacetBox * pila);
		void readList(FacetBox * pila);
		int  linePointIntersection(const Vector3D& r, const Vector3D& a, const Vector3D& b);
		void auxFun(int a, int b, int aa, int bb, int cc, int i, const Vector3D& A0, const Vector3D& B0, const Vector3D& C0, FacetBox * pila );
		void checkR (const Vector3D& r0, const Vector3D& r1, const Vector3D& A0, const Vector3D& B0, const Vector3D& C0, FacetBox * pila);

};


#endif

