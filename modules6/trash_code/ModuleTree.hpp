#ifndef MODULETREE_H
#define MODULETREE_H

using namespace std;

#include <iostream>
#include "FacetBox.hpp"
#include "ModuleSpaces.hpp"
#include <math.h>

class ModuleTree {
       	//private:	
	public:	
		FacetBox box;
		string s;
		
		ModuleTree * left;
		ModuleTree * right;
	//public:
		ModuleTree() {};	
		ModuleTree(const FacetBox& box, string S);
		FacetBox getBox() const;
		int getN() const;
		string getString() const;
		ModuleTree * getLeft() const;
		ModuleTree * getRight() const;
		void growBranch(const Vector3D& a, const Vector3D& b);
		FacetBox search(int l, int * S, ModuleTree * tree);
		void crunch(double a);
};
#endif
