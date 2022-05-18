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
		
		ModuleTree * left;
		ModuleTree * right;
	//public:
		ModuleTree() {};	
		ModuleTree(const FacetBox& box);
		FacetBox getBox() const;
		ModuleTree * getLeft() const;
		ModuleTree * getRight() const;
		void growBranch(const Vector3D& a, const Vector3D& b);
};
#endif