using namespace std;

#include "ModuleTree.hpp"
#include <cmath>
#include <iomanip>
#include <fstream>

ModuleTree::ModuleTree(const FacetBox& B) {

	//DATA
	box.copy(B);

	//BRANCHES
	left = NULL;
	right = NULL;

}

FacetBox ModuleTree::getBox() const {
	return box;
}

void ModuleTree::growBranch(const Vector3D& a, const Vector3D& b) {
	
	ModuleSpaces mod = ModuleSpaces(a, b, box);	
	
	this->left = new ModuleTree;
	ModuleTree * node0 = new ModuleTree(mod[0]);
	this->left = node0;

	this->right = new ModuleTree;
        ModuleTree * node1 = new ModuleTree(mod[1]);
        this->right = node1;
}

