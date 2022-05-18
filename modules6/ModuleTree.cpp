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

int ModuleTree::getN() const {
	return box.getN();
}

void ModuleTree::growBranch(const Vector3D& a, const Vector3D& b) {
	
	ModuleSpaces mod = ModuleSpaces(a, b, box);	
	
	if (mod[0].getN() != 0) {	
		this->left = new ModuleTree;
		ModuleTree * node0 = new ModuleTree(mod[0]);
		this->left = node0;
	} else {
		//this->left = new ModuleTree;
                //ModuleTree * node0 = new ModuleTree(box);
                //this->left = node0;
	}

	if (mod[1].getN() != 0) {
		this->right = new ModuleTree;
        	ModuleTree * node1 = new ModuleTree(mod[1]);
        	this->right = node1;
	} else {
		//this->right = new ModuleTree;
                //ModuleTree * node1 = new ModuleTree(box);
                //this->right = node1;
	}
}

//FacetBox ModuleTree::search(int l, int * S, ModuleTree * tree) {	
//
//	ModuleTree * user = tree;
//	ModuleTree * tmp;	
//
//	if (user->)
//
//}

