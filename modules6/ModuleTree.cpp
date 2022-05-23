using namespace std;

#include "ModuleTree.hpp"
#include <cmath>
#include <iomanip>
#include <fstream>

ModuleTree::ModuleTree(const FacetBox& B, string S) {

	//DATA
	box.copy(B);
	s = S;	
	//BRANCHES
	left = NULL;
	right = NULL;

}

FacetBox ModuleTree::getBox() const {
	return box;
}

void ModuleTree::crunch(double a) {
	box.crunch(a, box.getCenter());
}

int ModuleTree::getN() const {
	return box.getN();
}

void ModuleTree::growBranch(const Vector3D& a, const Vector3D& b) {
	
		
	ModuleSpaces mod = ModuleSpaces(a, b, box);	
	
	if (mod.getCut() != 0)	{

		this->left = new ModuleTree;
		ModuleTree * node0 = new ModuleTree(mod[0], "left");
		this->left = node0;

		this->right = new ModuleTree;
        	ModuleTree * node1 = new ModuleTree(mod[1], "right");
        	this->right = node1;
	} else {
		
	} 

	mod.restart();
}

//FacetBox ModuleTree::search(int l, int * S, ModuleTree * tree) {	
//
//	ModuleTree * user = tree;
//	ModuleTree * tmp;	
//
//	if (user->)
//
//}

string ModuleTree::getString() const {return s;}
