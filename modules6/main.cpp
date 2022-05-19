#include <vector>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <GL/glut.h>
#include "Vector3D.cpp"
#include "Matrix3D.cpp"
#include "Vector4D.cpp"
#include "Matrix4D.cpp"
#include "Quaternion.cpp"
#include "Facet.cpp"
#include "Octahedron.cpp"
#include "PlaneQuaternion.cpp"
#include "FacetBox.cpp"
#include "FacetGash.cpp"
#include "Vector3DBox.cpp"
#include "QuaternionBox.cpp"
#include "QuaternionBoxBox.cpp"
#include "Dodecahedron.cpp"
#include "STL.cpp"
#include "Hypercube.cpp"
#include "ModuleSpaces.cpp"
#include "ModuleTree.cpp"

//////////////////////////////////////
//                                  
//                                  
//        VARIABLES GLOBALES PARA EL TECLADO        
//                                  
//                                  
//////////////////////////////////////

/*variables*/ 
int ciclo = 0;
int cicloSegund = 0;
int color = 0;
double count = 0.25 * M_PI;
double angle = 0.0;
double count2 = 0.25 * 3.14159265358979;
double count3 = 1e-2;
double rotSpeed = 0.0;
double rotAxe = 0.0;
double rad = 5.0;
double rot = 0.0;
int iter0 = 0;
int iter = 9;
int iter1 = 0;
int iter2 = 0;
int pass00 = 0;
int pass0 = 0;
int pass = 0;
int pass1 = 0;
int pass2 = 0;
int ITT = 35;
int stlP = 8;


//////////////////////////////////////
//                                  
//                                  
//        BASE EUCLIDEANA Y CENTRO        
//                                  
//                                  
//////////////////////////////////////



Vector3D origen = Vector3D(0.0, 0.0, 0.0);
Vector3D I = Vector3D(1.0, 0.0, 0.0);
Vector3D J = Vector3D(0.0, 1.0, 0.0);
Vector3D K = Vector3D(0.0, 0.0, 1.0);
int faces[1440][5];


//////////////////////////////////////
//                                  
//                                  
//        FUNCIONES PARA DIBUJAR OBJETOS EN OPENGL        
//                                  
//                                  
//////////////////////////////////////



Vector3D piecewise(double t, const Vector3D& a, const Vector3D& b);
void drawLine(const Vector3D& a, const Vector3D& b);
void drawFacet(const Facet& f, int R, int G, int B);
void drawFacet2(const Vector3D& a, const Vector3D& b, const Vector3D& c, int R, int G, int B);
void drawOctahedron(const Octahedron& octa, int R, int G, int B);
void drawPlaneQuaternion(const PlaneQuaternion& plane, int R, int G, int B);
void drawFacetBox(const FacetBox& box, int R, int G, int B);
void drawFacetBoxSTL(const FacetBox& box, string fname);
void drawFacetGash(FacetGash gash, int R, int G, int B);
void axis();


void drawHypercube(int p, const Hypercube& cube) {
	
		p = p%9;
	        if (p == 0) drawFacetBox(cube[0], 255,   0,   0);
                if (p == 1) drawFacetBox(cube[1],   0, 255,   0);
                if (p == 2) drawFacetBox(cube[2],   0,   0, 255);
                if (p == 3) drawFacetBox(cube[3], 255, 255,   0);
                if (p == 4) drawFacetBox(cube[4], 255,   0, 255);
                if (p == 5) drawFacetBox(cube[5],   0, 255, 255);
                if (p == 6) drawFacetBox(cube[6], 255,   0, 255);
                if (p == 7) drawFacetBox(cube[7], 255,   0,   0);
		if (p == 8) {
			drawFacetBox(cube[0], 255,   0,   0);
                	drawFacetBox(cube[1],   0, 255,   0);
                	drawFacetBox(cube[2],   0,   0, 255);
                	drawFacetBox(cube[3], 255, 255,   0);
                	drawFacetBox(cube[4], 255,   0, 255);
                	drawFacetBox(cube[5],   0, 255, 255);
                	drawFacetBox(cube[6], 255,   0, 255);
                	//drawFacetBox(cube[7], 255,   0,   0);
		}
}


bool ternaryCat(int i, int j, int k, int l){

                        if (i<0) i*= -1;
                        if (j<0) j*= -1;
                        if (k<0) k*= -1;
                        if (l<0) l*= -1;
                        while ((i>0 && j>0) || (j>0 && k>0) || (k>0 && i >0) ||
                               (i>0 && l>0) || (l>0 && j>0) || (l>0 && k>0)) {
                                if ((i%3==1 && j%3==1) || (j%3==1 && k%3==1) || (k%3==1 && i%3==1) ||
                                    (i%3==1 && l%3==1) || (l%3==1 && j%3==1) || (l%3==1 && k%3==1))
                                        return true;
                                i/=3;
                                j/=3;
                                k/=3;
                                l/=3;
                           }

                        return false;
        }


void drawLattice(int p, int n, int ternary, double Rs, const Matrix4D& M) {
	
	p %= n;	
	double ds = Rs/(double) n;
	double rs = 0.5 * ds;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				for (int l = p; l < p+1; l++) {
				
					int ttt = ternaryCat(i, j, k, l);
					if (ttt == ternary) {			
						Vector4D a0 = Vector4D(i * ds, j * ds, k * ds, l * ds);
						a0 = a0 - 0.5*(n-1) * Vector4D(ds, ds, ds, ds);
						Hypercube cube = Hypercube(rs, a0, M, 1.75 * Rs);
						drawHypercube(8, cube);
					}
				}
}


Vector3D spherical(double rs, double phi, double teta) {
	return rs * Vector3D(
			cos(phi) * sin(teta),
			sin(phi) * sin(teta),
			cos(teta));
}


void drawTree(int l, int * S, ModuleTree * tree, int R, int G, int B) {

        ModuleTree * user = tree;
        int i = 0;
        while (i < l) {
                if (S[i] == 0) {
				drawFacetBox(user-> left->getBox(), R, G, B);	
				user = user-> left;
		} else {
                                drawFacetBox(user->right->getBox(), R, G, B);  
                                user = user->right;
			}
                i++;
        }
}


int COUNT = 10;
int treenode = 0; 
void escTree(ModuleTree * root, int space, string s) {
	
	if (root == NULL)
		return;
	
	space += COUNT;
	escTree(root->left, space, "left");
	
	//cout << endl;
	for (int i = COUNT; i < space; i++)
		cout << " ";

	cout << root->getString() << "\n";
		treenode++;
	escTree(root->right, space, "right");
}

bool isLeaf(ModuleTree * node) {
	return (node->left == nullptr && node->right == nullptr);
}

void printRootToLeafPaths(ModuleTree * node, vector<string> &path, int io) {
    // base case
    if (node == nullptr) {
        return;
    }
 
    // include the current node to the path
    path.push_back(node->getString());
 
    // if a leaf node is found, print the path
    if (isLeaf(node)) {

	cout << "drawFacetBoxSTL(";    
        for (string data: path) {
            cout << data << "->";
        }
	cout << "getBox(), \"chok.stl\");";
        cout << endl;
    }
 
    // recur for the left and right subtree
    printRootToLeafPaths(node->left, path, io);
    printRootToLeafPaths(node->right, path, io);
 
    // backtrack: remove the current node after the left, and right subtree are done
    path.pop_back();
}

void printRootToLeafPaths(ModuleTree * node, int io) {
    // vector to store root-to-leaf path
    vector<string> path;
 
    printRootToLeafPaths(node, path, io);
}
//////////////////////////////////////
//
//
//        PIPELINE RENDERIZAR
//
//	  *Setup() := lo puedes ocupar como una función que ejecuta una sola vez (en el primer frame) dentro del render
//	  *Update() := Se puede utilizar para actualizar algún objeto
//	  *Draw() := lo puedes ocupar para todo y también dibujar solo que con cuidado por que se ejecuta en cada frame
//
//
//////////////////////////////////////



/*Funciones para dibujar sin pensar en OpenGL*/
void Setup();
void Draw();
void updateProcessingProto();
void ProcessingProto();
void interface();

//////////////////////////////////////
//                                  
//                                  
//        La variable ciclo es el número de FRAMES que lleva el sistema
//        Inicia con 0        
//                                  
//                                  
//////////////////////////////////////

Dodecahedron D = Dodecahedron(1, origen);
FacetBox box0;
ModuleSpaces mod;
double tetaa = 0.25 * M_PI;
double phii = 0.25 * M_PI;

ModuleTree * tree;

int * bi = (int *) malloc (5 * sizeof(int));
void Setup() {

	if (ciclo == 0) {




		Vector3D w[12] = {
		Vector3D(       0.0,     0.5617,  0.917587),
    		Vector3D(       0.0,    -0.5617,  0.917587),
    		Vector3D(       0.0,     0.5617, -0.917587),
    		Vector3D(       0.0,    -0.5617, -0.917587),
    		Vector3D(  0.917587,        0.0,    0.5671),
    		Vector3D( -0.917587,        0.0,    0.5671),
    		Vector3D(  0.917587,        0.0,   -0.5671),
    		Vector3D( -0.917587,        0.0,   -0.5671),
    		Vector3D(    0.5671,   0.917587,       0.0),
    		Vector3D(   -0.5671,   0.917587,       0.0),
    		Vector3D(    0.5671,  -0.917587,       0.0),
    		Vector3D(   -0.5671,  -0.917587,       0.0)
		};

		
		for (int i = 0; i < 36; i++)
			box0.push(D[i]);
	
		double gold = 0.5 * (1 + sqrt(5));
		
		/*0*/
		tree = new ModuleTree(box0, "tree");	

		/*1*/
		tree->growBranch(10.0 * w[0], (0.1*gold) * w[0]); 
		
		/*2*/
		tree->left ->growBranch(10.0 * w[1], (0.1*gold) * w[1]); 	
		tree->right ->growBranch(10.0 * w[1], (0.1*gold) * w[1]);
		
		/*3*/
		tree->left ->left ->growBranch(10.0 * w[2], (0.1*gold) * w[2]);
		tree->left ->right->growBranch(10.0 * w[2], (0.1*gold) * w[2]);
		tree->right->left ->growBranch(10.0 * w[2], (0.1*gold) * w[2]);
                tree->right->right->growBranch(10.0 * w[2], (0.1*gold) * w[2]);

		/*4*/

		
		tree->left ->left ->growBranch(10.0 * w[3], (0.1*gold) * w[3]);
		tree->left ->right->left ->growBranch(10.0 * w[3], (0.1*gold) * w[3]);
		tree->left ->right->right->growBranch(10.0 * w[3], (0.1*gold) * w[3]);
		
		tree->right->left ->growBranch(10.0 * w[3], (0.1*gold) * w[3]);
		tree->right->right->left ->growBranch(10.0 * w[3], (0.1*gold) * w[3]);
		tree->right->right->right->growBranch(10.0 * w[3], (0.1*gold) * w[3]);

	
		/*5*/
		tree->left->left->growBranch(10.0 * w[4], (0.1*gold) * w[4]);
		tree->left->right->left->growBranch(10.0 * w[4], (0.1*gold) * w[4]);
		tree->left->right->right->growBranch(10.0 * w[4], (0.1*gold) * w[4]);
		tree->right->left->left->growBranch(10.0 * w[4], (0.1*gold) * w[4]);
		tree->right->left->right->growBranch(10.0 * w[4], (0.1*gold) * w[4]);
		tree->right->right->left->left->growBranch(10.0 * w[4], (0.1*gold) * w[4]);
		tree->right->right->left->right->growBranch(10.0 * w[4], (0.1*gold) * w[4]);
		tree->right->right->right->left->growBranch(10.0 * w[4], (0.1*gold) * w[4]);
		tree->right->right->right->right->growBranch(10.0 * w[4], (0.1*gold) * w[4]);

		/*6*/
		tree->left->left->left->growBranch(10.0 * w[5], (0.1*gold) * w[5]);
		tree->left->left->right->growBranch(10.0 * w[5], (0.1*gold) * w[5]);
		tree->left->right->left->left->growBranch(10.0 * w[5], (0.1*gold) * w[5]);
		tree->left->right->left->right->growBranch(10.0 * w[5], (0.1*gold) * w[5]);
		tree->left->right->right->left->growBranch(10.0 * w[5], (0.1*gold) * w[5]);
		tree->left->right->right->right->growBranch(10.0 * w[5], (0.1*gold) * w[5]);
		tree->right->left->left->left->growBranch(10.0 * w[5], (0.1*gold) * w[5]);
		tree->right->left->left->right->growBranch(10.0 * w[5], (0.1*gold) * w[5]);
		tree->right->left->right->left->growBranch(10.0 * w[5], (0.1*gold) * w[5]);
		tree->right->left->right->right->growBranch(10.0 * w[5], (0.1*gold) * w[5]);
		tree->right->right->left->left->left->growBranch(10.0 * w[5], (0.1*gold) * w[5]);
		tree->right->right->left->left->right->growBranch(10.0 * w[5], (0.1*gold) * w[5]);
		tree->right->right->left->right->left->growBranch(10.0 * w[5], (0.1*gold) * w[5]);
		tree->right->right->left->right->right->growBranch(10.0 * w[5], (0.1*gold) * w[5]);
		tree->right->right->right->left->left->growBranch(10.0 * w[5], (0.1*gold) * w[5]);
		tree->right->right->right->left->right->growBranch(10.0 * w[5], (0.1*gold) * w[5]);
		tree->right->right->right->right->left->growBranch(10.0 * w[5], (0.1*gold) * w[5]);
		tree->right->right->right->right->right->growBranch(10.0 * w[5], (0.1*gold) * w[5]);

		/*6*/
		tree->left->left->left->left->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->left->left->left->right->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->left->left->right->left->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->left->left->right->right->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->left->right->left->left->left->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->left->right->left->left->right->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->left->right->left->right->left->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->left->right->left->right->right->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->left->right->right->left->left->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->left->right->right->left->right->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->left->right->right->right->left->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->left->right->right->right->right->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->right->left->left->left->left->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->right->left->left->left->right->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->right->left->left->right->left->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->right->left->left->right->right->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->right->left->right->left->left->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->right->left->right->left->right->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->right->left->right->right->left->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->right->left->right->right->right->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->right->right->left->left->left->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->right->right->left->left->right->left->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->right->right->left->left->right->right->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->right->right->left->right->left->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->right->right->left->right->right->left->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->right->right->left->right->right->right->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->right->right->right->left->left->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->right->right->right->left->right->left->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->right->right->right->left->right->right->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->right->right->right->right->left->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->right->right->right->right->right->left->growBranch(10.0 * w[6], (0.1*gold) * w[6]);
		tree->right->right->right->right->right->right->growBranch(10.0 * w[6], (0.1*gold) * w[6]);


		/*7*/
		tree->left->left->left->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->left->left->left->right->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->left->left->left->right->right->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->left->left->right->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->left->left->right->right->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->left->right->left->left->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->left->right->left->left->right->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->left->right->left->left->right->right->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->left->right->left->right->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->left->right->left->right->right->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->left->right->left->right->right->right->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->left->right->right->left->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->left->right->right->left->right->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->left->right->right->left->right->right->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->left->right->right->right->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->left->right->right->right->right->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->left->right->right->right->right->right->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->left->left->left->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->left->left->left->right->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->left->left->left->right->right->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->left->left->right->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->left->left->right->right->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->left->left->right->right->right->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->left->right->left->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->left->right->left->right->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->left->right->left->right->right->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->left->right->right->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->left->right->right->right->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->left->right->right->right->right->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->right->left->left->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->right->left->left->right->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->right->left->left->right->right->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->right->left->left->right->right->right->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->right->left->right->left->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->right->left->right->left->right->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->right->left->right->right->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->right->left->right->right->right->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->right->left->right->right->right->right->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->right->right->left->left->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->right->right->left->left->right->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->right->right->left->right->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->right->right->left->right->right->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->right->right->left->right->right->right->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->right->right->right->left->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->right->right->right->left->right->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->right->right->right->right->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);
		tree->right->right->right->right->right->right->left->growBranch(10.0 * w[7], (0.1*gold) * w[7]);

		/*8*/
		tree->left->left->left->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->left->left->left->right->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->left->left->left->right->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->left->left->right->left->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->left->left->right->left->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->left->left->right->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->left->right->left->left->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->left->right->left->left->right->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->left->right->left->left->right->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->left->right->left->right->left->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->left->right->left->right->left->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->left->right->left->right->right->left->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->left->right->left->right->right->left->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->left->right->left->right->right->right->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->left->right->left->right->right->right->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->left->right->right->left->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->left->right->right->left->right->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->left->right->right->left->right->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->left->right->right->right->left->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->left->right->right->right->left->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->left->right->right->right->right->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->left->right->right->right->right->right->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->left->right->right->right->right->right->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->left->left->left->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->left->left->left->right->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->left->left->left->right->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->left->left->right->left->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->left->left->right->left->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->left->left->right->right->left->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->left->left->right->right->left->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->left->left->right->right->right->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->left->left->right->right->right->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->left->right->left->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->left->right->left->right->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->left->right->left->right->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->left->right->right->left->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->left->right->right->left->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->left->right->right->right->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->left->right->right->right->right->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->left->right->right->right->right->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->right->left->left->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->right->left->left->right->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->right->left->left->right->right->left->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->right->left->left->right->right->left->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->right->left->left->right->right->right->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->right->left->left->right->right->right->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->right->left->right->left->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->right->left->right->left->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->right->left->right->right->left->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->right->left->right->right->left->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->right->left->right->right->right->left->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->right->left->right->right->right->left->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->right->left->right->right->right->right->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->right->left->right->right->right->right->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->right->right->left->left->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->right->right->left->left->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->right->right->left->right->left->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->right->right->left->right->left->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->right->right->left->right->right->left->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->right->right->left->right->right->left->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->right->right->left->right->right->right->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->right->right->left->right->right->right->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->right->right->right->left->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->right->right->right->left->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->right->right->right->right->left->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->right->right->right->right->left->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->right->right->right->right->right->left->growBranch(10.0 * w[8], (0.1*gold) * w[8]);
		tree->right->right->right->right->right->right->right->growBranch(10.0 * w[8], (0.1*gold) * w[8]);

		/*9*/
tree->left->left->left->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->left->left->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->left->left->right->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->left->left->right->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->left->left->right->right->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->left->left->right->right->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->left->right->left->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->left->right->left->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->left->right->left->right->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->left->right->left->right->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->left->right->right->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->left->right->right->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->right->left->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->right->left->left->right->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->right->left->left->right->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->right->left->right->left->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->right->left->right->left->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->right->left->right->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->right->left->right->right->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->right->left->right->right->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->right->left->right->right->right->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->right->left->right->right->right->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->right->right->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->right->right->left->right->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->right->right->left->right->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->right->right->left->right->right->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->right->right->left->right->right->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->right->right->right->left->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->right->right->right->left->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->right->right->right->left->right->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->right->right->right->left->right->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->right->right->right->right->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->right->right->right->right->right->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->right->right->right->right->right->right->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->left->right->right->right->right->right->right->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->left->left->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->left->left->left->right->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->left->left->left->right->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->left->left->left->right->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->left->left->right->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->left->left->right->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->left->left->right->right->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->left->left->right->right->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->left->left->right->right->right->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->left->left->right->right->right->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->left->right->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->left->right->left->right->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->left->right->left->right->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->left->right->left->right->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->left->right->right->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->left->right->right->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->left->right->right->right->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->left->right->right->right->right->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->left->right->right->right->right->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->left->left->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->left->left->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->left->left->right->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->left->left->right->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->left->left->right->right->left->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->left->left->right->right->left->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->left->left->right->right->left->right->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->left->left->right->right->left->right->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->left->left->right->right->right->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->left->left->right->right->right->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->left->left->right->right->right->right->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->left->left->right->right->right->right->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->left->right->left->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->left->right->left->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->left->right->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->left->right->right->left->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->left->right->right->left->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->left->right->right->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->left->right->right->right->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->left->right->right->right->left->right->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->left->right->right->right->left->right->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->left->right->right->right->right->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->left->right->right->right->right->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->left->right->right->right->right->right->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->left->right->right->right->right->right->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->right->left->left->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->right->left->left->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->right->left->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->right->left->right->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->right->left->right->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->right->left->right->right->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->right->left->right->right->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->right->left->right->right->right->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->right->left->right->right->right->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->right->right->left->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->right->right->left->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->right->right->left->right->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->right->right->left->right->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->right->right->right->left->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->right->right->right->left->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->right->right->right->left->right->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->right->right->right->left->right->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->right->right->right->right->left->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->right->right->right->right->left->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->right->right->right->right->right->left->growBranch(10.0 * w[9], (0.1*gold) * w[9]);
tree->right->right->right->right->right->right->right->right->growBranch(10.0 * w[9], (0.1*gold) * w[9]);



/*10*/

tree->left->left->left->left->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->left->left->left->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->left->left->left->right->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->left->left->left->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->left->left->right->left->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->left->left->right->left->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->left->left->right->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->left->left->right->right->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->left->left->right->right->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->left->left->right->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->left->right->left->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->left->right->left->left->right->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->left->right->left->left->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->left->right->left->right->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->left->right->left->right->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->left->right->left->right->right->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->left->right->left->right->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->left->right->right->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->left->right->right->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->left->right->right->right->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->left->right->right->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->left->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->left->left->right->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->left->left->right->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->left->left->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->left->right->left->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->left->right->left->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->left->right->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->left->right->right->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->left->right->right->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->left->right->right->right->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->left->right->right->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->right->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->right->left->right->left->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->right->left->right->left->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->right->left->right->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->right->left->right->right->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->right->left->right->right->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->right->left->right->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->right->right->left->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->right->right->left->left->right->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->right->right->left->left->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->right->right->left->right->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->right->right->left->right->right->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->right->right->left->right->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->right->right->right->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->right->right->right->right->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->right->right->right->right->right->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->right->right->right->right->right->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->right->right->right->right->right->right->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->left->right->right->right->right->right->right->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->left->left->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->left->left->left->right->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->left->left->left->right->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->left->left->left->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->left->left->right->left->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->left->left->right->left->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->left->left->right->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->left->left->right->right->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->left->left->right->right->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->left->left->right->right->right->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->left->left->right->right->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->left->right->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->left->right->left->right->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->left->right->left->right->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->left->right->left->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->left->right->right->left->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->left->right->right->left->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->left->right->right->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->left->right->right->right->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->left->right->right->right->right->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->left->right->right->right->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->left->left->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->left->left->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->left->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->left->right->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->left->right->left->right->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->left->right->left->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->left->right->right->left->left->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->left->right->right->left->left->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->left->right->right->left->left->right->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->left->right->right->left->left->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->left->right->right->left->right->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->left->right->right->left->right->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->left->right->right->left->right->right->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->left->right->right->left->right->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->left->right->right->right->left->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->left->right->right->right->left->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->left->right->right->right->left->right->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->left->right->right->right->left->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->left->right->right->right->right->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->left->right->right->right->right->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->left->right->right->right->right->right->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->left->right->right->right->right->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->right->left->left->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->right->left->left->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->right->left->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->right->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->right->right->left->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->right->right->left->left->right->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->right->right->left->left->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->right->right->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->right->right->right->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->right->right->right->left->right->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->right->right->right->left->right->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->right->right->right->left->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->right->right->right->right->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->right->right->right->right->left->right->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->right->right->right->right->left->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->right->right->right->right->right->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->right->right->right->right->right->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->right->right->right->right->right->right->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->left->right->right->right->right->right->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->left->left->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->left->left->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->left->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->left->right->left->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->left->right->left->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->left->right->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->left->right->right->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->left->right->right->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->left->right->right->right->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->left->right->right->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->right->left->left->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->right->left->left->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->right->left->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->right->left->right->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->right->left->right->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->right->left->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->right->right->left->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->right->right->left->left->right->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->right->right->left->left->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->right->right->left->right->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->right->right->left->right->right->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->right->right->left->right->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->right->right->right->left->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->right->right->right->left->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->right->right->right->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->right->right->right->right->left->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->right->right->right->right->left->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->right->right->right->right->right->left->growBranch(10.0 * w[10], (0.1*gold) * w[10]);
tree->right->right->right->right->right->right->right->right->right->growBranch(10.0 * w[10], (0.1*gold) * w[10]);


tree->left->left->left->left->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->left->left->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->left->left->left->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->left->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->left->left->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->left->left->right->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->left->right->left->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->left->right->left->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->left->right->left->left->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->left->right->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->left->right->left->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->left->right->right->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->left->right->right->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->left->right->right->left->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->left->right->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->left->right->right->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->right->left->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->right->left->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->right->left->left->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->right->left->left->right->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->right->left->right->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->right->left->right->left->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->right->left->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->right->left->right->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->right->left->right->right->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->right->right->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->right->right->left->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->right->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->right->right->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->left->right->right->right->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->left->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->left->left->right->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->left->left->right->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->left->left->right->left->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->left->left->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->left->right->left->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->left->right->left->left->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->left->right->left->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->left->right->right->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->left->right->right->left->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->left->right->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->left->right->right->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->right->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->right->left->right->left->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->right->left->right->left->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->right->left->right->left->left->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->right->left->right->left->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->right->left->right->right->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->right->left->right->right->left->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->right->left->right->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->right->right->left->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->right->right->left->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->right->right->left->left->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->right->right->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->right->right->left->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->right->right->left->right->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->right->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->right->right->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->right->right->right->right->right->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->right->right->right->right->right->left->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->right->right->right->right->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->left->right->right->right->right->right->right->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->left->left->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->left->left->left->right->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->left->left->left->right->left->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->left->left->left->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->left->left->right->left->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->left->left->right->left->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->left->left->right->left->left->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->left->left->right->left->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->left->left->right->right->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->left->left->right->right->left->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->left->left->right->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->left->left->right->right->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->left->right->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->left->right->left->right->left->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->left->right->left->right->left->left->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->left->right->left->right->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->left->right->left->right->left->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->left->right->left->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->left->right->left->right->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->left->right->right->left->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->left->right->right->left->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->left->right->right->left->left->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->left->right->right->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->left->right->right->left->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->left->right->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->left->right->right->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->left->right->right->right->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->left->right->right->right->right->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->left->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->left->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->left->left->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->left->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->right->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->right->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->right->left->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->right->left->right->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->right->right->left->left->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->right->right->left->left->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->right->right->left->left->left->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->right->right->left->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->right->right->left->left->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->right->right->left->left->right->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->right->right->left->right->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->right->right->left->right->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->right->right->left->right->left->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->right->right->left->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->right->right->left->right->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->right->right->left->right->right->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->right->right->right->left->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->right->right->right->left->left->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->right->right->right->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->right->right->right->left->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->right->right->right->left->right->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->right->right->right->right->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->right->right->right->right->left->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->right->right->right->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->right->right->right->right->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->left->right->right->right->right->right->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->right->left->left->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->right->left->left->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->right->left->left->left->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->right->left->left->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->right->left->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->right->right->left->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->right->right->left->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->right->right->left->left->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->right->right->left->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->right->right->right->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->right->right->right->left->right->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->right->right->right->left->right->left->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->right->right->right->left->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->right->right->right->right->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->right->right->right->right->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->right->right->right->right->left->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->right->right->right->right->right->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->right->right->right->right->right->left->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->right->right->right->right->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->left->right->right->right->right->right->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->left->left->left->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->left->left->left->left->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->left->left->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->left->left->left->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->left->left->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->left->right->left->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->left->right->left->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->left->right->left->left->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->left->right->left->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->left->right->right->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->left->right->right->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->left->right->right->left->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->left->right->right->right->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->left->right->right->right->left->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->left->right->right->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->left->right->right->right->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->right->left->left->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->right->left->left->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->right->left->left->left->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->right->left->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->right->left->left->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->right->left->right->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->right->left->right->left->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->right->left->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->right->left->right->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->right->right->left->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->right->right->left->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->right->right->left->left->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->right->right->left->left->right->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->right->right->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->right->right->left->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->right->right->left->right->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->right->right->left->right->right->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->right->right->right->left->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->right->right->right->left->left->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->right->right->right->left->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->right->right->right->left->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->right->right->right->right->left->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->right->right->right->right->left->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->right->right->right->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->right->right->right->right->right->right->left->growBranch(10.0 * w[11], (0.1*gold) * w[11]);
tree->right->right->right->right->right->right->right->right->right->right->growBranch(10.0 * w[11], (0.1*gold) * w[11]);


		escTree(tree, 0, "root");
		cout << "\n\n\n";
		printRootToLeafPaths(tree, 11);
		





drawFacetBoxSTL(tree->left->left->left->left->left->left->getBox(), "chok0.stl");
drawFacetBoxSTL(tree->left->left->left->left->left->right->left->getBox(), "chok1.stl");
drawFacetBoxSTL(tree->left->left->left->left->left->right->right->getBox(), "chok2.stl");
drawFacetBoxSTL(tree->left->left->left->left->right->left->left->getBox(), "chok3.stl");
drawFacetBoxSTL(tree->left->left->left->left->right->left->right->getBox(), "chok4.stl");
drawFacetBoxSTL(tree->left->left->left->left->right->right->left->left->getBox(), "chok5.stl");
drawFacetBoxSTL(tree->left->left->left->left->right->right->left->right->getBox(), "chok6.stl");
drawFacetBoxSTL(tree->left->left->left->left->right->right->right->left->getBox(), "chok7.stl");
drawFacetBoxSTL(tree->left->left->left->left->right->right->right->right->getBox(), "chok8.stl");
drawFacetBoxSTL(tree->left->left->left->right->left->left->left->getBox(), "chok9.stl");
drawFacetBoxSTL(tree->left->left->left->right->left->left->right->left->getBox(), "chok10.stl");
drawFacetBoxSTL(tree->left->left->left->right->left->left->right->right->getBox(), "chok11.stl");
drawFacetBoxSTL(tree->left->left->left->right->left->right->left->left->getBox(), "chok12.stl");
drawFacetBoxSTL(tree->left->left->left->right->left->right->left->right->getBox(), "chok13.stl");
drawFacetBoxSTL(tree->left->left->left->right->left->right->right->getBox(), "chok14.stl");
drawFacetBoxSTL(tree->left->left->left->right->right->left->left->getBox(), "chok15.stl");
drawFacetBoxSTL(tree->left->left->left->right->right->left->right->left->getBox(), "chok16.stl");
drawFacetBoxSTL(tree->left->left->left->right->right->left->right->right->getBox(), "chok17.stl");
drawFacetBoxSTL(tree->left->left->left->right->right->right->left->left->getBox(), "chok18.stl");
drawFacetBoxSTL(tree->left->left->left->right->right->right->left->right->getBox(), "chok19.stl");
drawFacetBoxSTL(tree->left->left->left->right->right->right->right->left->getBox(), "chok20.stl");
drawFacetBoxSTL(tree->left->left->left->right->right->right->right->right->getBox(), "chok21.stl");
drawFacetBoxSTL(tree->left->left->right->left->left->left->getBox(), "chok22.stl");
drawFacetBoxSTL(tree->left->left->right->left->left->right->left->left->getBox(), "chok23.stl");
drawFacetBoxSTL(tree->left->left->right->left->left->right->left->right->getBox(), "chok24.stl");
drawFacetBoxSTL(tree->left->left->right->left->left->right->right->left->getBox(), "chok25.stl");
drawFacetBoxSTL(tree->left->left->right->left->left->right->right->right->left->getBox(), "chok26.stl");
drawFacetBoxSTL(tree->left->left->right->left->left->right->right->right->right->getBox(), "chok27.stl");
drawFacetBoxSTL(tree->left->left->right->left->right->left->left->getBox(), "chok28.stl");
drawFacetBoxSTL(tree->left->left->right->left->right->left->right->getBox(), "chok29.stl");
drawFacetBoxSTL(tree->left->left->right->left->right->right->left->left->getBox(), "chok30.stl");
drawFacetBoxSTL(tree->left->left->right->left->right->right->left->right->getBox(), "chok31.stl");
drawFacetBoxSTL(tree->left->left->right->left->right->right->right->left->left->getBox(), "chok32.stl");
drawFacetBoxSTL(tree->left->left->right->left->right->right->right->left->right->getBox(), "chok33.stl");
drawFacetBoxSTL(tree->left->left->right->left->right->right->right->right->left->getBox(), "chok34.stl");
drawFacetBoxSTL(tree->left->left->right->left->right->right->right->right->right->getBox(), "chok35.stl");
drawFacetBoxSTL(tree->left->left->right->right->left->left->getBox(), "chok36.stl");
drawFacetBoxSTL(tree->left->left->right->right->left->right->getBox(), "chok37.stl");
drawFacetBoxSTL(tree->left->left->right->right->right->left->getBox(), "chok38.stl");
drawFacetBoxSTL(tree->left->left->right->right->right->right->left->left->getBox(), "chok39.stl");
drawFacetBoxSTL(tree->left->left->right->right->right->right->left->right->getBox(), "chok40.stl");
drawFacetBoxSTL(tree->left->left->right->right->right->right->right->left->getBox(), "chok41.stl");
drawFacetBoxSTL(tree->left->left->right->right->right->right->right->right->getBox(), "chok42.stl");
drawFacetBoxSTL(tree->left->right->left->left->left->getBox(), "chok43.stl");
drawFacetBoxSTL(tree->left->right->left->left->right->left->left->getBox(), "chok44.stl");
drawFacetBoxSTL(tree->left->right->left->left->right->left->right->left->getBox(), "chok45.stl");
drawFacetBoxSTL(tree->left->right->left->left->right->left->right->right->getBox(), "chok46.stl");
drawFacetBoxSTL(tree->left->right->left->left->right->right->getBox(), "chok47.stl");
drawFacetBoxSTL(tree->left->right->left->right->left->left->left->getBox(), "chok48.stl");
drawFacetBoxSTL(tree->left->right->left->right->left->left->right->left->getBox(), "chok49.stl");
drawFacetBoxSTL(tree->left->right->left->right->left->left->right->right->getBox(), "chok50.stl");
drawFacetBoxSTL(tree->left->right->left->right->left->right->getBox(), "chok51.stl");
drawFacetBoxSTL(tree->left->right->left->right->right->left->left->getBox(), "chok52.stl");
drawFacetBoxSTL(tree->left->right->left->right->right->left->right->getBox(), "chok53.stl");
drawFacetBoxSTL(tree->left->right->left->right->right->right->left->getBox(), "chok54.stl");
drawFacetBoxSTL(tree->left->right->left->right->right->right->right->getBox(), "chok55.stl");
drawFacetBoxSTL(tree->left->right->right->left->left->getBox(), "chok56.stl");
drawFacetBoxSTL(tree->left->right->right->left->right->left->left->left->getBox(), "chok57.stl");
drawFacetBoxSTL(tree->left->right->right->left->right->left->left->right->left->getBox(), "chok58.stl");
drawFacetBoxSTL(tree->left->right->right->left->right->left->left->right->right->getBox(), "chok59.stl");
drawFacetBoxSTL(tree->left->right->right->left->right->left->right->getBox(), "chok60.stl");
drawFacetBoxSTL(tree->left->right->right->left->right->right->left->left->getBox(), "chok61.stl");
drawFacetBoxSTL(tree->left->right->right->left->right->right->left->right->getBox(), "chok62.stl");
drawFacetBoxSTL(tree->left->right->right->left->right->right->right->getBox(), "chok63.stl");
drawFacetBoxSTL(tree->left->right->right->right->left->left->left->getBox(), "chok64.stl");
drawFacetBoxSTL(tree->left->right->right->right->left->left->right->left->left->getBox(), "chok65.stl");
drawFacetBoxSTL(tree->left->right->right->right->left->left->right->left->right->getBox(), "chok66.stl");
drawFacetBoxSTL(tree->left->right->right->right->left->left->right->right->getBox(), "chok67.stl");
drawFacetBoxSTL(tree->left->right->right->right->left->right->left->getBox(), "chok68.stl");
drawFacetBoxSTL(tree->left->right->right->right->left->right->right->left->getBox(), "chok69.stl");
drawFacetBoxSTL(tree->left->right->right->right->left->right->right->right->getBox(), "chok70.stl");
drawFacetBoxSTL(tree->left->right->right->right->right->left->getBox(), "chok71.stl");
drawFacetBoxSTL(tree->left->right->right->right->right->right->left->getBox(), "chok72.stl");
drawFacetBoxSTL(tree->left->right->right->right->right->right->right->left->left->getBox(), "chok73.stl");
drawFacetBoxSTL(tree->left->right->right->right->right->right->right->left->right->getBox(), "chok74.stl");
drawFacetBoxSTL(tree->left->right->right->right->right->right->right->right->left->getBox(), "chok75.stl");
drawFacetBoxSTL(tree->left->right->right->right->right->right->right->right->right->getBox(), "chok76.stl");
drawFacetBoxSTL(tree->right->left->left->left->left->getBox(), "chok77.stl");
drawFacetBoxSTL(tree->right->left->left->left->right->left->left->getBox(), "chok78.stl");
drawFacetBoxSTL(tree->right->left->left->left->right->left->right->left->getBox(), "chok79.stl");
drawFacetBoxSTL(tree->right->left->left->left->right->left->right->right->getBox(), "chok80.stl");
drawFacetBoxSTL(tree->right->left->left->left->right->right->getBox(), "chok81.stl");
drawFacetBoxSTL(tree->right->left->left->right->left->left->left->getBox(), "chok82.stl");
drawFacetBoxSTL(tree->right->left->left->right->left->left->right->left->getBox(), "chok83.stl");
drawFacetBoxSTL(tree->right->left->left->right->left->left->right->right->getBox(), "chok84.stl");
drawFacetBoxSTL(tree->right->left->left->right->left->right->getBox(), "chok85.stl");
drawFacetBoxSTL(tree->right->left->left->right->right->left->left->getBox(), "chok86.stl");
drawFacetBoxSTL(tree->right->left->left->right->right->left->right->getBox(), "chok87.stl");
drawFacetBoxSTL(tree->right->left->left->right->right->right->left->getBox(), "chok88.stl");
drawFacetBoxSTL(tree->right->left->left->right->right->right->right->getBox(), "chok89.stl");
drawFacetBoxSTL(tree->right->left->right->left->left->getBox(), "chok90.stl");
drawFacetBoxSTL(tree->right->left->right->left->right->left->left->left->getBox(), "chok91.stl");
drawFacetBoxSTL(tree->right->left->right->left->right->left->left->right->getBox(), "chok92.stl");
drawFacetBoxSTL(tree->right->left->right->left->right->left->right->left->left->getBox(), "chok93.stl");
drawFacetBoxSTL(tree->right->left->right->left->right->left->right->left->right->getBox(), "chok94.stl");
drawFacetBoxSTL(tree->right->left->right->left->right->left->right->right->getBox(), "chok95.stl");
drawFacetBoxSTL(tree->right->left->right->left->right->right->left->left->getBox(), "chok96.stl");
drawFacetBoxSTL(tree->right->left->right->left->right->right->left->right->getBox(), "chok97.stl");
drawFacetBoxSTL(tree->right->left->right->left->right->right->right->getBox(), "chok98.stl");
drawFacetBoxSTL(tree->right->left->right->right->left->left->left->left->getBox(), "chok99.stl");
drawFacetBoxSTL(tree->right->left->right->right->left->left->left->right->getBox(), "chok100.stl");
drawFacetBoxSTL(tree->right->left->right->right->left->left->right->left->getBox(), "chok101.stl");
drawFacetBoxSTL(tree->right->left->right->right->left->left->right->right->left->getBox(), "chok102.stl");
drawFacetBoxSTL(tree->right->left->right->right->left->left->right->right->right->getBox(), "chok103.stl");
drawFacetBoxSTL(tree->right->left->right->right->left->right->left->getBox(), "chok104.stl");
drawFacetBoxSTL(tree->right->left->right->right->left->right->right->left->getBox(), "chok105.stl");
drawFacetBoxSTL(tree->right->left->right->right->left->right->right->right->getBox(), "chok106.stl");
drawFacetBoxSTL(tree->right->left->right->right->right->left->getBox(), "chok107.stl");
drawFacetBoxSTL(tree->right->left->right->right->right->right->left->getBox(), "chok108.stl");
drawFacetBoxSTL(tree->right->left->right->right->right->right->right->left->left->getBox(), "chok109.stl");
drawFacetBoxSTL(tree->right->left->right->right->right->right->right->left->right->getBox(), "chok110.stl");
drawFacetBoxSTL(tree->right->left->right->right->right->right->right->right->left->getBox(), "chok111.stl");
drawFacetBoxSTL(tree->right->left->right->right->right->right->right->right->right->getBox(), "chok112.stl");
drawFacetBoxSTL(tree->right->right->left->left->left->left->left->getBox(), "chok113.stl");
drawFacetBoxSTL(tree->right->right->left->left->left->left->right->left->getBox(), "chok114.stl");
drawFacetBoxSTL(tree->right->right->left->left->left->left->right->right->getBox(), "chok115.stl");
drawFacetBoxSTL(tree->right->right->left->left->left->right->left->left->getBox(), "chok116.stl");
drawFacetBoxSTL(tree->right->right->left->left->left->right->left->right->getBox(), "chok117.stl");
drawFacetBoxSTL(tree->right->right->left->left->left->right->right->getBox(), "chok118.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->left->left->getBox(), "chok119.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->left->right->left->left->getBox(), "chok120.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->left->right->left->right->getBox(), "chok121.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->left->right->right->left->getBox(), "chok122.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->left->right->right->right->left->getBox(), "chok123.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->left->right->right->right->right->getBox(), "chok124.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->left->left->left->left->getBox(), "chok125.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->left->left->left->right->left->getBox(), "chok126.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->left->left->left->right->right->getBox(), "chok127.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->left->left->right->left->left->getBox(), "chok128.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->left->left->right->left->right->getBox(), "chok129.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->left->left->right->right->left->left->getBox(), "chok130.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->left->left->right->right->left->right->getBox(), "chok131.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->left->left->right->right->right->left->getBox(), "chok132.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->left->left->right->right->right->right->getBox(), "chok133.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->left->right->left->left->getBox(), "chok134.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->left->right->left->right->left->getBox(), "chok135.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->left->right->left->right->right->getBox(), "chok136.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->left->right->right->left->getBox(), "chok137.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->left->right->right->right->left->left->getBox(), "chok138.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->left->right->right->right->left->right->getBox(), "chok139.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->left->right->right->right->right->left->getBox(), "chok140.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->left->right->right->right->right->right->getBox(), "chok141.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->right->left->left->left->getBox(), "chok142.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->right->left->left->right->getBox(), "chok143.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->right->left->right->left->left->getBox(), "chok144.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->right->left->right->left->right->getBox(), "chok145.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->right->left->right->right->left->left->getBox(), "chok146.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->right->left->right->right->left->right->getBox(), "chok147.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->right->left->right->right->right->left->getBox(), "chok148.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->right->left->right->right->right->right->getBox(), "chok159.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->right->right->left->left->getBox(), "chok160.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->right->right->left->right->getBox(), "chok161.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->right->right->right->left->getBox(), "chok162.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->right->right->right->right->left->left->getBox(), "chok163.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->right->right->right->right->left->right->getBox(), "chok164.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->right->right->right->right->right->left->getBox(), "chok165.stl");
drawFacetBoxSTL(tree->right->right->left->left->right->right->right->right->right->right->right->right->getBox(), "chok166.stl");
drawFacetBoxSTL(tree->right->right->left->right->left->left->left->left->getBox(), "chok167.stl");
drawFacetBoxSTL(tree->right->right->left->right->left->left->left->right->left->getBox(), "chok168.stl");
drawFacetBoxSTL(tree->right->right->left->right->left->left->left->right->right->getBox(), "chok169.stl");
drawFacetBoxSTL(tree->right->right->left->right->left->left->right->getBox(), "chok170.stl");
drawFacetBoxSTL(tree->right->right->left->right->left->right->getBox(), "chok171.stl");
drawFacetBoxSTL(tree->right->right->left->right->right->left->left->left->getBox(), "chok172.stl");
drawFacetBoxSTL(tree->right->right->left->right->right->left->left->right->left->left->getBox(), "chok173.stl");
drawFacetBoxSTL(tree->right->right->left->right->right->left->left->right->left->right->getBox(), "chok174.stl");
drawFacetBoxSTL(tree->right->right->left->right->right->left->left->right->right->left->getBox(), "chok175.stl");
drawFacetBoxSTL(tree->right->right->left->right->right->left->left->right->right->right->getBox(), "chok176.stl");
drawFacetBoxSTL(tree->right->right->left->right->right->left->right->getBox(), "chok177.stl");
drawFacetBoxSTL(tree->right->right->left->right->right->right->left->left->getBox(), "chok178.stl");
drawFacetBoxSTL(tree->right->right->left->right->right->right->left->right->left->left->getBox(), "chok179.stl");
drawFacetBoxSTL(tree->right->right->left->right->right->right->left->right->left->right->getBox(), "chok180.stl");
drawFacetBoxSTL(tree->right->right->left->right->right->right->left->right->right->getBox(), "chok181.stl");
drawFacetBoxSTL(tree->right->right->left->right->right->right->right->left->left->getBox(), "chok182.stl");
drawFacetBoxSTL(tree->right->right->left->right->right->right->right->left->right->left->getBox(), "chok183.stl");
drawFacetBoxSTL(tree->right->right->left->right->right->right->right->left->right->right->getBox(), "chok184.stl");
drawFacetBoxSTL(tree->right->right->left->right->right->right->right->right->left->left->getBox(), "chok185.stl");
drawFacetBoxSTL(tree->right->right->left->right->right->right->right->right->left->right->getBox(), "chok186.stl");
drawFacetBoxSTL(tree->right->right->left->right->right->right->right->right->right->left->getBox(), "chok187.stl");
drawFacetBoxSTL(tree->right->right->left->right->right->right->right->right->right->right->getBox(), "chok188.stl");
drawFacetBoxSTL(tree->right->right->right->left->left->left->left->left->getBox(), "chok189.stl");
drawFacetBoxSTL(tree->right->right->right->left->left->left->left->right->getBox(), "chok190.stl");
drawFacetBoxSTL(tree->right->right->right->left->left->left->right->left->left->getBox(), "chok191.stl");
drawFacetBoxSTL(tree->right->right->right->left->left->left->right->left->right->getBox(), "chok192.stl");
drawFacetBoxSTL(tree->right->right->right->left->left->left->right->right->getBox(), "chok193.stl");
drawFacetBoxSTL(tree->right->right->right->left->left->right->getBox(), "chok194.stl");
drawFacetBoxSTL(tree->right->right->right->left->right->left->left->left->getBox(), "chok195.stl");
drawFacetBoxSTL(tree->right->right->right->left->right->left->left->right->left->getBox(), "chok196.stl");
drawFacetBoxSTL(tree->right->right->right->left->right->left->left->right->right->left->getBox(), "chok197.stl");
drawFacetBoxSTL(tree->right->right->right->left->right->left->left->right->right->right->getBox(), "chok198.stl");
drawFacetBoxSTL(tree->right->right->right->left->right->left->right->getBox(), "chok199.stl");
drawFacetBoxSTL(tree->right->right->right->left->right->right->left->left->getBox(), "chok200.stl");
drawFacetBoxSTL(tree->right->right->right->left->right->right->left->right->left->left->getBox(), "chok201.stl");
drawFacetBoxSTL(tree->right->right->right->left->right->right->left->right->left->right->getBox(), "chok202.stl");
drawFacetBoxSTL(tree->right->right->right->left->right->right->left->right->right->getBox(), "chok203.stl");
drawFacetBoxSTL(tree->right->right->right->left->right->right->right->left->left->getBox(), "chok204.stl");
drawFacetBoxSTL(tree->right->right->right->left->right->right->right->left->right->left->getBox(), "chok205.stl");
drawFacetBoxSTL(tree->right->right->right->left->right->right->right->left->right->right->getBox(), "chok206.stl");
drawFacetBoxSTL(tree->right->right->right->left->right->right->right->right->left->left->getBox(), "chok207.stl");
drawFacetBoxSTL(tree->right->right->right->left->right->right->right->right->left->right->getBox(), "chok208.stl");
drawFacetBoxSTL(tree->right->right->right->left->right->right->right->right->right->left->getBox(), "chok209.stl");
drawFacetBoxSTL(tree->right->right->right->left->right->right->right->right->right->right->getBox(), "chok210.stl");
drawFacetBoxSTL(tree->right->right->right->right->left->left->left->left->getBox(), "chok211.stl");
drawFacetBoxSTL(tree->right->right->right->right->left->left->left->right->left->getBox(), "chok212.stl");
drawFacetBoxSTL(tree->right->right->right->right->left->left->left->right->right->getBox(), "chok213.stl");
drawFacetBoxSTL(tree->right->right->right->right->left->left->right->left->left->getBox(), "chok214.stl");
drawFacetBoxSTL(tree->right->right->right->right->left->left->right->left->right->getBox(), "chok215.stl");
drawFacetBoxSTL(tree->right->right->right->right->left->left->right->right->getBox(), "chok216.stl");
drawFacetBoxSTL(tree->right->right->right->right->left->right->left->left->getBox(), "chok217.stl");
drawFacetBoxSTL(tree->right->right->right->right->left->right->left->right->getBox(), "chok218.stl");
drawFacetBoxSTL(tree->right->right->right->right->left->right->right->left->left->getBox(), "chok219.stl");
drawFacetBoxSTL(tree->right->right->right->right->left->right->right->left->right->getBox(), "chok220.stl");
drawFacetBoxSTL(tree->right->right->right->right->left->right->right->right->getBox(), "chok221.stl");
drawFacetBoxSTL(tree->right->right->right->right->right->left->left->left->getBox(), "chok222.stl");
drawFacetBoxSTL(tree->right->right->right->right->right->left->left->right->left->left->getBox(), "chok223.stl");
drawFacetBoxSTL(tree->right->right->right->right->right->left->left->right->left->right->getBox(), "chok224.stl");
drawFacetBoxSTL(tree->right->right->right->right->right->left->left->right->right->left->getBox(), "chok225.stl");
drawFacetBoxSTL(tree->right->right->right->right->right->left->left->right->right->right->left->getBox(), "chok226.stl");
drawFacetBoxSTL(tree->right->right->right->right->right->left->left->right->right->right->right->getBox(), "chok227.stl");
drawFacetBoxSTL(tree->right->right->right->right->right->left->right->left->getBox(), "chok228.stl");
drawFacetBoxSTL(tree->right->right->right->right->right->left->right->right->left->getBox(), "chok229.stl");
drawFacetBoxSTL(tree->right->right->right->right->right->left->right->right->right->left->getBox(), "chok230.stl");
drawFacetBoxSTL(tree->right->right->right->right->right->left->right->right->right->right->left->getBox(), "chok231.stl");
drawFacetBoxSTL(tree->right->right->right->right->right->left->right->right->right->right->right->getBox(), "chok232.stl");
drawFacetBoxSTL(tree->right->right->right->right->right->right->left->left->left->getBox(), "chok233.stl");
drawFacetBoxSTL(tree->right->right->right->right->right->right->left->left->right->getBox(), "chok234.stl");
drawFacetBoxSTL(tree->right->right->right->right->right->right->left->right->left->left->getBox(), "chok235.stl");
drawFacetBoxSTL(tree->right->right->right->right->right->right->left->right->left->right->getBox(), "chok236.stl");
drawFacetBoxSTL(tree->right->right->right->right->right->right->left->right->right->getBox(), "chok237.stl");
drawFacetBoxSTL(tree->right->right->right->right->right->right->right->left->left->getBox(), "chok238.stl");
drawFacetBoxSTL(tree->right->right->right->right->right->right->right->left->right->getBox(), "chok239.stl");
drawFacetBoxSTL(tree->right->right->right->right->right->right->right->right->left->getBox(), "chok240.stl");
drawFacetBoxSTL(tree->right->right->right->right->right->right->right->right->right->left->left->getBox(), "chok241.stl");
drawFacetBoxSTL(tree->right->right->right->right->right->right->right->right->right->left->right->getBox(), "chok242.stl");
drawFacetBoxSTL(tree->right->right->right->right->right->right->right->right->right->right->left->getBox(), "chok243.stl");
drawFacetBoxSTL(tree->right->right->right->right->right->right->right->right->right->right->right->getBox(), "chok245.stl");

	}

}

///////////////////     DRAW       ///////////////////////
void Draw() {

	if (ciclo > 0) {
	

		axis();
		drawFacetBox(tree->right->right->getBox(), 255, 0, 0);	
	}
}


void ProcessingProto() {

	Setup();
	Draw();
}









Vector3D piecewise(double t, const Vector3D& a, const Vector3D& b) {

        return (t*(a-b)) + b;
}

void drawLine(const Vector3D& a, const Vector3D& b) {

        glColor3ub(0, 0, 0);
	glLineWidth(2.0);
        glBegin(GL_LINES);
        glVertex3f(a.x(), a.y(), a.z());
        glVertex3f(b.x(), b.y(), b.z());
        glEnd();
}

void axis() {

	drawLine(Vector3D(20,  0,  0), origen);
	drawLine(Vector3D( 0, 20,  0), origen);
	drawLine(Vector3D( 0,  0, 20), origen);
}

void drawFacet(const Facet& f, int R, int G, int B) {

        Vector3D n = f[3];
        Vector3D a = f[0];
        Vector3D b = f[1];
        Vector3D c = f[2];
	
	glColor4f(R, G, B, 0.75 );
	glBegin(GL_TRIANGLES);
        glNormal3f( n.x(), n.y(), n.z());
        glVertex3f( a.x(), a.y(), a.z());
        glVertex3f( b.x(), b.y(), b.z());
        glVertex3f( c.x(), c.y(), c.z());
        glEnd();

	drawLine(a, b);
	drawLine(b, c);
}

void drawFacet2(const Vector3D& a, const Vector3D& b, const Vector3D& c, int R, int G, int B) {

        Vector3D n = unit( (b-a) % (c-a));
        glColor3ub(R, G, B);
        glBegin(GL_TRIANGLES);
        glNormal3f( n.x(), n.y(), n.z());
        glVertex3f( a.x(), a.y(), a.z());
        glVertex3f( b.x(), b.y(), b.z());
        glVertex3f( c.x(), c.y(), c.z());
        glEnd();
}

void drawOctahedron(const Octahedron& octa, int R, int G, int B) {

        for (int i = 0; i < 8; i++)
                drawFacet(octa[i], R, G, B);
}

void drawPlaneQuaternion(const PlaneQuaternion& plane, int R, int G, int B) {

        double l = abs(plane[2]-plane[3]);
        drawOctahedron(Octahedron(0.05 * l, plane[0]), R, G, B);
        drawLine(plane[0], plane[1]);
        drawFacet2(piecewise(2.0, plane[2], plane[1]), piecewise(2.0, plane[3], plane[1]), piecewise(2.0, plane[4], plane[1]), R, G, B);
        drawFacet2(piecewise(2.0, plane[4], plane[1]), piecewise(2.0, plane[5], plane[1]), piecewise(2.0, plane[2], plane[1]), R, G, B);
}

void drawFacetBox(const FacetBox& box, int R, int G, int B) {


        if (box.getN() != 0)
                for (int i = 0; i < box.getN(); i++) {
                        drawFacet(box[i], R, G, B);
                }
}

void drawFacetBoxSTL(const FacetBox& box, string fname) {

        STL stl = STL(fname);

        if (box.getN() != 0)
                for (int i = 0; i < box.getN(); i++) {
                        stl.facet(box[i][0], box[i][1], box[i][2]);
                }
        stl.closeSTL();
}

void drawFacetGash(FacetGash gash, int R, int G, int B) {

        drawPlaneQuaternion(gash.getBlade(), 50, 100, 255);
        //drawFacetBox(gash.getFacets(), R, G, B);

        for (int i = 0; i < gash.getMM().getM(); i++)
                for (int j = 0; j < gash.getMM()[i].getN(); j++)
                        drawOctahedron(Octahedron(0.1, gash.getMM()[i][j]), 255, 0, 0);


}






/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/**/

//////////////////////////////////////
//                                  //
//                                  //
//       OPENGL AS BACKGROUND       //
//                                  //
//                                  //
//////////////////////////////////////

/*Posición y color de luz*/
GLfloat light_diffuse[] = {1.0, 1.0, 1.0, 1.0};
GLfloat light_position[] = {20.0, 20.0, 20.25, 0.0};

/*Funciones de OpenGL*/
void display(void);
void init(double theta);
void TimerFunction(int value);
void keyboard(unsigned char key, int x, int y);
void ProcessMenu(int value);

int main(int argc, char **argv)
{
 	srand(time(NULL)); 
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(720, 720);
	glutCreateWindow(" JAZ 4D	U.U ");
	ProcessMenu(1);
	init(count);

	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutTimerFunc(20, TimerFunction, 1);

	glutMainLoop();
	return 0;             /* ANSI C requires main to return int. */
}

void display(void) {
  
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(1.0, 1.0, 1.0, 1.0);
	ProcessingProto();
	glutSwapBuffers();
}

void init(double theta)
{
  /* Setup data. */
	GLfloat ambientLight[] = { 0.3f, 0.3f, 0.3f, 1.0f };
	GLfloat diffuseLight[] = { 0.7f, 0.7f, 0.7f, 1.0f };
	GLfloat specular[] = { 1.0f, 1.0f, 1.0f, 1.0f};
	GLfloat specref[] = { 1.0f, 1.0f, 1.0f, 0.25f };

  /* Enable a single OpenGL light. */
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);

  /* Use depth buffering for hidden surface elimination. */
	glEnable(GL_DEPTH_TEST);
	glFrontFace(GL_CCW);
  	//glEnable(GL_CULL_FACE);

  /*Enable color tracking*/
	glEnable(GL_COLOR_MATERIAL);

  /* Set material properties to follow glColor values*/
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

  /*All materials have high shine*/
	glMaterialfv(GL_FRONT, GL_AMBIENT, specref);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, specref);
	
	glMateriali(GL_FRONT, GL_SHININESS, 128);
	

  /* Setup the view of the cube. */
	glMatrixMode(GL_PROJECTION);
	gluPerspective( /* field of view in degree */ 50.0,
			/* aspect ratio */ 1.0,
			/* Z near */ 0.5, 
			/* Z far */ 10000.0);
	glMatrixMode(GL_MODELVIEW);
  /* Adjust Board position to be asthetic angle. */
  //glTranslatef(0.0, 0.15, -0.0);
	glRotatef(90, 0.0, 0.0, 1.0);

	glEnable(GL_NORMALIZE);
}

void TimerFunction(int value) {

	count += 0.0;
	rotSpeed += 0.00;
	ciclo += 1;
	angle += 0.006283;

	if (count > 2 * M_PI) count = 0;
	if (ciclo > 100) ciclo = 1; //CICLO NUNCA ES CERO JAJAJA
	if (angle > 2 * M_PI) angle = 0;
	
	glLoadIdentity();



	///////////////////////////////
	//
	//
	//	CAMERA CONTROL
	//
	//
	//////////////////////////////


	gluLookAt( 
		rad * cos(count)*sin(count2), rad * sin(count) * sin(count2), rad * cos(count2),      	/* eye is at (0,0,5) */
              	0.0, 0.0, 1.0,      									/* center is at (0,0,0) */
              	0.0, 0.0, 1.0);      									/* up is in positive Y direction */

	glutPostRedisplay();
	glutTimerFunc(20, TimerFunction, 1);
}

void keyboard(unsigned char key, int x, int y) {
	GLint params[2];

	switch (key) {

		case 'b': 
			rotSpeed += 0.05;
      		break;

    		
		case 'B':
      			rotSpeed -= 0.05;
      		break;

    
		case 'C':
      			count3 += 0.0001;
      		break;

    		
		case 'c':
      			count3 -= 0.0001;
      		break;

    		
		case 'E':
      			count2 += 0.05;
      		break;

    
		case 'e':
      			count2 -= 0.05;
      		break;

    
		case 'r':
      			count += 0.05;
      		break;

    		
		case 'R':
      			count -= 0.05;
      		break;


    		case 'm':
			rotAxe += 0.05;
      		break;

    		
		case 'M':
      			rotAxe -= 0.05;
      		break;

    		
		case 'f':
      			rad += 0.2;
      		break;

    
		case 'F':
      			rad -= 0.2;
      		break;

    
		case 'v':
      			rot += 0.02;
      		break;

    
		case 'V':
      			rot -= 0.02;
      		break;

     
		case 'p':
      			stlP += 1;
		break;

		case 'P':
                        stlP -= 1;
                break;

		case 'i':
                        iter = 9;
                break;

                case 'I':
                        iter = 27;
                break;

	}

	glutPostRedisplay();
}

void ProcessMenu(int value) {
	switch(value) {
    		case 1:
      			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      			glEnable(GL_BLEND);
      			glEnable(GL_POINT_SMOOTH);
      			glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
      			glEnable(GL_LINE_SMOOTH);
      			glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
      			glEnable(GL_POLYGON_SMOOTH);
      			glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
      		break;

    		case 2:
      			glDisable(GL_BLEND);
      			glDisable(GL_LINE_SMOOTH);
      			glDisable(GL_POINT_SMOOTH);
      			glDisable(GL_POLYGON_SMOOTH);
      		break;
    
    		default:
      		break;
  	}

	glutPostRedisplay();
}
