//                              _____ _       _     _   
// ___  ___ _ ____   _____ _ __| ____(_) __ _| |__ | |_ 
/// __|/ _ \ '__\ \ / / _ \ '__|  _| | |/ _` | '_ \| __|
//\__ \  __/ |   \ V /  __/ |  | |___| | (_| | | | | |_ 
//|___/\___|_|    \_/ \___|_|  |_____|_|\__, |_| |_|\__|
//                                      |___/           

#include <cstdio>
#include <stdio.h>
#include <string>
#include <iostream>
#include <cmath>
#include "PovRayWriter.cpp"
#include "Vector3D.cpp"
#include "Matrix3D.cpp"
#include "Vector4D.cpp"
#include "Quaternion.cpp"
#include "Matrix4D.cpp"
#include <mpi.h>
#include "STL.cpp"
//#include "Intersection.cpp"

using namespace std;

Vector3D Origin =  Vector3D( 0.0, 0.0, 0.0 );
Vector3D center = Vector3D(0, 0, 0);
Vector3D I = Vector3D(1, 0, 0);
Vector3D J = Vector3D(0, 1, 0);
Vector3D K = Vector3D(0, 0, 1);
double PI = 3.1415926535;
void vectorPrint(Vector3D vv);

Vector3D line(double t, const Vector3D& a, const Vector3D& b) {

	return (t*(a-b)) + b;
}

void drawModuleTree(ModuleTree * root, POVRayWriter * pov, int R, int G, int B, int count) {
	
	int r0, g0, b0, r1, g1, b1;
	r0 = R%8; r1 = (R+4)%8;
	g0 = G%8; g1 = (G+4)%8;
	b0 = B%8; b1 = (B+4)%8;
	
	string TITLE1 = "0ac_cut0" + to_string(count) + ".stl";
        STLWriter stl (TITLE1);	
	for (int i = 0; i < root->getN0(); i++) {

            		Quaternion Aa = root->getFacet0(i).getA();
                        Quaternion Ca = root->getFacet0(i).getC();
                        Quaternion Ba = root->getFacet0(i).getB();
			

                      pov->drawFacet(Aa, Ca, Ba, r0, g0, b0, 0.0);
                      stl.facetSTL(Aa.V(), Ca.V(), Ba.V());

       }

       stl.closeSTLWriter();


       string TITLE2 = "0ac_cut0" + to_string(count+1) + ".stl";
        STLWriter stl0 (TITLE2);
       for (int i = 0; i < root->getN1(); i++) {

                        Quaternion Aa = root->getFacet1(i).getA();
                        Quaternion Ca = root->getFacet1(i).getC();
                        Quaternion Ba = root->getFacet1(i).getB();
			

                        pov->drawFacet(Aa, Ca, Ba, r1, g1, b1, 0.0);
                        stl0.facetSTL(Aa.V(), Ca.V(), Ba.V());

      }

       stl0.closeSTLWriter();
}

void render(int count, int core) {

	//FILE INFO
		string TITLE = "hyperbolic0" + to_string(count) + ".pov";
		//cout << "--->Making Picture " << TITLE << endl; 
   		POVRayWriter pov (TITLE);
	
	//Floor	
		double elevation = -1.0;
		double trans = 0.0;
		double iC = (double) count;
		pov.TILEFLOOR3(160, 160, 2.0, elevation, trans);
	//Ambient Light
		double inten0 = 0.5;//pov.piecewiseNumber(0.0, 1.0, iC, 160);
		pov.ambientLight(inten0);
		double piece;
		
		cout << " Core := " << core << "	photo(" << count << ")" << endl;
		
	//CHAMBER/.
   		//pov.cube(10.0, center, 0.2, 0.0, 0, center, 1.0, 7, 7, 7, 0.75);
		
		Vector3D centerSphere = Vector3D(0.0, 0.0, 0.0);
		double transparency = 0.0;	
		
	//Camera Position in Spherical Coordinates
		double radius = 0.08;//2.0 * pov.piecewiseNumber(0.25, 1.0, 96, 96);
		double angle = pov.piecewiseNumber(2*PI, 0.0, count, 5*96);
   		double teta =  0.2 * PI;
	        double phii =  (0.1 * PI);
		Vector3D position1 = Vector3D(radius*sin(teta)*cos(phii), radius*cos(teta), radius*sin(teta)*sin(phii));	
		
		radius = 0.6 * pov.piecewiseNumber(0.75, 0.25, 96, 96);
		Vector3D position2 = Vector3D(radius*sin(teta)*sin(phii), radius*cos(teta), radius*sin(teta)*cos(phii));
	//Light Positioning
		double inten = 0.5;
		pov.createLight(position1, inten, inten, inten);
		//pov.spotLight(position1, center);
		Vector3D auxVec = Vector3D(0.0, 0.25, 0.0);
		Vector3D vv = Vector3D(0.0, 0.5, 0.0);
		double r = 1.0;
		int n = 4;

		double ana = pov.piecewiseNumber(2*PI, 0.0, iC+160, 160);
		Vector3D centerBall = Vector3D(0.0, 0.2, 0.0);
	//TEXT
	//
		//pov.texReaderBeta("0101Ram.txt",  position1, Vector3D( 0.3, 0.3, 0), "Blue",
                //                96-iC, 96, 0.00064, 0.001, center + Vector3D(0, 0.2, 0), 7, 0, 0);
		

	//dodecahedron
	double gold= (1 + sqrt(5))/2;
     double g1= 1/gold;
     double g2= 1/(gold*gold);
     double s = 0.8;
     r = 0.25;

		pov.ultra_wide_Camera(position1, center, 100.0);	
		
		//Dodecahedron definition
		Dodecahedron D = Dodecahedron(0.025, center);
		
		//Rotate dodecahedron
		double ppso = pov.piecewiseNumber(2 * PI, 0.0, 1, 200);
		D.rotate(ppso, Vector3D(1, 1, 1));
		pov.drawDodecahedron(D, 2, 7, 0, position1);
		

		//Build Facet list
		FacetBox pila = FacetBox(D.getFacet(0));
		for (int i = 0; i < 36; i++)
			pila.pushFacet(D.getFacet(i));
	

		//Define Tree
		ModuleTree * root = NULL;
		root = root->newModuleTree(D.getQ(14), D.getQ(3), (1.0/gold), pila);

		root->updateLeft (D.getQ(4), D.getQ(0), (1.0/gold), root->module2.getFacetBox0());
		root->updateRight(D.getQ(4), D.getQ(0), (1.0/gold), root->module2.getFacetBox1());

		root->left ->updateLeft (D.getQ(19), D.getQ(16), 0.0, root-> left->module2.getFacetBox0());	
		root->left ->updateRight(D.getQ(19), D.getQ(16), 0.0, root-> left->module2.getFacetBox1());

		root->right->updateLeft (D.getQ(19), D.getQ(16), 0.0, root-> right->module2.getFacetBox0());
                root->right->updateRight(D.getQ(19), D.getQ(16), 0.0, root-> right->module2.getFacetBox1());


		//drawModuleTree(root, &pov, &stl, 0, 0, 7);
		cout << "\nleft////////////////\n" << abs(D.getQ(0));
		drawModuleTree(root->left->left,  &pov,  0, 0, 7, count);
		
		cout << "\nRight////////////////\n";
		drawModuleTree(root->left->right, &pov, 7, 0, 0, count+2);

		drawModuleTree(root->right->left, &pov, 7, 7, 0, count+4);
		drawModuleTree(root->right->right, &pov, 0, 7, 0, count+6);












		pov.flecha(0.5 * I, center, 0.001, 6, 0, 0);
		pov.flecha(0.5 * K, center, 0.001, 0, 6, 0);
		//pov.flecha(0.5 * J, center, 0.001, 0, 0, 6);
		
		pov.closePOVRayWriter();
}


int main(int argc, char **argv) {


		srand(time(NULL));
		int numOfScenes = 0;
		int count = 1;
		int secondsScene1 = 2;
		int thread, core, coreS;
		thread = MPI_Init(&argc, &argv);
                thread = MPI_Comm_rank(MPI_COMM_WORLD, &core);
                thread = MPI_Comm_size(MPI_COMM_WORLD, &coreS);
		
	while (count < secondsScene1) {
		
		if (core == 0)	render(count + 0, core);
		if (core == 1)	render(count + 1, core);
		if (core == 2)	render(count + 2, core);
		if (core == 3)	render(count + 3, core);
		

		count += 4;
	}










      	cout << "\n\n\n"; 
	thread = MPI_Finalize();
	return 0;
}






void vectorPrint(Vector3D vv) {
  cout << "\n\n(" << vv.x() << ", " << vv.y() << ", " << vv.z() << ")\n\n";
}

//mpic++ tutorial0.cpp -o go -pthread; time mpirun -np 4 ./go file.txt; for x in hyperbolic0*.pov; do povray $x -D +W600 +H450; done
//mpic++ tutorial0.cpp -o go -pthread; time mpirun -np 1 ./go file.txt; povray hyperbolic00.pov +W600 +H450; eog hyperbolic00.eog
//ffmpeg -r 32 -f image2 -s 1080x720 -i hyperbolic0%1d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p ball0.mp4
//ffmpeg -i escena3.mp4 -ss 02:20 -to 02:40 -c:v libx264 -crf 32 escena3edit.mp4
//for x in hyperbolic0*.pov; do povray $x +W540 +H720 -D +Q11; done; eog hyperbolic00.png
// mpic++ Algebra7.cpp -o go -pthread; time mpirun -np 1 ./go file.txt; povray hyperbolic00.pov +W600 +H550 -D +Q11;  eog hyperbolic00.png
//povray hyperbolic00.pov +W540 +H720 -D +Q11; eog hyperbolic00.png
//mpic++ Algebra7.cpp -o go -pthread; time mpirun -np 4 ./go file.txt; for x in hyperbolic0*.pov; do povray $x -D +W300 +H300 +Q6; done; eog hyperbolic00.png
