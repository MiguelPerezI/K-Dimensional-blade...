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
   		double teta =  0.4 * PI;
	        double phii =  0.1 * PI;
		Vector3D position1 = Vector3D(radius*sin(teta)*sin(phii), radius*cos(teta), radius*sin(teta)*cos(phii));	
		
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

     Vector3D w[20] = {
     
    Vector3D(  g2*r,    0.0*r, 1.0*r),
    Vector3D(   -g2*r,  0.0*r, 1.0*r),//1
    Vector3D(   -g1*r,   g1*r,  g1*r),//2
    Vector3D(   0.0*r,  1.0*r,  g2*r),//3
    Vector3D(    g1*r,   g1*r,  g1*r),//4
    Vector3D(   0.0*r, -1.0*r,  g2*r),//5
    Vector3D(    g1*r,  -g1*r,  g1*r),//7
    Vector3D(   -g1*r,  -g1*r,  g1*r),
    Vector3D(    g2*r,  0.0*r,-1.0*r),//9
    Vector3D(   -g2*r,  0.0*r,-1.0*r),
    Vector3D(   -g1*r,  -g1*r, -g1*r),//11
    Vector3D(   0.0*r, -1.0*r, -g2*r),
    Vector3D(    g1*r,  -g1*r, -g1*r),//13
    Vector3D(    g1*r,   g1*r, -g1*r),
    Vector3D(   0.0*r,  1.0*r, -g2*r),//15
    Vector3D(   -g1*r,   g1*r, -g1*r),
    Vector3D(   1.0*r,  -g2*r, 0.0*r),//17
    Vector3D(  -1.0*r,   g2*r, 0.0*r),
    Vector3D(  -1.0*r,  -g2*r, 0.0*r),
    Vector3D(   1.0*r,   g2*r, 0.0*r)
     };
	//QUATERNIONS
		double Radius = 0.248;
		double gCir = 0.25;
		GreatCircle circle[20];
	      	double ss = PI;
		
		double piece2 = pov.piecewiseNumber(PI, 0.0, count, 2*96);
		double piece3 = pov.piecewiseNumber(4*PI, 0.0, count, 96);
		
		Quaternion Q = Qan(0.2*PI, J);
		Vector3D head = Vector3D(0.037294, 0.0246646, 0.0246646);
		Vector3D tail = Vector3D(0, 0, 0);
		Vector3D auxx = Vector3D(0.05, 0, 0);

		Quaternion H0 = Quaternion(0, head);
		Quaternion H1 = Quaternion(0, tail);
		Quaternion H2 = Quaternion(0, auxx);

		H0 = Q.conjugate() * H0 * Q;
		H1 = Q.conjugate() * H1 * Q;
		H2 = Q.conjugate() * H2 * Q;

		//INTERSECTION
		
		Facet facet = Facet(head, tail, auxx);
		Dodecahedron D = Dodecahedron(0.025, center);
		//pov.drawDodecahedron(D, 2, 7, 0, position1);
		pov.ultra_wide_Camera(position1, D.getCenter(), 100.0);	
		
		double cutt = 0.9;




		FacetBox pila = FacetBox(D.getFacet(0));

		Intersection In = Intersection(D.getQ(0), cutt * D.getQ(0));
		In.updateOrientation();
		for (int i = 0; i < 36; i++) {
			In.intersectFacet(D.getFacet(i));
			In.readListC(D.getFacet(i), &pila);
		}

		In.readList(&pila);
		In.restart();
		
		
		int sa, sf;	
		
		In.updateIntersection(D.getQ(1), cutt * D.getQ(1));	
		In.updateOrientation();

		sa = 1;
		sf = pila.getN();
		for (int i = sa; i < sf; i++) {
                        In.intersectFacet(pila.getFacet(i));
                        In.readListC(pila.getFacet(i), &pila);
                }

		In.readList(&pila);
		In.restart();
		
		
		
		for (int j = 2; j < 20; j++) {
		In.updateIntersection(D.getQ(j), cutt * D.getQ(j));
                In.updateOrientation();

		sa = sf;
                sf = pila.getN();
                for (int i = sa; i < sf; i++) {
                        In.intersectFacet(pila.getFacet(i));
                        In.readListC(pila.getFacet(i), &pila);
                }

                In.readList(&pila);
                In.restart();
		}



      	STLWriter stl ("0ac_cut.stl"); 
		for (int i = sf; i < pila.getN(); i++) {
                        pov.drawFacet( pila.getFacet(i).getA(), pila.getFacet(i).getC(), pila.getFacet(i).getB(), 0, 0, 7, 0.0);
			stl.facetSTL(  pila.getFacet(i).getA().V(), pila.getFacet(i).getC().V(), pila.getFacet(i).getB().V());
		}











		FacetBox pila1 = FacetBox(D.getFacet(0));

                In.updateIntersection(D.getQ(0), cutt * D.getQ(0));
                for (int i = 0; i < 36; i++) {
                        In.intersectFacet(D.getFacet(i));
                        In.readListC(D.getFacet(i), &pila1);
                }

                In.readList(&pila1);
                In.restart();



                In.updateIntersection(D.getQ(1), cutt * D.getQ(1));

                sa = 1;
                sf = pila1.getN();
                for (int i = 0; i < 36; i++) {
                        In.intersectFacet(D.getFacet(i));
                        In.readListC(D.getFacet(i), &pila1);
                }

                In.readList(&pila1);
                In.restart();



                for (int j = 2; j < 20; j++) {
                In.updateIntersection(D.getQ(j), cutt * D.getQ(j));

                sa = sf;
                sf = pila1.getN();
                for (int i = 0; i < 36; i++) {
                        In.intersectFacet(D.getFacet(i));
                        In.readListC(D.getFacet(i), &pila1);
                }

                In.readList(&pila1);
                In.restart();
                }



                for (int i = 1; i < pila1.getN(); i++) {
                        

			pov.drawFacet(
					pila1.getFacet(i).getA(), 
					pila1.getFacet(i).getC(), 
					pila1.getFacet(i).getB(), 7, 0, 0, 0.0);

                        stl.facetSTL(pila1.getFacet(i).getA().V(), pila1.getFacet(i).getC().V(), pila1.getFacet(i).getB().V());
                }











		stl.closeSTLWriter();
		pov.flecha(0.5 * I, center, 0.001, 6, 0, 0);
		pov.flecha(0.5 * K, center, 0.001, 0, 6, 0);
		pov.flecha(0.5 * J, center, 0.001, 0, 0, 6);
		
		pov.closePOVRayWriter();
}


int main(int argc, char **argv) {


		srand(time(NULL));
		int numOfScenes = 0;
		int count = 80;
		int secondsScene1 = 81;
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
