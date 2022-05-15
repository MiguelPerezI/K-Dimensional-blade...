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
//////////////////////////////////////
//                                  //
//                                  //
//        VARIABLES GLOBALES        //
//                                  //
//                                  //
//////////////////////////////////////

/*
g++ ant.cpp Arrow.cpp dodecahedron.cpp Extra_Operators.hpp geometry.cpp MainOpenGL.cpp matrix.cpp PovRayWriter.cpp simplex.cpp Turtle.cpp VectorND.cpp -lm -lGL -lGLU -lglut

*/

/*variables*/ 
int ciclo = 0;
int cicloSegund = 0;
int color = 0;
double count = 0.0;
double angle = 0.0;
double count2 = 0.5 * 3.14159265358979;
double count3 = 0.25 * M_PI;
double rotSpeed = 0.0;
double rotAxe = 0.0;
double rad = 10.0;
double rot = 0.0;

int iter0 = 0;
int iter = 0;
int iter1 = 0;
int iter2 = 0;

int pass00 = 0;
int pass0 = 0;
int pass = 0;
int pass1 = 0;
int pass2 = 0;
int ITT = 35;
Vector3D origen = Vector3D(0.0, 0.0, 0.0);
Vector3D I = Vector3D(1.0, 0.0, 0.0);
Vector3D J = Vector3D(0.0, 1.0, 0.0);
Vector3D K = Vector3D(0.0, 0.0, 1.0);
int faces[1440][5];

Vector3D piecewise(double t, const Vector3D& a, const Vector3D& b) {
	
	return (t*(a-b)) + b;
}

void drawLine(const Vector3D& a, const Vector3D& b) {

	glColor3ub(0, 0, 0);
        glBegin(GL_LINES);
        glVertex3f(a.x(), a.y(), a.z());
        glVertex3f(b.x(), b.y(), b.z());
        glEnd();
}

void drawFacet(const Facet& f, int R, int G, int B) {

        Vector3D n = f[3];
	Vector3D a = f[0];
	Vector3D b = f[1];
	Vector3D c = f[2];
        glColor3ub(R, G, B);
        glBegin(GL_TRIANGLES);
        glNormal3f( n.x(), n.y(), n.z());
        glVertex3f( a.x(), a.y(), a.z());
        glVertex3f( b.x(), b.y(), b.z());
        glVertex3f( c.x(), c.y(), c.z());
        glEnd();
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
		for (int i = 0; i < box.getN(); i++)
			drawFacet(box[i], B, G, B);

	//drawOctahedron(Octahedron(0.1 * abs(box.getCenter()-box[0][0]), box.getCenter()), R, G, B);
}

void drawFacetGash(FacetGash gash, int R, int G, int B) {
	
	drawPlaneQuaternion(gash.getBlade(), 50, 100, 255);
	drawFacetBox(gash.getFacets(), R, G, B);
	
	for (int i = 0; i < gash.getMM().getM(); i++)
		for (int j = 0; j < gash.getMM()[i].getN(); j++)
			drawOctahedron(Octahedron(0.1, gash.getMM()[i][j]), 255, 0, 0);
	

}

/*Funciones para dibujar sin pensar en OpenGL*/
void Setup();
void Draw();
void updateProcessingProto();
void ProcessingProto();
void interface();



//////////////////////////////////////
//                                  //
//                                  //
//        Processing Prototype      //
//                                  //
//                                  //
//////////////////////////////////////

/*Here we build our memory space and filled it with data using initObject methods corresponding to each class.*/
/*initObjects methods are functions that should build memory space and fill it with data*/

Octahedron octa = Octahedron(1.0, origen);
double phii = 0.5 * M_PI;
double tetaa = 0.25 * M_PI;

//FacetGash In = FacetGash(
//		2*Vector3D(cos(phii) * sin(tetaa), sin(phii) * sin(tetaa), cos(tetaa)),
//		Vector3D(0.33, 0.33, 0.33));
//
///////////////////     SETUP       ///////////////////////
void Setup() {

  if (ciclo == 0) {
	
  	//for (int i = 0; i < 8; i++) { 
	//	In.cutFacet(octa[i]);
	//	//In.readListC();
	//}

	//	
	////In.restart();
	//cout << In;
	////In.updateOrientation();
  }
}

//////////////////    UPDATE AUXILIARY FUNCTION ///////////////////

/*In this function we call any method that updates an object in a class.*/
/*Our goal is to define our memory space with initial values.*/
/*Having a memory space filled with initial values we are now able to update these initial values.*/

void updateProcessingProto() {

	if (ciclo > 0) {
	
		/*For example here we are updating our matrix rotation system.*/


	}
}

///////////////////     DRAW       ///////////////////////
FacetBox pila;
/*Everything is made up of triangles and each class of geometrical objects have triangle drawing methods.*/
void Draw() {

  if (ciclo > 0) {

    	
	/*Draw Here*/

	/*Prepare Blade*/
  	FacetGash In = FacetGash(
              2 *   Vector3D(cos(count3) * sin(tetaa), sin(count3) * sin(tetaa), cos(tetaa)),
              rot * Vector3D(cos(count3) * sin(tetaa), sin(count3) * sin(tetaa), cos(tetaa))
	);
	  
	/*Divide the polytope in two*/
	for (int i = 0; i < 8; i++) {
                In.cutFacet(octa[i]);
                //In.readListC(octa[i], &pila);
        	int aa = In.checkFacet(octa[i]);
		if (aa == 1)
			pila.push(octa[i]);
	}

	/*Draw one of the first half not touching the boundary*/
	drawFacetBox(pila, 0, 255, 0);

	drawFacetGash(In, 255, 0, 255);
	
      	for (int i = 0; i < 8; i++) {
	
		for (int j = 0; j < 3; j++) {
		
			if (In.checkPoint(Quaternion(octa[i][j]), Quaternion(K)) == 1) 
				drawOctahedron(Octahedron(0.1, octa[i][j]), 0, 0, 255);
		}
	}

	In.restart();
	pila.empty();
  
  }
}


void ProcessingProto() {

  Setup();
  updateProcessingProto();
  Draw();
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
GLfloat light_position[] = {1.0, 1.0, 0.25, 0.0};

/*Funciones de OpenGL*/
void display(void);
void init(double theta);
void TimerFunction(int value);
void keyboard(unsigned char key, int x, int y);
void ProcessMenu(int value);

int main(int argc, char **argv)
{
  
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(1080, 720);
  glutCreateWindow(" ------- 120 - cell ------- ");
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
  GLfloat specref[] = { 1.0f, 1.0f, 1.0f, 1.0f };

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
  glMaterialfv(GL_FRONT, GL_SPECULAR, specref);
  glMateriali(GL_FRONT, GL_SHININESS, 128);

  /* Setup the view of the cube. */
  glMatrixMode(GL_PROJECTION);
  gluPerspective( /* field of view in degree */ 40.0,
                              /* aspect ratio */ 1.0,
                                    /* Z near */ 0.5, 
                                    /* Z far */ 10000.0);
  glMatrixMode(GL_MODELVIEW);
//  gluLookAt( 4.01, 4.01, 9.0,      /* eye is at (0,0,5) */
//              0.0, 0.0, 1.0,      /* center is at (0,0,0) */
//             0.0, 0.0, 1.0);      /* up is in positive Y direction */
//
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
  if (ciclo > 100) ciclo = 1;
  if (angle > 2 * M_PI) angle = 0;
	
  glLoadIdentity();
  gluLookAt( rad * cos(count)*sin(count2), rad * sin(count) * sin(count2), rad * cos(count2),      /* eye is at (0,0,5) */
              0.0, 0.0, 1.0,      /* center is at (0,0,0) */
              0.0, 0.0, 1.0);      /* up is in positive Y direction */

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
      count3 += 0.05;
      break;

    case 'c':
      count3 -= 0.05;
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
      rad += 0.05;
      break;

    case 'F':
      rad -= 0.05;
      break;

    case 'v':
      rot += 0.005;
      break;

    case 'V':
      rot -= 0.005;
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
