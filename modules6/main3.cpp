#include <vector>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <GL/glut.h>
#include "Vector3D.cpp"
#include "Quaternion.cpp"

//////////////////////////////////////
//                                  
//                                  
//  VARIABLES GLOBALES PARA EL TECLADO        
//                                  
//                                  
//////////////////////////////////////

// ciclo: Frame number.
int ciclo = 0;

// count: Is used to measure rotating angle length.
double count = 0.25 * M_PI;
double count2 = 0.25 * 3.14159265358979;
double count3 = 1e-2;

// rad: Radius between camaera and center.
double rad = 5.0;

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

//        La variable ciclo es el número de FRAMES que lleva el sistema
//        Inicia con 0        


void Setup() {

	if (ciclo == 0) {
        
        cout << "\n\n";
        cout << "|———————————————————————————————————|\n";
        cout << "|          _                        |\n"; 
        cout << "|          \\`*-.                    |\n"; 
        cout << "|           )  _`-.                 |\n"; 
        cout << "|          .  : `. .                |\n"; 
        cout << "|          : _   '  \\               |\n"; 
        cout << "|          ; *` _.   `*-._          |\n"; 
        cout << "|          `-.-'          `-.       |\n"; 
        cout << "|            ;       `       `.     |\n"; 
        cout << "|            :.       .        \\    |\n"; 
        cout << "|            . \\  .   :   .-'   .   |\n"; 
        cout << "|            '  `+.;  ;  '      :   |\n"; 
        cout << "|            :  '  |    ;       ;-. |\n"; 
        cout << "|            ; '   : :`-:     _.`* ;|\n";
        cout << "|[K-Blade] .*' /  .*' ; .*`- +'  `*'|\n";
        cout << "|         `*-*   `*-*  `*-*'        |\n";
        cout << "|———————————————————————————————————|\n";
        cout << "\n\n";
        cout << "———————————————————————————————————————————————————————————————————————\n";
        cout << "|- Testing Vector3D Class——————————————————————————————————————————————\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n";
        
        // --- Construction ---------------------------------------------------
        cout << "- List of constructors ————————————————————————————————————————————————\n";
        Vector3D p;                       // Default constructor
        cout << "Default constructor:       Vector3D p;              → p = " << p << "\n";
        
        Vector3D q(1.0, 2.0, 3.0);        // Constructor with components
        cout << "Component constructor:     Vector3D q(1.0, 2.0, 3.0) → q = " << q << "\n";
        
        Vector3D r = q;                   // Copy constructor
        cout << "Copy constructor:          Vector3D r = q;           → r = " << r << "\n";
        
        Vector3D r1{0.1, 0.1, 0.12};      // Brace-initialized constructor
        cout << "Brace constructor:         Vector3D r1{0.1,0.1,0.12}          → r1 = " << r1 << "\n";
       
        Vector3D sup = Vector3D(r1);  // Copy constructor (explicit expression)
        cout << "Copy constructor (explicit):   Vector3D sup = Vector3D(r1)  → sup = " << sup << "\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n"; 

        // --- Access & Modification -----------------------------------------
        cout << "- Access & Modification Examples ——————————————————————————————————————\n";
        cout << "";
        q[0] = 4.0;                       // Write to q[0]
        double zy = q[2];                 // Read from q[2]
        cout << "After q[0] = 4.0:          q = " << q << "\n";
        cout << "Reading q[2]:              q[2] = " << zy << "\n";
        
        r += Vector3D(2, 2, 2);           // Compound addition
        cout << "After r += Vector3D(2,2,2): r = " << r << "\n";
        
        r -= Vector3D(0.5, 0.5, 0.5);     // Compound subtraction
        cout << "After r -= Vector3D(0.5,...): r = " << r << "\n";
        
        r1 /= 0.5;                        // Compound scalar division
        cout << "After r1 /= 0.5: r1 = " << r1 << "\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n";

        // --- Arithmetic -----------------------------------------------------
        cout << "- Arithmetic ——————————————————————————————————————————————————————————\n";
        Vector3D u = q + r;               // Addition
        cout << "Vector addition:           u = q + r                → u = " << u << "\n";
        
        Vector3D v = r - q;               // Subtraction
        cout << "Vector subtraction:        v = r - q                → v = " << v << "\n";
        
        Vector3D w = 2.5 * u;             // Scalar multiplication
        cout << "Scalar multiplication:     w = 2.5 * u              → w = " << w << "\n";
        
        Vector3D n = unit(w);             // Normalization
        cout << "Normalized vector:         n = unit(w)              → n = " << n << "\n";
        
        double d  = u * r;                // Dot product
        cout << "Dot product:               d = u · r                → d = " << d << "\n";
        
        Vector3D c = u % r;               // Cross product
        cout << "Cross product:             c = u × r                → c = " << c << "\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n";
        
        // --- Geometry helpers ----------------------------------------------
        cout << "- Geometry helpers ————————————————————————————————————————————————————\n";
        Vector3D mid = line(0.5, q, r);   // Midpoint
        cout << "Midpoint of q and r:       mid = line(0.5, q, r)    → mid = " << mid << "\n";
        
        // --- Summary Output ------------------------------------------------
        cout << "\n--- Summary ---\n";
        cout << "q = " << q << "\n";
        cout << "r = " << r << "\n";
        cout << "u = " << u << "\n";
        cout << "v = " << v << "\n";
        cout << "w = " << w << "\n";
        cout << "n = " << n << "\n";
        cout << "u·r = " << d << "\n";
        cout << "u×r = " << c << "\n";
        cout << "midpoint = " << mid << "\n";
    
        cout << "QUATERNION Test -------------------------\n";

    }

}

///////////////////     DRAW       ///////////////////////
void Draw() {

	if (ciclo > 0) {
        /*Draw here with OpenGL*/	
	
	}
}


void ProcessingProto() {
	Setup();
	Draw();
}



void drawLine(const Vector3D& a, const Vector3D& b) {

        glColor3ub(0, 0, 0);
	glLineWidth(2.0);
        glBegin(GL_LINES);
        glVertex3f(a.x(), a.y(), a.z());
        glVertex3f(b.x(), b.y(), b.z());
        glEnd();
}

void drawLineColor(const Vector3D& a, const Vector3D& b, int R, int G, int B) {

        glColor3ub(R, G, B);
        glLineWidth(2.0);
        glBegin(GL_LINES);
        glVertex3f(a.x(), a.y(), a.z());
        glVertex3f(b.x(), b.y(), b.z());
        glEnd();
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
	glMaterialfv(GL_FRONT, GL_SPECULAR, specref);

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
	ciclo += 1;

	if (count > 2 * M_PI) count = 0;
	if (ciclo > 100) ciclo = 1; //CICLO NUNCA ES CERO JAJAJA
	
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

		case 'C':
      			count3 += 0.01;
      		break;

    		
		case 'c':
      			count3 -= 0.01;
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

    		
		case 'f':
      			rad += 0.2;
      		break;

    
		case 'F':
      			rad -= 0.2;
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
