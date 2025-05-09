#include <vector>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <GL/glut.h>
#include "Vector3D.cpp"
#include "Vector4D.cpp"
#include "Quaternion.cpp"
#include "Facet.cpp"
#include "FacetBox.hpp"
#include "Dodecahedron.hpp"

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

// OpenGL functions
void drawFacet( const Facet& f, 
                int R, int G, int Bi, 
                float alpha);

void drawSphere(const Vector3D& pos, float radius, int slices, int stacks);
void drawAxes(float length);
void drawText3D(const Vector3D& pos, const char* text);

Vector3D d_ui(1.0, 0.0, 0.0);
Vector3D e_ui(0.0, 1.0, 0.0);
Vector3D f_ui(0.0, 0.0, 1.0);
Facet f_1(d_ui, e_ui, f_ui);

Dodecahedron dodeca(2.0, Vector3D{0,0,0});;
FacetBox box_out_0;
FacetBox box_out_1;
FacetBox box_out_2;
FacetBox box_out_3;


// Helper: convert HSV→RGB (all in [0,1])                                                  
struct Color { float r, g, b; };
Color hsv2rgb(float h, float s, float v) {
    h = fmodf(h, 360.0f) / 60.0f;
    int i = int(floor(h));
    float f = h - i;
    float p = v * (1 - s);
    float q = v * (1 - s * f);
    float t = v * (1 - s * (1 - f));
    switch(i) {
      case 0: return {v, t, p};
      case 1: return {q, v, p};
      case 2: return {p, v, t};
      case 3: return {p, q, v};
      case 4: return {t, p, v};
      default:return {v, p, q};
    }
}


// Subdivide each triangle in `initial` n times, returning only the final mesh.
FacetBox refine(const FacetBox& initial, int n) {
    FacetBox curr = initial;
    FacetBox next;

    for (int pass = 0; pass < n; ++pass) {
        next.clear();  // throw away last level
        for (size_t i = 0; i < curr.size(); ++i) {
            // fromFacet() returns 3 new sub-triangles around the centroid
            FacetBox tiny = FacetBox::fromFacet(curr[i]);
            next += tiny;  // append those three
        }
        std::swap(curr, next);
    }

    return curr;  // only the final, fully-refined mesh
}


void Setup() {

	if (ciclo == 0) {

        cout << "\n———————————————————————————————————————————————————————————————————————\n";
        cout <<   "|- Hyperbolic Dodecahedron ———————————————————————————————————————————\n";
        cout <<   "———————————————————————————————————————————————————————————————————————\n\n";

        cout << "      _----------_,\n";
        cout << "    ,\"__         _-:, \n";
        cout << "   /    \"\"--_--\"\"...:\\ \n";
        cout << "  /         |.........\\ \n";
        cout << " /          |..........\\ \n";
        cout << "/,         _'_........./. \n";
        cout << "! -,    _-\"   \"-_... ,;;: \n";
        cout << "\\   -_-\"         \"-_/;;;. \n";
        cout << " \\   \\             /;;;. \n";
        cout << "  \\   \\           /;;;. \n";
        cout << "   '.  \\         /;;;' \n";
        cout << "     \"-_\\_______/;;'        ~[Hyperbolic Dodecahedron Class]\n\n";

        // 1) Build your initial mesh (e.g. from a dodecahedron) 
        Dodecahedron d1(0.93, Vector3D(0, 0, 0), FaceMode::Pentagons);
        FacetBox box_sum_1(d1.getFacets());
        int levels = 5;  // or however many you like
        box_out_1 = box_sum_1.refine(levels, FacetBox::SubdivisionMode::Midpoint4);

        // 2) Mutate in-place: project every triangle into hyperbolic space
        box_out_1.applyHyperboloid();
        // Now 'box_out_1' contains hyperolic-mapped facets

        // 3) Or, produce a new hyperbolic mesh without touching the original:
        //FacetBox hyperBox = box_out_1.hyperboloid();

    }

}

///////////////////     DRAW       ///////////////////////
void Draw() {
    
    extern void Setup(); // assume you define this elsewhere
    static int ciclo = 1;  // or however you manage visibility
	if (ciclo > 0) {
        /*Draw here with OpenGL*/	
        //drawFacet(f_1, 200, 10, 40, 0.75f);
        //drawSphere(f_1[0], 0.1f, 6, 6);
        //drawSphere(f_1[1], 0.1f, 6, 6);
        //drawSphere(f_1[2], 0.1f, 6, 6);


        size_t total = box_out_1.size();      // e.g. 36*3
        for(size_t i = 0; i < total; ++i) {
            // pick hue from 0°→360° across the range
            float hue = float(i) / float(total) * 360.0f;
            Color c = hsv2rgb(hue, 0.8f, 1.0f);   // 80% saturation, full value
            // convert to 0–255 ints
            int R = int(c.r * 255), G = int(c.g * 255), B = int(c.b * 255);
            drawFacet(box_out_1[i], R, G, B, 1.0f);
        }


	}
}


void ProcessingProto() {
	//extern void Setup();  // your existing setup
    //Setup();
	Draw();
}


//-----------------------------------------------------------------------------
// Utility: draw a small sphere (GLUT) at a 3D position
//----------------------------------------------------------------------------- 
void drawSphere(const Vector3D& pos, float radius = 0.05f, int slices = 12, int stacks = 12) {
    glPushMatrix();
    glTranslatef(pos.x(), pos.y(), pos.z());
    // You can use GLUT's sphere or a display list for better performance
    glutSolidSphere(radius, slices, stacks);
    glPopMatrix();
}

//-----------------------------------------------------------------------------
// Draw coordinate axes (X=red, Y=green, Z=blue)
//-----------------------------------------------------------------------------
void drawAxes(float length) {
    glLineWidth(2.0f);
    glBegin(GL_LINES);
      // X axis
      glColor3f(1,0,0);
      glVertex3f(-length,0,0);
      glVertex3f(length,0,0);
      // Y axis
      glColor3f(0,1,0);
      glVertex3f(0,-length,0);
      glVertex3f(0,length,0);
      // Z axis
      glColor3f(0,0,1);
      glVertex3f(0,0,-length);
      glVertex3f(0,0,length);
    glEnd();
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

// Draws a filled, semi-transparent triangle plus its outline.
//
//  f     – your Facet
//  r,g,b – color components in [0..255]
//  alpha – [0..1] opacity (default 0.75)
inline void drawFacet(const Facet& f,
                      int r, int g, int b,
                      float alpha)
{
    // 1) Fetch geometry
    Vector3D normal = f.getNormal();   // (x,y,z)
    Vector3D A = f[0], B = f[1], C = f[2];

    // 2) Convert color to floats
    const float Rf = r / 255.0f;
    const float Gf = g / 255.0f;
    const float Bf = b / 255.0f;

    // 3) Preserve polygon & color state
    glPushAttrib(GL_COLOR_BUFFER_BIT | GL_POLYGON_BIT);

    // 4) Draw filled triangle
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0f, 1.0f);
    glColor4f(Rf, Gf, Bf, alpha);
    glBegin(GL_TRIANGLES);
        glNormal3f(normal.x(), normal.y(), normal.z());
        glVertex3f(A.x(), A.y(), A.z());
        glVertex3f(B.x(), B.y(), B.z());
        glVertex3f(C.x(), C.y(), C.z());
    glEnd();
    glDisable(GL_POLYGON_OFFSET_FILL);

    // 5) Draw outline
    glColor4f(0.0f, 0.0f, 0.0f, alpha);   // black lines
    glLineWidth(0.25f);
    glBegin(GL_LINE_LOOP);
        glVertex3f(A.x(), A.y(), A.z());
        glVertex3f(B.x(), B.y(), B.z());
        glVertex3f(C.x(), C.y(), C.z());
    glEnd();

    // 6) Restore state
    glPopAttrib();
}

//---------------------------------------------------------------------------
// Utility: draw text label at a 3D world coordinate
//---------------------------------------------------------------------------
void drawText3D(const Vector3D& pos, const char* text) {
    // Project 3D point to window coordinates
    GLint viewport[4];
    GLdouble modelview[16], projection[16];
    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
    glGetDoublev(GL_PROJECTION_MATRIX, projection);
    GLdouble winX, winY, winZ;
    gluProject(pos.x(), pos.y(), pos.z(), modelview, projection, viewport, &winX, &winY, &winZ);
    // Flip Y for raster position
    winY = viewport[3] - winY;

    // Prepare for 2D overlay
    glMatrixMode(GL_PROJECTION); glPushMatrix();
    glLoadIdentity();
    glOrtho(0, viewport[2], 0, viewport[3], -1, 1);
    glMatrixMode(GL_MODELVIEW); glPushMatrix();
    glLoadIdentity();

    // Draw text
    glRasterPos2i((int)winX, (int)winY);
    for (const char* c = text; *c; ++c) {
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *c);
    }

    // Restore matrices
    glPopMatrix();
    glMatrixMode(GL_PROJECTION); glPopMatrix();
    glMatrixMode(GL_MODELVIEW);

    // Restore state
    glPopAttrib();
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

//-----------------------------------------------------------------------------
// Global state for interactive camera
//----------------------------------------------------------------------------- 
static float g_angleX = 20.0f, g_angleY = -30.0f; // view angles (degrees)
static float g_zoom   = 1.0f;                    // zoom factor
static float g_panX   = 0.0f, g_panY = 0.0f;      // camera pan offsets
static int   g_lastX  = 0, g_lastY = 0;          // last mouse coords
static bool  g_leftDown   = false;               // rotating
static bool  g_middleDown = false;               // panning
static bool  g_rightDown  = false;               // zooming
static bool  g_showHelp  = true;                    // HUD toggle

//-----------------------------------------------------------------------------
// Forward declarations
//----------------------------------------------------------------------------- 
void initGL();
void reshape(int w, int h);
void display();
void mouseButton(int button, int state, int x, int y);
void mouseMotion(int x, int y);
void ProcessMenu(int value);
void reshape(int w, int h);
void MenuHandler(int choice);
void createUI();
void drawHUD();

//-----------------------------------------------------------------------------
// Main
//-----------------------------------------------------------------------------
int main(int argc, char** argv)
{
    srand((unsigned)time(nullptr));
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(720, 720);
    glutCreateWindow(" JAZ 4D   U.U ");

    // Enable smoothing & blending by default
    ProcessMenu(1);
    initGL();
    createUI();
    
    // Objects setup
    Setup();

    // Register callbacks
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouseButton);
    glutMotionFunc(mouseMotion);

    // Enter main loop
    glutMainLoop();
    return 0;
}

//-----------------------------------------------------------------------------
// Setup OpenGL lights, materials, and default projection
//-----------------------------------------------------------------------------
void initGL()
{
    // Lighting
    GLfloat ambient[]  = {0.3f, 0.3f, 0.3f, 1.0f};
    GLfloat diffuse[]  = {0.7f, 0.7f, 0.7f, 1.0f};
    GLfloat specular[] = {1, 1, 1, 1};
    GLfloat position[] = {20, 20, 20.25f, 0};

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glLightfv(GL_LIGHT0, GL_AMBIENT,  ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE,  diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
    glLightfv(GL_LIGHT0, GL_POSITION, position);

    // Material
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
    GLfloat shine[] = {128.0f};
    glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, shine);

    // Depth test
    glEnable(GL_DEPTH_TEST);
    glFrontFace(GL_CCW);

    // Normalize normals for scaled geometry
    glEnable(GL_NORMALIZE);

    // Clear color
    glClearColor(1,1,1,1);
}

//-----------------------------------------------------------------------------
// Handle window size changes
//-----------------------------------------------------------------------------
void reshape(int w, int h)
{
    glViewport(0,0,w,h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(50.0, (double)w/h, 0.5, 1000.0);
    glMatrixMode(GL_MODELVIEW);
}

//-----------------------------------------------------------------------------
// Draw help overlay
//-----------------------------------------------------------------------------
void drawHUD() {
    if (!g_showHelp) return;
    const char* lines[] = {
        "L-drag: Rotate",
        "M-drag: Pan",
        "R-drag/Wheel: Zoom",
        "[H]: Toggle Help",
        "Right-click: UI Menu"
    };
    glMatrixMode(GL_PROJECTION); glPushMatrix();
    glLoadIdentity(); glOrtho(0,1,0,1,-1,1);
    glMatrixMode(GL_MODELVIEW); glPushMatrix(); glLoadIdentity();

    glColor3f(0,0,0);
    float y = 0.95f;
    for(auto &ln: lines) {
        glRasterPos2f(0.02f,y);
        for(const char* c=ln;*c;++c)
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,*c);
        y -= 0.05f;
    }

    glMatrixMode(GL_MODELVIEW); glPopMatrix();
    glMatrixMode(GL_PROJECTION); glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
}

//-----------------------------------------------------------------------------
// Main display: apply interactive camera, then draw
//-----------------------------------------------------------------------------
void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    /* - Camera: pan, zoom, rotate*/
    // Move back and zoom
    glTranslatef(0, 0, -5.0f * g_zoom);

    // Apply rotations
    glRotatef(g_angleX, 1, 0, 0);
    glRotatef(g_angleY, 0, 1, 0);

    // Draw axes at origin
    drawAxes(10.0f);


    ProcessingProto();   // calls Setup() then Draw()
    drawHUD();
    
    glutSwapBuffers();
}



//-----------------------------------------------------------------------------
// Mouse button: track left/middle/right for rotate/pan/zoom, handle wheel
//----------------------------------------------------------------------------- 
void mouseButton(int button, int state, int x, int y)
{
    if (button == GLUT_LEFT_BUTTON) {
        g_leftDown   = (state == GLUT_DOWN);
    }
    else if (button == GLUT_MIDDLE_BUTTON) {
        g_middleDown = (state == GLUT_DOWN);
    }
    else if (button == GLUT_RIGHT_BUTTON) {
        g_rightDown  = (state == GLUT_DOWN);
    }
    else if (button == 3) {           // wheel up
        g_zoom *= 1.05f;
        glutPostRedisplay();
    }
    else if (button == 4) {           // wheel down
        g_zoom /= 1.05f;
        glutPostRedisplay();
    }
    g_lastX = x; g_lastY = y;
}


//-----------------------------------------------------------------------------
// Mouse drag: update angles/pan/zoom
//-----------------------------------------------------------------------------
void mouseMotion(int x, int y)
{
    int dx = x - g_lastX;
    int dy = y - g_lastY;

    if (g_leftDown) {
        g_angleY += dx * 0.5f;
        g_angleX += dy * 0.5f;
        g_angleX = fmaxf(-90.0f, fminf(90.0f, g_angleX));
    }
    else if (g_middleDown) {
        // pan: move camera laterally
        g_panX += dx * 0.01f * g_zoom;
        g_panY -= dy * 0.01f * g_zoom;
    }
    else if (g_rightDown) {
        g_zoom *= 1.0f - dy * 0.005f;
        g_zoom = fmaxf(0.1f, fminf(10.0f, g_zoom));
    }

    g_lastX = x; g_lastY = y;
    glutPostRedisplay();
}

//-----------------------------------------------------------------------------
// UI menu callback
//----------------------------------------------------------------------------- 
void MenuHandler(int choice) {
    switch(choice) {
        case 1: g_showHelp = !g_showHelp; break;       // toggle HUD
        case 2: g_angleX=20; g_angleY=-30; g_zoom=1; g_panX=g_panY=0; break; // reset
        case 3: exit(0); break;                      // quit
    }
    glutPostRedisplay();
}

//-----------------------------------------------------------------------------
// Create right-click menu
//-----------------------------------------------------------------------------
void createUI() {
    int menu = glutCreateMenu(MenuHandler);
    glutAddMenuEntry("Toggle Help", 1);
    glutAddMenuEntry("Reset View",   2);
    glutAddMenuEntry("Quit",         3);
    glutAttachMenu(GLUT_RIGHT_BUTTON);
}

//-----------------------------------------------------------------------------
// Your existing menu callback for smoothing/blending
//-----------------------------------------------------------------------------
void ProcessMenu(int value)
{
    if (value == 1) {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_POINT_SMOOTH);   glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
        glEnable(GL_LINE_SMOOTH);    glHint(GL_LINE_SMOOTH_HINT,  GL_NICEST);
        glEnable(GL_POLYGON_SMOOTH); glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
    } else {
        glDisable(GL_BLEND);
        glDisable(GL_POINT_SMOOTH);
        glDisable(GL_LINE_SMOOTH);
        glDisable(GL_POLYGON_SMOOTH);
    }
    glutPostRedisplay();
}

