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
#include "Cube.hpp"

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

// Glass transparency toggle
bool g_glassMode = false;

// Fisheye camera parameters
static bool g_fisheyeMode = false;
static float g_fisheyeStrength = 1.0f;  // 0.0 = normal, 2.0 = extreme fisheye

// Surface proximity and collision detection
static float g_minDistanceToSurface = 0.1f;
static float g_surfaceZoomFactor = 1.0f;
static bool g_enableSurfaceZoom = true;

// UI state for plane selection
int g_currentOrientation = 0;  // 0=XY, 1=YZ, 2=XZ
int g_currentLayer = 4;        // Current layer within the selected orientation
int g_maxLayers = 8;           // Maximum layers available (matches cube2/cube3 subdivision)

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

// Main.c compatible triangle rendering
void drawFacetMainCStyle(const Facet& f, int index);

void drawSphere(const Vector3D& pos, float radius, int slices, int stacks);
void drawAxes(float length);
void drawText3D(const Vector3D& pos, const char* text);

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

//==============================================================================
// Fisheye Camera Functions
//==============================================================================

//-----------------------------------------------------------------------------
// Apply fisheye distortion projection
//-----------------------------------------------------------------------------
void applyFisheyeProjection(int w, int h) {
    if (!g_fisheyeMode) {
        // Standard perspective projection
        gluPerspective(50.0, (double)w/h, 0.1, 1000.0);
        return;
    }
    
    // Custom fisheye projection matrix
    glLoadIdentity();
    
    float aspect = (float)w / h;
    float baseFov = 50.0f;
    
    // Fisheye parameters - increase FOV dramatically for fisheye effect
    float fisheyeFactor = 1.0f + g_fisheyeStrength * 2.5f;  // 1.0 to 3.5
    float fisheyeFov = baseFov * fisheyeFactor;
    
    // Clamp to reasonable values (but allow wide angles)
    fisheyeFov = fminf(fisheyeFov, 170.0f);
    
    // Apply perspective with fisheye FOV
    gluPerspective(fisheyeFov, aspect, 0.05, 1000.0);
    
    // Apply additional barrel distortion matrix
    if (g_fisheyeStrength > 0.1f) {
        GLfloat distortionMatrix[16] = {
            1.0f - g_fisheyeStrength * 0.2f, 0, 0, 0,
            0, 1.0f - g_fisheyeStrength * 0.2f, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1
        };
        glMultMatrixf(distortionMatrix);
    }
}

/*---------------------------------------
Global setup parameters 
-----------------------------------------*/
int cube_dim = 2.0;
int N = 37 ;
Cube cube(cube_dim, Vector3D{0,0.0001,0}, N);
FacetBox plane_subcells;

double scal = 0.0;
double radi = 0.0;

//==============================================================================
// Simple Cube Reflection Function
//==============================================================================
    
/**
 * @brief Reflect a subdivided cube using checkboard pattern (i+j+k)%2
 * @param cube Target cube to modify and update
 * @param center Center of the reflection sphere 
 * @param radius Radius of the reflection sphere
 * @returns void
 */
void applySigmaTransformationToCube(Cube& cube, const Vector3D& center, double radius) {
    if (!cube.hasSubdivision()) {
        return;
    }
    
    int n = cube.getSubdivisionLevels();
    int centerIdx = n / 2;
    // Apply transformation to all vertices
    for (int i = -centerIdx; i < centerIdx+1; i++) {
        for (int j = -centerIdx; j < centerIdx+1; j++) {
            for (int k = -centerIdx; k < centerIdx+1; k++) {
                for (int l = 0; l < 8; l++) {
                    Vector3D v_p = cube.getSubCell(i, j, k).vertices[l];
                    Vector3D new_pos = sigma(v_p, center, radius);
                    cube.updateSubCellVertex(i, j, k, l, new_pos);
                }
            }
        }
    }
    
    // Refresh triangulation
    cube.refreshTriangulation();
    cout << "  Triangulation refreshed for cube\n";                                                                                                                                     
}


void Setup() {

	if (ciclo == 0) {

        cout << "\n———————————————————————————————————————————————————————————————————————\n";
        cout <<   "|- CUBEs                    ———————————————————————————————————————————\n";
        cout <<   "———————————————————————————————————————————————————————————————————————\n\n";
        cout << "⠀⠀⠀⠀⠀⠀⠀⠀⠀⠄⠀⡀⠀⠀⠀⠀⠀⠀⢀⠀⠀⡀⠀⢀⠀⢀⡀⣤⡢⣤⡤⡀⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀\n";
        cout << "⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⡄⡄⠐⡀⠈⣀⠀⡠⡠⠀⣢⣆⢌⡾⢙⠺⢽⠾⡋⣻⡷⡫⢵⣭⢦⣴⠦⠀⢠⠀⠀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀\n";
        cout << "⠀⠀⠀⠀⠀⠀⠀⠀⢠⣤⣽⣥⡈⠧⣂⢧⢾⠕⠞⠡⠊⠁⣐⠉⠀⠉⢍⠀⠉⠌⡉⠀⠂⠁⠱⠉⠁⢝⠻⠎⣬⢌⡌⣬⣡⣀⣢⣄⡄⠀⡀⠀⠀⠀⠀⠀⠀⠀\n";
        cout << "⠀⠀⠀⠀⠀⠀⠀⢀⢸⣿⣿⢿⣾⣯⣑⢄⡂⠀⠄⠂⠀⠀⢀⠀⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⠄⠐⠀⠀⠀⠀⣄⠭⠂⠈⠜⣩⣿⢝⠃⠀⠁⠀⠀⠀⠀⠀⠀⠀\n";
        cout << "⠀⠀⠀⠀⠀⠀⠀⢀⣻⡟⠏⠀⠚⠈⠚⡉⡝⢶⣱⢤⣅⠈⠀⠄⠀⠀⠀⠀⠀⠠⠀⠀⡂⠐⣤⢕⡪⢼⣈⡹⡇⠏⠏⠋⠅⢃⣪⡏⡇⡍⠀⠀⠀⠀⠀⠀⠀⠀\n";
        cout << "⠀⠀⠀⠀⠀⠀⠀⠀⠺⣻⡄⠀⠀⠀⢠⠌⠃⠐⠉⢡⠱⠧⠝⡯⣮⢶⣴⣤⡆⢐⣣⢅⣮⡟⠦⠍⠉⠀⠁⠐⠀⠀⠀⠄⠐⠡⣽⡸⣎⢁⠀⠀⠀⠀⠀⠀⠀⠀\n";
        cout << "⠀⠀⠀⠀⠀⠀⠀⢈⡻⣧⠀⠁⠐⠀⠀⠀⠀⠀⠀⠊⠀⠕⢀⡉⠈⡫⠽⡿⡟⠿⠟⠁⠀⠀⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠬⠥⣋⡯⠀⠀⠀⠀⠀⠀⠀⠀⠀\n";
        cout << "⠀⠀⠀⠀⠀⠀⠀⡀⣾⡍⠕⡀⠀⠀⠀⠄⠠⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠥⣤⢌⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠀⠀⠄⢀⠀⢝⢞⣫⡆⡄⠀⠀⠀⠀⠀⠀⠀⠀\n";
        cout << "⠀⠀⠀⠀⠀⠀⠀⠀⣽⡶⡄⠐⡀⠀⠀⠀⠀⠀⠀⢀⠀⠄⠀⠀⠀⠄⠁⠇⣷⡆⠀⠀⠀⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡸⢝⣮⠍⠀⠀⠀⠀⠀⠀⠀⠀⠀\n";
        cout << "⠀⠀⠀⠀⠀⠀⢀⠀⢾⣷⠀⠠⡀⠀⠀⠀⠀⢀⠀⠀⠀⠀⠀⠁⡁⠀⠀⣾⡥⠖⠀⠀⠀⠂⠀⠀⠀⠀⠀⠁⠀⡀⠁⠀⠀⠻⢳⣻⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀\n";
        cout << "⠀⠀⠀⠀⠀⠀⠀⠀⣞⡙⠨⣀⠠⠄⠀⠂⠀⠀⠀⠈⢀⠀⠀⠀⠀⠀⠤⢚⢢⣟⠀⠀⠀⠀⡐⠀⠀⡀⠀⠀⠀⠀⠁⠈⠌⠊⣯⣮⡏⠡⠂⠀⠀⠀⠀⠀⠀⠀\n";
        cout << "⠀⠀⠀⠀⠀⠀⠀⠀⣻⡟⡄⡡⣄⠀⠠⠀⠀⡅⠀⠐⠀⡀⠀⡀⠀⠄⠈⠃⠳⠪⠤⠀⠀⠀⠀⡀⠀⠂⠀⠀⠀⠁⠈⢠⣠⠒⠻⣻⡧⠀⠀⠀⠀⠀⠀⠀⠀⠀\n";
        cout << "⠀⠀⠀⠀⠀⠀⠀⠀⠪⡎⠠⢌⠑⡀⠂⠀⠄⠠⠀⠠⠀⠁⡀⠠⠠⡀⣀⠜⢏⡅⠀⠀⡀⠁⠀⠀⠁⠁⠐⠄⡀⢀⠀⠀⠄⢑⣿⣿⣿⡀⠀⠀⠀⠀⠀⠀⠀⠀\n";
        cout << "⠀⠀⠀⠀⠀⠀⠀⠀⠼⣻⠧⣣⣀⠐⠨⠁⠕⢈⢀⢀⡁⠀⠈⠠⢀⠀⠐⠜⣽⡗⡤⠀⠂⠀⠠⠀⢂⠠⠀⠁⠀⠀⠔⠀⠑⣨⣿⢯⠋⡅⠀⠀⠀⠀⠀⠀⠀⠀\n";
        cout << "⠀⠀⠀⠀⠀⠀⠀⠀⡚⣷⣭⠎⢃⡗⠄⡄⢀⠁⠀⠅⢀⢅⡀⠠⠀⢠⡀⡩⠷⢇⠀⡀⠄⡀⠄⠂⠀⠀⠄⠉⡠⠃⠴⠀⠈⢁⣿⡛⡯⠀⠀⠀⠀⠀⠀⠀⠀⠀\n";
        cout << "⠀⠀⠀⠀⠀⠀⠀⠘⡬⡿⣿⡏⡻⡯⠌⢁⢛⠠⠓⠐⠐⠐⠌⠃⠋⠂⡠⢰⣈⢏⣠⠂⠈⠀⠠⠒⠡⠄⠢⠤⠨⠢⡬⠆⠿⢷⢿⡽⡧⠉⠊⠀⠀⠀⠀⠀⠀⠀\n";
        cout << "⠀⠀⠀⠀⠀⠀⠀⠀⠺⣷⣺⣗⣿⡶⡎⡅⣣⢎⠠⡅⣢⡖⠴⠬⡈⠂⡨⢡⠾⣣⣢⠀⠀⡹⠄⡄⠄⡇⣰⡖⡊⠔⢹⣄⣿⣭⣵⣿⢷⡀⠀⠀⠀⠀⠀⠀⠀⠀\n";
        cout << "⠀⠀⠀⠀⠀⠀⠀⠀⠩⣿⣿⣲⣿⣷⣟⣼⠟⣬⢉⡠⣪⢜⣂⣁⠥⠓⠚⡁⢶⣷⣠⠂⡄⡢⣀⡐⠧⢆⣒⡲⡳⡫⢟⡃⢪⡧⣟⡟⣯⠐⠀⠀⠀⠀⠀⠀⠀⠀\n";
        cout << "⠀⠀⠀⠀⠀⠀⠀⠀⢺⠟⢿⢟⢻⡗⡮⡿⣲⢷⣆⣏⣇⡧⣄⢖⠾⡷⣿⣤⢳⢷⣣⣦⡜⠗⣭⢂⠩⣹⢿⡲⢎⡧⣕⣖⣓⣽⡿⡖⡿⠀⠀⠀⠀⠀⠀⠀⠀⠀\n";
        cout << "⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠉⠂⠂⠏⠿⢻⣥⡪⢽⣳⣳⣥⡶⣫⣍⢐⣥⣻⣾⡻⣅⢭⡴⢭⣿⠕⣧⡭⣞⣻⣣⣻⢿⠟⠛⠙⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀\n";
        cout << "⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠄⠋⠫⠯⣍⢻⣿⣿⣷⣕⣵⣹⣽⣿⣷⣇⡏⣿⡿⣍⡝⠵⠯⠁⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀\n";
        cout << "⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⠠⠁⠋⢣⠓⡍⣫⠹⣿⣿⣷⡿⠯⠺⠁⠁⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀\n";
        cout << "⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠋⢀⠋⢈⡿⠿⠁⠉⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀\n";
        cout << "                              ~[Cube Class Example - Enhanced Camera]\n\n";
                                                                                                                                                                                        
        // ==================== CUBE SETUP ====================


        // ==================== ADVANCED SUBDIVISION DEMOS - Reflection of the subdivision in Spheres ====================
        applySigmaTransformationToCube(cube, Vector3D{0,0,0}, 0.5); 
        applySigmaTransformationToCube(cube, cube.getSubcellCenter(1,1,1), cube.getSubcellRadius(1,1,1));
        cube.writeSTL_s("/home/mike666/Downloads/mesh_output.stl", "MyCube", "checkerboard", 2, 2, 9);
    }
}

///////////////////     DRAW       ///////////////////////
int ii = 4;
int jj = 5;
int kk = 9;
void Draw() {
    
    extern void Setup(); // assume you define this elsewhere
    static int ciclo = 1;  // or however you manage visibility
	if (ciclo > 0) {
        /*Draw here with OpenGL*/
        // In your draw‐all loop:
        plane_subcells = cube.getCheckerboardFacets(ii, kk, jj);
        for (size_t i = 0; i < plane_subcells.size(); ++i) {
            drawFacetMainCStyle(plane_subcells[i], (int)(i + 500));
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

    if (g_glassMode) {
        // Glass transparency effect
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glDepthMask(GL_FALSE);  // Don't write to depth buffer for transparency

        // Glass material properties
        glColor4f(0.7f, 0.9f, 1.0f, 0.3f);  // Light blue tint with transparency

        // Enhanced material for glass effect
        GLfloat glass_ambient[] = {0.2f, 0.3f, 0.4f, 0.3f};
        GLfloat glass_diffuse[] = {0.3f, 0.5f, 0.7f, 0.3f};
        GLfloat glass_specular[] = {1.0f, 1.0f, 1.0f, 0.3f};
        GLfloat glass_shininess[] = {100.0f};

        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, glass_ambient);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, glass_diffuse);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, glass_specular);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, glass_shininess);

        glutSolidSphere(radius, slices, stacks);

        glDepthMask(GL_TRUE);  // Re-enable depth writing
    } else {
        // Regular solid sphere
        glColor3f(0.8f, 0.2f, 0.2f);  // Red color for solid mode
        glutSolidSphere(radius, slices, stacks);
    }

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
// Main.c compatible triangle drawing function - uses exact same color scheme
void drawFacetMainCStyle(const Facet& f, int index)
{
    // Get triangle geometry
    Vector3D normal = f.getNormal();
    Vector3D A = f[0], B = f[1], C = f[2];

    // Preserve polygon & color state
    glPushAttrib(GL_COLOR_BUFFER_BIT | GL_POLYGON_BIT);

    // Draw filled triangle with polygon offset to prevent z-fighting
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0f, 1.0f);

    // Use Main.c color scheme - cycle through 4 colors based on triangle index
    //if (index % 4 == 0) glColor3ub(200, 200, 200);  // Light gray
    //if (index % 4 == 1) glColor3ub(200, 0, 0);      // Red
    //if (index % 4 == 2) glColor3ub(0, 200, 0);      // Green
    //if (index % 4 == 3) glColor3ub(0, 0, 200);      // Blue
    glColor3ub(200, 200, 200);

    // Render triangle with proper normal - same style as Main.c geometry.c
    glBegin(GL_TRIANGLES);
        glNormal3f(normal.x(), normal.y(), normal.z());
        glVertex3f(A.x(), A.y(), A.z());
        glVertex3f(B.x(), B.y(), B.z());
        glVertex3f(C.x(), C.y(), C.z());
    glEnd();
    glDisable(GL_POLYGON_OFFSET_FILL);

    // Draw triangle border lines
    glColor3ub(0, 0, 0);  // Black border
    glLineWidth(1.0f);
    glBegin(GL_LINE_LOOP);
        glVertex3f(A.x(), A.y(), A.z());
        glVertex3f(B.x(), B.y(), B.z());
        glVertex3f(C.x(), C.y(), C.z());
    glEnd();

    // Restore state
    glPopAttrib();
}

// Keep the old function for backwards compatibility
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
void keyboard(unsigned char key, int x, int y);

//-----------------------------------------------------------------------------
// Main
//-----------------------------------------------------------------------------
int main(int argc, char** argv)
{
    srand((unsigned)time(nullptr));
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(720, 720);
    glutCreateWindow(" JAZ 4D Enhanced Camera U.U ");

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
    glutKeyboardFunc(keyboard);

    // Enter main loop
    glutMainLoop();
    return 0;
}

//-----------------------------------------------------------------------------
// Setup OpenGL lights, materials, and default projection
//-----------------------------------------------------------------------------
void initGL()
{
    // Lighting - matching Main.c exactly
    GLfloat ambient[]  = {0.3f, 0.3f, 0.3f, 1.0f};
    GLfloat diffuse[]  = {0.7f, 0.7f, 0.7f, 1.0f};
    GLfloat specular[] = {1.0f, 1.0f, 1.0f, 1.0f};
    GLfloat position[] = {1.0f, 1.0f, 0.25f, 0.0f};  // Match Main.c light position

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glLightfv(GL_LIGHT0, GL_AMBIENT,  ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE,  diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
    glLightfv(GL_LIGHT0, GL_POSITION, position);

    // Material properties - matching Main.c exactly
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);  // Same as Main.c
    GLfloat specref[] = {1.0f, 1.0f, 1.0f, 1.0f};
    GLfloat shine[] = {128.0f};  // Same shininess as Main.c
    glMaterialfv(GL_FRONT, GL_SPECULAR, specref);
    glMaterialfv(GL_FRONT, GL_SHININESS, shine);

    // Depth test
    glEnable(GL_DEPTH_TEST);
    glFrontFace(GL_CCW);

    // Normalize normals for scaled geometry
    glEnable(GL_NORMALIZE);

    // Clear color - match Main.c white background
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
}

//-----------------------------------------------------------------------------
// Handle window size changes with fisheye support
//-----------------------------------------------------------------------------
void reshape(int w, int h)
{
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    // Apply fisheye projection
    applyFisheyeProjection(w, h);
    
    glMatrixMode(GL_MODELVIEW);
}

//-----------------------------------------------------------------------------
// Draw help overlay with new controls
//-----------------------------------------------------------------------------
void drawHUD() {
    if (!g_showHelp) return;
    const char* lines[] = {
        "L-drag: Rotate",
        "M-drag: Pan", 
        "R-drag/Wheel: Zoom",
        "[H]: Toggle Help",
        "[G]: Toggle Glass Mode",
        "[F]: Toggle Fisheye Mode",
        "[[]: Fisheye Strength -/+",
        "[O]: Cycle Orientation (XY/YZ/XZ)",
        "[+/-]: Change Layer",
        "[R]: Toggle Reflection",
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

    // Display current plane selection status
    y -= 0.05f;
    char status[100];
    const char* orientName = g_currentOrientation == 0 ? "XY" :
                            g_currentOrientation == 1 ? "YZ" : "XZ";
    sprintf(status, "Current: %s plane, Layer %d/%d", orientName, g_currentLayer, g_maxLayers-1);
    glColor3f(0.8f, 0.2f, 0.2f);  // Red color for status
    glRasterPos2f(0.02f, y);
    for(const char* c=status; *c; ++c)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);

    // Add fisheye status
    if (g_fisheyeMode) {
        y -= 0.05f;
        char fisheyeStatus[100];
        sprintf(fisheyeStatus, "Fisheye: ON (Strength: %.1f)", g_fisheyeStrength);
        glColor3f(0.2f, 0.2f, 0.8f);  // Blue color for fisheye
        glRasterPos2f(0.02f, y);
        for(const char* c=fisheyeStatus; *c; ++c)
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);
    }
    
    glMatrixMode(GL_MODELVIEW); glPopMatrix();
    glMatrixMode(GL_PROJECTION); glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
}

//-----------------------------------------------------------------------------
// Main display with surface zoom and fisheye
//-----------------------------------------------------------------------------
void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    // Apply surface zoom to the base zoom
    float effectiveZoom = g_zoom / g_surfaceZoomFactor;
    
    // Move back and zoom (with surface zoom effect)
    glTranslatef(g_panX, g_panY, -5.0f * effectiveZoom);

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
// Keyboard handler with fisheye and surface zoom controls
//-----------------------------------------------------------------------------
void keyboard(unsigned char key, int x, int y) {
    switch (key) {
        case 'g':
        case 'G':
            g_glassMode = !g_glassMode;
            printf("Glass mode: %s\n", g_glassMode ? "ON" : "OFF");
            glutPostRedisplay();
            break;
            
        case 'h':
        case 'H':
            g_showHelp = !g_showHelp;
            glutPostRedisplay();
            break;

        // Fisheye controls
        case 'f':
        case 'F':
            g_fisheyeMode = !g_fisheyeMode;
            printf("Fisheye mode: %s\n", g_fisheyeMode ? "ON" : "OFF");
            reshape(glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
            glutPostRedisplay();
            break;
            
        case '[':
            g_fisheyeStrength = fmaxf(0.0f, g_fisheyeStrength - 0.1f);
            printf("Fisheye strength: %.1f\n", g_fisheyeStrength);
            if (g_fisheyeMode) {
                reshape(glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
            }
            glutPostRedisplay();
            break;
            
        case ']':
            g_fisheyeStrength = fminf(2.0f, g_fisheyeStrength + 0.1f);
            printf("Fisheye strength: %.1f\n", g_fisheyeStrength);
            if (g_fisheyeMode) {
                reshape(glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
            }
            glutPostRedisplay();
            break;
            
        // Plane orientation controls
        case 'o':
            kk += 1;
            glutPostRedisplay();
            break;
        case 'O':
            kk -= 1;
            glutPostRedisplay();
            break;

        // Layer controls
        case '+':
            ii += 1;
            glutPostRedisplay();
            break;
        case '-':
            ii -= 1;
            glutPostRedisplay();
            break;
        case 27: // ESC key
            exit(0);
            break;
    }
}

//-----------------------------------------------------------------------------
// UI menu callback
//----------------------------------------------------------------------------- 
void MenuHandler(int choice) {
    switch(choice) {
        case 1: g_showHelp = !g_showHelp; break;       // toggle HUD
        case 2: g_angleX=20; g_angleY=-30; g_zoom=1; g_panX=g_panY=0; 
                g_fisheyeMode=false; g_fisheyeStrength=1.0f; 
                g_enableSurfaceZoom=true; break; // reset
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
