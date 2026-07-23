// main19_truncated_cuboctahedron.cpp
//
// Minimal FreeGLUT demo that draws exactly ONE truncated cuboctahedron built
// from a basic (non-subdivided) Cube via Cube::getFacetsTruncatedCuboctahedron().
//
// The truncated cuboctahedron (great rhombicuboctahedron) is the Archimedean solid
// with 6 octagonal faces (one per cube FACE) + 8 hexagonal faces (one per cube
// CORNER) + 12 square faces (one per cube EDGE) = 26 faces, 48 vertices, 72 edges,
// vertex configuration 4.6.8, triangulated as 92 triangles. It is built on the fly
// from the cube's 8 corner vertices (trilinear interpolation at 48 fixed normalized
// positions), so it needs no extra storage and deforms with the lattice under
// inversion exactly like the plain cube, truncated cube and truncated octahedron.
//
// It is regular only at the fixed a/b encoded in Cube.hpp (a = 1/(1+2√2),
// b = (1+√2)/(1+2√2)); there is no shape morph. The interactive knobs are hollow
// on/off and the inset ratio:
//   hollow = false -> solid 92-tri mesh (6 octagons + 8 hexagons + 12 squares)
//   hollow = true  -> 288-tri frame: each face inset toward its own centroid by
//                     `inset`, inner face skipped as a hole.
//
// Controls:
//   L-drag        rotate    |  R-drag / wheel  zoom
//   ,  /  .       decrease / increase the hollow inset (border thickness, step 0.05)
//   o             toggle hollow (each face: inset toward center, inner face skipped)
//   a             toggle auto-sweep of the inset (thin -> thick -> thin)
//   c             write one STL of the current shape to ~/Downloads
//   h             toggle HUD      |  Esc  quit
//
// Build:
//   g++ -std=c++17 main19_truncated_cuboctahedron.cpp -o main19_truncated_cuboctahedron \
//       -lGL -lGLU -lglut -lm
//   ./main19_truncated_cuboctahedron          # xvfb-run ./... if headless

#include <vector>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string>
#include <GL/glut.h>
#include "Vector3D.cpp"
#include "Quaternion.cpp"
#include "Facet.cpp"
#include "FacetBox.hpp"
#include "Cube.hpp"

using namespace std;

// ----------------------------------------------------------------------------
// Scene: one basic cube -> one truncated cuboctahedron (92 solid / 288 hollow).
// ----------------------------------------------------------------------------
Cube   cube(1.0, Vector3D{0, 0, 0});   // basic, non-subdivided; verts_ populated
bool   g_hollow     = true;            // hollow out each face (inset frame; inner face skipped)
double g_inset      = 0.5;             // hollow border inset ratio (0,1); bigger = bigger hole
bool   g_autoSweep  = false;           // auto-animate the inset
double g_sweepPhase = 0.0;             // auto-sweep phase (radians)

// ----------------------------------------------------------------------------
// Interactive camera state (mirrors main9.cpp / main18).
// ----------------------------------------------------------------------------
static float g_angleX = 20.0f, g_angleY = -30.0f;  // view angles (degrees)
static float g_zoom   = 1.0f;                       // zoom factor
static float g_panX   = 0.0f, g_panY = 0.0f;        // camera pan offsets
static int   g_lastX  = 0, g_lastY = 0;             // last mouse coords
static bool  g_leftDown   = false;                  // rotating
static bool  g_middleDown = false;                  // panning
static bool  g_rightDown  = false;                  // zooming
static bool  g_showHelp   = true;                   // HUD toggle

// Forward declarations.
void initGL();
void reshape(int w, int h);
void display();
void Draw();
void drawHUD();
void drawAxes(float length);
void drawFacetMainCStyle(const Facet& f, int index);
void mouseButton(int button, int state, int x, int y);
void mouseMotion(int x, int y);
void keyboard(unsigned char key, int x, int y);
void sweepTimer(int value);
void Setup();

// ----------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------
int main(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(720, 720);
    glutCreateWindow(" Truncated Cuboctahedron  (hollow inset)  U.U ");

    initGL();
    Setup();

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouseButton);
    glutMotionFunc(mouseMotion);
    glutKeyboardFunc(keyboard);

    glutMainLoop();
    return 0;
}

// ----------------------------------------------------------------------------
// One-time info.
// ----------------------------------------------------------------------------
void Setup()
{
    cout << "\n--- Truncated Cuboctahedron demo ---\n";
    cout << "6 octagons (cube faces) + 8 hexagons (cube corners) + 12 squares (cube edges) "
            "= 26 faces; 92 solid / 288 hollow tris.\n";
    cout << "hollow = " << (g_hollow ? "ON" : "OFF") << "  inset = " << g_inset
         << "  (each face inset toward its center; inner face skipped)\n";
    cout << "Controls:  , . inset   o hollow   a auto-sweep   c STL   h HUD   Esc quit\n\n";
}

// ----------------------------------------------------------------------------
// OpenGL lights, material, depth, smoothing.
// ----------------------------------------------------------------------------
void initGL()
{
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

    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
    GLfloat shine[] = {128.0f};
    glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, shine);

    glEnable(GL_DEPTH_TEST);
    glFrontFace(GL_CCW);
    glEnable(GL_NORMALIZE);              // normals stay unit under scaling

    // Smoothing / blending (for clean outlines).
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

    glClearColor(1, 1, 1, 1);
}

// ----------------------------------------------------------------------------
void reshape(int w, int h)
{
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(50.0, (double)w / (double)h, 0.5, 1000.0);
    glMatrixMode(GL_MODELVIEW);
}

// ----------------------------------------------------------------------------
// Per-frame draw: rebuild the truncated cuboctahedron at the current
// hollow/inset and render it. (92 solid / 288 hollow facets/frame is trivial.)
// ----------------------------------------------------------------------------
void Draw()
{
    FacetBox fb = cube.getFacetsTruncatedCuboctahedron(g_hollow, g_inset);
    size_t total = fb.size();                 // 92 solid / 288 hollow
    for (size_t i = 0; i < total; ++i)
        drawFacetMainCStyle(fb[i], (int)i);
}

// ----------------------------------------------------------------------------
void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    glTranslatef(g_panX, g_panY, -5.0f * g_zoom);
    glRotatef(g_angleX, 1, 0, 0);
    glRotatef(g_angleY, 0, 1, 0);

    drawAxes(2.0f);
    Draw();
    drawHUD();

    glutSwapBuffers();
}

// ----------------------------------------------------------------------------
// Gray lit fill + black outline (copied from main16.cpp:drawFacetMainCStyle).
// ----------------------------------------------------------------------------
void drawFacetMainCStyle(const Facet& f, int /*index*/)
{
    Vector3D normal = f.getNormal();
    Vector3D A = f[0], B = f[1], C = f[2];

    glPushAttrib(GL_COLOR_BUFFER_BIT | GL_POLYGON_BIT);

    // Filled triangle (gray), offset so the outline does not z-fight.
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0f, 1.0f);
    glColor3ub(200, 200, 200);
    glBegin(GL_TRIANGLES);
        glNormal3f(normal.x(), normal.y(), normal.z());
        glVertex3f(A.x(), A.y(), A.z());
        glVertex3f(B.x(), B.y(), B.z());
        glVertex3f(C.x(), C.y(), C.z());
    glEnd();
    glDisable(GL_POLYGON_OFFSET_FILL);

    // Black outline.
    glColor3ub(0, 0, 0);
    glLineWidth(1.0f);
    glBegin(GL_LINE_LOOP);
        glVertex3f(A.x(), A.y(), A.z());
        glVertex3f(B.x(), B.y(), B.z());
        glVertex3f(C.x(), C.y(), C.z());
    glEnd();

    glPopAttrib();
}

// ----------------------------------------------------------------------------
void drawAxes(float length)
{
    glLineWidth(2.0f);
    glBegin(GL_LINES);
        glColor3f(1, 0, 0); glVertex3f(-length, 0, 0); glVertex3f(length, 0, 0);
        glColor3f(0, 1, 0); glVertex3f(0, -length, 0); glVertex3f(0, length, 0);
        glColor3f(0, 0, 1); glVertex3f(0, 0, -length); glVertex3f(0, 0, length);
    glEnd();
}

// ----------------------------------------------------------------------------
// HUD overlay: hollow/inset + controls.
// ----------------------------------------------------------------------------
void drawHUD()
{
    if (!g_showHelp) return;

    glMatrixMode(GL_PROJECTION); glPushMatrix();
    glLoadIdentity(); glOrtho(0, 1, 0, 1, -1, 1);
    glMatrixMode(GL_MODELVIEW); glPushMatrix(); glLoadIdentity();

    glDisable(GL_LIGHTING);
    glColor3f(0, 0, 0);

    char buf[160];
    snprintf(buf, sizeof(buf), "truncated cuboctahedron  4.6.8  (6 oct + 8 hex + 12 sq)");
    glRasterPos2f(0.02f, 0.95f);
    for (const char* c = buf; *c; ++c) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *c);

    snprintf(buf, sizeof(buf), "hollow: %s   inset = %.2f   (%s)   %s tris",
             g_hollow ? "ON (frame)" : "OFF (solid)", g_inset,
             g_hollow ? "face inset toward center; inner face skipped"
                      : "solid 92-tri mesh",
             g_hollow ? "288" : "92");
    glRasterPos2f(0.02f, 0.91f);
    for (const char* c = buf; *c; ++c) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *c);

    const char* lines[] = {
        ",  .  scrub inset (border)",
        "o    toggle hollow",
        "a    auto-sweep inset",
        "c    write STL",
        "h    toggle HUD",
        "Esc  quit"
    };
    float y = 0.86f;
    for (auto& ln : lines) {
        glRasterPos2f(0.02f, y);
        for (const char* c = ln; *c; ++c) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);
        y -= 0.04f;
    }

    glMatrixMode(GL_MODELVIEW); glPopMatrix();
    glMatrixMode(GL_PROJECTION); glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glEnable(GL_LIGHTING);
}

// ----------------------------------------------------------------------------
// Mouse: orbit / pan / zoom (mirrors main9.cpp / main18).
// ----------------------------------------------------------------------------
void mouseButton(int button, int state, int x, int y)
{
    if (button == GLUT_LEFT_BUTTON)        g_leftDown   = (state == GLUT_DOWN);
    else if (button == GLUT_MIDDLE_BUTTON) g_middleDown = (state == GLUT_DOWN);
    else if (button == GLUT_RIGHT_BUTTON)  g_rightDown  = (state == GLUT_DOWN);
    else if (button == 3) { g_zoom *= 1.05f; glutPostRedisplay(); }   // wheel up
    else if (button == 4) { g_zoom /= 1.05f; glutPostRedisplay(); }   // wheel down
    g_lastX = x; g_lastY = y;
}

void mouseMotion(int x, int y)
{
    int dx = x - g_lastX;
    int dy = y - g_lastY;

    if (g_leftDown) {
        g_angleY += dx * 0.5f;
        g_angleX += dy * 0.5f;
        g_angleX = fmaxf(-90.0f, fminf(90.0f, g_angleX));
    } else if (g_middleDown) {
        g_panX += dx * 0.01f * g_zoom;
        g_panY -= dy * 0.01f * g_zoom;
    } else if (g_rightDown) {
        g_zoom *= 1.0f - dy * 0.005f;
        g_zoom = fmaxf(0.1f, fminf(10.0f, g_zoom));
    }

    g_lastX = x; g_lastY = y;
    glutPostRedisplay();
}

// ----------------------------------------------------------------------------
// Keyboard: inset scrub, hollow toggle, auto-sweep, STL, HUD, quit.
// ----------------------------------------------------------------------------
void keyboard(unsigned char key, int /*x*/, int /*y*/)
{
    switch (key) {
        case 27: exit(0); break;                       // Esc

        case '.': case '>':                            // increase inset (bigger hole / thinner border)
            g_inset += 0.05; if (g_inset > 0.95) g_inset = 0.95;
            cout << "inset = " << g_inset << "\n";
            glutPostRedisplay(); break;

        case ',': case '<':                            // decrease inset (thicker border / smaller hole)
            g_inset -= 0.05; if (g_inset < 0.05) g_inset = 0.05;
            cout << "inset = " << g_inset << "\n";
            glutPostRedisplay(); break;

        case 'o': case 'O':                            // toggle hollow (inset frame vs solid)
            g_hollow = !g_hollow;
            cout << "hollow " << (g_hollow ? "ON (frame)" : "OFF (solid)") << "\n";
            glutPostRedisplay(); break;

        case 'a': case 'A':                            // toggle auto-sweep
            g_autoSweep = !g_autoSweep;
            cout << "auto-sweep " << (g_autoSweep ? "ON" : "OFF") << "\n";
            if (g_autoSweep) glutTimerFunc(33, sweepTimer, 0);
            glutPostRedisplay(); break;

        case 'c': case 'C': {                          // write STL at current hollow/inset
            const string path = "/home/mike666/Downloads/truncated_cuboctahedron.stl";
            cube.writeSTL_s_truncated_cuboctahedron(path, "TruncatedCuboctahedron",
                                                    "full", 0, 9, 2, g_hollow, g_inset);
            cout << "wrote STL: " << path << "  (hollow=" << (g_hollow ? "ON" : "off")
                 << ", inset=" << g_inset << ")\n";
            break;
        }

        case 'h': case 'H':
            g_showHelp = !g_showHelp; glutPostRedisplay(); break;
    }
}

// ----------------------------------------------------------------------------
// Auto-sweep the inset over [0.10, 0.90] (thin <-> thick border). Triangle wave
// so the motion reverses smoothly at both ends. ~30 Hz.
// ----------------------------------------------------------------------------
void sweepTimer(int /*value*/)
{
    if (!g_autoSweep) return;
    g_sweepPhase += 0.03;
    double p = fmod(g_sweepPhase, 2.0 * M_PI);          // 0..2PI
    double tri = (p < M_PI) ? (p / M_PI) : (2.0 - p / M_PI);   // 0..1..0
    g_inset = 0.10 + 0.80 * tri;                        // 0.10 .. 0.90 .. 0.10
    glutPostRedisplay();
    glutTimerFunc(33, sweepTimer, 0);
}