// main_cuboctahedron_lattice.cpp
//
// Octagon-connected cuboctahedron lattice with sphere inversion.
//
// A truncated cuboctahedron is placed at EVERY subcell of the subdivided cube.
// Each cuboctahedron's 6 octagons sit on its 6 subcell faces, so two
// face-adjacent cuboctahedra share a coincident octagon on the shared face — the
// octagons become the "windows" that connect the lattice. The hexagon + square
// frames bound the octahedral voids (at lattice vertices) and cubic voids (at
// lattice edges) left open between the cells. This is the cantitruncated honeycomb
// skeleton with the filler cells (truncated octahedra + cubes) removed — only the
// cuboctahedra, stuck together through their octagons.
//
// To avoid z-fighting / wasted triangles on the coincident shared octagons, each
// shared octagon is drawn ONCE (Cube::getCuboctahedronLatticeFacets culls a cell's
// -x/-y/-z octagon when the neighbour on that side exists; that neighbour draws it
// via its + face). Hexagon/square frames are never coincident (distinct void-
// bounding planes) so they are always emitted.
//
// Inversion: a single sphere inversion (Mobius) at a mouse-driven center deforms
// every subcell's 8 corner vertices (captured once as a pristine identity lattice,
// recomputed from it each change so there is no frame-over-frame compounding).
// Because each cell is rebuilt from its (deformed) subcell corners, the whole
// lattice deforms together with no special handling.
//
// Controls:
//   L-drag (plain)  move inversion center in XY on the current z-plane
//   Alt + L-drag    move inversion center in Z
//   Shift + L-drag  pan camera
//   M-drag / wheel  orbit camera          |  R-drag  zoom          |  wheel  zoom
//   r / R           decrease / increase inversion radius
//   0               reset inversion center + radius
//   i               toggle inversion on/off (off = identity lattice, no deformation)
//   o               toggle hollow (on = octagonal windows; off = solid cuboctahedra)
//   ,  /  .         decrease / increase hollow inset (border thickness)
//   c               write an STL of the current (deformed) lattice to ~/Downloads
//   h               toggle HUD      |  Esc  quit
//
// Build:
//   g++ -std=c++17 main_cuboctahedron_lattice.cpp -o main_cuboctahedron_lattice \
//       -lGL -lGLU -lglut -lm
//   ./main_cuboctahedron_lattice          # xvfb-run ./... if headless

#include <vector>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string>
#include <fstream>
#include <GL/glut.h>
#include "Vector3D.cpp"
#include "Quaternion.cpp"
#include "Facet.cpp"
#include "FacetBox.hpp"
#include "Cube.hpp"

using namespace std;

// ----------------------------------------------------------------------------
// Scene: one subdivided cube carrying a cuboctahedron at every subcell.
// ----------------------------------------------------------------------------
int   cube_dim = 2.0;
int   N        = 6;                                   // cuboctahedra per axis (N^3 cells)
Cube  cube(cube_dim, Vector3D{0, 0.0001, 0}, N);

// Inversion state (single mouse-driven sphere inversion).
static Vector3D g_invCenter(0, 0, 0);
static double   g_invRadius = 1.0;
static bool     g_inversionOn = true;     // i: toggle -> off renders the identity lattice
static bool     g_geomDirty = true;       // recompute lattice when center/radius/on-off changes
static unsigned g_geomVersion = 0;       // cache key for the batched mesh
static std::vector<float> g_identityPositions;   // pristine identity lattice (captured once)

// Lattice knobs.
static bool   g_hollow = true;            // on = octagonal windows; off = solid cuboctahedra
static double g_inset  = 0.5;             // hollow border inset ratio (0,1)

// ----------------------------------------------------------------------------
// Interactive camera state.
// ----------------------------------------------------------------------------
static float g_angleX = 20.0f, g_angleY = -30.0f;
static float g_zoom   = 1.0f;
static float g_panX   = 0.0f, g_panY = 0.0f;
static int   g_lastX  = 0, g_lastY = 0;
static bool  g_leftDown = false, g_middleDown = false, g_rightDown = false;
static bool  g_showHelp = true;
static GLdouble g_modelview[16]  = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
static GLdouble g_projection[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
static GLint    g_viewport[4]    = {0,0,720,720};

// Forward declarations.
void initGL();
void reshape(int w, int h);
void display();
void Draw();
void drawHUD();
void drawAxes(float length);
void mouseButton(int button, int state, int x, int y);
void mouseMotion(int x, int y);
void keyboard(unsigned char key, int x, int y);
void Setup();

// ----------------------------------------------------------------------------
// Single sphere inversion at g_invCenter, regularized so a vertex almost on the
// center cannot "explode" to infinity. A plain inversion sends points very near
// the center off to ~r^2/|d| -> huge (a subcell vertex 1e-4 from the center would
// jump ~1e4 units and drag its facets into a spike). Clamp the effective squared
// distance to minDist^2 = (r^2 / g_maxInvDist)^2 so NO deformed point lands
// farther than g_maxInvDist from the center. Points beyond minDist invert exactly
// as before; only the near-singular core is tamed (it stays near the center
// instead of shooting away). The clamp is continuous with the normal inversion
// at the boundary (|sigma| = r^2/|d| = g_maxInvDist there), so it is smooth.
// ----------------------------------------------------------------------------
static double g_maxInvDist = 2.0;   // deformed points stay within this of the center (set in Setup to cube_dim)

static inline Vector3D sigmaCenter(const Vector3D& p) {
    Vector3D d = p - g_invCenter;
    double ds = d * d;
    if (ds < 1e-12) return p;
    const double minDist = (g_invRadius * g_invRadius) / g_maxInvDist;   // saturate closer than this
    const double minDs2  = minDist * minDist;
    const double effDs   = (ds > minDs2) ? ds : minDs2;
    const double f = (g_invRadius * g_invRadius) / effDs;
    return g_invCenter + f * d;
}

// Capture the pristine identity lattice ONCE, before any deformation.
static void captureIdentityLattice() {
    cube.fillVertexLattice(g_identityPositions);   // n^3*8*3 floats
    cout << "[anim] identity lattice captured: " << (g_identityPositions.size() / 3)
         << " verts (" << (g_identityPositions.size() * sizeof(float) / (1024.0 * 1024.0))
         << " MB)\n";
}

// Push the single mouse-driven inversion of the identity lattice back into the
// Cube's subcells (same iteration order fillVertexLattice uses: physical
// [0,n)^3, vertex 0..7). Recomputed only when g_geomDirty is set; when inversion
// is off (g_inversionOn=false) the pristine identity positions are written back.
static void updateAnimatedGeometry() {
    const int n = cube.getSubdivisionLevels();
    const int center = n / 2;
    const float* id = g_identityPositions.data();
    size_t idx = 0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k) {
                const int lx = i - center, ly = j - center, lz = k - center;
                for (int l = 0; l < 8; ++l) {
                    Vector3D p(id[idx], id[idx + 1], id[idx + 2]);
                    cube.updateSubCellVertex(lx, ly, lz, l, g_inversionOn ? sigmaCenter(p) : p);
                    idx += 3;
                }
            }
    ++g_geomVersion;   // invalidate the batched-mesh cache
}

// ----------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------
int main(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(720, 720);
    glutCreateWindow(" Octagon-connected Cuboctahedron Lattice  (hollow + inversion)  U.U ");

    initGL();
    Setup();

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouseButton);
    glutMotionFunc(mouseMotion);
    glutKeyboardFunc(keyboard);

    glutPostRedisplay();
    glutMainLoop();
    return 0;
}

// ----------------------------------------------------------------------------
// One-time info + capture/deform + an initial STL of the identity lattice.
// ----------------------------------------------------------------------------
void Setup()
{
    cout << "\n--- Octagon-connected cuboctahedron lattice ---\n";
    cout << "A truncated cuboctahedron at every one of the " << N*N*N
         << " subcells; neighbours share an octagonal face (window).\n";
    cout << "Octahedral voids at lattice vertices + cubic voids at lattice edges left open.\n";
    cout << "hollow=" << (g_hollow ? "ON (octagonal windows)" : "OFF (solid)")
         << "  inset=" << g_inset << "\n";
    cout << "Controls: L-drag center XY | Alt+L Z | Shift+L pan | M-drag orbit | wheel zoom\n";
    cout << "          r/R radius | 0 reset | i inversion on/off | o hollow | ,. inset | c STL | h HUD | Esc quit\n\n";

    // Bound the deformed lattice to within the cube radius of the inversion center
    // so a vertex almost on the center can't explode to infinity (see sigmaCenter).
    g_maxInvDist = cube_dim;

    captureIdentityLattice();
    updateAnimatedGeometry();

    // Write an STL of the (current, t=0) lattice so there is immediate output.
    const string path = "/home/mike666/Downloads/cuboctahedron_lattice.stl";
    cube.writeSTL_s_cuboctahedron_lattice(path, "CuboctahedronLattice", g_hollow, g_inset);
    cout << "wrote initial STL: " << path << "\n";
}

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
    glEnable(GL_NORMALIZE);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

    glClearColor(1, 1, 1, 1);
    printf("GL renderer: %s\n", (const char*)glGetString(GL_RENDERER));
}

// ----------------------------------------------------------------------------
void reshape(int w, int h)
{
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(50.0, (double)w / (double)h, 0.1, 1000.0);
    glMatrixMode(GL_MODELVIEW);
}

// ----------------------------------------------------------------------------
void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    glTranslatef(g_panX, g_panY, -5.0f * g_zoom);
    glRotatef(g_angleX, 1, 0, 0);
    glRotatef(g_angleY, 0, 1, 0);

    // Cache camera matrices for mouse->world unprojection (left-drag center).
    glGetDoublev(GL_MODELVIEW_MATRIX, g_modelview);
    glGetDoublev(GL_PROJECTION_MATRIX, g_projection);
    glGetIntegerv(GL_VIEWPORT, g_viewport);

    // Deform the lattice to the current inversion center/radius only when changed.
    if (g_geomDirty) { updateAnimatedGeometry(); g_geomDirty = false; }

    drawAxes(4.0f);
    Draw();
    drawHUD();

    glutSwapBuffers();
}

// ----------------------------------------------------------------------------
// Batched render: cache vertex/normal arrays, rebuild only when the knobs
// (hollow/inset) or the geometry version (inversion) change — camera-only frames
// reuse them. One glDrawArrays for fills + one for outlines.
// ----------------------------------------------------------------------------
void Draw()
{
    static std::vector<float> s_verts, s_norms;
    static size_t   s_N = 0;
    static unsigned s_geom = 0xFFFFFFFFu;
    static double   s_inset = -1.0;
    static bool     s_hollow = false;

    if (s_hollow != g_hollow || s_inset != g_inset || s_geom != g_geomVersion) {
        FacetBox fb = cube.getCuboctahedronLatticeFacets(g_hollow, g_inset);
        s_N = fb.size();
        s_verts.resize(s_N * 9);
        s_norms.resize(s_N * 9);
        for (size_t i = 0; i < s_N; ++i) {
            const Facet& f = fb[i];
            Vector3D n = f.getNormal();
            Vector3D A = f[0], B = f[1], C = f[2];
            float* v  = &s_verts[i * 9];
            float* nv = &s_norms[i * 9];
            v[0]=(float)A.x(); v[1]=(float)A.y(); v[2]=(float)A.z();
            v[3]=(float)B.x(); v[4]=(float)B.y(); v[5]=(float)B.z();
            v[6]=(float)C.x(); v[7]=(float)C.y(); v[8]=(float)C.z();
            for (int k = 0; k < 3; ++k) {
                nv[k*3]   = (float)n.x();
                nv[k*3+1] = (float)n.y();
                nv[k*3+2] = (float)n.z();
            }
        }
        s_hollow = g_hollow; s_inset = g_inset; s_geom = g_geomVersion;
    }

    if (s_N > 0) {
        glPushAttrib(GL_COLOR_BUFFER_BIT | GL_POLYGON_BIT | GL_ENABLE_BIT);
        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_NORMAL_ARRAY);
        glVertexPointer(3, GL_FLOAT, 0, s_verts.data());
        glNormalPointer(GL_FLOAT, 0, s_norms.data());
        // Filled gray triangles (polygon offset so outlines don't z-fight).
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(1.0f, 1.0f);
        glColor3ub(200, 200, 200);
        glDrawArrays(GL_TRIANGLES, 0, (GLsizei)(s_N * 3));
        glDisable(GL_POLYGON_OFFSET_FILL);
        // Black outlines: one draw call via line polygon mode.
        glColor3ub(0, 0, 0);
        glLineWidth(1.0f);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glDrawArrays(GL_TRIANGLES, 0, (GLsizei)(s_N * 3));
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glDisableClientState(GL_NORMAL_ARRAY);
        glDisableClientState(GL_VERTEX_ARRAY);
        glPopAttrib();
    }
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
void drawHUD()
{
    if (!g_showHelp) return;

    glMatrixMode(GL_PROJECTION); glPushMatrix();
    glLoadIdentity(); glOrtho(0, 1, 0, 1, -1, 1);
    glMatrixMode(GL_MODELVIEW); glPushMatrix(); glLoadIdentity();

    glDisable(GL_LIGHTING);
    glColor3f(0, 0, 0);

    char buf[220];
    snprintf(buf, sizeof(buf), "octagon-connected cuboctahedron lattice  (%dx%dx%d cells, voids open)",
             N, N, N);
    glRasterPos2f(0.02f, 0.95f);
    for (const char* c = buf; *c; ++c) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *c);

    snprintf(buf, sizeof(buf), "inversion: %s   InvCenter: (%.4f, %.4f, %.4f)  r=%.3f",
             g_inversionOn ? "ON" : "OFF", g_invCenter.x(), g_invCenter.y(), g_invCenter.z(), g_invRadius);
    glColor3f(0.0f, 0.6f, 0.0f);
    glRasterPos2f(0.02f, 0.91f);
    for (const char* c = buf; *c; ++c) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);

    snprintf(buf, sizeof(buf), "hollow: %s  inset=%.2f  (%s)",
             g_hollow ? "ON (octagon windows)" : "OFF (solid)", g_inset,
             g_hollow ? "see through the lattice along the octagons"
                      : "cuboctahedra stuck face-to-face via octagons");
    glColor3f(0.2f, 0.6f, 0.6f);
    glRasterPos2f(0.02f, 0.87f);
    for (const char* c = buf; *c; ++c) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);

    const char* lines[] = {
        "L-drag: center XY | Alt+L: Z | Shift+L: pan | M-drag: orbit | wheel: zoom",
        "r/R: radius -/+  |  0: reset  |  i: inversion on/off  |  o: hollow  |  ,. : inset",
        "c: write STL  |  h: toggle HUD  |  Esc: quit"
    };
    float y = 0.83f;
    glColor3f(0.5f, 0.5f, 0.5f);
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
// Unproject the screen cursor to the horizontal plane z = g_invCenter.z() and set
// the center's XY there. y is flipped (GLUT top-down -> GL bottom-up).
// ----------------------------------------------------------------------------
static void placeInvCenterXY(int x, int y) {
    GLdouble winY = (GLdouble)(g_viewport[3] - y);
    GLdouble ox, oy, oz, fx, fy, fz;
    if (!gluUnProject((GLdouble)x, winY, 0.0, g_modelview, g_projection, g_viewport, &ox, &oy, &oz)) return;
    if (!gluUnProject((GLdouble)x, winY, 1.0, g_modelview, g_projection, g_viewport, &fx, &fy, &fz)) return;
    Vector3D O(ox, oy, oz), F(fx, fy, fz);
    Vector3D D = F - O;
    if (fabs(D.z()) < 1e-9) return;                 // ray parallel to the plane: no hit
    double t = (g_invCenter.z() - O.z()) / D.z();
    Vector3D P = O + D * t;
    g_invCenter = Vector3D(P.x(), P.y(), g_invCenter.z());
    g_geomDirty = true;
}

// ----------------------------------------------------------------------------
void mouseButton(int button, int state, int x, int y)
{
    int mods = glutGetModifiers();
    if (button == GLUT_LEFT_BUTTON) {
        g_leftDown = (state == GLUT_DOWN);
        if (state == GLUT_DOWN && !(mods & GLUT_ACTIVE_SHIFT) && !(mods & GLUT_ACTIVE_ALT))
            placeInvCenterXY(x, y);                  // plain left press: snap center to cursor
    }
    else if (button == GLUT_MIDDLE_BUTTON) g_middleDown = (state == GLUT_DOWN);
    else if (button == GLUT_RIGHT_BUTTON)  g_rightDown  = (state == GLUT_DOWN);
    else if (button == 3) { g_zoom *= 1.05f; glutPostRedisplay(); }    // wheel up
    else if (button == 4) { g_zoom /= 1.05f; glutPostRedisplay(); }   // wheel down
    g_lastX = x; g_lastY = y;
}

void mouseMotion(int x, int y)
{
    int dx = x - g_lastX;
    int dy = y - g_lastY;
    int mods = glutGetModifiers();

    if (g_middleDown) {                              // wheel-click + hold = orbit
        g_angleY += dx * 0.5f;
        g_angleX += dy * 0.5f;
        g_angleX = fmaxf(-90.0f, fminf(90.0f, g_angleX));
    }
    else if (g_leftDown) {
        if (mods & GLUT_ACTIVE_ALT) {                // Alt + left: move center in Z
            double nz = g_invCenter.z() - dy * 0.01 * (double)g_zoom;
            g_invCenter = Vector3D(g_invCenter.x(), g_invCenter.y(), nz);
            g_geomDirty = true;
        }
        else if (mods & GLUT_ACTIVE_SHIFT) {        // Shift + left: pan
            g_panX += dx * 0.01f * g_zoom;
            g_panY -= dy * 0.01f * g_zoom;
        }
        else {                                       // plain left: center XY
            placeInvCenterXY(x, y);
        }
    }
    else if (g_rightDown) {                          // right: zoom
        g_zoom *= 1.0f - dy * 0.005f;
        g_zoom = fmaxf(0.1f, fminf(10.0f, g_zoom));
    }

    g_lastX = x; g_lastY = y;
    glutPostRedisplay();
}

// ----------------------------------------------------------------------------
void keyboard(unsigned char key, int /*x*/, int /*y*/)
{
    switch (key) {
        case 27: exit(0); break;                     // Esc

        case 'r': g_invRadius -= 0.05; if (g_invRadius < 0.05) g_invRadius = 0.05;
                  g_geomDirty = true; cout << "inv radius = " << g_invRadius << "\n";
                  glutPostRedisplay(); break;
        case 'R': g_invRadius += 0.05; if (g_invRadius > 4.0) g_invRadius = 4.0;
                  g_geomDirty = true; cout << "inv radius = " << g_invRadius << "\n";
                  glutPostRedisplay(); break;

        case '0': g_invCenter = Vector3D(0.0, 0.0, 1e-9); g_invRadius = 0.5;
                  g_geomDirty = true; cout << "inversion reset\n";
                  glutPostRedisplay(); break;

        case 'i': case 'I': g_inversionOn = !g_inversionOn; g_geomDirty = true;
                  cout << "inversion " << (g_inversionOn ? "ON" : "OFF (identity lattice)") << "\n";
                  glutPostRedisplay(); break;

        case 'o': case 'O': g_hollow = !g_hollow;
                  cout << "hollow " << (g_hollow ? "ON (octagon windows)" : "OFF (solid)") << "\n";
                  glutPostRedisplay(); break;

        case '.': case '>': g_inset += 0.05; if (g_inset > 0.95) g_inset = 0.95;
                  cout << "inset = " << g_inset << "\n"; glutPostRedisplay(); break;
        case ',': case '<': g_inset -= 0.05; if (g_inset < 0.05) g_inset = 0.05;
                  cout << "inset = " << g_inset << "\n"; glutPostRedisplay(); break;

        case 'c': case 'C': {
            const string path = "/home/mike666/Downloads/cuboctahedron_lattice.stl";
            cube.writeSTL_s_cuboctahedron_lattice(path, "CuboctahedronLattice", g_hollow, g_inset);
            cout << "wrote STL: " << path << "  (hollow=" << (g_hollow ? "ON" : "off")
                 << ", inset=" << g_inset << ")\n";
            break;
        }

        case 'h': case 'H': g_showHelp = !g_showHelp; glutPostRedisplay(); break;
    }
}
