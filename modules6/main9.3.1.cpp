#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <set>
#include <unordered_map>
#include <stdexcept>
#include <stdio.h>
#include <math.h>
#include <GL/glut.h>
#include "Vector3D.cpp"
#include "Vector4D.cpp"
#include "Quaternion.cpp"
#include "Facet.cpp"
#include "FacetBox.hpp"

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

// Blender-style vertex selection (forward declarations — defined later)
const FacetBox& currentActiveBox();
void buildUniqueVertexList();
void drawVertexMarkers();
int  pickVertex(int mouseX, int mouseY, double pixelThreshold);
void handleVertexPick(int mouseX, int mouseY, bool shift);

Vector3D d_ui(1.0, 0.0, 0.0);
Vector3D e_ui(0.0, 1.0, 0.0);
Vector3D f_ui(0.0, 0.0, 1.0);
Facet f_1(d_ui, e_ui, f_ui);

// Loaded STL mesh, centered+scaled into the open unit ball (Klein disk model).
FacetBox g_boxKlein;
// Same mesh after applyHyperboloid() — projected to the Poincaré disk.
FacetBox g_boxHyper;
// Round-trip: inverse map applied to g_boxHyper, back to the Klein disk.
FacetBox g_boxBack;
// Active view: 0 = Original (Klein), 1 = Hyperbolic (Poincaré), 2 = Round-trip.
int g_view = 0;

//////////////////////////////////////
//                                  //
//  BLENDER-STYLE VERTEX SELECTION  //
//                                  //
//////////////////////////////////////
// A unique mesh vertex is identified by ONE representative facet/corner; the
// same logical vertex is shared by several facets, so we keep a single ref to
// store each mesh point once. Selection is by vertex identity, therefore it
// survives switching between the Klein / Poincaré / round-trip views.
struct VertRef { int facetIdx; int cornerIdx; };
std::vector<VertRef> g_uniqueVerts;    // one representative per unique vertex (built from g_boxKlein)
std::set<int>        g_selectedVerts;  // indices into g_uniqueVerts currently highlighted
const double g_pickRadiusPx = 10.0;    // screen-space click tolerance (pixels)
bool         g_showVerts    = false;   // [V] toggle the faint all-vertex dot overlay

// Click-vs-drag bookkeeping so a left click picks a vertex while a left drag
// still orbits the camera (Blender treats a near-stationary press as a click).
int  g_clickStartX = 0, g_clickStartY = 0;
bool g_maybeClick   = false;   // left button went down, not yet a drag
bool g_dragged      = false;   // movement exceeded the threshold → orbit, not pick
const int g_clickThreshold = 5;


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


//-----------------------------------------------------------------------------
// ASCII STL loader → FacetBox.
// Parses every "vertex x y z" line; groups each 3 consecutive vertices into a
// Facet. Facet's ctor recomputes the normal, so the STL "facet normal" lines
// are ignored. Throws std::runtime_error on a missing file or zero facets.
//-----------------------------------------------------------------------------
FacetBox loadSTL_ascii(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("loadSTL_ascii: cannot open " + path);

    FacetBox box;
    std::string line;
    std::vector<Vector3D> verts;   // accumulate 3 vertices per facet
    while (std::getline(in, line)) {
        std::istringstream iss(line);
        std::string tok;
        iss >> tok;
        if (tok == "vertex") {
            double x, y, z;
            if (!(iss >> x >> y >> z))
                throw std::runtime_error("loadSTL_ascii: malformed vertex line");
            verts.emplace_back(x, y, z);
            if (verts.size() == 3) {
                box.push(verts[0], verts[1], verts[2]);   // FacetBox.hpp:159
                verts.clear();
            }
        }
    }
    if (box.size() == 0)
        throw std::runtime_error("loadSTL_ascii: no facets parsed from " + path);
    return box;
}

//-----------------------------------------------------------------------------
// Center a mesh on its vertex centroid (FacetBox::center(), FacetBox.hpp:222)
// and uniformly scale it so the farthest vertex lies at distance `target` (<1)
// from the origin. This places every vertex strictly inside the open unit ball,
// which is required by Quaternion::toHyperboloid() (it throws when ‖x‖ ≥ 1).
//-----------------------------------------------------------------------------
FacetBox centerAndScale(const FacetBox& src, double target = 0.9) {
    Vector3D c = src.center();

    // Pass 1: max distance from the centroid.
    double maxDist = 0.0;
    for (size_t i = 0; i < src.size(); ++i) {
        const Facet& f = src[i];
        for (int k = 0; k < 3; ++k) {
            double dd = abs(f[k] - c);          // abs(Vector3D) — Vector3D.hpp:56
            if (dd > maxDist) maxDist = dd;
        }
    }
    if (maxDist < 1e-12) maxDist = 1.0;         // guard degenerate/empty input

    // Pass 2: build a new box, each vertex = (v - c) * (target / maxDist).
    double s = target / maxDist;
    FacetBox out;
    for (size_t i = 0; i < src.size(); ++i) {
        const Facet& f = src[i];
        Vector3D a = s * (f[0] - c);             // double * Vector3D — Vector3D.hpp:48
        Vector3D b = s * (f[1] - c);
        Vector3D d = s * (f[2] - c);
        out.push(a, b, d);
    }
    return out;
}

//-----------------------------------------------------------------------------
// Inverse of Facet::applyHyperboloid() (Poincaré disk → Klein disk).
// toHyperboloid() maps x ↦ y = x / (1 + √(1−‖x‖²)); its inverse is
// x = 2y / (1 + ‖y‖²). For ‖y‖ < 1 the result stays ‖x‖ < 1, so the round-trip
// never leaves the unit ball.
//-----------------------------------------------------------------------------
Vector3D invHyperVertex(const Vector3D& y) {
    double n2 = y * y;                  // dot product = ‖y‖² — Vector3D.hpp:50
    double f = 2.0 / (1.0 + n2);
    return f * y;                       // double * Vector3D — Vector3D.hpp:48
}

FacetBox inverseHyperboloid(const FacetBox& p) {
    FacetBox out;
    for (size_t i = 0; i < p.size(); ++i) {
        const Facet& f = p[i];
        out.push(invHyperVertex(f[0]),
                 invHyperVertex(f[1]),
                 invHyperVertex(f[2]));
    }
    return out;
}


void Setup() {

	if (ciclo == 0) {

        cout << "\n———————————————————————————————————————————————————————————————————————\n";
        cout <<   "|- Hyperbolic STL Mesh (cubo.stl)———————————————————————————————————————————\n";
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
        cout << "     \"-_\\_______/;;'        ~[Hyperbolic STL Mesh]\n\n";

        // 1) Load the external mesh and fit it inside the open unit ball.
        //    toHyperboloid() throws when ‖x‖ ≥ 1, so we center + scale first.
        try {
            g_boxKlein = centerAndScale(
                loadSTL_ascii("/home/mike666/Downloads/cubo.stl"), 0.9);

            // Remesh before projecting: one Midpoint4 pass (×4 triangles) then
            // one Centroid3 pass (×3). Finer triangles make the Poincaré disk
            // curve the edges smoothly instead of bending a few large chords.
            // New vertices are convex combos of vertices already inside the
            // unit ball, so the mesh stays inside it (still safe to project).
            g_boxKlein = g_boxKlein.refine(1, FacetBox::SubdivisionMode::Midpoint4);
            g_boxKlein = g_boxKlein.refine(1, FacetBox::SubdivisionMode::Centroid3);

            // 2) Send it into hyperbolic space: Klein disk → Poincaré disk.
            g_boxHyper = g_boxKlein.hyperboloid();   // FacetBox.hpp:257

            // 3) Bring it back: apply the inverse map (Poincaré → Klein).
            g_boxBack = inverseHyperboloid(g_boxHyper);

            std::cout << "Loaded cubo.stl: "
                      << g_boxKlein.size() << " facets (Klein), "
                      << g_boxHyper.size() << " (Poincare), "
                      << g_boxBack.size()  << " (round-trip)\n";

            // 4) Index every unique mesh vertex so the user can click-select
            //    them Blender-style. Built from the Klein mesh (the canonical
            //    geometry); the same vertex ids map onto the other two views.
            buildUniqueVertexList();
        } catch (const std::exception& e) {
            std::cerr << "Setup error: " << e.what() << "\n";
            // leave boxes empty → Draw() renders nothing (size()==0 loops are no-ops)
        }

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


        const FacetBox& active = (g_view == 0) ? g_boxKlein
                              : (g_view == 1) ? g_boxHyper
                                             : g_boxBack;
        size_t total = active.size();
        for(size_t i = 0; i < total; ++i) {
            // pick hue from 0°→360° across the range
            float hue = float(i) / float(total) * 360.0f;
            Color c = hsv2rgb(hue, 0.8f, 1.0f);   // 80% saturation, full value
            // convert to 0–255 ints
            int R = int(c.r * 255), G = int(c.g * 255), B = int(c.b * 255);
            drawFacet(active[i], R, G, B, 1.0f);
        }

        // Overlay vertex dots: gray = unselected, bright-orange = selected.
        drawVertexMarkers();


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
void keyboard(unsigned char key, int x, int y);

//=============================================================================
//  BLENDER-STYLE VERTEX SELECTION
//=============================================================================
// Which of the three meshes is currently on screen.
const FacetBox& currentActiveBox() {
    return (g_view == 0) ? g_boxKlein
         : (g_view == 1) ? g_boxHyper
                         : g_boxBack;
}

// Walk every facet corner of the Klein mesh and keep one representative per
// unique spatial vertex. Shared corners (the same point used by several
// facets) collapse to a single entry, so the selection list holds each mesh
// point once. A spatial-hash grid makes this O(n) instead of O(n²) — essential
// here because the refined mesh has ~10^5 corners. Because the Poincaré and
// round-trip meshes are pointwise images of the Klein mesh, the grouping is
// valid for all three views: a selected vertex stays selected when you switch
// between [1]/[2]/[3].
void buildUniqueVertexList() {
    g_uniqueVerts.clear();
    g_selectedVerts.clear();
    const FacetBox& box = g_boxKlein;
    if (box.size() == 0) { std::cout << "Unique vertices indexed: 0\n"; return; }

    const double tol = 1e-6;       // merge tolerance (shared verts are bit-identical)
    const double cs  = 1e-3;       // grid cell: >> tol, << smallest vertex spacing
    const long long OFF = (1LL << 20);
    auto pack = [&](long long rx, long long ry, long long rz) -> long long {
        return (rx + OFF) | ((ry + OFF) << 21) | ((rz + OFF) << 42);
    };
    // cell key -> indices into g_uniqueVerts that live in that cell
    std::unordered_map<long long, std::vector<int>> grid;
    grid.reserve(1 << 15);

    for (size_t fi = 0; fi < box.size(); ++fi) {
        const Facet& f = box[fi];
        for (int ci = 0; ci < 3; ++ci) {
            const Vector3D& p = f[ci];
            long long rx = (long long)floor(p.x() / cs);
            long long ry = (long long)floor(p.y() / cs);
            long long rz = (long long)floor(p.z() / cs);

            // A duplicate, if any, sits in this cell or one of its 26 neighbours
            // (covers the tolerance even when a point straddles a cell boundary).
            bool dup = false;
            for (int dz = -1; dz <= 1 && !dup; ++dz)
            for (int dy = -1; dy <= 1 && !dup; ++dy)
            for (int dx = -1; dx <= 1 && !dup; ++dx) {
                auto it = grid.find(pack(rx + dx, ry + dy, rz + dz));
                if (it == grid.end()) continue;
                for (int j : it->second) {
                    const VertRef& r = g_uniqueVerts[j];
                    if (abs(p - box[r.facetIdx][r.cornerIdx]) < tol) { dup = true; break; }
                }
            }
            if (!dup) {
                int newIdx = (int)g_uniqueVerts.size();
                g_uniqueVerts.push_back({(int)fi, ci});
                grid[pack(rx, ry, rz)].push_back(newIdx);
            }
        }
    }
    std::cout << "Unique vertices indexed: " << g_uniqueVerts.size() << "\n";
}

// Draw the selected vertices as bright-orange dots with a thin black outline,
// always on top (depth test off) so they stay visible through the mesh — like
// Blender's selected verts. The optional faint-gray overlay of *every* vertex
// (Blender's full edit-mode dot field) is gated behind [V].
void drawVertexMarkers() {
    if (g_uniqueVerts.empty()) return;
    const FacetBox& active = currentActiveBox();

    glPushAttrib(GL_LIGHTING_BIT | GL_DEPTH_BUFFER_BIT |
                 GL_COLOR_BUFFER_BIT | GL_POINT_BIT);
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);            // selected verts always visible (X-ray)
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_POINT_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

    // Optional overlay of EVERY vertex as faint gray dots (Blender edit-mode
    // look). Toggled with [V]; off by default so it doesn't clutter the mesh.
    if (g_showVerts) {
        glPointSize(3.0f);
        glColor4ub(60, 60, 60, 120);
        glBegin(GL_POINTS);
        for (size_t i = 0; i < g_uniqueVerts.size(); ++i) {
            if (g_selectedVerts.count((int)i)) continue;
            const VertRef& r = g_uniqueVerts[i];
            const Vector3D& p = active[r.facetIdx][r.cornerIdx];
            glVertex3f(p.x(), p.y(), p.z());
        }
        glEnd();
    }

    // Selected vertices: black outline + bright-orange core.
    if (!g_selectedVerts.empty()) {
        glPointSize(14.0f);
        glColor3ub(0, 0, 0);
        glBegin(GL_POINTS);
        for (int idx : g_selectedVerts) {
            const VertRef& r = g_uniqueVerts[idx];
            const Vector3D& p = active[r.facetIdx][r.cornerIdx];
            glVertex3f(p.x(), p.y(), p.z());
        }
        glEnd();

        glPointSize(9.0f);
        glColor3ub(255, 140, 0);     // bright orange
        glBegin(GL_POINTS);
        for (int idx : g_selectedVerts) {
            const VertRef& r = g_uniqueVerts[idx];
            const Vector3D& p = active[r.facetIdx][r.cornerIdx];
            glVertex3f(p.x(), p.y(), p.z());
        }
        glEnd();
    }

    glPopAttrib();
}

// Project every unique vertex to screen space and return the index of the one
// nearest the cursor within pixelThreshold — front-most (smallest depth) wins
// ties, matching Blender's "pick the closest to the view". Returns -1 when no
// vertex is close enough (i.e. the click landed on empty space).
int pickVertex(int mouseX, int mouseY, double pixelThreshold) {
    if (g_uniqueVerts.empty()) return -1;
    const FacetBox& active = currentActiveBox();

    GLint    viewport[4];
    GLdouble modelview[16], projection[16];

    // Rebuild the exact projection used by reshape()...
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glGetIntegerv(GL_VIEWPORT, viewport);
    double h = viewport[3] ? (double)viewport[3] : 1.0;
    gluPerspective(100.0, (double)viewport[2] / h, 0.001, 1000.0);
    glGetDoublev(GL_PROJECTION_MATRIX, projection);
    glPopMatrix();

    // ...and the exact camera modelview used by display().
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glTranslatef(0, 0, -5.0f * g_zoom);
    glRotatef(g_angleX, 1, 0, 0);
    glRotatef(g_angleY, 0, 1, 0);
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
    glPopMatrix();

    int    bestIdx  = -1;
    double bestWinZ = 1e9;
    double bestDist = 1e9;
    for (size_t i = 0; i < g_uniqueVerts.size(); ++i) {
        const VertRef& r = g_uniqueVerts[i];
        const Vector3D& p = active[r.facetIdx][r.cornerIdx];
        GLdouble winX, winY, winZ;
        gluProject(p.x(), p.y(), p.z(),
                   modelview, projection, viewport,
                   &winX, &winY, &winZ);
        // GLUT mouse Y is top-down; GL window Y is bottom-up — flip to compare.
        double dx = winX - (double)mouseX;
        double dy = winY - (double)(viewport[3] - mouseY);
        double dist = sqrt(dx * dx + dy * dy);
        if (dist <= pixelThreshold) {
            if (winZ < bestWinZ - 1e-9 ||
                (fabs(winZ - bestWinZ) <= 1e-9 && dist < bestDist)) {
                bestWinZ = winZ;
                bestDist = dist;
                bestIdx  = (int)i;
            }
        }
    }
    return bestIdx;
}

// Apply Blender selection semantics to a left click at (mouseX, mouseY):
//   * click a vertex        -> select only that one (replace the selection)
//   * shift+click a vertex  -> toggle it in/out of the selection (multi-select)
//   * click empty space     -> forget the whole selection
//   * shift+click empty     -> keep the current selection
void handleVertexPick(int mouseX, int mouseY, bool shift) {
    int picked = pickVertex(mouseX, mouseY, g_pickRadiusPx);
    if (picked >= 0) {
        if (shift) {
            if (g_selectedVerts.count(picked)) g_selectedVerts.erase(picked);
            else                               g_selectedVerts.insert(picked);
        } else {
            g_selectedVerts.clear();
            g_selectedVerts.insert(picked);
        }
        std::cout << "Vertex " << picked << " -> selected total: "
                  << g_selectedVerts.size() << "\n";
    } else {
        if (!shift && !g_selectedVerts.empty()) {
            std::cout << "Cleared " << g_selectedVerts.size()
                      << " selected vertices\n";
            g_selectedVerts.clear();
        }
    }
    glutPostRedisplay();
}

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
    gluPerspective(100.0, (double)w/h, 0.001, 1000.0);
    glMatrixMode(GL_MODELVIEW);
}

//-----------------------------------------------------------------------------
// Draw help overlay
//-----------------------------------------------------------------------------
void drawHUD() {
    if (!g_showHelp) return;
    const char* lines[] = {
        "L-click: Select vertex (orange)",
        "Shift+L-click: Add/Remove vertex",
        "L-click empty: Clear selection",
        "L-drag: Rotate",
        "M-drag: Pan",
        "R-drag/Wheel: Zoom",
        "[1]: View: Original (Klein)",
        "[2]: View: Hyperbolic (Poincare)",
        "[3]: View: Round-trip Back",
        "[V]: Toggle vertex dots",
        "[C]: Clear selection",
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

    // Current view mode label
    const char* modeName = (g_view == 0) ? "Mode: Original (Klein)"
                         : (g_view == 1) ? "Mode: Hyperbolic (Poincare)"
                                         : "Mode: Round-trip Back";
    glRasterPos2f(0.02f, y);
    for(const char* c=modeName; *c; ++c)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);

    // Live selection count
    y -= 0.05f;
    char selBuf[64];
    snprintf(selBuf, sizeof(selBuf), "Selected: %lu vertices",
             (unsigned long)g_selectedVerts.size());
    glRasterPos2f(0.02f, y);
    for(const char* c=selBuf; *c; ++c)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);

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
    //drawAxes(10.0f);


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
        if (state == GLUT_DOWN) {
            g_leftDown    = true;
            g_clickStartX = x; g_clickStartY = y;
            g_maybeClick  = true;
            g_dragged     = false;
        } else {                      // release
            g_leftDown = false;
            // A click = press + release with almost no movement. If the cursor
            // traveled past the threshold it was a camera orbit, not a pick.
            if (g_maybeClick && !g_dragged) {
                int mods = glutGetModifiers();
                handleVertexPick(x, y, (mods & GLUT_ACTIVE_SHIFT) != 0);
            }
            g_maybeClick = false;
            g_dragged    = false;
        }
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
        // Only treat this as a drag (camera orbit) once the cursor has moved
        // past a small threshold — otherwise a stationary click would jitter
        // the view instead of selecting a vertex.
        if (!g_dragged) {
            int mvx = x - g_clickStartX; if (mvx < 0) mvx = -mvx;
            int mvy = y - g_clickStartY; if (mvy < 0) mvy = -mvy;
            if (mvx > g_clickThreshold || mvy > g_clickThreshold)
                g_dragged = true;
        }
        if (g_dragged) {
            g_angleY += dx * 0.5f;
            g_angleX += dy * 0.5f;
            g_angleX = fmaxf(-90.0f, fminf(90.0f, g_angleX));
        }
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
// Keyboard handler: switch view (Klein / Poincaré / round-trip), help, quit
//-----------------------------------------------------------------------------
void keyboard(unsigned char key, int /*x*/, int /*y*/) {
    switch (key) {
        case '1': g_view = 0; std::cout << "View: Original (Klein)\n";      break;
        case '2': g_view = 1; std::cout << "View: Hyperbolic (Poincare)\n"; break;
        case '3': g_view = 2; std::cout << "View: Round-trip Back\n";       break;
        case 'h': case 'H': g_showHelp = !g_showHelp; break;
        case 'v': case 'V': g_showVerts = !g_showVerts; break;   // toggle vertex dots
        case 'c': case 'C':                           // clear the selection
            if (!g_selectedVerts.empty())
                std::cout << "Cleared " << g_selectedVerts.size() << " selected vertices\n";
            g_selectedVerts.clear();
            break;
        case 'q': case 'Q': case 27: exit(0); break;   // 27 = ESC
    }
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
        case 4: g_view = 0; break;                   // view: Original (Klein)
        case 5: g_view = 1; break;                   // view: Hyperbolic (Poincare)
        case 6: g_view = 2; break;                   // view: Round-trip Back
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
    glutAddMenuEntry("View: Original (Klein)",      4);
    glutAddMenuEntry("View: Hyperbolic (Poincare)", 5);
    glutAddMenuEntry("View: Round-trip Back",       6);
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

