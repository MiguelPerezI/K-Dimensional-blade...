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
#include "Cube.hpp"            // truncated cuboctahedron source (replaces the STL file)

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

// Spherical-inversion spawning (doc/inversion-center-radius.md) — defined later.
Vector3D computeCircumcenter3D(const Vector3D& A, const Vector3D& B, const Vector3D& C);
bool     computeInversionFromSelected(const std::vector<Vector3D>& verts,
                                      const Vector3D& meshCentroid,
                                      Vector3D& outCenter, double& outRadius);
void     ensureOutwardNormals(FacetBox& box);
void     spawnInversion();
void     undoSpawn();
void     redoSpawn();
void     globalInvert();       // SHIFT+ENTER: spawn inverted copies of the whole scene (keep originals)

// Unified vertex index (loaded mesh + every spawned child). A VertRef addresses
// ONE mesh: meshId = -1 -> loaded mesh (resolve via currentActiveBox()); >=0 ->
// g_spawned[meshId].mesh. Defined with the selection globals below; forward-declared
// here so meshOf() can be declared before its first use.
struct VertRef;
const FacetBox& meshOf(const VertRef& r);
void            appendUniqueVerts(const FacetBox& box, int meshId);
void            rebuildSpawnedVerts();
void     rebuildMesh();        // (re)build the 3 view meshes + vertex index from the truncated cuboctahedron
void     scrubInset(double delta);  // scrub hollow inset (clamped), rebuild if hollow
void     toggleHollowShape();       // solid <-> hollow frame, then rebuild
void     writeShapeSTL();           // export current shape as ASCII STL

Vector3D d_ui(1.0, 0.0, 0.0);
Vector3D e_ui(0.0, 1.0, 0.0);
Vector3D f_ui(0.0, 0.0, 1.0);
Facet f_1(d_ui, e_ui, f_ui);

// Loaded mesh, centered+scaled into the open unit ball (Klein disk model).
FacetBox g_boxKlein;
// Same mesh after applyHyperboloid() — projected to the Poincaré disk.
FacetBox g_boxHyper;
// Round-trip: inverse map applied to g_boxHyper, back to the Klein disk.
FacetBox g_boxBack;
// Active view: 0 = Original (Klein), 1 = Hyperbolic (Poincaré), 2 = Round-trip.
int g_view = 0;

// Truncated cuboctahedron (great rhombicuboctahedron, vertex config 4.6.8) is the
// mesh source here — built on demand from a basic Cube via
// getFacetsTruncatedCuboctahedron(hollow, inset), replacing the STL file used in
// main9.3.2. Rebuilt whenever hollow/inset changes (rebuildMesh()), which also
// re-derives the three views and the selectable-vertex index.
Cube   g_cube(1.0, Vector3D(0, 0, 0));   // basic, non-subdivided; 8 corners populated
bool   g_hollow = false;                 // false = solid 92-tri mesh; true = hollow 288-tri frame
double g_inset  = 0.5;                   // hollow border inset ratio, scrubbed with , / . (0.05..0.95)

//////////////////////////////////////
//                                  //
//  BLENDER-STYLE VERTEX SELECTION  //
//                                  //
//////////////////////////////////////
// A unique mesh vertex is identified by ONE representative facet/corner; the
// same logical vertex is shared by several facets, so we keep a single ref to
// store each mesh point once. Selection is by vertex identity, therefore it
// survives switching between the Klein / Poincaré / round-trip views.
struct VertRef { int meshId; int facetIdx; int cornerIdx; };  // meshId: -1=loaded, >=0=g_spawned[idx]
std::vector<VertRef> g_uniqueVerts;    // one representative per unique vertex (loaded first, then spawned)
std::set<int>        g_selectedVerts;  // indices into g_uniqueVerts currently highlighted
size_t g_loadedVertCount = 0;          // loaded-mesh unique verts occupy [0, g_loadedVertCount); spawned after
const double g_pickRadiusPx = 10.0;    // screen-space click tolerance (pixels)
bool         g_showVerts    = false;   // [V] toggle the faint all-vertex dot overlay

// Click-vs-drag bookkeeping so a left click picks a vertex while a left drag
// still orbits the camera (Blender treats a near-stationary press as a click).
int  g_clickStartX = 0, g_clickStartY = 0;
bool g_maybeClick   = false;   // left button went down, not yet a drag
bool g_dragged      = false;   // movement exceeded the threshold → orbit, not pick
const int g_clickThreshold = 5;

//////////////////////////////////////
//                                  //
//   SPHERICAL-INVERSION SPAWNING   //
//                                  //
//////////////////////////////////////
// From a set of selected mesh vertices that define a face (≥3) or an edge
// (exactly 2), compute the inversion sphere per doc/inversion-center-radius.md
// (center = circumcenter + √(R²−ρ²)·n̂, radius R = φ·ρ) and spawn an inverted
// copy of the active mesh via FacetBox::sigma(). ENTER spawns; Ctrl+Z / Ctrl+Y
// step the spawn history back / forward (LIFO undo, redo stack).
const double g_phi          = (1.0 + std::sqrt(5.0)) / 2.0; // golden ratio φ (doc §1)
const double g_minInvRadius = 0.05;                         // doc §2 clamps
const double g_maxInvRadius = 25.0;
struct SpawnRecord { FacetBox mesh; Vector3D center; double radius; int view; int batch; };
std::vector<SpawnRecord> g_spawned;   // visible inverted children (back = newest)
std::vector<SpawnRecord> g_undone;    // redo stack (undone spawns, restore with Ctrl+Y)
bool g_showInvSpheres = false;        // [B] toggle inversion-sphere wireframe overlay

// Batching: one Enter or one Shift+Enter can push SEVERAL children at once (a mixed
// selection spawns a copy of every touched mesh; Shift+Enter spawns a copy of the
// whole scene). All meshes pushed by a single action share a batch id, so a single
// Ctrl+Z removes the whole action together and Ctrl+Y restores it. Monotonic; a new
// action takes ++g_spawnBatch.
int g_spawnBatch = 0;


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


//-----------------------------------------------------------------------------
// Rebuild the three view meshes + the unique-vertex selection index from the
// truncated cuboctahedron at the current hollow/inset. Called once at startup
// and whenever the user scrubs inset / toggles hollow: the geometry (and hence
// the selectable vertex set) changes, so the selection and any spawned
// inversions are cleared. centerAndScale fits the mesh in the open unit ball so
// toHyperboloid() (which throws at ‖x‖ ≥ 1) is safe.
//-----------------------------------------------------------------------------
void rebuildMesh() {
    FacetBox raw = g_cube.getFacetsTruncatedCuboctahedron(g_hollow, g_inset);
    g_boxKlein = centerAndScale(raw, 0.9);          // fit inside the open unit ball
    g_boxHyper = g_boxKlein.hyperboloid();          // Klein → Poincaré
    g_boxBack  = inverseHyperboloid(g_boxHyper);    // Poincaré → Klein (round-trip)
    buildUniqueVertexList();                         // reindex selectable vertices
    g_selectedVerts.clear();                         // prior selection no longer valid
    g_spawned.clear();                               // spawned children referenced the
    g_undone.clear();                                //   old geometry — drop them too
    std::cout << "Rebuilt truncated cuboctahedron: "
              << (g_hollow ? "hollow" : "solid") << ", inset=" << g_inset
              << " -> " << g_boxKlein.size() << " tris, "
              << g_uniqueVerts.size() << " unique verts (view " << g_view << ")\n";
}

// Scrub the hollow border inset by `delta` (clamped to [0.05, 0.95]). A solid
// mesh ignores inset, so we only rebuild (and reindex vertices) when hollow.
void scrubInset(double delta) {
    g_inset += delta;
    if (g_inset > 0.95) g_inset = 0.95;
    if (g_inset < 0.05) g_inset = 0.05;
    std::cout << "inset = " << g_inset
              << (g_hollow ? "  -> rebuilding" : "  (no effect unless hollow: 'o')") << "\n";
    if (g_hollow) rebuildMesh();
}

// Toggle solid (92-tri) <-> hollow (288-tri frame) and rebuild the mesh.
void toggleHollowShape() {
    g_hollow = !g_hollow;
    std::cout << "hollow " << (g_hollow ? "ON (frame)" : "OFF (solid)") << "\n";
    rebuildMesh();
}

// Export the current shape (at the active hollow/inset) as an ASCII STL.
void writeShapeSTL() {
    const std::string path = "/home/mike666/Downloads/truncated_cuboctahedron.stl";
    g_cube.writeSTL_s_truncated_cuboctahedron(path, "TruncatedCuboctahedron",
                                             "full", 0, 9, 2, g_hollow, g_inset);
    std::cout << "wrote STL: " << path << "  (hollow=" << (g_hollow ? "ON" : "off")
              << ", inset=" << g_inset << ")\n";
}


void Setup() {

	if (ciclo == 0) {

        cout << "\n———————————————————————————————————————————————————————————————————————\n";
        cout <<   "|- Truncated Cuboctahedron (4.6.8)——————————————————————————————————————————\n";
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
        cout << "     \"-_\\_______/;;'        ~[Truncated Cuboctahedron]\n\n";

        // 1) Build the truncated cuboctahedron from the cube and derive the three
        //    view meshes + the selectable-vertex index. rebuildMesh() is also
        //    called whenever the user scrubs inset / toggles hollow.
        try {
            rebuildMesh();
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

        // Spawned spherically-inverted meshes (ENTER to spawn; Ctrl+Z/Y history).
        // Opaque (alpha 1.0) with the SAME per-facet rainbow as the active mesh, so
        // a child reads as a solid copy of the original — same color & texture —
        // rather than a translucent single-hue blob.
        for (size_t s = 0; s < g_spawned.size(); ++s) {
            const FacetBox& m = g_spawned[s].mesh;
            size_t n = m.size();
            for (size_t i = 0; i < n; ++i) {
                float hue = float(i) / float(n) * 360.0f;   // same rainbow mapping as the active mesh
                Color c = hsv2rgb(hue, 0.8f, 1.0f);
                int R = int(c.r * 255), G = int(c.g * 255), B = int(c.b * 255);
                drawFacet(m[i], R, G, B, 1.0f);             // opaque — same texture as the original
            }
        }

        // Optional wireframe of every inversion sphere (toggle with [B]) so the
        // doc's construction (center + radius the sphere passes the verts on) is
        // visible. Drawn unlit & translucent so it never hides the mesh.
        if (g_showInvSpheres) {
            glPushAttrib(GL_LIGHTING_BIT | GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
            glDisable(GL_LIGHTING);
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            for (size_t s = 0; s < g_spawned.size(); ++s) {
                const Vector3D& ctr = g_spawned[s].center;
                glPushMatrix();
                glTranslatef(ctr.x(), ctr.y(), ctr.z());
                glColor4ub(60, 120, 255, 110);
                glutWireSphere((double)g_spawned[s].radius, 24, 16);
                glPopMatrix();
            }
            glPopAttrib();
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

// Resolve a VertRef to the FacetBox that holds its vertex. meshId = -1 -> the
// loaded mesh in the CURRENT view (the Poincaré / round-trip meshes are pointwise
// images, so the facet/corner layout — and thus the index — is valid in all three
// views); meshId >= 0 -> that spawned child (stored in its spawn-view coords).
const FacetBox& meshOf(const VertRef& r) {
    return (r.meshId < 0) ? currentActiveBox() : g_spawned[(size_t)r.meshId].mesh;
}

// Append one representative per unique spatial vertex of `box` to g_uniqueVerts,
// tagged with `meshId`. The spatial-hash grid is LOCAL to this call, so dedup is
// WITHIN this mesh only — a spawned child's glued vert (coincident with the
// parent's) stays a separate, selectable point on the child. The dup-check
// resolves through `box` itself (not currentActiveBox), so the loaded mesh dedups
// in Klein coords regardless of g_view. O(n) per mesh.
void appendUniqueVerts(const FacetBox& box, int meshId) {
    if (box.size() == 0) return;

    const double tol = 1e-6;       // merge tolerance (shared verts are bit-identical)
    const double cs  = 1e-3;       // grid cell: >> tol, << smallest vertex spacing
    const long long OFF = (1LL << 20);
    auto pack = [&](long long rx, long long ry, long long rz) -> long long {
        return (rx + OFF) | ((ry + OFF) << 21) | ((rz + OFF) << 42);
    };
    std::unordered_map<long long, std::vector<int>> grid;   // local to this mesh
    grid.reserve(1 << 15);

    for (size_t fi = 0; fi < box.size(); ++fi) {
        const Facet& f = box[fi];
        for (int ci = 0; ci < 3; ++ci) {
            const Vector3D& p = f[ci];
            long long rx = (long long)floor(p.x() / cs);
            long long ry = (long long)floor(p.y() / cs);
            long long rz = (long long)floor(p.z() / cs);

            // A duplicate sits in this cell or one of its 26 neighbours.
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
                g_uniqueVerts.push_back({meshId, (int)fi, ci});
                grid[pack(rx, ry, rz)].push_back(newIdx);
            }
        }
    }
}

// Rebuild the spawned-vertex tail of g_uniqueVerts from the current g_spawned,
// keeping the loaded-mesh head [0, g_loadedVertCount) untouched. Called after every
// spawn / undo / redo. Spawns are LIFO on the back of g_spawned, so spawned-vert
// indices are prefix-stable and surviving selections stay valid; any selected
// index that no longer exists (a removed spawn) is pruned.
void rebuildSpawnedVerts() {
    g_uniqueVerts.resize(g_loadedVertCount);
    for (size_t s = 0; s < g_spawned.size(); ++s)
        appendUniqueVerts(g_spawned[s].mesh, (int)s);
    for (auto it = g_selectedVerts.begin(); it != g_selectedVerts.end(); ) {
        if (*it >= (int)g_uniqueVerts.size()) it = g_selectedVerts.erase(it);
        else                                  ++it;
    }
}

// Index the loaded mesh's unique vertices from the Klein mesh, then append every
// spawned child's vertices so they are clickable too. After rebuildMesh g_spawned
// is empty, so rebuildSpawnedVerts() is a no-op there.
void buildUniqueVertexList() {
    g_uniqueVerts.clear();
    g_selectedVerts.clear();
    appendUniqueVerts(g_boxKlein, -1);
    g_loadedVertCount = g_uniqueVerts.size();
    rebuildSpawnedVerts();
    std::cout << "Unique vertices indexed: " << g_uniqueVerts.size()
              << " (loaded " << g_loadedVertCount << ", spawned "
              << (g_uniqueVerts.size() - g_loadedVertCount) << ")\n";
}

// Draw the selected vertices as bright-orange dots with a thin black outline,
// always on top (depth test off) so they stay visible through the mesh — like
// Blender's selected verts. The optional faint-gray overlay of *every* vertex
// (Blender's full edit-mode dot field) is gated behind [V].
void drawVertexMarkers() {
    if (g_uniqueVerts.empty()) return;

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
            const Vector3D& p = meshOf(r)[r.facetIdx][r.cornerIdx];
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
            const Vector3D& p = meshOf(r)[r.facetIdx][r.cornerIdx];
            glVertex3f(p.x(), p.y(), p.z());
        }
        glEnd();

        glPointSize(9.0f);
        glColor3ub(255, 140, 0);     // bright orange
        glBegin(GL_POINTS);
        for (int idx : g_selectedVerts) {
            const VertRef& r = g_uniqueVerts[idx];
            const Vector3D& p = meshOf(r)[r.facetIdx][r.cornerIdx];
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
        const Vector3D& p = meshOf(r)[r.facetIdx][r.cornerIdx];
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

//=============================================================================
//  SPHERICAL-INVERSION SPAWNING  (doc/inversion-center-radius.md)
//=============================================================================
// 3-D circumcenter of triangle A,B,C — the center of its circumscribed circle,
// lying in the triangle's plane. Any 3 of a face's vertices give the same point
// when the face is concyclic (doc §3). Throws if the triangle is degenerate.
Vector3D computeCircumcenter3D(const Vector3D& A, const Vector3D& B, const Vector3D& C) {
    Vector3D u = B - A;                 // edge AB
    Vector3D v = C - A;                 // edge AC
    double uu = u * u;                  // |u|²  (dot — Vector3D.cpp:107)
    double vv = v * v;                  // |v|²
    double uv = u * v;                  // u·v
    double det = uu * vv - uv * uv;     // = |u × v|² (squared, = 4·area²)
    if (det < 1e-18)
        throw std::runtime_error("computeCircumcenter3D: degenerate (collinear) triangle");
    // Solve  [uu uv ; uv vv]·[α;β] = [uu/2 ; vv/2]  →  O = A + α·u + β·v.
    double alpha = vv * (uu - uv) / (2.0 * det);
    double beta  = uu * (vv - uv) / (2.0 * det);
    return A + alpha * u + beta * v;    // double*Vector3D — Vector3D.cpp:88
}

// Turn the selected vertex positions into an inversion (center, radius) using
// the doc's §1 construction. ≥3 non-collinear verts define a *face*; exactly 2
// (or an all-collinear set) define an *edge* via the doc-spirited extension.
// Returns false only if the selection has fewer than 2 distinct points.
bool computeInversionFromSelected(const std::vector<Vector3D>& verts,
                                  const Vector3D& meshCentroid,
                                  Vector3D& outCenter, double& outRadius) {
    if (verts.size() < 2) return false;

    // --- find any non-collinear triple for the face construction ---
    int i0 = -1, i1 = -1, i2 = -1;
    for (size_t i = 0; i < verts.size() && i0 < 0; ++i)
        for (size_t j = i + 1; j < verts.size() && i0 < 0; ++j)
            for (size_t k = j + 1; k < verts.size() && i0 < 0; ++k)
                if (!areColinear(verts[i], verts[j], verts[k]))   // Vector3D.cpp:165
                    { i0 = (int)i; i1 = (int)j; i2 = (int)k; }

    Vector3D base;     // point on the face/edge we lift the center from
    Vector3D normal;   // outward unit direction we lift along
    double rho;        // circumradius ρ (face) or half-edge (edge)

    if (i0 >= 0) {
        // ----- face construction (doc §1) -----
        const Vector3D& A = verts[i0];
        const Vector3D& B = verts[i1];
        const Vector3D& C = verts[i2];
        base   = computeCircumcenter3D(A, B, C);
        rho    = abs(base - A);                       // face circumradius
        normal = unit((B - A) % (C - A));             // face normal
        if (normal * (base - meshCentroid) < 0.0)     // flip to point outward,
            normal = -normal;                         // away from the body (doc §8)
    } else {
        // ----- edge construction (exactly-2 / all-collinear fallback) -----
        // Sphere through the two farthest selected points, centered on the ray
        // from the mesh centroid through their midpoint; R = φ·(half-edge). The
        // child is glued along the shared edge.
        int pa = 0, pb = 1; double best = -1.0;
        for (size_t i = 0; i < verts.size(); ++i)
            for (size_t j = i + 1; j < verts.size(); ++j) {
                double d = abs(verts[i] - verts[j]);
                if (d > best) { best = d; pa = (int)i; pb = (int)j; }
            }
        if (best < 1e-12) return false;               // all points coincident
        const Vector3D& p = verts[pa];
        const Vector3D& q = verts[pb];
        base   = (p + q) * 0.5;                       // edge midpoint
        rho    = abs(p - q) * 0.5;                    // half-edge = ρ
        normal = unit(base - meshCentroid);           // outward from body
        if (abs(normal) < 1e-9) normal = Vector3D(0, 0, 1);  // guard midpoint≈centroid
    }

    // Both paths finish identically: R = φ·ρ (clamped per doc §2), center lifted
    // along the outward normal by √(R²−ρ²); the sphere passes through every
    // selected vertex by construction (Pythagoras on the normal lift).
    double R = fmax(g_minInvRadius, fmin(g_maxInvRadius, g_phi * rho));
    double offset = sqrt(fmax(0.0, R * R - rho * rho));
    outCenter = base + offset * normal;
    outRadius = abs(outCenter - verts[0]);            // == R
    return true;
}

// Spherical inversion reverses orientation, so a child whose source was wound
// CCW-outward comes back wound inward → it would render inside-out and the next
// construction's outward normal would be wrong. Per face, flip the winding of
// any facet whose normal points toward the cell centroid (doc §8).
void ensureOutwardNormals(FacetBox& box) {
    if (box.size() == 0) return;
    Vector3D c = box.center();                        // FacetBox.hpp:222
    for (size_t i = 0; i < box.size(); ++i) {
        Vector3D A = box[i][0], B = box[i][1], C = box[i][2];
        Vector3D fc = (A + B + C) / 3.0;              // face centroid
        Vector3D n  = unit((B - A) % (C - A));        // current-winding normal
        if (n * (fc - c) < 0.0)                       // points inward → reverse
            box.replace(i, A, C, B);                  // swap B,C — FacetBox.hpp:191
    }
}

// ENTER: spawn a spherically-inverted copy of EVERY mesh that has a selected
// vertex, through the ONE sphere defined by the full selection. Verts may come
// from any mesh (a mixed selection is NOT refused — they just define the sphere);
// each touched mesh gets a new child glued to the sphere, and the originals stay.
// A single-mesh selection spawns one child (as before). All children of one Enter
// share a batch, so a single Ctrl+Z removes them together. New children's verts
// are appended to the clickable index (rebuildSpawnedVerts) so they can be selected
// for the next spawn (honeycomb: a child of a child). Coherent when all selected
// verts share one view; children spawned in a different view invert best-effort.
void spawnInversion() {
    if (g_selectedVerts.size() < 2) {
        std::cout << "Spawn: select >=2 vertices defining a face/edge first.\n";
        return;
    }

    // Resolve selected vertex positions (each in its own mesh's space) and collect
    // the distinct meshes touched by the selection.
    std::vector<Vector3D> verts;
    verts.reserve(g_selectedVerts.size());
    std::set<int> touched;                              // meshId of each touched mesh (-1 = loaded)
    for (int idx : g_selectedVerts) {
        const VertRef& r = g_uniqueVerts.at(idx);
        const FacetBox& m = meshOf(r);
        verts.push_back(m[r.facetIdx][r.cornerIdx]);
        touched.insert(r.meshId);
    }

    // ONE sphere from the full selection (coherent when all verts share one view).
    Vector3D center; double radius;
    if (!computeInversionFromSelected(verts, currentActiveBox().center(), center, radius)) {
        std::cout << "Spawn: selected vertices are degenerate (coincident).\n";
        return;
    }

    // Invert each touched mesh through that sphere, each as a new child (originals
    // untouched). Reserve so existing child refs in g_spawned stay valid while we push.
    int batch = ++g_spawnBatch;
    size_t origCount = g_spawned.size();
    g_spawned.reserve(origCount + touched.size());

    bool crossView = false;
    int spawnedCount = 0;
    for (int meshId : touched) {                        // -1 (loaded) sorts first
        if (meshId < 0) {
            const FacetBox& src = currentActiveBox();
            if (src.size() == 0) continue;
            FacetBox child = src.sigma(center, radius);  // FacetBox.hpp:321
            ensureOutwardNormals(child);
            g_spawned.push_back({child, center, radius, g_view, batch});
            ++spawnedCount;
        } else {
            size_t s = (size_t)meshId;
            if (s >= origCount) continue;               // safety (shouldn't happen)
            const FacetBox& src = g_spawned[s].mesh;    // stable: reserved
            if (src.size() == 0) continue;
            FacetBox child = src.sigma(center, radius);
            ensureOutwardNormals(child);
            int srcView = g_spawned[s].view;
            g_spawned.push_back({child, center, radius, srcView, batch});
            if (srcView != g_view) crossView = true;
            ++spawnedCount;
        }
    }

    if (spawnedCount == 0) {
        std::cout << "Spawn: no non-empty source mesh among the selected verts.\n";
        return;
    }
    g_undone.clear();                                   // a new action clears redo
    rebuildSpawnedVerts();                              // make the new children's verts clickable

    std::cout << "Spawned " << spawnedCount << " inverted mesh"
              << (spawnedCount == 1 ? "" : "es") << " -> " << g_spawned.size()
              << " total | center " << center << "radius " << radius
              << " | phi*rho = " << (radius / g_phi)
              << " | from " << verts.size() << " verts on " << touched.size()
              << " mesh" << (touched.size() == 1 ? "" : "es")
              << " (view " << g_view << ")\n";
    if (crossView)
        std::cout << "  note: some children were spawned in a different view; their"
                     " copies invert best-effort (spawn in one view for coherence).\n";
    glutPostRedisplay();
}

// Ctrl+Z: revert the most recent spawn BATCH — every mesh sharing the top record's
// batch id (one Enter = 1 mesh; one Shift+Enter = the whole spawned scene). They are
// contiguous at the back, so we pop until the batch changes; the batch goes onto the
// redo stack as a unit so Ctrl+Y restores it. g_spawned is LIFO, so spawned-vert
// indices stay prefix-stable and only the removed batch's selection is pruned.
void undoSpawn() {
    if (g_spawned.empty()) {
        std::cout << "Undo: nothing to revert (no spawns).\n";
        return;
    }
    int batch = g_spawned.back().batch;
    size_t before = g_spawned.size();
    while (!g_spawned.empty() && g_spawned.back().batch == batch) {
        g_undone.push_back(g_spawned.back());
        g_spawned.pop_back();
    }
    rebuildSpawnedVerts();                              // drop the undone batch's verts, prune selection
    std::cout << "Undid spawn batch of " << (before - g_spawned.size())
              << " -> " << g_spawned.size() << " visible, " << g_undone.size()
              << " on redo stack.\n";
    glutPostRedisplay();
}

// Ctrl+Y: re-apply the most recently undone spawn BATCH (all meshes sharing the top
// redo record's batch id, restored in their original order).
void redoSpawn() {
    if (g_undone.empty()) {
        std::cout << "Redo: nothing to re-apply.\n";
        return;
    }
    int batch = g_undone.back().batch;
    while (!g_undone.empty() && g_undone.back().batch == batch) {
        g_spawned.push_back(g_undone.back());
        g_undone.pop_back();
    }
    rebuildSpawnedVerts();                              // re-append the redone batch's verts
    std::cout << "Redid spawn batch -> " << g_spawned.size()
              << " visible, " << g_undone.size() << " on redo stack.\n";
    glutPostRedisplay();
}

// SHIFT+ENTER: spawn an inverted copy of the WHOLE scene through the sphere
// defined by the selected vertices — the active mesh (current view) AND every
// spawned child — each added as a NEW child. The originals are kept, so after the
// press you have the original scene PLUS a full inverted copy of it. Selected verts
// lie on the sphere, so each inverted copy is glued to the verts it shares with the
// sphere. A mixed selection is fine (the verts just define the sphere). All copies
// share one batch, so a single Ctrl+Z removes the whole spawned scene. Coherent
// when the children are in the current view; children spawned in a different view
// invert best-effort (spawn in one view, Klein recommended, for coherence).
void globalInvert() {
    if (g_selectedVerts.size() < 2) {
        std::cout << "Global invert: select >=2 verts defining the sphere first.\n";
        return;
    }
    const FacetBox& active = currentActiveBox();
    if (active.size() == 0) {
        std::cout << "Global invert: no mesh loaded.\n";
        return;
    }

    // ONE sphere from the full selection (verts from any mesh allowed).
    std::vector<Vector3D> verts;
    verts.reserve(g_selectedVerts.size());
    for (int idx : g_selectedVerts) {
        const VertRef& r = g_uniqueVerts.at(idx);
        const FacetBox& m = meshOf(r);
        verts.push_back(m[r.facetIdx][r.cornerIdx]);
    }
    Vector3D center; double radius;
    if (!computeInversionFromSelected(verts, active.center(), center, radius)) {
        std::cout << "Global invert: selected verts are degenerate (coincident).\n";
        return;
    }

    // Invert the active mesh + every existing child through the same sphere, each as
    // a new child (originals untouched). Reserve so existing child refs stay valid.
    int batch = ++g_spawnBatch;
    size_t origCount = g_spawned.size();
    g_spawned.reserve(origCount + 1 + origCount);

    {   // active mesh (current view)
        FacetBox child = active.sigma(center, radius);  // FacetBox.hpp:321
        ensureOutwardNormals(child);
        g_spawned.push_back({child, center, radius, g_view, batch});
    }
    bool crossView = false;
    for (size_t s = 0; s < origCount; ++s) {            // stable: reserved, indices < origCount
        const FacetBox& src = g_spawned[s].mesh;
        if (src.size() == 0) continue;
        FacetBox child = src.sigma(center, radius);
        ensureOutwardNormals(child);
        int srcView = g_spawned[s].view;
        g_spawned.push_back({child, center, radius, srcView, batch});
        if (srcView != g_view) crossView = true;
    }

    g_undone.clear();                                   // a new action clears redo
    rebuildSpawnedVerts();                               // make the new copies' verts clickable

    std::cout << "Global invert: spawned " << (g_spawned.size() - origCount)
              << " inverted meshes (active + children) -> " << g_spawned.size()
              << " total | center " << center << "radius " << radius
              << " | " << g_selectedVerts.size() << " verts targeted (view " << g_view
              << "). Ctrl+Z reverts the batch.\n";
    if (crossView)
        std::cout << "  note: some children were spawned in a different view; their"
                     " copies invert best-effort (spawn in one view for coherence).\n";
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
        "Enter: Spawn from selected meshes",
        "Shift+Enter: Spawn inverted whole scene",
        "Ctrl+Z: Undo batch | Ctrl+Y: Redo batch",
        "Spawned-mesh verts are clickable too (honeycomb)",
        "[B]: Toggle inversion spheres",
        "[,]/[.]: Decrease/Increase inset",
        "[O]: Toggle hollow (solid/frame)",
        "[W]: Write STL of current shape",
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

    // Live spawn history count
    y -= 0.05f;
    char spawnBuf[64];
    snprintf(spawnBuf, sizeof(spawnBuf), "Spawned: %lu (redo: %lu)",
             (unsigned long)g_spawned.size(), (unsigned long)g_undone.size());
    glRasterPos2f(0.02f, y);
    for(const char* c=spawnBuf; *c; ++c)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);

    // Current shape (truncated cuboctahedron hollow/inset)
    y -= 0.05f;
    char shapeBuf[96];
    snprintf(shapeBuf, sizeof(shapeBuf), "Shape: %s  inset=%.2f  (%s tris)",
             g_hollow ? "hollow (frame)" : "solid", g_inset,
             g_hollow ? "288" : "92");
    glRasterPos2f(0.02f, y);
    for(const char* c=shapeBuf; *c; ++c)
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
        case 13: case 10: {                            // ENTER (GLUT '\r' / '\n')
            int mods = glutGetModifiers();
            if (mods & GLUT_ACTIVE_SHIFT) globalInvert();   // SHIFT+ENTER: invert everything
            else                          spawnInversion(); // ENTER: spawn one child
            break;
        }
        case 26:          undoSpawn();    break;      // Ctrl+Z  (control-char code)
        case 25:          redoSpawn();    break;      // Ctrl+Y
        case 'b': case 'B': g_showInvSpheres = !g_showInvSpheres; break; // toggle spheres
        case '.': case '>': scrubInset(+0.05); break;   // increase hollow inset
        case ',': case '<': scrubInset(-0.05); break;   // decrease hollow inset
        case 'o': case 'O': toggleHollowShape(); break; // solid <-> hollow frame
        case 'w': case 'W': writeShapeSTL();   break;   // write STL of current shape
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
        case 4: g_view = 0; break;   // view: Original (Klein)
        case 5: g_view = 1; break;   // view: Hyperbolic (Poincare)
        case 6: g_view = 2; break;   // view: Round-trip Back
        case 7: spawnInversion(); break;             // spawn inverted mesh (Enter)
        case 8: undoSpawn();    break;               // undo spawn (Ctrl+Z)
        case 9: redoSpawn();    break;               // redo spawn (Ctrl+Y)
        case 10: g_showInvSpheres = !g_showInvSpheres; break; // toggle spheres (B)
        case 11: scrubInset(+0.05); break;        // increase inset (.)
        case 12: scrubInset(-0.05); break;        // decrease inset (,)
        case 13: toggleHollowShape(); break;      // toggle hollow (O)
        case 14: writeShapeSTL(); break;          // write STL (W)
        case 15: globalInvert(); break;           // spawn inverted copy of whole scene (Shift+Enter)
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
    glutAddMenuEntry("Spawn Inverted Mesh (Enter)", 7);
    glutAddMenuEntry("Undo Spawn (Ctrl+Z)",         8);
    glutAddMenuEntry("Redo Spawn (Ctrl+Y)",         9);
    glutAddMenuEntry("Toggle Inversion Spheres (B)",10);
    glutAddMenuEntry("Increase Inset (.)",           11);
    glutAddMenuEntry("Decrease Inset (,)",           12);
    glutAddMenuEntry("Toggle Hollow (O)",            13);
    glutAddMenuEntry("Write STL (W)",                14);
    glutAddMenuEntry("Global Invert: Spawn Copy (Shift+Enter)", 15);
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

