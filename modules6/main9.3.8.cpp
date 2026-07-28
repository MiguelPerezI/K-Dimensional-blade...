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
#include <charconv>               // std::to_chars — shortest round-trip decimal (greyBucket)
#include <filesystem>             // std::filesystem::directory_iterator — scan g_stlDir for the Load STL submenu
#include <algorithm>              // std::sort — keep the .stl list ordered across rescans
#include <cctype>                 // tolower — case-insensitive .stl extension match
#include <csignal>                // std::signal/SIGPIPE/SIG_IGN — survive a dead ffmpeg pipe
#include <ctime>                  // time/localtime/strftime — timestamped render output filename
#define GL_GLEXT_PROTOTYPES 1     // expose FBO prototypes (glGenFramebuffers, ...) in glext.h
#include <GL/glut.h>
#include <GL/glext.h>             // offscreen framebuffer object -> render-to-video at 1280x720
#include "Vector3D.cpp"
#include "Vector4D.cpp"
#include "Quaternion.cpp"
#include "Facet.cpp"
#include "FacetBox.hpp"
#include "Cube.hpp"            // polyhedron factory (cube / truncated cube / octahedron / cuboctahedron)
#include "Dodecahedron.hpp"     // separate class; header-only (do NOT link Dodecahedron.cpp — it is stale)

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
void handleVertexFocus(int mouseX, int mouseY);   // right-click -> camera pivot

// Camera: current pivot (focused vertex, or origin) and the shared modelview
// used by BOTH display() and pickVertex() so picking stays accurate after a focus.
Vector3D focusTarget();
void     applyCameraModelview();

// Soft animated background gradient (very light pale blue <-> pale rose, a
// desaturated "jet" feel) + the timer that keeps it breathing when the user
// isn't interacting. Defined just before display().
void drawBackgroundGradient();
void backgroundTimer(int value);

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
void     remesh(FacetBox::SubdivisionMode mode);  // R/M: subdivide every facet (centroid3 / midpoint4)

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

// Externally-loaded ASCII STL mesh (the raw, un-centered source for Shape::STL).
// rebuildMesh() centers+scales it into the unit ball exactly like the procedural
// shapes, so arbitrary STL geometry drops in with no special-case handling.
FacetBox g_stlMesh;
std::string g_stlPath;                 // last loaded STL path (L key reloads it; menu marks it active)
std::vector<std::string> g_stlFiles;   // .stl filenames found in g_stlDir (drives the Load STL submenu)
std::string g_stlDir = "/home/mike666/Downloads";   // folder scanned for the Load STL submenu

// Render-to-video state. renderVideo() drives the keyframe playback deterministically
// (advancePlayback(1/fps) per frame) and pipes each rendered back buffer to ffmpeg as
// raw rgb24. g_rendering gates drawHUD()/drawVertexMarkers() so the captured frames
// contain only the mesh + breathing background (no on-screen UI), and switches the
// background phase to video time so it breathes smoothly in the output. Declared up
// here (before Draw()) because g_rendering is read at the vertex-marker gate in Draw().
static bool   g_rendering = false;       // true only during a render-to-video pass
static FILE*  g_ffPipe = nullptr;        // ffmpeg stdin (raw rgb24 frames)
static int    g_renderFps = 24;          // selected capture/output fps
static int    g_renderW = 0, g_renderH = 0;   // output dims (fixed 1280x720 16:9); set in renderVideo()
static std::vector<unsigned char> g_frameBuf; // g_renderW*g_renderH*3 scratch
static int    g_renderFrame = 0, g_renderTotal = 0;  // progress
static bool   g_renderAbort = false;     // set by the capture hook if ffmpeg's pipe breaks
static int    g_lastRenderFps = 0;       // marks the "* active" entry in the Render submenu
static bool   g_cliRender = false;       // -render <fps> given on the command line
static int    g_cliRenderFps = 0;

// Offscreen render target for render-to-video: a 1280x720 (16:9) framebuffer object the
// video is drawn into, independent of the 720x720 on-screen window (no resize / event
// pump). Created lazily by ensureRenderFBO() and reused across renders.
static GLuint g_rFbo = 0, g_rColorTex = 0, g_rDepthRb = 0;
static const int g_renderOutW = 1280, g_renderOutH = 720;   // fixed YouTube HD output

// Active view: 0 = Original (Klein), 1 = Hyperbolic (Poincaré), 2 = Round-trip.
int g_view = 0;
// Jet-palette phase, advanced by every remesh() press so each refinement pass
// recolors the whole scene with a rotated jet distribution (distinct bands each
// press, long cycle before repeating). Draw() samples jet(i/total + g_colorPhase).
float g_colorPhase = 0.0f;
// Per-facet coloring mode. 0 = jet (rotated by g_colorPhase). 1 = special 5-bucket
// grayscale: each facet is black→grey→white by greyBucket(facet). Toggled with G.
int g_colorMode = 0;

// Active polyhedron source for rebuildMesh(). Switched at runtime from the
// middle-click "Polyhedron" submenu (direct) or the S key / "Next Polyhedron"
// menu entry (cycle). Default preserves the original truncated-cuboctahedron
// behavior. g_hollow/g_inset affect every shape: hollow turns each flat face
// into an inset picture-frame border (Cube 48 tris, Dodecahedron 120 tris,
// Truncated Cube octagonal holes; TO/TC as before), and inset scrubs its
// thickness. Solid meshes are unchanged.
enum class Shape {
    Cube,
    Dodecahedron,
    TruncatedCube,
    TruncatedOctahedron,
    TruncatedCuboctahedron,
    STL                  // an externally-loaded ASCII STL mesh (g_stlMesh); ignores hollow/inset
};

const char* shapeName(Shape s) {
    switch (s) {
        case Shape::Cube:                   return "Cube";
        case Shape::Dodecahedron:           return "Dodecahedron";
        case Shape::TruncatedCube:          return "Truncated Cube";
        case Shape::TruncatedOctahedron:    return "Truncated Octahedron";
        case Shape::TruncatedCuboctahedron: return "Truncated Cuboctahedron";
        case Shape::STL:                     return "STL";
    }
    return "Unknown";
}

// Mesh source — built on demand from a basic Cube (or a Dodecahedron) inside
// rebuildMesh(), which also re-derives the three views and the selectable-vertex
// index. Rebuilt whenever the shape, hollow, or inset changes.
Cube   g_cube(1.0, Vector3D(0, 0, 0));   // basic, non-subdivided; 8 corners populated
bool   g_hollow = false;                 // false = solid mesh; true = hollow frame (truncated family)
double g_inset  = 0.5;                   // hollow border inset ratio, scrubbed with , / . (0.05..0.95)
double g_bgLightness = 0.28;             // background wash lightness (0=color, 1=white); set in main(), not in any menu/key
Shape  g_shape = Shape::TruncatedCuboctahedron;  // default keeps prior behavior

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
std::set<int>        g_selectedVerts;  // indices into g_uniqueVerts currently highlighted (left-click, orange)
int g_focusVert = -1;                  // index into g_uniqueVerts the camera orbits/zooms around
                                       // (right-click, blue). -1 = world origin. Cleared on rebuild/reset/menu.
size_t g_loadedVertCount = 0;          // loaded-mesh unique verts occupy [0, g_loadedVertCount); spawned after
const double g_pickRadiusPx = 10.0;    // screen-space click tolerance (pixels)
bool         g_showVerts    = false;   // [V] toggle the faint all-vertex dot overlay

// Click-vs-drag bookkeeping so a left click picks a vertex while a left drag
// still orbits the camera (Blender treats a near-stationary press as a click).
int  g_clickStartX = 0, g_clickStartY = 0;
bool g_maybeClick   = false;   // left button went down, not yet a drag
bool g_dragged      = false;   // movement exceeded the threshold → orbit, not pick
const int g_clickThreshold = 5;

// Same click-vs-drag bookkeeping for the right button, so a stationary right-click
// focuses a vertex (sets the camera pivot) while a right-drag still zooms.
int  g_rightClickStartX = 0, g_rightClickStartY = 0;
bool g_rightMaybeClick  = false;
bool g_rightDragged     = false;

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


// Classic MATLAB "jet" colormap: dark blue → blue → cyan → yellow → red → dark red.
// Piecewise-linear across six stops (same shape as matplotlib's jet). Input t may
// be any real; it is wrapped into [0,1) so adding g_colorPhase just rotates the band
// instead of running off the end. Returns RGB in [0,1].
Color jet(float t) {
    t = fmodf(t, 1.0f); if (t < 0.0f) t += 1.0f;     // wrap into [0,1)
    const float pos[6]   = {0.000f, 0.125f, 0.375f, 0.625f, 0.875f, 1.000f};
    const float col[6][3] = {
        {0.00f, 0.00f, 0.50f},  // 0.000 dark blue
        {0.00f, 0.00f, 1.00f},  // 0.125 blue
        {0.00f, 1.00f, 1.00f},  // 0.375 cyan
        {1.00f, 1.00f, 0.00f},  // 0.625 yellow
        {1.00f, 0.00f, 0.00f},  // 0.875 red
        {0.50f, 0.00f, 0.00f},  // 1.000 dark red
    };
    int k = 1;
    while (k < 5 && t > pos[k]) ++k;                 // find segment [pos[k-1], pos[k]]
    float f = (t - pos[k-1]) / (pos[k] - pos[k-1]);   // 0..1 within that segment
    return { col[k-1][0] + f*(col[k][0]-col[k-1][0]),
             col[k-1][1] + f*(col[k][1]-col[k-1][1]),
             col[k-1][2] + f*(col[k][2]-col[k-1][2]) };
}


// Special grayscale criterion (g_colorMode == 1). For a facet with vertices A,B,C:
//   3_sum   = A + B + C                       (vector sum of the three vertices)
//   scalar  = 3_sum.x + 3_sum.y + 3_sum.z      (a double, e.g. 12.123234)
//   take the LAST DIGIT of scalar, then mod 5  -> one of {0,1,2,3,4}
// The 5 buckets map to a black→white ramp: 0 = black, 4 = white, 1..3 = the grey
// gradients in between. "Last digit" means the rightmost digit of the number as
// you would write it (12.123234 -> '4'), so we render the value with std::to_chars'
// shortest round-trip decimal and read back from the end. That avoids picking any
// arbitrary fixed precision and ignores FP noise; the sign is ignored (the last
// digit of -12.123234 is also '4'). If to_chars ever emits exponent form, we read
// the last digit of the mantissa (before the 'e'), which is the same significant
// digit. Falls back to bucket 0 if no digit is found.
inline int greyBucket(const Facet& f) {
    Vector3D s3 = f[0] + f[1] + f[2];            // 3_sum = A + B + C
    double scalar = s3.x() + s3.y() + s3.z();    // e.g. 12.123234
    char buf[64];
    auto res = std::to_chars(buf, buf + sizeof(buf), scalar);   // shortest decimal
    const char* end = res.ptr;                   // one past the last char written
    // if exponent form ("1.23e-05"), stop scanning at the 'e' so we read the
    // mantissa's last digit, not the exponent's
    for (const char* p = buf; p < end; ++p)
        if (*p == 'e' || *p == 'E') { end = p; break; }
    int last = 0;
    for (const char* p = end - 1; p >= buf; --p)            // scan back for last digit
        if (*p >= '0' && *p <= '9') { last = *p - '0'; break; }
    return last % 5;                             // mod 5 -> {0,1,2,3,4}
}
// Map a 0..4 bucket to an even 0..255 grey byte: 0,64,128,191,255
// (black, 3 grey gradients, white).
inline int greyByte(int bucket) {
    return (int)std::lround(bucket * 255.0 / 4.0);
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
    FacetBox raw;
    switch (g_shape) {
        case Shape::Cube:
            raw = g_cube.getFacetsFrame(g_hollow, g_inset);  // solid 12 tris / hollow 48-tri frame
            break;
        case Shape::Dodecahedron: {
            // Separate class; build fresh each rebuild. FaceMode::Triangles = 36 tris solid.
            // hollow -> 120-tri 12-pentagon inset frame; inset scrubs the border thickness.
            Dodecahedron dodec(1.0, Vector3D(0, 0, 0), FaceMode::Triangles);
            raw = dodec.getFacetsFrame(g_hollow, g_inset);
            break;
        }
        case Shape::TruncatedCube:
            raw = g_cube.getFacetsOctagonal(g_hollow, g_inset);  // inset scrubs the octagonal hole
            break;
        case Shape::TruncatedOctahedron:
            raw = g_cube.getFacetsTruncatedOctahedron(0.5, g_hollow, g_inset); // s=0.5 regular TO
            break;
        case Shape::TruncatedCuboctahedron:
            raw = g_cube.getFacetsTruncatedCuboctahedron(g_hollow, g_inset);
            break;
        case Shape::STL:
            raw = g_stlMesh;                 // externally-loaded STL; hollow/inset do not apply
            break;
    }
    g_boxKlein = centerAndScale(raw, 0.9);          // fit inside the open unit ball
    g_boxHyper = g_boxKlein.hyperboloid();          // Klein → Poincaré
    g_boxBack  = inverseHyperboloid(g_boxHyper);    // Poincaré → Klein (round-trip)
    buildUniqueVertexList();                         // reindex selectable vertices
    g_selectedVerts.clear();                         // prior selection no longer valid
    g_spawned.clear();                               // spawned children referenced the
    g_undone.clear();                                //   old geometry — drop them too
    if (g_shape == Shape::STL) {
        std::cout << "Loaded STL " << g_stlPath << ": "
                  << g_boxKlein.size() << " tris, "
                  << g_uniqueVerts.size() << " unique verts (view " << g_view << ")\n";
    } else {
        std::cout << "Rebuilt " << shapeName(g_shape) << ": "
                  << (g_hollow ? "hollow" : "solid") << ", inset=" << g_inset
                  << " -> " << g_boxKlein.size() << " tris, "
                  << g_uniqueVerts.size() << " unique verts (view " << g_view << ")\n";
    }
}

// Switch the active polyhedron source and rebuild. Called from the middle-click
// "Polyhedron" submenu (direct selection, IDs 20-24).
void setShape(Shape s) {
    g_shape = s;
    rebuildMesh();
    std::cout << "Shape: " << shapeName(g_shape) << "\n";
}

// Cycle to the next polyhedron in submenu order (S key / "Next Polyhedron" entry).
void nextShape() {
    switch (g_shape) {
        case Shape::Cube:                   g_shape = Shape::Dodecahedron;           break;
        case Shape::Dodecahedron:           g_shape = Shape::TruncatedCube;          break;
        case Shape::TruncatedCube:          g_shape = Shape::TruncatedOctahedron;    break;
        case Shape::TruncatedOctahedron:    g_shape = Shape::TruncatedCuboctahedron; break;
        case Shape::TruncatedCuboctahedron: g_shape = Shape::Cube;                   break;
        case Shape::STL:                     g_shape = Shape::Cube;                   break;  // S exits a loaded STL back to Cube
    }
    rebuildMesh();
    std::cout << "Shape: " << shapeName(g_shape) << "\n";
}

// Scrub the hollow border inset by `delta` (clamped to [0.05, 0.95]). A solid
// mesh ignores inset, so we only rebuild (and reindex vertices) when hollow.
void scrubInset(double delta) {
    if (g_shape == Shape::STL) { std::cout << "STL mesh: hollow/inset do not apply\n"; return; }
    g_inset += delta;
    if (g_inset > 0.95) g_inset = 0.95;
    if (g_inset < 0.05) g_inset = 0.05;
    std::cout << "inset = " << g_inset
              << (g_hollow ? "  -> rebuilding" : "  (no effect unless hollow: 'o')") << "\n";
    if (g_hollow) rebuildMesh();
}

// Toggle solid (92-tri) <-> hollow (288-tri frame) and rebuild the mesh.
void toggleHollowShape() {
    if (g_shape == Shape::STL) { std::cout << "STL mesh: hollow/inset do not apply\n"; return; }
    g_hollow = !g_hollow;
    std::cout << "hollow " << (g_hollow ? "ON (frame)" : "OFF (solid)") << "\n";
    rebuildMesh();
}

// Export the current shape (at the active hollow/inset) as an ASCII STL.
void writeShapeSTL() {
    if (g_shape != Shape::TruncatedCuboctahedron) {
        std::cout << "STL export is only implemented for the Truncated Cuboctahedron "
                  << "(current shape: " << shapeName(g_shape) << "). Not written.\n";
        return;
    }
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
            int R, G, B;
            if (g_colorMode == 1) {
                // special 5-bucket grayscale criterion (black → grey → white)
                int gb = greyByte(greyBucket(active[i]));
                R = G = B = gb;
            } else {
                // jet colormap across the facets, rotated by g_colorPhase so every
                // remesh() press gives a fresh jet distribution
                float t = float(i) / float(total) + g_colorPhase;
                Color c = jet(t);
                R = int(c.r * 255); G = int(c.g * 255); B = int(c.b * 255);
            }
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
                int R, G, B;
                if (g_colorMode == 1) {
                    int gb = greyByte(greyBucket(m[i]));   // same grayscale criterion as the active mesh
                    R = G = B = gb;
                } else {
                    float t = float(i) / float(n) + g_colorPhase;   // same jet mapping as the active mesh
                    Color c = jet(t);
                    R = int(c.r * 255); G = int(c.g * 255); B = int(c.b * 255);
                }
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
        if (!g_rendering) drawVertexMarkers();   // hide interaction dots in rendered video


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
    //glBegin(GL_LINE_LOOP);
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
// Camera keyframe timeline (Blender-style)
//-----------------------------------------------------------------------------
// Snapshot the live orbit camera into ordered keyframes (K); during playback (P)
// the idle timer tweens between consecutive keyframes and writes the result back
// into the live globals, so any manual orbit/zoom/pan done BETWEEN saves is
// ignored: the path A->B depends only on A, B, and the easing — not on how you
// wandered. Loops forever. Euler + shortest-yaw rotation, smoothstep easing, lerp
// on the linear channels and the world-space pivot point.
struct CamKey {
    float    angleX, angleY;        // orbit pitch/yaw (degrees)
    float    zoom;                  // dolly distance factor
    float    panX, panY;            // screen-space pan
    Vector3D target;                // focusTarget() resolved at save (world space)
    int      focusVert;             // g_focusVert at save (restored on stop)
};
static std::vector<CamKey> g_timeline;       // ordered keyframes = the timeline
static bool     g_playing      = false;       // playback active (clock advances)?
static double   g_playSegmentT = 0.0;         // seconds into the current segment
static int      g_playSeg      = 0;           // tweens key[seg] -> key[(seg+1)%n]
static Vector3D g_playTarget;                // tweened pivot (world space)
static double   g_keyDur       = 2.0;         // seconds per segment (tune here)

// In-between keyframing: a cursor set by "Go To Keyframe" / "Go To Last Keyframe".
// -1 = append mode (K adds at the end, as before); >=0 = "Insert Keyframe Here" (I)
// drops a new keyframe right after keyframe g_insertAfter, i.e. between it and the
// next one — so you can refine any transition, not just add at the tail.
static int g_insertAfter = -1;

// GLUT menu identifiers so createUI() can destroy + rebuild the menus at runtime
// (the "Go To Keyframe" submenu must list one entry per keyframe, and keyframes are
// added/removed while the app runs). Zero until createUI() first runs.
static int g_menuPoly = 0, g_menuGoto = 0, g_menuMain = 0, g_menuStl = 0, g_menuRender = 0;

// Set by any timeline mutation (save/insert/drop/clear/load) to request a menu
// rebuild. Serviced by backgroundTimer(), NOT synchronously — destroying a GLUT
// menu from inside its own callback can be use-after-free, so the rebuild is
// deferred to the next timer tick (a safe context). ~30ms latency, imperceptible.
static bool g_menuDirty = false;

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
void scanSTLFolder();                         // refresh g_stlFiles from g_stlDir (Load STL submenu)
void loadSTLFile(const std::string& path);     // load an ASCII STL as the active mesh
void renderVideo(int fps);                     // render one timeline loop to ~/Downloads MP4 via ffmpeg
void renderCliTimer(int value);                // one-shot timer: fire renderVideo() once the window is realized
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

// World-space camera pivot: the right-clicked (focused) vertex, or the origin when
// none is focused. Resolved in the CURRENT view, so the pivot auto-follows the same
// logical vertex when the user switches Klein / Poincaré / round-trip (1/2/3).
Vector3D focusTarget() {
    if (g_focusVert >= 0 && g_focusVert < (int)g_uniqueVerts.size()) {
        const VertRef& r = g_uniqueVerts[g_focusVert];
        return meshOf(r)[r.facetIdx][r.cornerIdx];
    }
    return Vector3D();   // (0,0,0) — Vector3D defaults to origin
}

//-----------------------------------------------------------------------------
// Camera keyframe timeline helpers (Blender-style key framing)
//-----------------------------------------------------------------------------

// smoothstep easing: 0->1 with zero derivative at both ends — the ease-in/ease-out
// "Bezier feel" without needing actual Bezier handle math.
static double smoothstepD(double x) {
    if (x < 0.0) x = 0.0; if (x > 1.0) x = 1.0;
    return x * x * (3.0 - 2.0 * x);
}

// Advance the playback clock by dt seconds and write the tweened camera into the
// live globals (g_angleX/Y, g_zoom, g_panX/Y) plus g_playTarget. Called every timer
// tick. Loops: when a segment ends, g_playSeg wraps modulo n so a 2-key timeline
// plays A->B->A->B... and the last segment tweens back to keyframe 0.
static void advancePlayback(double dt) {
    const int n = (int)g_timeline.size();
    if (n == 0) { g_playing = false; return; }      // nothing to play
    if (n == 1) {                                   // single keyframe: hold it (no motion)
        const CamKey& k = g_timeline[0];
        g_angleX = k.angleX; g_angleY = k.angleY; g_zoom = k.zoom;
        g_panX = k.panX;   g_panY = k.panY;   g_playTarget = k.target;
        return;
    }
    g_playSegmentT += dt;
    while (g_playSegmentT >= g_keyDur) {
        g_playSegmentT -= g_keyDur;
        g_playSeg = (g_playSeg + 1) % n;             // loop wraparound
        if (!g_rendering)                            // side mission: announce the current
            std::cout << "Playback: now at keyframe " // keyframe during interactive playback
                      << (g_playSeg + 1) << "/" << n << "\n"; // (kept quiet during a video render)
    }
    const CamKey& a = g_timeline[g_playSeg];
    const CamKey& b = g_timeline[(g_playSeg + 1) % n];
    const double p = g_playSegmentT / g_keyDur;       // 0..1 within this segment
    const double e = smoothstepD(p);
    // linear channels
    g_zoom   = (float)(a.zoom   + e * (b.zoom   - a.zoom));
    g_panX   = (float)(a.panX   + e * (b.panX   - a.panX));
    g_panY   = (float)(a.panY   + e * (b.panY   - a.panY));
    g_angleX = (float)(a.angleX + e * (b.angleX - a.angleX));
    // yaw: take the SHORTER angular path so the camera never spins the long way
    double dy = b.angleY - a.angleY;
    while (dy >  180.0) dy -= 360.0;
    while (dy < -180.0) dy += 360.0;
    g_angleY = (float)(a.angleY + e * dy);
    // pivot: lerp world-space target points (smooth focus handoff between keys)
    g_playTarget = a.target + (e * (b.target - a.target));
}

// Snapshot the current live camera into the timeline (K).
static void saveKeyframe() {
    CamKey k;
    k.angleX = g_angleX; k.angleY = g_angleY; k.zoom = g_zoom;
    k.panX = g_panX;     k.panY = g_panY;
    k.target = focusTarget();
    k.focusVert = g_focusVert;
    g_timeline.push_back(k);
    std::cout << "Keyframe " << g_timeline.size() << " saved"
              << " (yaw=" << k.angleY << " pitch=" << k.angleX
              << " zoom=" << k.zoom << ")\n";
    g_menuDirty = true;                // rebuild "Go To Keyframe" submenu (count changed)
}

// Play/Pause (P). On start, jump to keyframe 0; on stop, snap to the nearest
// keyframe so focusTarget() matches the stored target (no pivot jump back to
// manual control).
static void togglePlayback() {
    if (g_timeline.empty()) {
        std::cout << "Timeline empty — press K to add keyframes\n";
        return;
    }
    if (!g_playing) {
        g_playing = true;
        g_playSeg = 0; g_playSegmentT = 0.0;
        const CamKey& k0 = g_timeline[0];            // init to key 0 so the first
        g_angleX = k0.angleX; g_angleY = k0.angleY; // rendered frame is correct
        g_zoom = k0.zoom; g_panX = k0.panX; g_panY = k0.panY; g_playTarget = k0.target;
        std::cout << "Playback: PLAY (" << g_timeline.size() << " keyframes, loop, "
                  << g_keyDur << "s/seg)\n";
    } else {
        g_playing = false;
        int n = (int)g_timeline.size();
        int snap = (n >= 2 && (g_playSegmentT / g_keyDur) > 0.5)
                 ? (g_playSeg + 1) % n : g_playSeg;
        const CamKey& k = g_timeline[snap];
        g_angleX = k.angleX; g_angleY = k.angleY; g_zoom = k.zoom;
        g_panX = k.panX; g_panY = k.panY; g_focusVert = k.focusVert;
        std::cout << "Playback: PAUSED — snapped to keyframe " << snap << "\n";
    }
}

// Clear the whole timeline (T).
static void clearTimeline() {
    if (g_playing) { g_playing = false; std::cout << "Playback stopped\n"; }
    size_t had = g_timeline.size();
    g_timeline.clear();
    g_playSeg = 0; g_playSegmentT = 0.0;
    g_insertAfter = -1;                // cursor invalid once the timeline is empty
    std::cout << "Timeline cleared (" << had << " keyframes removed)\n";
    g_menuDirty = true;                // rebuild "Go To Keyframe" submenu
}

// Remove the most recent keyframe (Backspace).
static void dropLastKeyframe() {
    if (g_timeline.empty()) { std::cout << "Timeline empty\n"; return; }
    g_timeline.pop_back();
    if ((int)g_timeline.size() <= g_playSeg) { g_playSeg = 0; g_playSegmentT = 0.0; }
    if (g_insertAfter >= (int)g_timeline.size()) g_insertAfter = -1;
    std::cout << "Removed last keyframe; " << g_timeline.size() << " remain\n";
    g_menuDirty = true;                // rebuild "Go To Keyframe" submenu (count changed)
}

// Jump the live camera to keyframe i (1-based in the menu, 0-based here): load its
// pose into the live globals, restore the focus vertex it was saved with, stop
// playback so the jumped pose holds (else the timer would clobber it next tick),
// and set the insert cursor to i so a subsequent "Insert Here" drops a new keyframe
// right after it (i.e. between keyframe i and the next). Called from the dynamic
// "Go To Keyframe" submenu (ids 1000+i) and from jumpToLastKeyframe().
static void goToKeyframe(int i) {
    const int n = (int)g_timeline.size();
    if (n == 0) { std::cout << "Timeline empty — press K to add keyframes\n"; return; }
    if (i < 0) i = 0; if (i >= n) i = n - 1;
    const CamKey& k = g_timeline[i];
    g_playing = false;                       // hold the jumped pose; don't let the timer tween over it
    g_angleX = k.angleX; g_angleY = k.angleY; g_zoom = k.zoom;
    g_panX = k.panX;     g_panY = k.panY;
    g_focusVert = k.focusVert;               // so focusTarget() == k.target (no pivot jump)
    g_playTarget = k.target;
    g_insertAfter = i;                        // next "Insert Here" goes right after this keyframe
    std::cout << "Jumped to keyframe " << (i + 1) << "/" << n
              << "  (insert-after cursor = " << (i + 1) << ")\n";
}

// Convenience: jump to the last keyframe so you can keep appending/inserting from
// the end of the timeline (menu "Go To Last Keyframe").
static void jumpToLastKeyframe() {
    goToKeyframe((int)g_timeline.size() - 1);
}

// Insert the current live camera as a NEW keyframe right after the insert cursor
// (set by the last "Go To Keyframe"), shifting the rest right — i.e. drop it
// BETWEEN that keyframe and the next. With no cursor (g_insertAfter < 0) or a cursor
// past the end, this is just an append (same as K). The cursor advances to the new
// keyframe so repeated inserts stack in order. Rebuilds the "Go To Keyframe"
// submenu so the new entry appears. (Keyboard 'I' / menu.)
static void insertKeyframeHere() {
    CamKey k;
    k.angleX = g_angleX; k.angleY = g_angleY; k.zoom = g_zoom;
    k.panX = g_panX;     k.panY = g_panY;
    k.target = focusTarget();
    k.focusVert = g_focusVert;
    const int n = (int)g_timeline.size();
    int pos;
    if (g_insertAfter < 0 || g_insertAfter >= n) {     // no cursor / past end -> append
        g_timeline.push_back(k);
        pos = n;
    } else {
        pos = g_insertAfter + 1;                       // insert right after the cursor
        g_timeline.insert(g_timeline.begin() + pos, k);
    }
    g_insertAfter = pos;                               // stack further inserts in order
    g_playing = false;
    g_menuDirty = true;                                // rebuild the "Go To Keyframe" submenu
    std::cout << "Inserted keyframe at position " << (pos + 1)
              << "  (" << g_timeline.size() << " total)\n";
}

//-----------------------------------------------------------------------------
// STL loading — pick any ASCII .stl in g_stlDir from the "Load STL..." submenu
// (nested in the Polyhedron submenu), or via the -stl CLI flag / L key. The loader
// itself (loadSTL_ascii) already existed but was unused; this wires it in.
//-----------------------------------------------------------------------------

// Refresh g_stlFiles with every .stl in g_stlDir (case-insensitive extension), sorted
// so the "Load STL..." submenu is stable across rescans. Called from createUI() on
// every (re)build, so newly-added files show up after a rebuild or "Rescan STL folder".
// A missing/unreadable folder just leaves the list empty (the submenu shows a
// placeholder) — no throw.
void scanSTLFolder() {
    g_stlFiles.clear();
    namespace fs = std::filesystem;
    std::error_code ec;
    if (!fs::is_directory(g_stlDir, ec)) return;          // missing folder -> empty list
    std::vector<std::string> names;
    for (auto& entry : fs::directory_iterator(g_stlDir, ec)) {
        std::string ext = entry.path().extension().string();
        for (auto& c : ext) c = (char)std::tolower((unsigned char)c);
        if (ext == ".stl") names.push_back(entry.path().filename().string());
    }
    std::sort(names.begin(), names.end());
    g_stlFiles = std::move(names);
}

// Load an ASCII STL as the active mesh. Parses via loadSTL_ascii (throws on a bad
// file); assigns g_stlMesh only on success so a failure leaves all state intact,
// switches g_shape to STL, rebuilds the three views, and flags a menu rebuild so
// the "Load STL..." submenu marks the active file. Called from the submenu, the
// L key, and the -stl CLI flag.
void loadSTLFile(const std::string& path) {
    try {
        FacetBox m = loadSTL_ascii(path);   // throws std::runtime_error on open/parse failure
        g_stlMesh = std::move(m);
        g_stlPath = path;
        g_shape = Shape::STL;
        rebuildMesh();                      // centers+scales, re-derives hyper/back, reindexes, clears spawns/selection
        g_menuDirty = true;                 // refresh the Load STL submenu's active marker (deferred to the timer)
        std::cout << "Loaded STL: " << path << "\n";
    } catch (const std::exception& e) {
        std::cerr << "STL load failed: " << e.what() << "\n";
    }
}

//-----------------------------------------------------------------------------
// Session save/load (.h3geo) — snapshot the whole working scene to disk and back.
// Text format, parsed with >>; doubles at precision 17 so coordinates round-trip.
// Saves: the loaded mesh (Klein view — the canonical; hyper/round-trip re-derived
// on load, exactly as rebuildMesh()/remesh() do), every spawned inversion (its
// mesh + center/radius/view/batch), the camera keyframes, the selection, and the
// scalar prefs (shape/hollow/inset/view/color/bg/focus/keyDur). Load via CLI
// "-load <path>" or the menu "Load Session" (which loads ./session.h3geo).
//-----------------------------------------------------------------------------
static void saveSession() {
    const std::string path = "session.h3geo";
    std::ofstream out(path);
    if (!out) { std::cout << "Save failed: cannot open " << path << " for writing\n"; return; }
    out.precision(17);                 // round-trip doubles exactly

    out << "H3GEO 1\n";
    out << (int)g_shape << ' ' << (g_hollow ? 1 : 0) << ' ' << g_inset << ' '
        << g_view << ' ' << g_colorMode << ' ' << g_colorPhase << ' '
        << g_bgLightness << ' ' << (g_showInvSpheres ? 1 : 0) << ' '
        << g_focusVert << ' ' << g_keyDur << "\n";

    // Loaded mesh (Klein) — the canonical view. Normals are NOT saved: Facet(A,B,C)
    // recomputes them from the cross product on construction.
    out << "klein " << g_boxKlein.size() << "\n";
    for (size_t i = 0; i < g_boxKlein.size(); ++i) {
        const Facet& f = g_boxKlein[i];
        Vector3D a = f[0], b = f[1], c = f[2];
        out << a.x() << ' ' << a.y() << ' ' << a.z() << ' '
            << b.x() << ' ' << b.y() << ' ' << b.z() << ' '
            << c.x() << ' ' << c.y() << ' ' << c.z() << "\n";
    }

    // Spawned inversions: each child lives in its spawn view's space and is rendered
    // as-is, so store its FacetBox directly (no re-projection on load).
    out << "spawns " << g_spawned.size() << "\n";
    for (auto const& s : g_spawned) {
        out << s.center.x() << ' ' << s.center.y() << ' ' << s.center.z() << ' '
            << s.radius << ' ' << s.view << ' ' << s.batch << ' '
            << s.mesh.size() << "\n";
        for (size_t i = 0; i < s.mesh.size(); ++i) {
            const Facet& f = s.mesh[i];
            Vector3D a = f[0], b = f[1], c = f[2];
            out << a.x() << ' ' << a.y() << ' ' << a.z() << ' '
                << b.x() << ' ' << b.y() << ' ' << b.z() << ' '
                << c.x() << ' ' << c.y() << ' ' << c.z() << "\n";
        }
    }

    // Camera keyframes.
    out << "keyframes " << g_timeline.size() << "\n";
    for (auto const& k : g_timeline) {
        out << k.angleX << ' ' << k.angleY << ' ' << k.zoom << ' '
            << k.panX << ' ' << k.panY << ' '
            << k.target.x() << ' ' << k.target.y() << ' ' << k.target.z() << ' '
            << k.focusVert << "\n";
    }

    // Selection (indices into g_uniqueVerts; valid again after buildUniqueVertexList()
    // restores the same vertex order).
    out << "selected " << g_selectedVerts.size() << "\n";
    for (int idx : g_selectedVerts) out << idx << ' ';
    out << "\n";

    out.close();
    std::cout << "Saved session: " << path
              << "  (shape " << shapeName(g_shape)
              << ", klein " << g_boxKlein.size() << " tris, "
              << g_spawned.size() << " spawns, "
              << g_timeline.size() << " keyframes)\n";
    std::cout << "  Load with: ./main9.3.7 -load " << path
              << "   (or menu: Load Session)\n";
}

// Load a .h3geo session. Parses into locals first; on ANY parse failure it prints an
// error and returns false WITHOUT touching live state. Only on full success does it
// commit. Selection/focus are restored AFTER buildUniqueVertexList (which clears
// both, since old indices are invalid until the vertex list is rebuilt).
static bool loadSession(const std::string& path) {
    std::ifstream in(path);
    if (!in) { std::cout << "Load failed: cannot open " << path << "\n"; return false; }
    std::string magic; int ver;
    in >> magic >> ver;
    if (magic != "H3GEO") {
        std::cout << "Load failed: not an H3GEO file (" << path << ")\n"; return false;
    }

    int shape, hollow, view, colorMode, showInvSpheres, focusVert;
    double inset, colorPhase, bgLightness, keyDur;
    in >> shape >> hollow >> inset >> view >> colorMode >> colorPhase
       >> bgLightness >> showInvSpheres >> focusVert >> keyDur;
    if (!in) { std::cout << "Load failed: truncated header\n"; return false; }

    // Loaded mesh (Klein)
    std::string tag; size_t nKlein;
    in >> tag >> nKlein;
    if (!in || tag != "klein") { std::cout << "Load failed: expected 'klein' section\n"; return false; }
    FacetBox klein;
    for (size_t i = 0; i < nKlein; ++i) {
        double ax,ay,az,bx,by,bz,cx,cy,cz;
        in >> ax >> ay >> az >> bx >> by >> bz >> cx >> cy >> cz;
        if (!in) { std::cout << "Load failed: truncated klein tris\n"; return false; }
        klein.push(Vector3D(ax,ay,az), Vector3D(bx,by,bz), Vector3D(cx,cy,cz));
    }

    // Spawns
    size_t nSpawns;
    in >> tag >> nSpawns;
    if (!in || tag != "spawns") { std::cout << "Load failed: expected 'spawns' section\n"; return false; }
    std::vector<SpawnRecord> spawned;
    for (size_t s = 0; s < nSpawns; ++s) {
        double ccx,ccy,ccz, radius; int sv, batch; size_t tk;
        in >> ccx >> ccy >> ccz >> radius >> sv >> batch >> tk;
        if (!in) { std::cout << "Load failed: truncated spawn header\n"; return false; }
        SpawnRecord rec;
        rec.center = Vector3D(ccx,ccy,ccz);
        rec.radius = radius; rec.view = sv; rec.batch = batch;
        for (size_t i = 0; i < tk; ++i) {
            double ax,ay,az,bx,by,bz,cx,cy,cz;
            in >> ax >> ay >> az >> bx >> by >> bz >> cx >> cy >> cz;
            if (!in) { std::cout << "Load failed: truncated spawn tris\n"; return false; }
            rec.mesh.push(Vector3D(ax,ay,az), Vector3D(bx,by,bz), Vector3D(cx,cy,cz));
        }
        spawned.push_back(std::move(rec));
    }

    // Keyframes
    size_t nKeys;
    in >> tag >> nKeys;
    if (!in || tag != "keyframes") { std::cout << "Load failed: expected 'keyframes' section\n"; return false; }
    std::vector<CamKey> timeline;
    for (size_t i = 0; i < nKeys; ++i) {
        CamKey k; double tx,ty,tz;
        in >> k.angleX >> k.angleY >> k.zoom >> k.panX >> k.panY
           >> tx >> ty >> tz >> k.focusVert;
        if (!in) { std::cout << "Load failed: truncated keyframes\n"; return false; }
        k.target = Vector3D(tx,ty,tz);
        timeline.push_back(k);
    }

    // Selection (best-effort: out-of-range indices are filtered out after the index
    // rebuild, since the file may have been saved under a different mesh).
    size_t nSel;
    in >> tag >> nSel;
    if (!in || tag != "selected") { std::cout << "Load failed: expected 'selected' section\n"; return false; }
    std::set<int> sel;
    for (size_t i = 0; i < nSel; ++i) {
        int idx; in >> idx;
        if (!in) break;
        sel.insert(idx);
    }

    // ── commit (all parsing succeeded) ──
    g_shape        = (Shape)shape;
    g_hollow       = (hollow != 0);
    g_inset        = inset;
    g_view         = view;
    g_colorMode    = colorMode;
    g_colorPhase   = (float)colorPhase;
    g_bgLightness  = bgLightness;
    g_showInvSpheres = (showInvSpheres != 0);
    g_keyDur       = keyDur;

    g_boxKlein = std::move(klein);
    g_boxHyper = g_boxKlein.hyperboloid();          // re-derive Poincaré view
    g_boxBack  = inverseHyperboloid(g_boxHyper);    // re-derive round-trip view
    g_spawned  = std::move(spawned);
    g_undone.clear();                               // redo stack not persisted

    g_timeline  = std::move(timeline);
    g_playing   = false;
    g_playSeg   = 0;
    g_playSegmentT = 0.0;
    g_insertAfter = -1;

    buildUniqueVertexList();                         // reindex base + spawns (identical
                                                     // content -> saved indices stay valid)

    g_selectedVerts.clear();
    for (int idx : sel)
        if (idx >= 0 && idx < (int)g_uniqueVerts.size())
            g_selectedVerts.insert(idx);
    g_focusVert = (focusVert >= 0 && focusVert < (int)g_uniqueVerts.size()) ? focusVert : -1;

    g_menuDirty = true;                              // rebuild "Go To Keyframe" submenu

    std::cout << "Loaded session: " << path
              << "  (shape " << shapeName(g_shape)
              << ", klein " << g_boxKlein.size() << " tris, "
              << g_spawned.size() << " spawns, "
              << g_timeline.size() << " keyframes)\n";
    glutPostRedisplay();
    return true;
}

// The exact camera modelview used by display() AND pickVertex(). Same as the old
// orbit-at-origin camera (push back, rotate X, rotate Y) with one added innermost
// translate of -focusTarget(): the focused vertex maps to the origin, so the scene
// rotates and zooms about that vertex instead of the world origin. With no focus
// (target = 0) the translate is a no-op and behavior is identical to main9.3.5.
void applyCameraModelview() {
    glLoadIdentity();
    // Camera-space translate (leftmost op): screen-space pan (g_panX/g_panY) plus
    // the fixed dolly push-back along -Z (5*g_zoom). g_panX/g_panY were updated on
    // middle-drag but never applied before — now they pan the view on screen.
    glTranslatef(g_panX, g_panY, -5.0f * g_zoom);
    glRotatef(g_angleX, 1.0f, 0.0f, 0.0f);
    glRotatef(g_angleY, 0.0f, 1.0f, 0.0f);
    // During playback, orbit about the tweened world-space target; otherwise about
    // the live focus vertex (or origin). With no focus and not playing, the
    // translate is a no-op and behavior is identical to before.
    Vector3D t = g_playing ? g_playTarget : focusTarget();
    glTranslatef(-(float)t.x(), -(float)t.y(), -(float)t.z());
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
    g_focusVert = -1;    // indices into the old g_uniqueVerts are now invalid
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

    // Focused (right-clicked) vertex: black outline + blue core, drawn on top of
    // the orange selection so a vertex that is both selected and focused reads blue.
    if (g_focusVert >= 0 && g_focusVert < (int)g_uniqueVerts.size()) {
        const VertRef& r = g_uniqueVerts[g_focusVert];
        const Vector3D& p = meshOf(r)[r.facetIdx][r.cornerIdx];

        glPointSize(14.0f);
        glColor3ub(0, 0, 0);                 // black outline
        glBegin(GL_POINTS);
        glVertex3f(p.x(), p.y(), p.z());
        glEnd();

        glPointSize(9.0f);
        glColor3ub(60, 160, 255);            // blue core
        glBegin(GL_POINTS);
        glVertex3f(p.x(), p.y(), p.z());
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
    applyCameraModelview();          // same transform as display() (incl. focus pivot)
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

// Right-click on a vertex -> make it the camera pivot (drawn blue) and look at it,
// so subsequent rotate/zoom happen about that vertex. Right-click on empty space
// does nothing (keeps the current focus/lookat), per the design choice.
void handleVertexFocus(int mouseX, int mouseY) {
    int picked = pickVertex(mouseX, mouseY, g_pickRadiusPx);
    if (picked >= 0) {
        g_focusVert = picked;
        std::cout << "Focus vertex " << picked << " -> camera lookat\n";
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

// R (centroid3) / M (midpoint4): subdivide EVERY facet once — the base mesh AND every
// spawned child — so the whole scene gets finer together. One pass per press; press
// again to refine deeper. The base is refined in Klein (the source of truth) and the 3
// views re-derived; each child is refined in place in its own view's coords. Refinement
// verts are convex combinations of existing verts, so the base stays inside the 0.9 ball
// (hyperboloid() never throws) and children stay within their hulls. fromFacet and
// subdivide4 both preserve CCW winding, so refined children stay outward (no
// ensureOutwardNormals needed). Not on the spawn undo stack; reset by rebuilding the
// shape (toggle hollow O, scrub inset ,/.).
void remesh(FacetBox::SubdivisionMode mode) {
    if (g_boxKlein.size() == 0) { std::cout << "Remesh: no mesh loaded.\n"; return; }
    g_boxKlein = g_boxKlein.refine(1, mode);       // FacetBox.hpp:275
    g_boxHyper = g_boxKlein.hyperboloid();         // FacetBox.hpp:257
    g_boxBack  = inverseHyperboloid(g_boxHyper);   // local (line 286)
    size_t kids = 0;
    for (size_t s = 0; s < g_spawned.size(); ++s) {
        if (g_spawned[s].mesh.size() == 0) continue;
        g_spawned[s].mesh = g_spawned[s].mesh.refine(1, mode);  // winding preserved -> outward stays
        ++kids;
    }
    buildUniqueVertexList();   // reindex base + all spawns; clears selection (verts changed)
    // Rotate the jet palette so this press gets a fresh color distribution:
    // golden-ratio conjugate keeps consecutive passes in distinct bands and gives a
    // long cycle before the palette repeats.
    g_colorPhase += 0.3819660f;
    if (g_colorPhase >= 1.0f) g_colorPhase -= 1.0f;
    const char* name = (mode == FacetBox::SubdivisionMode::Centroid3) ? "centroid3" : "midpoint4";
    std::cout << "Remesh (" << name << "): base -> " << g_boxKlein.size() << " tris, "
              << kids << " child" << (kids == 1 ? "" : "ren") << " refined -> "
              << g_uniqueVerts.size() << " unique verts (view " << g_view
              << ", jet phase " << g_colorPhase << ")\n";
    glutPostRedisplay();
}

//-----------------------------------------------------------------------------
// Main
//-----------------------------------------------------------------------------
int main(int argc, char** argv)
{
    // ── Background wash lightness — tune HERE, not in any menu/key ──────────
    // 0.0 = saturated blue/red, 1.0 = white. 0.78 = the very light pastel wash
    // that won't fight the objects; 0.48 (current) = noticeably more color.
    // Lower = more color, higher = paler.
    g_bgLightness = 0.48;

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

    // Optional: load a saved session (.h3geo) and/or an ASCII STL at startup.
    //   ./main9.3.7 -load session.h3geo
    //   ./main9.3.7 -stl path.stl
    // A missing/invalid file prints an error and continues with the default scene.
    // Both flags may be given; the later one determines the active mesh.
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-load" && i + 1 < argc)      loadSession(argv[++i]);
        else if (a == "-stl"  && i + 1 < argc) loadSTLFile(argv[++i]);
        else if (a == "-render" && i + 1 < argc) { g_cliRender = true; g_cliRenderFps = atoi(argv[++i]); }
    }

    // Register callbacks
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouseButton);
    glutMotionFunc(mouseMotion);
    glutKeyboardFunc(keyboard);

    // Drive the soft breathing background gradient (~30 fps) even while idle.
    glutTimerFunc(30, backgroundTimer, 0);

    // CLI -render: fire one render pass once the window is realized, then exit.
    if (g_cliRender) glutTimerFunc(250, renderCliTimer, 0);

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
        "R-click: Focus vertex (blue) + lookat",
        "L-drag: Rotate",
        "M-click: UI Menu",
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
        "[R]: Remesh centroid3 (finer)",
        "[M]: Remesh midpoint4 (finer)",
        "[G]: Color mode: jet / grey 5-bucket",
        "[,]/[.]: Decrease/Increase inset",
        "[O]: Toggle hollow (solid/frame)",
        "[W]: Write STL of current shape",
        "[L]: Reload last STL",
        "[S]: Next polyhedron",
        "[K]: Save camera keyframe",
        "[I]: Insert keyframe (after Go-To)",
        "[P]: Play/Pause camera timeline (loop)",
        "[T]: Clear camera timeline",
        "[Backspace]: Remove last keyframe",
        "Menu: Go To/Insert keyframe, Save/Load .h3geo",
        "Menu: Polyhedron > Load STL (from ~/Downloads)",
        "Menu: Render Video > fps -> ~/Downloads MP4 1280x720 (no HUD/cursor)",
        "[H]: Toggle Help",
        "Menu (M-click): Reset LookAt to Origin"
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

    // Current shape (polyhedron + hollow/inset)
    y -= 0.05f;
    char shapeBuf[128];
    snprintf(shapeBuf, sizeof(shapeBuf), "Shape: %s [%s]  inset=%.2f  (%lu tris)",
             shapeName(g_shape), g_hollow ? "hollow (frame)" : "solid", g_inset,
             (unsigned long)g_boxKlein.size());
    glRasterPos2f(0.02f, y);
    for(const char* c=shapeBuf; *c; ++c)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);

    // Current per-facet coloring mode
    y -= 0.05f;
    const char* colorName = (g_colorMode == 0) ? "Color: jet (rotated per remesh)"
                                              : "Color: grayscale 5-bucket (black->white)";
    glRasterPos2f(0.02f, y);
    for(const char* c=colorName; *c; ++c)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);

    // Camera keyframe timeline status (with the in-between insert cursor)
    y -= 0.05f;
    char camBuf[128];
    char insBuf[24];
    if (g_insertAfter < 0)
        snprintf(insBuf, sizeof(insBuf), "insert@end");
    else
        snprintf(insBuf, sizeof(insBuf), "insert@after %d", g_insertAfter + 1);
    if (g_timeline.empty())
        snprintf(camBuf, sizeof(camBuf), "Cam timeline: empty  %s", insBuf);
    else
        snprintf(camBuf, sizeof(camBuf), "Cam timeline: %lu keyframes  %s  seg %d/%lu  %s",
                 (unsigned long)g_timeline.size(),
                 g_playing ? "PLAY (loop)" : "idle",
                 g_playSeg, (unsigned long)g_timeline.size(), insBuf);
    glRasterPos2f(0.02f, y);
    for(const char* c=camBuf; *c; ++c)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);

    glMatrixMode(GL_MODELVIEW); glPopMatrix();
    glMatrixMode(GL_PROJECTION); glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
}

//-----------------------------------------------------------------------------
// Soft animated background: a very light wash that breathes between pale blue
// and pale rose (a desaturated "jet" feel) and back. Drawn first in display(),
// in screen space with lighting + depth off, so it always sits behind the 3D
// objects and never competes with their colors. Both endpoints ease toward each
// other over ~20 s, so the cool end and the warm end slowly swap: light blue ->
// light red -> light blue. All matrix + enable state is saved/restored, so the
// 3D perspective projection and lighting are untouched when the scene draws next.
//-----------------------------------------------------------------------------
void drawBackgroundGradient() {
    // Saturated base hues: the "most color" the wash gets at g_bgLightness = 0.
    // Chosen so g_bgLightness = 0.78 reproduces the original very-light pastels —
    // pale azure (0.80,0.88,0.95) and pale rose (0.96,0.85,0.82). Raising
    // g_bgLightness blends both toward white (paler); lowering it brings out more
    // color. g_bgLightness is set in main(), not in any menu/key, so this is a
    // compile-time knob, not an interactive control.
    const float blueR = 0.09f, blueG = 0.45f, blueB = 0.77f;   // base azure (g_bgLightness = 0)
    const float redR  = 0.82f, redG  = 0.32f, redB  = 0.18f;   // base rose (g_bgLightness = 0)

    // Blend each base hue toward white by g_bgLightness, clamped to [0,1].
    const double Ld = (g_bgLightness < 0.0) ? 0.0
                    : (g_bgLightness > 1.0) ? 1.0 : g_bgLightness;
    const float L = (float)Ld;
    float bR = blueR + L * (1.0f - blueR);
    float bG = blueG + L * (1.0f - blueG);
    float bB = blueB + L * (1.0f - blueB);
    float rR = redR  + L * (1.0f - redR);
    float rG = redG  + L * (1.0f - redG);
    float rB = redB  + L * (1.0f - redB);

    // Wall-clock phase so the breathing is smooth regardless of frame rate; the
    // timer below only triggers redraws, the speed comes from elapsed time.
    const double periodMs = 20000.0;                          // ~20 s per blue->red->blue
    // During a video render the capture loop is blocking, so wall-clock would advance by
    // the uneven real encode time per frame and the wash would jitter. Use the video frame
    // time instead so the breathing stays smooth in the output.
    double t = g_rendering ? ((double)g_renderFrame / (double)g_renderFps * 1000.0 / periodMs)
                           : ((double)glutGet(GLUT_ELAPSED_TIME) / periodMs);   // cycles
    double p = 0.5 + 0.5 * sin(t * 2.0 * M_PI);               // [0,1], smooth (triangle-free)

    // Top eases blue->red while bottom eases red->blue, so the vertical wash
    // inverts over time (p=0: top blue / bottom red; p=1: top red / bottom blue;
    // p=0.5: a uniform pale lavender).
    float topR = (float)(bR + p * (rR - bR));
    float topG = (float)(bG + p * (rG - bG));
    float topB = (float)(bB + p * (rB - bB));
    float botR = (float)(rR + p * (bR - rR));
    float botG = (float)(rG + p * (bG - rG));
    float botB = (float)(rB + p * (bB - rB));

    glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);

    glMatrixMode(GL_PROJECTION); glPushMatrix(); glLoadIdentity();
    glOrtho(0.0, 1.0, 0.0, 1.0, -1.0, 1.0);
    glMatrixMode(GL_MODELVIEW);  glPushMatrix(); glLoadIdentity();

    glBegin(GL_QUADS);
        glColor3f(botR, botG, botB); glVertex2f(0.0f, 0.0f);   // bottom-left
        glColor3f(botR, botG, botB); glVertex2f(1.0f, 0.0f);   // bottom-right
        glColor3f(topR, topG, topB); glVertex2f(1.0f, 1.0f);   // top-right
        glColor3f(topR, topG, topB); glVertex2f(0.0f, 1.0f);   // top-left
    glEnd();

    glMatrixMode(GL_MODELVIEW);  glPopMatrix();
    glMatrixMode(GL_PROJECTION); glPopMatrix();
    glMatrixMode(GL_MODELVIEW);

    glPopAttrib();
}

//-----------------------------------------------------------------------------
// Idle timer that keeps the window redrawing (~30 fps) so the breathing
// background gradient animates even when the user isn't interacting. Re-arms
// itself; the only per-frame cost beyond the normal scene is one fullscreen
// gradient quad, so this is cheap.
//-----------------------------------------------------------------------------
void backgroundTimer(int /*value*/) {
    // Wall-clock delta since the last tick, so the breathing background AND the
    // camera keyframe playback advance in real seconds regardless of frame rate.
    static double lastMs = (double)glutGet(GLUT_ELAPSED_TIME);
    double nowMs = (double)glutGet(GLUT_ELAPSED_TIME);
    double dt = (nowMs - lastMs) * 0.001;            // seconds
    if (dt < 0.0) dt = 0.0;                           // clock-oddity guard
    if (dt > 0.1) dt = 0.1;                            // clamp stalls -> no teleport
    lastMs = nowMs;

    if (g_playing) advancePlayback(dt);              // tween the camera timeline

    // Deferred menu rebuild (set by a timeline mutation inside a menu callback —
    // rebuilding there could be use-after-free, so it happens here, safely).
    if (g_menuDirty) { g_menuDirty = false; createUI(); }

    glutPostRedisplay();
    glutTimerFunc(30, backgroundTimer, 0);
}

//-----------------------------------------------------------------------------
// Main display: apply interactive camera, then draw
//-----------------------------------------------------------------------------
void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Soft breathing background wash (pale blue <-> pale rose). Drawn before the
    // camera/scene, in screen space, restoring all state — see drawBackgroundGradient().
    drawBackgroundGradient();

    /* - Camera: pan, zoom, rotate (orbit about the focused vertex) */
    applyCameraModelview();   // push back, rotate X/Y, then translate -focusTarget()

    // Draw axes at origin
    //drawAxes(10.0f);


    ProcessingProto();   // calls Setup() then Draw()
    if (!g_rendering) drawHUD();   // hide on-screen UI from the captured video frames
    
    // Render-to-video: read the just-drawn 1280x720 offscreen framebuffer (no HUD, no
    // vertex dots) and pipe it to ffmpeg as one raw rgb24 frame. The FBO is bound by
    // renderVideo(); read its color attachment here.
    if (g_rendering && g_ffPipe) {
        glReadBuffer(GL_COLOR_ATTACHMENT0);
        glPixelStorei(GL_PACK_ALIGNMENT, 1);
        glReadPixels(0, 0, g_renderW, g_renderH, GL_RGB, GL_UNSIGNED_BYTE, g_frameBuf.data());
        size_t expect = (size_t)g_renderW * g_renderH * 3;
        if (fwrite(g_frameBuf.data(), 1, expect, g_ffPipe) != expect) {
            std::cerr << "\nRender: pipe write failed (ffmpeg died?) - aborting\n";
            g_renderAbort = true;     // renderVideo()'s capture loop watches this and stops
        }
    }

    // During a render the FBO (not the window) holds the frame; skip the window swap so
    // the on-screen window simply stays frozen on its last interactive frame.
    if (!g_rendering) glutSwapBuffers();
}

//-----------------------------------------------------------------------------
// (Re)create the offscreen 1280x720 framebuffer object used by renderVideo(), if it
// doesn't already exist, and report whether it is framebuffer-complete. Rendering into
// an FBO (instead of the on-screen window) decouples the video resolution from the
// interactive viewport: the window stays 720x720 while the video is exactly 16:9, with
// no window resize / event-loop pump. The FBO is reused across renders; GL reclaims it
// at process exit.
//-----------------------------------------------------------------------------
static bool ensureRenderFBO() {
    if (g_rFbo) return true;                 // already created at the fixed 1280x720

    glGenFramebuffers(1, &g_rFbo);
    glBindFramebuffer(GL_FRAMEBUFFER, g_rFbo);

    // Color: an RGB8 texture we glReadPixels from (rgb24 matches ffmpeg's raw input).
    glGenTextures(1, &g_rColorTex);
    glBindTexture(GL_TEXTURE_2D, g_rColorTex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, g_renderOutW, g_renderOutH, 0,
                 GL_RGB, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                           GL_TEXTURE_2D, g_rColorTex, 0);

    // Depth: the scene uses depth testing, so the FBO needs a depth renderbuffer.
    glGenRenderbuffers(1, &g_rDepthRb);
    glBindRenderbuffer(GL_RENDERBUFFER, g_rDepthRb);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24,
                          g_renderOutW, g_renderOutH);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                              GL_RENDERBUFFER, g_rDepthRb);

    GLenum st = glCheckFramebufferStatus(GL_FRAMEBUFFER);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    if (st != GL_FRAMEBUFFER_COMPLETE) {
        std::cerr << "Render: offscreen framebuffer incomplete - cannot render at "
                  << g_renderOutW << "x" << g_renderOutH << " (GL "
                  << (const char*)glGetString(GL_VERSION) << ")\n";
        g_rFbo = 0;                          // don't pretend it's ready on the next call
        return false;
    }
    return true;
}

//-----------------------------------------------------------------------------
// Render one full loop of the camera keyframe timeline to an MP4 via ffmpeg.
//
// Drives the playback deterministically: instead of the wall-clock backgroundTimer, we
// step advancePlayback(1/fps) per frame so every output frame is reproducible. Each
// frame is drawn (HUD and vertex-marker overlays gated off by g_rendering) into a
// 1280x720 16:9 offscreen framebuffer, read back with glReadPixels, and piped to ffmpeg
// as raw rgb24; -vf vflip undoes OpenGL's bottom-up row order. A carriage-return progress
// bar prints to stderr. This blocks the GLUT loop for the render duration (the window
// is intentionally unresponsive meanwhile) — safe because the GL context stays current
// and display() does no event processing.
//-----------------------------------------------------------------------------
void renderVideo(int fps) {
    const int n = (int)g_timeline.size();
    if (n < 2) {
        std::cerr << "Render: need >= 2 keyframes to animate (have " << n << ")\n";
        return;
    }

    // A dead ffmpeg mid-render would raise SIGPIPE on fwrite and kill the app with no
    // cleanup; ignore it so a broken pipe just yields a short write (handled below).
    std::signal(SIGPIPE, SIG_IGN);

    // Confirm ffmpeg is on PATH before opening the capture pipe.
    FILE* ver = popen("ffmpeg -version", "r");
    if (!ver || pclose(ver) != 0) {
        std::cerr << "Render: ffmpeg not found on PATH (install ffmpeg to render video)\n";
        return;
    }

    // Fixed output size: 1280x720 (16:9, standard YouTube HD). The video is rendered into
    // an offscreen framebuffer at this size, independent of the 720x720 on-screen window,
    // so the result is exactly 16:9 with no window-resize/event-pump jitter.
    g_renderW = g_renderOutW;   // 1280
    g_renderH = g_renderOutH;   // 720
    if (!ensureRenderFBO())     // (re)creates g_rFbo at 1280x720; prints + returns on failure
        return;

    // Timestamped output so successive renders don't clobber one another.
    char path[256];
    time_t now = time(nullptr);
    std::tm* lt = std::localtime(&now);
    strftime(path, sizeof(path),
             "/home/mike666/Downloads/render_%Y%m%d_%H%M%S.mp4", lt);

    // Pipe raw rgb24 frames into ffmpeg -> H.264 MP4 (yuv420p for broad compatibility).
    char cmd[512];
    snprintf(cmd, sizeof(cmd),
        "ffmpeg -y -f rawvideo -pixel_format rgb24 -video_size %dx%d -framerate %d -i - "
        "-vf vflip -c:v libx264 -preset fast -crf 18 -pix_fmt yuv420p \"%s\"",
        g_renderW, g_renderH, fps, path);
    g_ffPipe = popen(cmd, "w");
    if (!g_ffPipe) {
        std::cerr << "Render: cannot start ffmpeg\n";
        return;
    }

    // Take over display: hide the OS cursor (polish — it isn't in the framebuffer
    // anyway), drive playback ourselves, allocate the frame scratch buffer.
    glutSetCursor(GLUT_CURSOR_NONE);
    g_renderFps   = fps;
    g_rendering   = true;
    g_playing     = true;               // so applyCameraModelview() follows g_playTarget
    g_renderAbort = false;
    g_frameBuf.assign((size_t)g_renderW * g_renderH * 3, 0);

    // Start exactly at keyframe 0 so frame 0 reads as key 0 (not one step into the path).
    g_playSeg = 0; g_playSegmentT = 0.0;
    const CamKey& k0 = g_timeline[0];
    g_angleX = k0.angleX; g_angleY = k0.angleY; g_zoom = k0.zoom;
    g_panX = k0.panX;     g_panY = k0.panY;     g_playTarget = k0.target;

    // Redirect all drawing into the 1280x720 offscreen FBO for the capture loop, with a
    // 16:9 perspective (matching reshape()'s fovy/aspect logic) so the scene frames
    // correctly instead of being stretched. Nothing in the render path resets the
    // viewport or PROJECTION matrix, so this persists across the per-frame display() calls.
    glBindFramebuffer(GL_FRAMEBUFFER, g_rFbo);
    glViewport(0, 0, g_renderW, g_renderH);
    glMatrixMode(GL_PROJECTION); glLoadIdentity();
    gluPerspective(100.0, (double)g_renderW / (double)g_renderH, 0.001, 1000.0);
    glMatrixMode(GL_MODELVIEW);
    glReadBuffer(GL_COLOR_ATTACHMENT0);   // glReadPixels target while the FBO is bound

    g_renderTotal = (int)((double)n * g_keyDur * fps + 0.5);   // one full loop, rounded
    g_renderFrame = 0;
    std::cerr << "Rendering " << g_renderTotal << " frames (" << ((double)n * g_keyDur)
              << "s) @ " << fps << " fps, " << g_renderOutW << "x" << g_renderOutH
              << " -> " << path << "\n";

    bool aborted = false;
    for (int i = 0; i < g_renderTotal; ++i) {
        g_renderFrame = i;
        display();                       // draws scene (no HUD/markers) + captures + swaps
        if (g_renderAbort) { aborted = true; break; }   // ffmpeg pipe broke mid-render

        // Carriage-return terminal progress bar.
        int pct   = (g_renderTotal > 0) ? (i * 100 / g_renderTotal) : 100;
        const int barW = 30;
        int filled = (g_renderTotal > 0) ? (i * barW / g_renderTotal) : barW;
        fprintf(stderr, "\rRendering [");
        for (int b = 0; b < barW; ++b) fputc(b < filled ? '#' : ' ', stderr);
        fprintf(stderr, "] %3d%% (%d/%d)", pct, i + 1, g_renderTotal);
        fflush(stderr);

        advancePlayback(1.0 / (double)fps);   // step to the next frame's camera state
    }
    fprintf(stderr, "\n");

    // Flush + close ffmpeg (pclose blocks until encoding finishes and returns its status).
    fflush(g_ffPipe);
    int rc = pclose(g_ffPipe);
    g_ffPipe = nullptr;

    // Return rendering to the on-screen window: unbind the offscreen FBO and restore the
    // window's viewport + projection (reshape() reapplies both) so the interactive view
    // is intact after the render.
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    reshape(glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));

    // Restore interactive state.
    g_rendering = false;
    g_playing   = false;
    glutSetCursor(GLUT_CURSOR_INHERIT);
    g_lastRenderFps = fps;               // for the "* active" marker in the Render submenu
    g_menuDirty = true;                  // defer a menu rebuild to refresh that marker

    if (aborted)
        std::cerr << "Render aborted (ffmpeg pipe closed early)\n";
    else if (rc != 0)
        std::cerr << "Render: ffmpeg exited with status " << rc
                  << " (video may be incomplete)\n";
    std::cerr << "Render done: " << path << " (" << g_renderTotal << " frames, "
              << ((double)g_renderTotal / fps) << "s @ " << fps << " fps)\n";
}

// One-shot timer callback for the CLI -render path: by the time this fires (~250 ms
// after glutMainLoop starts) the window is realized, so glReadPixels in renderVideo()
// has a valid framebuffer.
void renderCliTimer(int /*value*/) {
    renderVideo(g_cliRenderFps);
    exit(0);
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
        if (state == GLUT_DOWN) {
            g_rightDown        = true;
            g_rightClickStartX = x; g_rightClickStartY = y;
            g_rightMaybeClick  = true;
            g_rightDragged     = false;
        } else {                      // release
            g_rightDown = false;
            // Stationary right-click (press + release, near-zero travel) focuses
            // the vertex under the cursor and sets the camera lookat; a right-
            // drag is left to zoom (see mouseMotion).
            if (g_rightMaybeClick && !g_rightDragged)
                handleVertexFocus(x, y);
            g_rightMaybeClick = false;
            g_rightDragged    = false;
        }
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
        // Only zoom once the cursor has moved past the click threshold — a
        // stationary right-click focuses a vertex, it shouldn't also zoom.
        if (!g_rightDragged) {
            int mvx = x - g_rightClickStartX; if (mvx < 0) mvx = -mvx;
            int mvy = y - g_rightClickStartY; if (mvy < 0) mvy = -mvy;
            if (mvx > g_clickThreshold || mvy > g_clickThreshold)
                g_rightDragged = true;
        }
        if (g_rightDragged) {
            g_zoom *= 1.0f - dy * 0.005f;
            g_zoom = fmaxf(0.1f, fminf(10.0f, g_zoom));
        }
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
        case 'l': case 'L':                            // reload the most recently loaded STL
            if (!g_stlPath.empty()) loadSTLFile(g_stlPath);
            else std::cout << "No STL loaded yet -- use Polyhedron > Load STL menu\n";
            break;
        case 'r': case 'R': remesh(FacetBox::SubdivisionMode::Centroid3); break; // remesh centroid3
        case 'm': case 'M': remesh(FacetBox::SubdivisionMode::Midpoint4);   break; // remesh midpoint4
        case 'g': case 'G': {                                              // toggle per-facet coloring
            g_colorMode = (g_colorMode + 1) % 2;
            std::cout << "Color mode: " << (g_colorMode == 0 ? "jet (rotated per remesh)"
                                                             : "grayscale 5-bucket criterion")
                      << "\n";
            break;
        }
        case 'k': case 'K': saveKeyframe(); break;      // save camera keyframe (timeline)
        case 'i': case 'I': insertKeyframeHere(); break; // insert keyframe after the Go-To cursor
        case 'p': case 'P': togglePlayback(); break;     // play/pause camera timeline (loop)
        case 't': case 'T': clearTimeline(); break;      // clear camera timeline
        case 8:           dropLastKeyframe(); break;     // BACKSPACE: remove last keyframe
        case 's': case 'S': nextShape(); break;   // cycle polyhedron (Cube/Dodec/TruncCube/TruncOct/TruncCuboct)
        case 'q': case 'Q': case 27: exit(0); break;   // 27 = ESC
    }
    glutPostRedisplay();
}

//-----------------------------------------------------------------------------
// UI menu callback
//-----------------------------------------------------------------------------
void MenuHandler(int choice) {
    // "Load STL..." submenu: ids 2000+i map to g_stlFiles[i]; id 3500 = Rescan STL
    // folder; id 3999 = "(no .stl)" placeholder (falls through every range below,
    // so it is ignored). Must be checked BEFORE the keyframe range, since the old
    // `choice >= 1000` guard would otherwise swallow STL ids (2000+) as keyframes.
    if (choice == 3500) {                       // Rescan STL folder
        g_menuDirty = true;                      // deferred rebuild re-scans g_stlDir
        glutPostRedisplay();
        return;
    }
    if (choice >= 2000 && choice < 3000) {       // STL submenu
        int idx = choice - 2000;
        if (idx >= 0 && idx < (int)g_stlFiles.size())
            loadSTLFile(g_stlDir + "/" + g_stlFiles[idx]);
        glutPostRedisplay();
        return;
    }
    // Dynamic "Go To Keyframe" submenu: one entry per keyframe, ids 1000+i. The
    // "(no keyframes)" placeholder uses id 9999 (idx huge -> ignored here).
    if (choice >= 1000 && choice < 2000) {
        int idx = choice - 1000;
        if (idx >= 0 && idx < (int)g_timeline.size()) goToKeyframe(idx);
        glutPostRedisplay();
        return;
    }
    // "Render Video..." submenu: ids 4000+fps (4008/4016/4024/4032/4064). Picking a
    // frame rate starts renderVideo() immediately (one blocking render pass).
    // g_menuDirty refreshes the "* active" marker afterward via the deferred rebuild.
    if (choice >= 4000 && choice < 4100) {
        int fps = choice - 4000;
        if (fps == 8 || fps == 16 || fps == 24 || fps == 32 || fps == 64)
            renderVideo(fps);
        glutPostRedisplay();
        return;
    }
    switch(choice) {
        case 1: g_showHelp = !g_showHelp; break;       // toggle HUD
        case 2: g_angleX=20; g_angleY=-30; g_zoom=1; g_panX=g_panY=0; g_focusVert=-1; break; // reset (also clear focus)
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
        case 16: remesh(FacetBox::SubdivisionMode::Centroid3); break; // remesh centroid3 (R)
        case 17: remesh(FacetBox::SubdivisionMode::Midpoint4);   break; // remesh midpoint4 (M)
        case 18: g_colorMode = (g_colorMode + 1) % 2;                 // toggle coloring (G)
            std::cout << "Color mode: " << (g_colorMode == 0 ? "jet (rotated per remesh)"
                                                             : "grayscale 5-bucket criterion")
                      << "\n";
            break;
        case 19: g_focusVert = -1;                 // clear focus -> pivot back to origin
            std::cout << "Camera lookat reset to origin\n";
            break;
        case 20: setShape(Shape::Cube);                   break;  // Polyhedron submenu
        case 21: setShape(Shape::Dodecahedron);           break;
        case 22: setShape(Shape::TruncatedCube);          break;
        case 23: setShape(Shape::TruncatedOctahedron);    break;
        case 24: setShape(Shape::TruncatedCuboctahedron); break;
        case 25: nextShape(); break;                      // Next Polyhedron (S)
        case 26: saveKeyframe();      break;               // Save camera keyframe (K)
        case 27: togglePlayback();    break;               // Play/Pause camera timeline (P)
        case 28: clearTimeline();     break;               // Clear camera timeline (T)
        case 29: dropLastKeyframe();  break;               // Remove last keyframe (Backspace)
        case 30: insertKeyframeHere(); break;             // Insert keyframe after Go-To cursor (I)
        case 31: jumpToLastKeyframe();   break;           // Go to last keyframe (continue adding)
        case 32: saveSession();          break;           // Save session (.h3geo)
        case 33: loadSession("session.h3geo"); break;     // Load session (.h3geo)
    }
    glutPostRedisplay();
}

//-----------------------------------------------------------------------------
// Create popup menu (middle-click; right-click is used to focus a vertex)
//-----------------------------------------------------------------------------
void createUI() {
    // Re-callable: destroy any previously-built menus first. The "Go To Keyframe"
    // submenu must reflect the live keyframe count, so createUI() is re-invoked
    // after every timeline mutation (save/insert/drop/clear/load). On the first call
    // from main() all ids are 0, so the destroy calls are skipped.
    if (g_menuMain) glutDestroyMenu(g_menuMain);
    if (g_menuGoto) glutDestroyMenu(g_menuGoto);
    if (g_menuPoly) glutDestroyMenu(g_menuPoly);
    if (g_menuStl) glutDestroyMenu(g_menuStl);
    if (g_menuRender) glutDestroyMenu(g_menuRender);

    // --- Load STL submenu (built first so it can nest inside g_menuPoly) ---
    // IDs 2000+i map to g_stlFiles[i] in MenuHandler; 3500 = Rescan, 3999 = the
    // "(no .stl)" placeholder. The active file is prefixed "* ".
    scanSTLFolder();
    g_menuStl = glutCreateMenu(MenuHandler);
    if (g_stlFiles.empty()) {
        glutAddMenuEntry("(no .stl in ~/Downloads)", 3999);
    } else {
        glutAddMenuEntry("Rescan STL folder", 3500);
        for (int i = 0; i < (int)g_stlFiles.size(); ++i) {
            const std::string& fn = g_stlFiles[i];
            bool active = (g_shape == Shape::STL) && (g_stlPath == g_stlDir + "/" + fn);
            std::string lbl = active ? ("* " + fn) : fn;
            glutAddMenuEntry(lbl.c_str(), 2000 + i);
        }
    }

    // --- Polyhedron submenu (direct shape selection; IDs 20-24) + nested Load STL ---
    g_menuPoly = glutCreateMenu(MenuHandler);
    glutAddMenuEntry("Cube",                     20);
    glutAddMenuEntry("Dodecahedron",             21);
    glutAddMenuEntry("Truncated Cube",           22);
    glutAddMenuEntry("Truncated Octahedron",     23);
    glutAddMenuEntry("Truncated Cuboctahedron",  24);
    glutAddSubMenu("Load STL...", g_menuStl);    // pick any .stl found in ~/Downloads

    // --- Go To Keyframe submenu (dynamic; one entry per saved keyframe) ---
    // IDs 1000+i map to keyframe i in MenuHandler. Empty -> a placeholder entry
    // (id 9999) so the submenu is never blank.
    g_menuGoto = glutCreateMenu(MenuHandler);
    if (g_timeline.empty()) {
        glutAddMenuEntry("(no keyframes - press K)", 9999);
    } else {
        for (int i = 0; i < (int)g_timeline.size(); ++i) {
            const CamKey& k = g_timeline[i];
            char lbl[96];
            snprintf(lbl, sizeof(lbl),
                     "Keyframe %d   (yaw %.0f  pitch %.0f  zoom %.2f)",
                     i + 1, (double)k.angleY, (double)k.angleX, (double)k.zoom);
            glutAddMenuEntry(lbl, 1000 + i);
        }
    }

    // --- Render Video submenu (pick a frame rate -> starts the render) ---
    // IDs 4008/4016/4024/4032/4064 = fps + 4000; MenuHandler maps the [4000,4100)
    // range to renderVideo(fps). The last-used fps is prefixed "* ".
    g_menuRender = glutCreateMenu(MenuHandler);
    {
        const int rates[] = {8, 16, 24, 32, 64};
        char lbl[32];
        for (int r : rates) {
            bool active = (g_lastRenderFps == r);
            snprintf(lbl, sizeof(lbl), "%s%d fps", active ? "* " : "", r);
            glutAddMenuEntry(lbl, 4000 + r);
        }
    }

    // --- Main menu ---
    g_menuMain = glutCreateMenu(MenuHandler);
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
    glutAddMenuEntry("Remesh: Centroid3 (R)", 16);
    glutAddMenuEntry("Remesh: Midpoint4 (M)",  17);
    glutAddMenuEntry("Toggle Coloring: Jet / Grey (G)", 18);
    glutAddMenuEntry("Reset Camera LookAt to Origin", 19);
    glutAddMenuEntry("Next Polyhedron (S)", 25);     // mirrors the S keyboard shortcut
    glutAddMenuEntry("Save Camera Keyframe (K)", 26);
    glutAddMenuEntry("Insert Keyframe Here (I)", 30);
    glutAddMenuEntry("Go To Last Keyframe", 31);
    glutAddMenuEntry("Play/Pause Camera Timeline (P)", 27);
    glutAddMenuEntry("Clear Camera Timeline (T)", 28);
    glutAddMenuEntry("Remove Last Keyframe (Backspace)", 29);
    glutAddMenuEntry("Save Session (.h3geo)", 32);
    glutAddMenuEntry("Load Session (.h3geo)", 33);
    glutAddSubMenu("Render Video...", g_menuRender); // pick fps -> render timeline to MP4
    glutAddSubMenu("Polyhedron", g_menuPoly);        // direct-selection submenu
    glutAddSubMenu("Go To Keyframe", g_menuGoto);   // one entry per saved keyframe
    glutAttachMenu(GLUT_MIDDLE_BUTTON);   // right button is reserved for vertex-focus
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

