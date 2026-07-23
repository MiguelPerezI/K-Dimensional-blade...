#include <vector>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cstring>      // strlen() for the GLSL shader info-log paths
// Enable GL extension function prototypes (glCreateShader, glShaderSource,
// glCompileShader, glUseProgram, glUniform*, ...) from <GL/glext.h> so we can
// compile/use GLSL shaders without a loader like GLEW. Must be defined before
// <GL/gl.h> is included (it is, transitively, via <GL/glut.h>).
#define GL_GLEXT_PROTOTYPES 1
#include <GL/glut.h>
#include "Vector3D.cpp"
#include "Vector4D.cpp"
#include "Quaternion.cpp"
#include "Facet.cpp"
#include "FacetBox.hpp"
#include "Dodecahedron.hpp"
#include "Cube.hpp"
#include <ctime>
#include <cstdlib>
#include <sys/stat.h>   // mkdir()
#include <fstream>      // std::ofstream for the snapshot writer
#include <cstdint>      // uint32_t/uint16_t for binary STL

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
static float g_fisheyeStrength = 0.5f;  // 0.0 = normal, 2.0 = extreme fisheye

// Surface proximity and collision detection
static float g_minDistanceToSurface = 0.00000001f;
static float g_surfaceZoomFactor = 1.0f;
static bool g_enableSurfaceZoom = true;

// UI state for plane selection
int g_currentOrientation = 0;  // 0=XY, 1=YZ, 2=XZ
int g_currentLayer = 4;        // Current layer within the selected orientation
int g_maxLayers = 8;           // Maximum layers available (matches cube2/cube3 subdivision)

// STL snapshot counter — each press of 'c' writes a uniquely-named file to ~/Downloads
int g_snapshotCount = 0;

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
int N = 12 ;                                   // cuboctahedron lattice: one cell per subcell (keep modest)
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

// ii, jj, kk : integers that represent the module action on the subecube lattice for cube.
int ii = 4;
int jj = 5;
int kk = 9;

//==============================================================================
// Realtime inversion animation — CPU port of main17_gpu's updateAnimatedGeometry.
//   pass 1 center = Vector3D(0, 0, 0.1*sin(t)), radius = 0.5   (oscillates in z)
//   pass 2 center/radius fixed at subcell (1,1,1) — as in main15 Setup()
// The cube's PRISTINE identity lattice is captured once; each frame the composed
// sigma is recomputed from it (never in-place) so frame-over-frame compounding is
// avoided. The CPU render path reads the Cube's subcell vertices directly
// (getCheckerboardFacets), so updateAnimatedGeometry() just writes the deformed
// vertices back into the subcells — no per-frame retriangulation. Triangulation
// is refreshed only on an STL snapshot ('c') so writeSTL() sees the current frame.
//==============================================================================
struct SigmaParams { Vector3D c1; double r1; Vector3D c2; double r2; Vector3D c3; double r3; };

static std::vector<float> g_identityPositions;  // pristine identity lattice (captured once)
static double g_animTime   = 0.0;               // animation clock (advanced by the timer)
static bool   g_animPaused = false;
static const double g_timeStep = 0.02;           // animTime advance per ~33ms tick (main17_gpu)

// --- Slow sequential inversion (new feature) ---
// The three inversions ramp up one at a time via per-pass blends alpha in [0,1],
// then ping-pong back (undo in reverse: 3->2->1). g_invProgress in [0,6) is the
// phase clock (6 ramps per cycle). g_invSpeed is a RATE: larger = faster.
// Per-tick advance = INV_RATE_PER_TICK * g_invSpeed; at ~30 Hz that's
// ~0.3*g_invSpeed/sec, so at speed=1.0 one inversion ramp (~1.0 unit) takes
// ~3.3 s, full cycle ~20 s. Default 0.5 => ~6.7 s per inversion (visibly slow).
static bool   g_slowMode      = true;    // true=slow sequential, false=original simultaneous
static double g_invSpeed      = 0.5;     // inversion advance rate (larger=faster); CLI/keyboard adjustable
static double g_invProgress   = 0.0;    // slow-mode phase clock in [0,6)
static const double INV_RATE_PER_TICK = 0.01;
static const double g_invSpeedStep    = 0.1;  // keyboard step for s/S
struct InvBlend { double a1, a2, a3; };

// --- Per-frame STL capture for one slow-mode cycle (new) ---
// --capture-cycle: write one STL per animation frame for a full 6-phase cycle,
// then exit. --outdir <dir>: where to write frame_XXXXXX.stl (default
// ~/Downloads/cycle_<timestamp>/). Captures the SAME checkerboard selection the
// render shows (cube.getCheckerboardFacetsCuboctahedronLattice(ii,kk,jj, g_hollow, g_inset)), matching the 'c' key.
static bool        g_captureCycle = false;  // --capture-cycle flag
static std::string g_captureDir;           // --outdir target directory
static bool        g_capturing    = false; // active capture state (true while capturing)
static int         g_frameIndex   = 0;    // current frame number (0-based)
static double      g_prevProgress = -1.0;  // progress at last captured frame (-1 = none yet)
static bool        g_animStepped  = false; // set by animTimer when it advanced, consumed in display
static bool        g_binarySTL    = false; // --binary: write binary STL (else ASCII)
static bool        g_hollow       = false; // --hollow / [O]: frame holes per face (inset + skip inner face)
static double      g_morphS       = 0.5;   // --morph: TO morph s (0.5=regular TO, 1.0=cuboctahedron, <0.5 stretched)
static double      g_inset        = 0.5;   // --inset: hollow border inset (0,1); bigger = bigger hole / thinner border

// --- Mouse-driven single inversion at the center (replaces the slow 3-pass) ---
static Vector3D g_invCenter(0.0, 1e-9, 0.0);   // inversion center, mouse-controlled; desk = XZ plane (y~0)
static double   g_invRadius = 0.5;             // inversion sphere radius; r/R to adjust
static double   g_invDragGain = 0.005;         // world units per pixel for left-drag on the XZ desk (lower = less sensitive); scaled by g_zoom
static bool     g_geomDirty = true;            // recompute the lattice when center/radius changes (true => first frame)
static unsigned g_geomVersion = 0;            // bumped each time the lattice is deformed (cache key for the batched mesh)
static float    g_lastDrawMs = 0.0f;          // last frame's draw work (ms), shown in the HUD
// Cached camera matrices for mouse->world unprojection (refreshed each display()).
static GLdouble g_modelview[16]   = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
static GLdouble g_projection[16]  = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
static GLint    g_viewport[4]     = {0,0,720,720};

// Same animated inversion parameters as main17_gpu::currentSigmaParams(t).
static SigmaParams currentSigmaParams(double t) {
    SigmaParams s;
    s.c1 = Vector3D(0.1*sin(2*t), 0.0, 0.1 * sin(t) * cos(t));   // inversion point oscillates in z
    s.r1 = 0.5;
    s.c2 = cube.getSubcellCenter(1, 1, 1);
    s.r2 = cube.getSubcellRadius(1, 1, 1);
    s.c3 = cube.getSubcellCenter(-1, -1, -1);
    s.r3 = cube.getSubcellRadius(-1, -1, -1);
    return s;
}
// Map the slow-mode phase clock g_invProgress in [0,6) to per-pass blends
// (a1,a2,a3). Forward: inv1 up, inv2 up, inv3 up (reaches full composition at
// p=3). Reverse undo: inv3 down, inv2 down, inv1 down (back to identity at
// p=6 ~ 0). Continuous at every boundary; each a is in [0,1] by construction.
static InvBlend currentAlphas(double progress) {
    progress = fmod(progress, 6.0);  if (progress < 0.0) progress += 6.0;
    if      (progress < 1.0) return {progress,        0.0,       0.0};       // inv1 up
    else if (progress < 2.0) return {1.0, progress - 1.0,       0.0};       // inv2 up
    else if (progress < 3.0) return {1.0, 1.0,        progress - 2.0};      // inv3 up
    else if (progress < 4.0) return {1.0, 1.0,        4.0 - progress};      // undo inv3
    else if (progress < 5.0) return {1.0, 5.0 - progress, 0.0};             // undo inv2
    else                     return {6.0 - progress, 0.0, 0.0};            // undo inv1
}
static const char* phaseName(double progress) {
    progress = fmod(progress, 6.0);  if (progress < 0.0) progress += 6.0;
    if (progress < 1.0) return "inv1->";
    if (progress < 2.0) return "inv2->";
    if (progress < 3.0) return "inv3->";
    if (progress < 4.0) return "<-inv3";
    if (progress < 5.0) return "<-inv2";
    return "<-inv1";
}

// Three-pass sphere inversion with a singular-center guard, generalized with a
// per-pass blend alpha in [0,1]: out_i = (1-a_i)*in_i + a_i*sigma_i(in_i).
// (The free sigma() in Vector3D.cpp throws on a singular center; here we leave
// the point untouched instead, so the animation never crashes.) At a=(1,1,1)
// this is identical to the original all-at-once composed inversion.
static inline Vector3D composeSigmaBlended(const Vector3D& p, const SigmaParams& s, const InvBlend& b) {
    // pass 1: out = (1-a1)*p + a1*sigma1(p)
    Vector3D d1 = p - s.c1;  double ds1 = d1 * d1;
    Vector3D s1v = (ds1 < 1e-12) ? p : s.c1 + (s.r1 * s.r1 / ds1) * d1;
    Vector3D q  = p + b.a1 * (s1v - p);
    // pass 2: out = (1-a2)*q + a2*sigma2(q)
    Vector3D d2 = q - s.c2;  double ds2 = d2 * d2;
    Vector3D s2v = (ds2 < 1e-12) ? q : s.c2 + (s.r2 * s.r2 / ds2) * d2;
    Vector3D r2v = q + b.a2 * (s2v - q);
    // pass 3: out = (1-a3)*r2v + a3*sigma3(r2v)
    Vector3D d3 = r2v - s.c3; double ds3 = d3 * d3;
    Vector3D s3v = (ds3 < 1e-12) ? r2v : s.c3 + (s.r3 * s.r3 / ds3) * d3;
    return r2v + b.a3 * (s3v - r2v);
}
// Original behavior: all three passes at full strength (slow mode disabled).
static inline Vector3D composeSigma(const Vector3D& p, const SigmaParams& s) {
    return composeSigmaBlended(p, s, {1.0, 1.0, 1.0});
}

// Single sphere inversion at the mouse-driven g_invCenter, regularized so a vertex
// almost on the center cannot "explode" to infinity (a subcell vertex ~1e-4 from
// the center would otherwise jump ~1e4 units and drag its facets into a spike).
// Clamp the effective squared distance to minDist^2 = (r^2 / g_maxInvDist)^2 so NO
// deformed point lands farther than g_maxInvDist from the center. Points beyond
// minDist invert exactly as before; only the near-singular core is tamed (it
// stays near the center). Continuous with the normal inversion at the boundary.
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
    cube.fillVertexLattice(g_identityPositions);  // n^3*8*3 floats — same call as main17_gpu
    cout << "[anim] identity lattice captured: " << (g_identityPositions.size() / 3)
         << " verts (" << (g_identityPositions.size() * sizeof(float) / (1024.0 * 1024.0))
         << " MB)\n";
}

// Push the single mouse-driven inversion of the identity lattice back into the
// Cube's subcells (same iteration order fillVertexLattice() uses: physical
// [0,n)^3, vertex 0..7). The lattice holds ONE full-strength sphere inversion at
// g_invCenter; recomputed only when g_geomDirty is set (mouse moved the center or
// r/R changed the radius).
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
                    cube.updateSubCellVertex(lx, ly, lz, l, sigmaCenter(p));
                    idx += 3;
                }
            }
    ++g_geomVersion;   // invalidate the batched-mesh cache (positions changed)
}

// ~30 Hz redraw clock. The slow-transition phase clock is removed: g_invProgress
// is frozen and g_animStepped is never set, so the --capture-cycle recorder never
// triggers (kept inert on purpose). The lattice is now driven by the mouse
// (g_invCenter/g_invRadius) and is recomputed on demand in display() when
// g_geomDirty is set. [Space] still toggles g_animPaused (inert HUD indicator).
static void animTimer(int /*value*/) {
    glutPostRedisplay();
    glutTimerFunc(33, animTimer, 0);
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
        cout << "modules <ii, jj, kk> = <" << ii << ", " << jj << ", " << kk << ">\n";
        // ==================== ADVANCED SUBDIVISION DEMOS - Reflection of the subdivision in Spheres ====================
        // Capture the pristine identity lattice once, then deform to the t=0 state.
        // The animation (animTimer advances g_animTime; display() deforms each frame)
        // recomputes this every frame; at t=0 it equals the two baked reflections
        // main15 used to apply here.
        // Bound the deformed lattice to within the cube radius of the inversion
        // center so a vertex almost on the center can't explode (see sigmaCenter).
        g_maxInvDist = cube_dim;

        captureIdentityLattice();
        updateAnimatedGeometry();
        cube.writeSTL_s_cuboctahedron_lattice_checkerboard("/home/mike666/Downloads/mesh_output_cuboctahedron_lattice_module.stl", "CuboctahedronLattice", ii, kk, jj, g_hollow, g_inset);
    }
}

///////////////////     DRAW       ///////////////////////
///////////////////////////////////////////////////////////////////////////////
// Gouraud shader (per-vertex ambient + diffuse + specular). Ported from
// matthewachan/opengl-shaders (gouraud_vertex.vs / gouraud_fragment.fs), retargeted
// to GLSL #version 130 so it runs on any GL 3.0+ compatibility context (Mesa's
// default) without requesting a specific context version. Scoped to the lattice
// FILLS only -- outlines / sphere / glass / axes / HUD stay fixed-function under
// glUseProgram(0). Toggle with [L]; if the context can't compile the GLSL the toggle
// is inert and a warning is printed (fixed-function stays usable). No perf gain on
// llvmpipe -- this is a visual/quality change.
//
// Why this slots in with almost no plumbing:
//  * The module already draws fills via glVertexPointer (conventional attrib 0) +
//    glNormalPointer (conventional attrib 2). We glBindAttribLocation position->0 and
//    vertexNormal->2, so the existing client arrays feed the shader unchanged.
//  * g_modelview / g_projection are already cached each frame (see display()); we
//    feed them straight to the modelViewMatrix / projectionMatrix uniforms.
//  * The modelview is rotation+translation only, so the shader's
//    normalize(modelViewMatrix * normal) is correct (matches GL_NORMALIZE).
//  * GL_LIGHT0 is directional with a constant eye-space direction (set at identity in
//    initGL), which maps 1:1 to directionalLight.
///////////////////////////////////////////////////////////////////////////////
static bool   g_useShader    = false;   // [L] toggle: Gouraud shading on the lattice fills
static GLuint g_gouraudProg  = 0;
static bool   g_shaderOK     = false;   // compiled+linked successfully?

// Uniform locations (queried once at link time; -1 = not found / optimized out).
static GLint g_u_mv = -1, g_u_proj = -1, g_u_ambient = -1, g_u_specPow = -1;
static GLint g_u_matColour = -1, g_u_matUseColour = -1, g_u_matReflect = -1;
static GLint g_u_dirColour = -1, g_u_dirDir = -1, g_u_dirIntensity = -1;
static GLint g_u_ptColour = -1, g_u_ptPos = -1, g_u_ptIntensity = -1;
static GLint g_u_ptAttC = -1, g_u_ptAttL = -1, g_u_ptAttE = -1;

// GLSL #version 130: no layout(location=) (bound via glBindAttribLocation before
// link), in/out (not attribute/varying), custom fragment output. Works on GL 3.0+.
static const char* kGouraudVS =
"#version 130\n"
"in vec3 position;\n"
"in vec3 vertexNormal;\n"
"out vec4 color;\n"
"uniform mat4 modelViewMatrix;\n"
"uniform mat4 projectionMatrix;\n"
"uniform vec3 ambientLight;\n"
"uniform float specularPower;\n"
"struct Material { vec3 colour; int useColour; float reflectance; };\n"
"struct Attenuation { float constant; float linear; float exponent; };\n"
"struct PointLight { vec3 colour; vec3 position; float intensity; Attenuation att; };\n"
"struct DirectionalLight { vec3 colour; vec3 direction; float intensity; };\n"
"uniform Material material;\n"
"uniform PointLight pointLight;\n"
"uniform DirectionalLight directionalLight;\n"
"vec4 calcLightColour(vec3 lc, float li, vec3 p, vec3 toLight, vec3 n){\n"
"  float diff = max(dot(n, toLight), 0.0);\n"
"  vec4 dc = vec4(lc, 1.0) * li * diff;\n"
"  vec3 cd = normalize(-p);\n"
"  vec3 refl = normalize(reflect(-toLight, n));\n"
"  float sf = pow(max(dot(cd, refl), 0.0), specularPower);\n"
"  vec4 sc = li * sf * material.reflectance * vec4(lc, 1.0);\n"
"  return dc + sc;\n"
"}\n"
"vec4 calcPointLight(PointLight L, vec3 p, vec3 n){\n"
"  vec3 d = L.position - p;\n"
"  vec4 c = calcLightColour(L.colour, L.intensity, p, normalize(d), n);\n"
"  float dist = length(d);\n"
"  float attInv = L.att.constant + L.att.linear * dist + L.att.exponent * dist * dist;\n"
"  return c / attInv;\n"
"}\n"
"void main(){\n"
"  vec4 mv = modelViewMatrix * vec4(position, 1.0);\n"
"  gl_Position = projectionMatrix * mv;\n"
"  vec3 mn = normalize(modelViewMatrix * vec4(vertexNormal, 0.0)).xyz;\n"
"  color = vec4(ambientLight, 1.0);\n"
"  color += calcLightColour(directionalLight.colour, directionalLight.intensity, mv.xyz, normalize(directionalLight.direction), mn);\n"
"  color += calcPointLight(pointLight, mv.xyz, mn);\n"
"}\n";

static const char* kGouraudFS =
"#version 130\n"
"struct Material { vec3 colour; int useColour; float reflectance; };\n"
"in vec4 color;\n"
"out vec4 fragColor;\n"
"uniform Material material;\n"
"void main(){ fragColor = color * vec4(material.colour, 1.0); }\n";

static GLuint gouraudCompileShader(GLenum type, const char* src) {
    GLuint s = glCreateShader(type);
    GLint len = (GLint)strlen(src);
    glShaderSource(s, 1, &src, &len);
    glCompileShader(s);
    GLint ok = GL_FALSE;
    glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
    if (ok != GL_TRUE) {
        GLint n = 0; glGetShaderiv(s, GL_INFO_LOG_LENGTH, &n);
        std::vector<char> log((size_t)(n > 1 ? n : 1));
        glGetShaderInfoLog(s, (GLsizei)log.size(), nullptr, log.data());
        printf("Gouraud %s shader compile failed:\n%s\n",
               type == GL_VERTEX_SHADER ? "vertex" : "fragment", log.data());
        glDeleteShader(s);
        return 0;
    }
    return s;
}

static void initGouraudShader() {
    GLuint vs = gouraudCompileShader(GL_VERTEX_SHADER, kGouraudVS);
    if (!vs) { g_shaderOK = false; return; }
    GLuint fs = gouraudCompileShader(GL_FRAGMENT_SHADER, kGouraudFS);
    if (!fs) { glDeleteShader(vs); g_shaderOK = false; return; }

    GLuint prog = glCreateProgram();
    // Bind the shader inputs to the conventional vertex/normal slots so the existing
    // glVertexPointer (attrib 0) and glNormalPointer (attrib 2) client arrays feed
    // the shader with no extra setup. Must be done before linking.
    glBindAttribLocation(prog, 0, "position");
    glBindAttribLocation(prog, 2, "vertexNormal");
    glAttachShader(prog, vs);
    glAttachShader(prog, fs);
    glLinkProgram(prog);
    glDeleteShader(vs);
    glDeleteShader(fs);

    GLint ok = GL_FALSE;
    glGetProgramiv(prog, GL_LINK_STATUS, &ok);
    if (ok != GL_TRUE) {
        GLint n = 0; glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &n);
        std::vector<char> log((size_t)(n > 1 ? n : 1));
        glGetProgramInfoLog(prog, (GLsizei)log.size(), nullptr, log.data());
        printf("Gouraud program link failed:\n%s\n", log.data());
        glDeleteProgram(prog);
        g_shaderOK = false;
        return;
    }
    g_gouraudProg = prog;
    g_shaderOK = true;

    g_u_mv           = glGetUniformLocation(prog, "modelViewMatrix");
    g_u_proj         = glGetUniformLocation(prog, "projectionMatrix");
    g_u_ambient      = glGetUniformLocation(prog, "ambientLight");
    g_u_specPow      = glGetUniformLocation(prog, "specularPower");
    g_u_matColour    = glGetUniformLocation(prog, "material.colour");
    g_u_matUseColour = glGetUniformLocation(prog, "material.useColour");
    g_u_matReflect   = glGetUniformLocation(prog, "material.reflectance");
    g_u_dirColour    = glGetUniformLocation(prog, "directionalLight.colour");
    g_u_dirDir       = glGetUniformLocation(prog, "directionalLight.direction");
    g_u_dirIntensity = glGetUniformLocation(prog, "directionalLight.intensity");
    g_u_ptColour     = glGetUniformLocation(prog, "pointLight.colour");
    g_u_ptPos        = glGetUniformLocation(prog, "pointLight.position");
    g_u_ptIntensity  = glGetUniformLocation(prog, "pointLight.intensity");
    g_u_ptAttC       = glGetUniformLocation(prog, "pointLight.att.constant");
    g_u_ptAttL       = glGetUniformLocation(prog, "pointLight.att.linear");
    g_u_ptAttE       = glGetUniformLocation(prog, "pointLight.att.exponent");
    printf("Gouraud shader loaded (toggle with [L]).\n");
}

// Set the Gouraud uniforms for this frame from the cached camera matrices and the
// fixed-function light/material values (matches initGL's GL_LIGHT0). Call with
// g_gouraudProg already active (glUseProgram), after the camera transforms so
// g_modelview / g_projection are fresh.
static void setupGouraudUniforms() {
    if (!g_shaderOK || g_gouraudProg == 0) return;

    // GL stores matrices column-major; glUniformMatrix4fv(GL_FALSE) reads column-major.
    GLfloat mv[16], pr[16];
    for (int i = 0; i < 16; ++i) { mv[i] = (GLfloat)g_modelview[i]; pr[i] = (GLfloat)g_projection[i]; }
    if (g_u_mv   != -1) glUniformMatrix4fv(g_u_mv,   1, GL_FALSE, mv);
    if (g_u_proj != -1) glUniformMatrix4fv(g_u_proj, 1, GL_FALSE, pr);

    // GL_LIGHT0 is directional: position {1,1,0.25,0} set at identity modelview in
    // initGL, so its eye-space direction is the CONSTANT {1,1,0.25} (it tracks the
    // camera). The shader lights in view space, so pass it directly (no per-frame
    // transform) -- this matches the fixed-function look.
    GLfloat dx = 1.0f, dy = 1.0f, dz = 0.25f;
    GLfloat invlen = 1.0f / sqrtf(dx*dx + dy*dy + dz*dz);
    dx *= invlen; dy *= invlen; dz *= invlen;
    if (g_u_ambient      != -1) glUniform3f(g_u_ambient,      0.3f, 0.3f, 0.3f);  // GL_AMBIENT
    if (g_u_specPow      != -1) glUniform1f(g_u_specPow,      128.0f);            // GL_SHININESS
    if (g_u_dirColour    != -1) glUniform3f(g_u_dirColour,    0.7f, 0.7f, 0.7f);  // GL_DIFFUSE
    if (g_u_dirDir       != -1) glUniform3f(g_u_dirDir,       dx, dy, dz);
    if (g_u_dirIntensity != -1) glUniform1f(g_u_dirIntensity, 1.0f);

    // Material: GL_COLOR_MATERIAL(AMBIENT_AND_DIFFUSE) + glColor3ub(200,200,200).
    const GLfloat c = 200.0f / 255.0f;
    if (g_u_matColour    != -1) glUniform3f(g_u_matColour,    c, c, c);
    if (g_u_matUseColour != -1) glUniform1i(g_u_matUseColour, 1);
    if (g_u_matReflect   != -1) glUniform1f(g_u_matReflect,   1.0f);              // GL_SPECULAR {1,1,1}

    // No point light in the scene: intensity 0 -> its term is 0 (att.constant=1 avoids
    // a divide-by-zero even though the result is 0).
    if (g_u_ptColour     != -1) glUniform3f(g_u_ptColour,     0.0f, 0.0f, 0.0f);
    if (g_u_ptPos        != -1) glUniform3f(g_u_ptPos,        0.0f, 0.0f, 0.0f);
    if (g_u_ptIntensity  != -1) glUniform1f(g_u_ptIntensity,  0.0f);
    if (g_u_ptAttC       != -1) glUniform1f(g_u_ptAttC, 1.0f);
    if (g_u_ptAttL       != -1) glUniform1f(g_u_ptAttL, 0.0f);
    if (g_u_ptAttE       != -1) glUniform1f(g_u_ptAttE, 0.0f);
}

// Background palette ([B] cycles). Light entries keep the black facet outlines
// visible; the dark ones give a dramatic lit look with the Gouraud shader on
// (outlines will be near-invisible on dark -- switch back to a light one for structure).
static const GLfloat g_bgPalette[][3] = {
    {0.86f, 0.89f, 0.93f},   // 0: light steel blue (default)
    {0.95f, 0.93f, 0.88f},   // 1: warm off-white
    {0.82f, 0.83f, 0.85f},   // 2: pale neutral gray
    {0.10f, 0.12f, 0.18f},   // 3: deep navy
    {0.18f, 0.20f, 0.23f},   // 4: charcoal
};
static const int g_bgPaletteCount = (int)(sizeof(g_bgPalette) / sizeof(g_bgPalette[0]));
static int g_bgIndex = 0;

void Draw() {
    
    extern void Setup(); // assume you define this elsewhere
    static int ciclo = 1;  // or however you manage visibility
	if (ciclo > 0) {
        /*Draw here with OpenGL*/
        // In your draw‐all loop:
        // Batched render: cache vertex/normal arrays, rebuild only when the lattice
        // knobs (inset/hollow) or geometry (g_geomVersion) change — camera-only
        // frames reuse them. One glDrawArrays for the fills + one for the outlines.
        static std::vector<float> s_verts, s_norms;
        static size_t   s_N = 0;
        static unsigned s_geom = 0xFFFFFFFFu;
        static double   s_inset = -1.0;
        static bool     s_hollow = false;
        static int      s_ii = -1, s_jj = -1, s_kk = -1;   // modular sublattice selection (i/j/k keys)
        if (s_inset != g_inset || s_hollow != g_hollow ||
            s_geom != g_geomVersion ||
            s_ii != ii || s_jj != jj || s_kk != kk) {
            FacetBox fb = cube.getCheckerboardFacetsCuboctahedronLattice(ii, kk, jj, g_hollow, g_inset);
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
            s_inset=g_inset; s_hollow=g_hollow; s_geom=g_geomVersion; s_ii=ii; s_jj=jj; s_kk=kk;
            plane_subcells = fb;  // keep the global in sync for any other consumer
        }
        if (s_N > 0) {
            glPushAttrib(GL_COLOR_BUFFER_BIT | GL_POLYGON_BIT | GL_ENABLE_BIT);
            glEnableClientState(GL_VERTEX_ARRAY);
            glEnableClientState(GL_NORMAL_ARRAY);
            glVertexPointer(3, GL_FLOAT, 0, s_verts.data());
            glNormalPointer(GL_FLOAT, 0, s_norms.data());
            // Filled triangles (polygon offset so outlines don't z-fight). With the
            // Gouraud shader ON, lighting is per-vertex in the shader (it reads attrib0=
            // position and attrib2=normal from the client arrays already enabled above),
            // so glColor3ub is skipped. With the shader OFF, fixed-function GL_LIGHTING +
            // GL_COLOR_MATERIAL uses glColor3ub as the ambient+diffuse material colour.
            glEnable(GL_POLYGON_OFFSET_FILL);
            glPolygonOffset(1.0f, 1.0f);
            const bool useShader = (g_useShader && g_shaderOK && g_gouraudProg != 0);
            if (useShader) {
                glUseProgram(g_gouraudProg);
                setupGouraudUniforms();
            } else {
                glColor3ub(200, 200, 200);
            }
            glDrawArrays(GL_TRIANGLES, 0, (GLsizei)(s_N * 3));
            if (useShader) glUseProgram(0);   // outlines + the rest stay fixed-function
            glDisable(GL_POLYGON_OFFSET_FILL);
            // Black outlines: every facet edge, one draw call via line polygon mode.
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
static float g_angleX = 80.0f, g_angleY = 0.0f; // top-down desk view: look down -Y at the XZ plane (orbit with middle-drag)
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
// Parse slow-inversion options: --speed <v>, --mode slow|original, --slow,
// --original, or a bare positional number (=> speed). Defaults: slow, 0.5.
// Called after glutInit so GLUT has already stripped its own flags from argv.
// Create a directory path (like mkdir -p); ignores already-exists errors.
static void ensureDir(const std::string& path) {
    std::string cur;
    for (char c : path) {
        cur += c;
        if (c == '/') mkdir(cur.c_str(), 0755);   // create each prefix ending in '/'
    }
    mkdir(path.c_str(), 0755);                     // final component
}

static void parseArgs(int argc, char** argv) {
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if      (a == "--speed" && i + 1 < argc)  g_invSpeed = atof(argv[++i]);
        else if (a == "--mode"  && i + 1 < argc)  g_slowMode = (std::string(argv[++i]) != "original");
        else if (a == "--slow")                   g_slowMode = true;
        else if (a == "--original")               g_slowMode = false;
        else if (a == "--capture-cycle" || a == "--capture")  g_captureCycle = true;
        else if (a == "--outdir" && i + 1 < argc)             g_captureDir = argv[++i];
        else if (a == "--binary" || a == "--binary-stl")      g_binarySTL = true;
        else if (a == "--hollow")                             g_hollow = true;   // frame holes per face (hollow)
        else if (a == "--morph" && i + 1 < argc)             g_morphS = atof(argv[++i]);   // TO morph s (0.5=regular TO)
        else if (a == "--inset" && i + 1 < argc)              g_inset  = atof(argv[++i]);   // hollow border inset (0,1)
        else { char* end = nullptr; double v = strtod(argv[i], &end);  // positional => speed
               if (end != argv[i] && *end == '\0') g_invSpeed = v; }
    }
    if (g_invSpeed < 1e-4) g_invSpeed = 1e-4;   // never stall
    if (g_invSpeed > 10.0) g_invSpeed = 10.0;
    if (g_morphS < 0.05)  g_morphS = 0.05;       // TO morph s in (0, 1]
    if (g_morphS > 1.0)   g_morphS = 1.0;
    if (g_inset  < 0.05)  g_inset  = 0.05;       // hollow inset in (0, 1)
    if (g_inset  > 0.95)  g_inset  = 0.95;

    // --capture-cycle: per-frame STL dump for one slow-mode cycle, then exit.
    if (g_captureCycle) {
        if (!g_slowMode) printf("[capture] --capture-cycle needs slow mode; enabling SLOW.\n");
        g_slowMode = true;                                  // the cycle is a slow-mode concept
        if (g_captureDir.empty()) {                          // default: grouped, timestamped subfolder
            const char* home = getenv("HOME");
            std::string base = home ? (std::string(home) + "/Downloads") : "/home/mike666/Downloads";
            time_t now = time(nullptr);
            struct tm* lt = localtime(&now);
            char ts[32];
            strftime(ts, sizeof(ts), "%Y%m%d_%H%M%S", lt);
            g_captureDir = base + "/cycle_truncated_octahedron_" + ts;
        }
        ensureDir(g_captureDir);
        g_capturing = true;                                  // start capturing on the first display
        long expect = (long)(6.0 / (INV_RATE_PER_TICK * g_invSpeed) + 0.5);
        printf("[capture] per-frame %s STL for one cycle -> %s\n",
               g_binarySTL ? "BINARY" : "ASCII", g_captureDir.c_str());
        printf("[capture] ~%ld frames expected (speed=%.3f); program exits when the cycle wraps.\n",
               expect, g_invSpeed);
    }

    printf("[inv] single inversion at mouse center; r=%.3f (r/R to adjust)  hollow=%s  s=%.2f  inset=%.2f\n",
           g_invRadius, g_hollow ? "ON" : "OFF", g_morphS, g_inset);
}

int main(int argc, char** argv)
{
    srand((unsigned)time(nullptr));
    glutInit(&argc, argv);
    parseArgs(argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(720, 720);
    glutCreateWindow(" JAZ 4D Enhanced Camera U.U  (Cuboctahedron Lattice Module)");

    // Enable smoothing & blending by default
    ProcessMenu(1);
    initGL();
    createUI();
    
    // Objects setup (captures the identity lattice, deforms to t=0)
    Setup();

    // Register callbacks
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouseButton);
    glutMotionFunc(mouseMotion);
    glutKeyboardFunc(keyboard);

    // No animation clock: the scene is static unless the mouse/keyboard changes
    // it, so we redraw only on input (orbit/pan/zoom/center-drag/keys/reshape).
    // Force one initial redraw so the inverted lattice shows at startup.
    glutPostRedisplay();

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
    GLfloat position[] = {5.0f, 5.0f, 5.0f, 0.0f};  // Match Main.c light position

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

    // Clear color - default light steel-blue (softer than pure white; the black facet
    // outlines stay visible). [B] cycles a runtime palette (see g_bgPalette / display()).
    glClearColor(0.86f, 0.89f, 0.93f, 1.0f);

    // Report the GL renderer so it's clear whether we're on hardware or software
    // (llvmpipe) GL — software rasterization of ~414k triangles is the remaining
    // cost once the immediate-mode overhead is gone.
    printf("GL renderer: %s\n", (const char*)glGetString(GL_RENDERER));
    printf("GL version:  %s\n", (const char*)glGetString(GL_VERSION));

    // Load the Gouraud shader. If the GL context can't compile it, g_shaderOK stays
    // false and the [L] toggle is inert (fixed-function stays fully usable).
    initGouraudShader();
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
        "[L]: Toggle Gouraud Shader (lattice shading)",
        "[B]: Cycle Background Colour",
        "[C]: Save STL Snapshot (~/Downloads)",
        "[Space]: Pause/Resume Anim",
        "[M]: Toggle Slow/Original Inversion",
        "[s/S]: Inversion speed +/-",
        "[-/=]: Inset -/+ (hollow border thickness; same as ,/.)",
        "[,/.]: Inset -/+ (hollow border thickness)",
        "[O]: Toggle Hollow (hex/square hole per face)",
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

    // Fisheye status + controls (f toggles, z/Z adjusts strength) — always shown.
    {
        y -= 0.05f;
        char fisheyeStatus[120];
        sprintf(fisheyeStatus, "Fisheye: %s (Strength: %.1f)  [f toggle, z/Z strength]",
                g_fisheyeMode ? "ON" : "OFF", g_fisheyeStrength);
        glColor3f(0.2f, 0.2f, 0.8f);  // Blue color for fisheye
        glRasterPos2f(0.02f, y);
        for(const char* c=fisheyeStatus; *c; ++c)
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);
    }

    // Inversion-center status: live center (x,y,z) and radius
    {
        y -= 0.05f;
        char invStatus[200];
        sprintf(invStatus, "InvCenter: (%.4f, %.4f, %.4f)  r=%.3f",
                g_invCenter.x(), g_invCenter.y(), g_invCenter.z(), g_invRadius);
        glColor3f(0.0f, 0.6f, 0.0f);  // green
        glRasterPos2f(0.02f, y);
        for(const char* c=invStatus; *c; ++c)
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);
    }

    // Input hints
    {
        y -= 0.05f;
        const char* hint = "L-drag: center XZ | Alt+L: Y | Shift+L: pan | wheel-click: orbit | wheel: zoom | r/R: radius | 0: reset | c: STL";
        glColor3f(0.5f, 0.5f, 0.5f);
        glRasterPos2f(0.02f, y);
        for(const char* c=hint; *c; ++c)
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);
    }

    // Frame draw time (performance readout). With input-only redraw this is the
    // last frame's work; orbit/pan reuse the cached mesh so it stays small.
    {
        y -= 0.05f;
        char fr[120];
        sprintf(fr, "Frame: %.1f ms  (~%.0f fps)", g_lastDrawMs,
                g_lastDrawMs > 0.0f ? 1000.0f / g_lastDrawMs : 0.0f);
        glColor3f(0.8f, 0.5f, 0.2f);
        glRasterPos2f(0.02f, y);
        for(const char* c=fr; *c; ++c)
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);
    }

    // Cuboctahedron lattice (octagon-connected) inset + hollow status
    {
        y -= 0.05f;
        char toStatus[160];
        sprintf(toStatus, "Cuboctahedron lattice (octagon-connected)  Inset=%.2f  Hollow: %s",
                g_inset, g_hollow ? "ON (octagon windows)" : "OFF (solid)");
        glColor3f(0.2f, 0.6f, 0.6f);  // teal
        glRasterPos2f(0.02f, y);
        for(const char* c=toStatus; *c; ++c)
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);
    }

    // Gouraud shader + background status
    {
        y -= 0.05f;
        char sh[120];
        sprintf(sh, "Gouraud: %s  |  Background: %d/%d",
                g_shaderOK ? (g_useShader ? "ON" : "OFF") : "unavailable",
                g_bgIndex + 1, g_bgPaletteCount);
        glColor3f(0.5f, 0.3f, 0.7f);  // purple
        glRasterPos2f(0.02f, y);
        for(const char* c=sh; *c; ++c)
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);
    }
    
    glMatrixMode(GL_MODELVIEW); glPopMatrix();
    glMatrixMode(GL_PROJECTION); glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
}

//-----------------------------------------------------------------------------
// Main display with surface zoom and fisheye
//-----------------------------------------------------------------------------
// Forward decl: defined after writeFacetBoxSTL (used here for per-frame cycle capture).
static void maybeCaptureFrame();

void display()
{
    int frameT0 = glutGet(GLUT_ELAPSED_TIME);   // measure this frame's draw work
    // Background from the current palette entry ([B] cycles); light entries keep the
    // black facet outlines visible.
    glClearColor(g_bgPalette[g_bgIndex][0], g_bgPalette[g_bgIndex][1],
                 g_bgPalette[g_bgIndex][2], 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    // Apply surface zoom to the base zoom
    float effectiveZoom = g_zoom / g_surfaceZoomFactor;
    
    // Move back and zoom (with surface zoom effect)
    glTranslatef(g_panX, g_panY, -5.0f * effectiveZoom);

    // Apply rotations
    glRotatef(g_angleX, 1, 0, 0);
    glRotatef(g_angleY, 0, 1, 0);

    // Cache the camera matrices for mouse->world unprojection (used by the
    // left-drag inversion-center placement). Read here, after the camera
    // translate+rotate, so the ray is in the same world space as the lattice
    // (the lattice is drawn in world space with no extra model transform).
    glGetDoublev(GL_MODELVIEW_MATRIX, g_modelview);
    glGetDoublev(GL_PROJECTION_MATRIX, g_projection);
    glGetIntegerv(GL_VIEWPORT, g_viewport);

    // Deform the lattice to the current inversion center/radius only when they
    // changed (mouse moved the center or r/R changed the radius). Setup() calls
    // this once before the first frame, so the lattice starts fully inverted.
    if (g_geomDirty) { updateAnimatedGeometry(); g_geomDirty = false; }

    // Per-frame STL capture for one cycle (if --capture-cycle). Done after the
    // mesh is deformed so the STL matches the rendered frame.
    if (g_capturing) maybeCaptureFrame();

    // Draw axes at origin
    drawAxes(10.0f);

    ProcessingProto();   // calls Setup() then Draw()

    g_lastDrawMs = (float)(glutGet(GLUT_ELAPSED_TIME) - frameT0);   // draw work done; shown in HUD
    drawHUD();
    
    glutSwapBuffers();
}

//-----------------------------------------------------------------------------
// Move the inversion center on the XZ "desk" plane by the mouse drag delta.
// Desk = XZ plane (Y up): dx -> X, dy -> Z. Low gain (g_invDragGain) makes it
// less sensitive than the old absolute cursor unproject. World-frame (independent
// of camera orientation); leaves Y unchanged. No snap on press -- you grab and
// push the center; lift and re-drag to continue.
//-----------------------------------------------------------------------------
static void dragInvCenterXZ(int dx, int dy) {
    const double k = g_invDragGain * (double)g_zoom;
    g_invCenter = Vector3D(g_invCenter.x() + dx * k,
                           g_invCenter.y(),
                           g_invCenter.z() + dy * k);   // drag down (dy>0) -> +Z
    g_geomDirty = true;
}

//-----------------------------------------------------------------------------
// Mouse button: wheel-click = orbit, left = inversion center (plain XZ / Alt Y),
// Shift+left = pan, right = zoom, wheel roll = zoom. Plain left only arms the
// drag (the center is moved by the drag delta in mouseMotion; no snap on press).
//-----------------------------------------------------------------------------
void mouseButton(int button, int state, int x, int y)
{
    if (button == GLUT_LEFT_BUTTON) {
        g_leftDown = (state == GLUT_DOWN);   // plain/Alt/Shift+left all arm the drag; motion handler does the work
    }
    else if (button == GLUT_MIDDLE_BUTTON) {
        g_middleDown = (state == GLUT_DOWN);   // wheel-click + hold = orbit camera
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
// Mouse drag: wheel-click = orbit, Alt+left = center Y (height above the XZ
// desk), Shift+left = pan, plain left = center XZ (drag delta), right = zoom.
//-----------------------------------------------------------------------------
void mouseMotion(int x, int y)
{
    int dx = x - g_lastX;
    int dy = y - g_lastY;
    int mods = glutGetModifiers();

    if (g_middleDown) {                       // wheel-click + hold = orbit (was left)
        g_angleY += dx * 0.5f;
        g_angleX += dy * 0.5f;
        g_angleX = fmaxf(-90.0f, fminf(90.0f, g_angleX));
    }
    else if (g_leftDown) {
        if (mods & GLUT_ACTIVE_ALT) {          // Alt + left: move center in Y (height above the XZ desk; XZ frozen)
            double ny = g_invCenter.y() - dy * 0.01 * (double)g_zoom;   // drag up = +Y
            g_invCenter = Vector3D(g_invCenter.x(), ny, g_invCenter.z());
            g_geomDirty = true;
        }
        else if (mods & GLUT_ACTIVE_SHIFT) {  // Shift + left: drag/pan camera (was middle)
            g_panX += dx * 0.01f * g_zoom;
            g_panY -= dy * 0.01f * g_zoom;
        }
        else {                                 // plain left: move center on the XZ desk by drag delta
            dragInvCenterXZ(dx, dy);
        }
    }
    else if (g_rightDown) {
        g_zoom *= 1.0f - dy * 0.005f;
        g_zoom = fmaxf(0.1f, fminf(10.0f, g_zoom));
    }

    g_lastX = x; g_lastY = y;
    glutPostRedisplay();
}

//-----------------------------------------------------------------------------
// Write an arbitrary FacetBox to an ASCII STL file (same format as Cube::writeSTL).
//-----------------------------------------------------------------------------
static void writeFacetBoxSTL(const FacetBox& facets, const std::string& path,
                             const char* name) {
    std::ofstream stl(path);
    if (!stl.is_open())
        throw std::runtime_error("captureSTLSnapshot: cannot create " + path);
    stl << "solid " << name << "\n";
    for (size_t i = 0; i < facets.size(); ++i) {
        const Facet& f = facets[i];
        Vector3D n = f.getNormal();
        Vector3D a = f[0], b = f[1], c = f[2];
        stl << "facet normal " << std::scientific
            << n.x() << " " << n.y() << " " << n.z() << "\n";
        stl << "\touter loop\n";
        stl << "\t\tvertex " << a.x() << " " << a.y() << " " << a.z() << "\n";
        stl << "\t\tvertex " << b.x() << " " << b.y() << " " << b.z() << "\n";
        stl << "\t\tvertex " << c.x() << " " << c.y() << " " << c.z() << "\n";
        stl << "\tendloop\n";
        stl << "endfacet\n";
    }
    stl << "endsolid " << name << "\n";
}

// Write the FacetBox as BINARY STL (little-endian, 50 bytes/triangle) to `path`.
// ~4.8x smaller than ASCII; Blender auto-detects binary via file size (80+4+50*N).
static void writeFacetBoxSTLBinary(const FacetBox& facets, const std::string& path,
                                   const char* /*name*/) {
    std::ofstream stl(path, std::ios::binary);
    if (!stl.is_open())
        throw std::runtime_error("writeFacetBoxSTLBinary: cannot create " + path);
    const char header[80] = "Binary STL - main15_slow_inversion";  // rest padded with '\0' -> 80 bytes
    stl.write(header, 80);
    const uint32_t n = (uint32_t)facets.size();
    stl.write(reinterpret_cast<const char*>(&n), 4);
    for (size_t i = 0; i < facets.size(); ++i) {
        const Facet& f = facets[i];
        Vector3D nv = f.getNormal();
        Vector3D a = f[0], b = f[1], c = f[2];
        float vals[12] = {
            (float)nv.x(), (float)nv.y(), (float)nv.z(),
            (float)a.x(),  (float)a.y(),  (float)a.z(),
            (float)b.x(),  (float)b.y(),  (float)b.z(),
            (float)c.x(),  (float)c.y(),  (float)c.z()
        };
        stl.write(reinterpret_cast<const char*>(vals), 48);
        const uint16_t attr = 0;
        stl.write(reinterpret_cast<const char*>(&attr), 2);
    }
}

// Dispatch to ASCII or binary STL based on the --binary flag.
static void writeFacetBox(const FacetBox& facets, const std::string& path, const char* name) {
    if (g_binarySTL) writeFacetBoxSTLBinary(facets, path, name);
    else             writeFacetBoxSTL(facets, path, name);
}

// Write the current on-screen mesh (checkerboard selection, same as Draw()/'c')
// as STL (ASCII or binary per --binary) to <g_captureDir>/frame_NNNNNN.stl.
static void captureFrameSTL(int frameIndex) {
    char fname[64];
    snprintf(fname, sizeof(fname), "frame_%06d.stl", frameIndex);
    std::string path = g_captureDir + "/" + fname;
    try {
        FacetBox sel = cube.getCheckerboardFacetsCuboctahedronLattice(ii, kk, jj, g_hollow, g_inset);  // matches Draw()/captureSTLSnapshot()
        writeFacetBox(sel, path, "CuboctahedronLattice");
    } catch (const std::exception& e) {
        fprintf(stderr, "[capture] FAILED frame %d -> %s: %s -- aborting capture.\n",
                frameIndex, path.c_str(), e.what());
        g_capturing = false;
        exit(1);
    }
}

// Drive per-frame capture for one cycle. Called once per display() after the
// mesh is deformed. Frame 0 is the identity (progress 0, captured on the first
// display). Each subsequent animation step captures one frame; when the phase
// clock wraps (g_invProgress decreases) one full cycle is complete and we exit.
static void maybeCaptureFrame() {
    if (g_prevProgress < 0.0) {
        // First frame of the run: mesh is at progress 0 (identity). Capture once.
        captureFrameSTL(g_frameIndex);
        printf("[capture] frame %06d  progress=%.3f  -> %s/frame_%06d.stl\n",
               g_frameIndex, g_invProgress, g_captureDir.c_str(), g_frameIndex);
        g_frameIndex++;
        g_prevProgress = 0.0;
        return;
    }
    if (!g_animStepped) return;  // skip non-animation displays (resize/expose)
    if (g_invProgress < g_prevProgress) {
        // Phase clock wrapped past 6 -> one full cycle complete. The wrapped frame
        // would be the identity again (duplicate of frame 0), so stop without it.
        printf("[capture] cycle complete: %d frames written to %s -- exiting.\n",
               g_frameIndex, g_captureDir.c_str());
        fflush(stdout);
        exit(0);
    }
    captureFrameSTL(g_frameIndex);
    if (g_frameIndex % 100 == 0)
        printf("[capture] frame %06d  progress=%.3f\n", g_frameIndex, g_invProgress);
    g_frameIndex++;
    g_prevProgress = g_invProgress;
    g_animStepped = false;
}

//-----------------------------------------------------------------------------
// Capture an STL snapshot of EXACTLY what is on screen to ~/Downloads.
// Bound to the 'c' key. The render (Draw()) shows the checkerboard selection
// cube.getCheckerboardFacetsCuboctahedronLattice(ii, kk, jj, g_hollow, g_inset) — predicate (i%ii==0 && k%jj==0) || j%kk==0,
// i.e. the cells selected by the current <ii,jj,kk> modular-algebra values. The
// snapshot builds that SAME FacetBox (same call, same current animated subcell
// vertices) and writes it, so the STL matches the OpenGL render triangle-for-
// triangle — NOT the full mesh.
// Each snapshot gets a unique timestamped + counted filename so repeated presses
// never overwrite previous captures.
//-----------------------------------------------------------------------------
void captureSTLSnapshot()
{
    // Resolve ~/Downloads from $HOME (fall back to the hardcoded path used elsewhere)
    const char* home = getenv("HOME");
    std::string dir = home ? (std::string(home) + "/Downloads") : "/home/mike666/Downloads";

    // Best-effort: ensure the directory exists (ignore "already exists")
    mkdir(dir.c_str(), 0755);

    // Build a unique filename: YYYYMMDD_HHMMSS + per-session counter
    time_t now = time(nullptr);
    struct tm* lt = localtime(&now);
    char ts[32];
    strftime(ts, sizeof(ts), "%Y%m%d_%H%M%S", lt);

    char fname[256];
    snprintf(fname, sizeof(fname), "mesh_snapshot_truncated_octahedron_%s_%03d.stl", ts, g_snapshotCount++);

    std::string path = dir + "/" + fname;

    try {
        // EXACTLY the render's selection: Draw() does getCheckerboardFacetsCuboctahedronLattice(ii, kk, jj, ...).
        // It reads the current (animated) subcell vertices directly, so this matches
        // the last rendered frame — no refreshTriangulation needed, and the STL
        // matches what OpenGL shows triangle-for-triangle.
        FacetBox sel = cube.getCheckerboardFacetsCuboctahedronLattice(ii, kk, jj, g_hollow, g_inset);
        writeFacetBox(sel, path, "CuboctahedronLattice");
        cout << "[snapshot] rendered selection saved -> " << path
             << "  (facets: " << sel.size() << ")"
             << "  (t=" << g_animTime << ", " << (g_animPaused ? "PAUSED" : "running") << ")\n";
    } catch (const std::exception& e) {
        cerr << "[snapshot] FAILED: " << e.what() << "\n";
    }
}

//-----------------------------------------------------------------------------
// Keyboard handler with fisheye and surface zoom controls
//-----------------------------------------------------------------------------
void keyboard(unsigned char key, int x, int y) {
    switch (key) {
        // Gouraud shader on the lattice fills ([L]). Scoped to the fills; outlines and
        // the rest stay fixed-function. Inert (with a hint) if the shader failed to
        // compile/link for this GL context.
        case 'l':
        case 'L':
            if (!g_shaderOK) {
                printf("Gouraud shader unavailable on this GL context (using fixed-function).\n");
            } else {
                g_useShader = !g_useShader;
                printf("Gouraud shading: %s\n", g_useShader ? "ON" : "OFF");
            }
            glutPostRedisplay();
            break;

        // Cycle the background colour ([B]). Light entries keep the black outlines visible.
        case 'b':
        case 'B':
            g_bgIndex = (g_bgIndex + 1) % g_bgPaletteCount;
            printf("Background: %d/%d\n", g_bgIndex + 1, g_bgPaletteCount);
            glutPostRedisplay();
            break;

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

        // Save an STL snapshot of the ENTIRE mesh to ~/Downloads
        case 'c':
        case 'C':
            captureSTLSnapshot();
            break;

        // Pause/resume the inversion animation (t holds while paused)
        case ' ':
            g_animPaused = !g_animPaused;
            printf("Animation: %s (t=%.3f)\n", g_animPaused ? "PAUSED" : "RUNNING", g_animTime);
            break;

        // Fisheye controls
        case 'f':
        case 'F':
            g_fisheyeMode = !g_fisheyeMode;
            printf("Fisheye mode: %s\n", g_fisheyeMode ? "ON" : "OFF");
            reshape(glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
            glutPostRedisplay();
            break;
            
        // Fisheye strength: 'z' increase (more fisheye), 'Z' decrease (less).
        case 'z':
            g_fisheyeStrength = fminf(2.0f, g_fisheyeStrength + 0.1f);
            printf("Fisheye strength: %.1f\n", g_fisheyeStrength);
            if (g_fisheyeMode) {
                reshape(glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
            }
            glutPostRedisplay();
            break;
            
        case 'Z':
            g_fisheyeStrength = fmaxf(0.0f, g_fisheyeStrength - 0.1f);
            printf("Fisheye strength: %.1f\n", g_fisheyeStrength);
            if (g_fisheyeMode) {
                reshape(glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
            }
            glutPostRedisplay();
            break;
            
        // Switch 'k' module
        case 'k':
            kk += 1;
            cout << "modules <ii, jj, kk> = <" << ii << ", " << jj << ", " << kk << ">\n";
            glutPostRedisplay();
            break;
        case 'K':
            kk -= 1; if (kk < 1) kk = 1;
            cout << "modules <ii, jj, kk> = <" << ii << ", " << jj << ", " << kk << ">\n";
            glutPostRedisplay();
            break;

        // Switch 'i' module
        case 'i':
            ii += 1;
            cout << "modules <ii, jj, kk> = <" << ii << ", " << jj << ", " << kk << ">\n";
            glutPostRedisplay();
            break;
        case 'I':
            ii -= 1; if (ii < 1) ii = 1;
            cout << "modules <ii, jj, kk> = <" << ii << ", " << jj << ", " << kk << ">\n";
            glutPostRedisplay();
            break;

        // Switch 'j' module
        case 'j':
            jj += 1;
            cout << "modules <ii, jj, kk> = <" << ii << ", " << jj << ", " << kk << ">\n";
            glutPostRedisplay();
            break;
        case 'J':
            jj -= 1; if (jj < 1) jj = 1;
            cout << "modules <ii, jj, kk> = <" << ii << ", " << jj << ", " << kk << ">\n";
            glutPostRedisplay();
            break;

        // Slow-transition mode toggle (inert now): the slow transition was removed;
        // there is a single full inversion at the mouse center. Kept for compatibility.
        case 'm':
        case 'M':
            if (g_capturing) {
                printf("Inversion mode locked to SLOW during cycle capture.\n");
                break;
            }
            g_slowMode = !g_slowMode;
            printf("Slow transition removed: single inversion at mouse center (r/R for radius).\n");
            glutPostRedisplay();
            break;

        // Inversion speed (inert now): the slow transition was removed; the
        // inversion center is mouse-driven. Kept for compatibility.
        case 's':
            g_invSpeed = fmin(10.0, g_invSpeed + g_invSpeedStep);
            printf("Slow transition removed: single inversion at mouse center (r/R for radius).\n");
            glutPostRedisplay();
            break;
        case 'S':
            g_invSpeed = fmax(1e-4, g_invSpeed - g_invSpeedStep);
            printf("Slow transition removed: single inversion at mouse center (r/R for radius).\n");
            glutPostRedisplay();
            break;

        // Inset scrub (live): '-' decrease (thicker border / smaller hole),
        // '='/'+' increase (bigger hole / thinner border). Same step/clamp as
        // ','/'.'. (The lattice has no TO morph, so these keys scrub inset here.)
        case '-':
            g_inset = fmax(0.05, g_inset - 0.05);
            printf("Inset = %.2f\n", g_inset);
            glutPostRedisplay();
            break;
        case '=':
        case '+':
            g_inset = fmin(0.95, g_inset + 0.05);
            printf("Inset = %.2f\n", g_inset);
            glutPostRedisplay();
            break;

        // Inset scrub: ',' decrease (thicker border / smaller hole), '.' increase (bigger hole).
        case ',':
            g_inset = fmax(0.05, g_inset - 0.05);
            printf("Inset = %.2f\n", g_inset);
            glutPostRedisplay();
            break;
        case '.':
        case '>':
            g_inset = fmin(0.95, g_inset + 0.05);
            printf("Inset = %.2f\n", g_inset);
            glutPostRedisplay();
            break;

        // Toggle hollow: each face inset toward its centroid, inner face skipped
        // => 144 frame tris/cell (hex/square holes) instead of 44 solid.
        case 'o':
        case 'O':
            g_hollow = !g_hollow;
            printf("Hollow (truncated octahedron): %s\n", g_hollow ? "ON (frame holes)" : "OFF (solid)");
            glutPostRedisplay();
            break;

        // Inversion radius: 'r' increase, 'R' decrease (0.05 step).
        case 'r':
            g_invRadius = fmin(5.0, g_invRadius + 0.05);
            g_geomDirty = true;
            printf("Inversion radius = %.3f\n", g_invRadius);
            glutPostRedisplay();
            break;
        case 'R':
            g_invRadius = fmax(0.05, g_invRadius - 0.05);
            g_geomDirty = true;
            printf("Inversion radius = %.3f\n", g_invRadius);
            glutPostRedisplay();
            break;

        // Reset inversion center to the origin and radius to 0.5.
        case '0':
            g_invCenter = Vector3D(0.0, 1e-9, 0.0);   // back to the XZ desk origin (y~0)
            g_invRadius = 0.5;
            g_geomDirty = true;
            printf("Inversion center reset to origin, r=0.5\n");
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
