#include <vector>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <cstring>        // memcpy (matMul)
#include <GL/glut.h>
#include <GL/glx.h>        // glXGetProcAddressARB — manual modern-GL loader (no GLEW)
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

//==============================================================================
// GPU rasterization (ported from main17_gpu.cpp).
// The per-vertex inversion math stays on the CPU (updateAnimatedGeometry mutates
// the cube's subcells, exactly as the CPU original); the octagonal checkerboard
// mesh is expanded into a flat VBO each frame and drawn with one indexed
// glDrawElements (fill + outline) through a GLSL Blinn-Phong shader — replacing
// the ~1.8M glBegin/glEnd batches the CPU Draw() issued per frame. The shader
// paints flat gray + black outline, matching the original drawFacetMainCStyle.
// No CUDA/thrust: the speedup is GPU-side rasterization, exactly like main17_gpu.
//==============================================================================

// Manual loader for post-GL-1.1 functions (no GLEW dependency). GL 1.x functions
// (glDrawElements, glClear, glPolygonMode, ...) link directly against libGL.
#define GLDECL(ret, name, args) typedef ret (*PFN_##name) args; static PFN_##name name
GLDECL(void,  glGenBuffers, (GLsizei, GLuint*));
GLDECL(void,  glDeleteBuffers, (GLsizei, const GLuint*));
GLDECL(void,  glBindBuffer, (GLenum, GLuint));
GLDECL(void,  glBufferData, (GLenum, GLsizeiptr, const void*, GLenum));
GLDECL(void,  glBufferSubData, (GLenum, GLintptr, GLsizeiptr, const void*));
GLDECL(GLuint, glCreateShader, (GLenum));
GLDECL(void,  glDeleteShader, (GLuint));
GLDECL(void,  glShaderSource, (GLuint, GLsizei, const GLchar* const*, const GLint*));
GLDECL(void,  glCompileShader, (GLuint));
GLDECL(GLuint, glCreateProgram, (void));
GLDECL(void,  glDeleteProgram, (GLuint));
GLDECL(void,  glAttachShader, (GLuint, GLuint));
GLDECL(void,  glBindAttribLocation, (GLuint, GLuint, const GLchar*));
GLDECL(void,  glLinkProgram, (GLuint));
GLDECL(void,  glUseProgram, (GLuint));
GLDECL(void,  glGetShaderiv, (GLuint, GLenum, GLint*));
GLDECL(void,  glGetShaderInfoLog, (GLuint, GLsizei, GLsizei*, GLchar*));
GLDECL(void,  glGetProgramiv, (GLuint, GLenum, GLint*));
GLDECL(void,  glGetProgramInfoLog, (GLuint, GLsizei, GLsizei*, GLchar*));
GLDECL(GLint,  glGetUniformLocation, (GLuint, const GLchar*));
GLDECL(void,  glUniformMatrix4fv, (GLint, GLsizei, GLboolean, const GLfloat*));
GLDECL(void,  glUniform3f, (GLint, GLfloat, GLfloat, GLfloat));
GLDECL(void,  glUniform1i, (GLint, GLint));
GLDECL(void,  glEnableVertexAttribArray, (GLuint));
GLDECL(void,  glDisableVertexAttribArray, (GLuint));
GLDECL(void,  glVertexAttribPointer, (GLuint, GLint, GLenum, GLboolean, GLsizei, const void*));
#undef GLDECL

static void loadGL() {
    #define L(name) name = (PFN_##name)glXGetProcAddressARB((const GLubyte*)#name)
    L(glGenBuffers); L(glDeleteBuffers); L(glBindBuffer); L(glBufferData); L(glBufferSubData);
    L(glCreateShader); L(glDeleteShader); L(glShaderSource); L(glCompileShader);
    L(glCreateProgram); L(glDeleteProgram); L(glAttachShader); L(glLinkProgram);
    L(glUseProgram); L(glGetShaderiv); L(glGetShaderInfoLog);
    L(glGetProgramiv); L(glGetProgramInfoLog); L(glGetUniformLocation); L(glBindAttribLocation);
    L(glUniformMatrix4fv); L(glUniform3f); L(glUniform1i);
    L(glEnableVertexAttribArray); L(glDisableVertexAttribArray); L(glVertexAttribPointer);
    #undef L
}

// Shaders: #version 460 core. VS passes world-space position (uMVP has view+proj);
// FS computes a flat per-triangle face normal from screen-space derivatives of the
// world position (dFdx/dFdy), then Blinn-Phong matching the fixed-function light
// from initGL(). uMode=1 => black outline. (Verbatim from main17_gpu.)
static const char* kVS =
    "#version 460 core\n"
    "in vec3 aPos;\n"
    "uniform mat4 uMVP;\n"
    "out vec3 vWorldPos;\n"
    "void main() {\n"
    "    vWorldPos = aPos;\n"
    "    gl_Position = uMVP * vec4(aPos, 1.0);\n"
    "}\n";

static const char* kFS =
    "#version 460 core\n"
    "in vec3 vWorldPos;\n"
    "uniform int   uMode;        // 0 = lit fill, 1 = black outline\n"
    "uniform vec3  uLightDir;    // world space (not necessarily normalized)\n"
    "uniform vec3  uCamPos;      // world space\n"
    "out vec4 fragColor;\n"
    "void main() {\n"
    "    if (uMode == 1) { fragColor = vec4(0.0, 0.0, 0.0, 1.0); return; }\n"
    "    vec3 N = normalize(cross(dFdx(vWorldPos), dFdy(vWorldPos)));\n"
    "    vec3 V = normalize(uCamPos - vWorldPos);\n"
    "    if (dot(N, V) < 0.0) N = -N;        // two-sided lighting\n"
    "    vec3 L = normalize(uLightDir);\n"
    "    vec3 H = normalize(L + V);\n"
    "    vec3 base = vec3(200.0 / 255.0);   // gray (matches glColor3ub(200,200,200))\n"
    "    float diff = max(dot(N, L), 0.0);\n"
    "    float spec = pow(max(dot(N, H), 0.0), 128.0);\n"
    "    vec3 col = 0.3 * base              // ambient  (GL_LIGHT0 ambient = 0.3)\n"
    "             + 0.7 * diff * base       // diffuse  (GL_LIGHT0 diffuse = 0.7)\n"
    "             + 1.0 * spec * vec3(1.0); // specular (GL_LIGHT0 specular = 1, white)\n"
    "    fragColor = vec4(col, 1.0);\n"
    "}\n";

static GLuint compileShader(GLenum type, const char* src) {
    GLuint s = glCreateShader(type);
    glShaderSource(s, 1, &src, nullptr);
    glCompileShader(s);
    GLint ok = 0; glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        char log[4096] = {0}; glGetShaderInfoLog(s, sizeof(log), nullptr, log);
        std::cerr << "[trunc-gpu] " << (type == GL_VERTEX_SHADER ? "VERTEX" : "FRAGMENT")
                  << " shader compile failed:\n" << log << std::endl;
        glDeleteShader(s);
        return 0;
    }
    return s;
}

static GLuint buildProgram() {
    GLuint vs = compileShader(GL_VERTEX_SHADER, kVS);
    GLuint fs = compileShader(GL_FRAGMENT_SHADER, kFS);
    if (!vs || !fs) return 0;
    GLuint prog = glCreateProgram();
    glAttachShader(prog, vs);
    glAttachShader(prog, fs);
    glBindAttribLocation(prog, 0, "aPos");   // attribute 0 = aPos
    glLinkProgram(prog);
    glDeleteShader(vs);
    glDeleteShader(fs);
    GLint ok = 0; glGetProgramiv(prog, GL_LINK_STATUS, &ok);
    if (!ok) {
        char log[4096] = {0}; glGetProgramInfoLog(prog, sizeof(log), nullptr, log);
        std::cerr << "[trunc-gpu] program link failed:\n" << log << std::endl;
        glDeleteProgram(prog);
        return 0;
    }
    return prog;
}

// Small column-major mat4 helpers (OpenGL convention). (Verbatim from main17_gpu.)
static void matMul(const float* A, const float* B, float* out) {
    float tmp[16];
    for (int c = 0; c < 4; ++c)
        for (int r = 0; r < 4; ++r) {
            float s = 0.0f;
            for (int k = 0; k < 4; ++k) s += A[k * 4 + r] * B[c * 4 + k];
            tmp[c * 4 + r] = s;
        }
    memcpy(out, tmp, sizeof(tmp));
}
static void xformPoint(const float* M, float x, float y, float z, float out[3]) {
    out[0] = M[0]*x + M[4]*y + M[8]*z  + M[12];
    out[1] = M[1]*x + M[5]*y + M[9]*z  + M[13];
    out[2] = M[2]*x + M[6]*y + M[10]*z + M[14];
}
// Transpose of the upper-left 3x3 of `view`, padded to 4x4 (translation = 0). Maps
// eye-space directions/points (light dir, camera origin) back to world space.
static void buildInvRotation(const float* view, float invR[16]) {
    memset(invR, 0, sizeof(float) * 16);
    for (int c = 0; c < 3; ++c)
        for (int r = 0; r < 3; ++r)
            invR[c * 4 + r] = view[r * 4 + c];  // transpose of 3x3
    invR[15] = 1.0f;
}

// GPU render state.
static GLuint g_vbo = 0;       // octagonal vertex positions (selected cells x 126 verts, dynamic)
static GLuint g_ibo = 0;       // triangle indices (rebuilt on ii/jj/kk or hollow change)
static GLuint g_program = 0;   // linked shader program
static GLint  g_locMVP = -1, g_locLightDir = -1, g_locCamPos = -1, g_locMode = -1;
static GLsizei g_indexCount = 0;
static bool g_selectionDirty = true;   // (re)build IBO + size VBO on first frame & on key change
static std::vector<float> g_deformedOctagonal;  // per-frame expanded vertex buffer, uploaded to g_vbo

// Forward declarations (defined after `cube` and the camera globals).
static void rebuildBuffers();
static void updateOctagonalGeometry();
static void drawGPU();

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
struct SigmaParams { Vector3D c1; double r1; Vector3D c2; double r2; };//Vector3D c3; double r3; };

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
// render shows (cube.getCheckerboardFacets(ii,kk,jj)), matching the 'c' key.
static bool        g_captureCycle = false;  // --capture-cycle flag
static std::string g_captureDir;           // --outdir target directory
static bool        g_capturing    = false; // active capture state (true while capturing)
static int         g_frameIndex   = 0;    // current frame number (0-based)
static double      g_prevProgress = -1.0;  // progress at last captured frame (-1 = none yet)
static bool        g_animStepped  = false; // set by animTimer when it advanced, consumed in display
static bool        g_binarySTL    = false; // --binary: write binary STL (else ASCII)
static bool        g_hollow       = false; // --hollow / [O]: omit center-fan tris (octagonal hole per face)

// Same animated inversion parameters as main17_gpu::currentSigmaParams(t).
static SigmaParams currentSigmaParams(double t) {
    SigmaParams s;
    s.c1 = Vector3D(0.1*sin(2*t), 0.0, 0.1 * sin(t) * cos(t));   // inversion point oscillates in z
    s.r1 = 0.5;
    s.c2 = cube.getSubcellCenter(1, 1, 1);
    s.r2 = cube.getSubcellRadius(1, 1, 1);
    //s.c3 = cube.getSubcellCenter(-1, -1, -1);
    //s.r3 = cube.getSubcellRadius(-1, -1, -1);
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
    return r2v;
    // pass 3: out = (1-a3)*r2v + a3*sigma3(r2v)
    //Vector3D d3 = r2v - s.c3; double ds3 = d3 * d3;
    //Vector3D s3v = (ds3 < 1e-12) ? r2v : s.c3 + (s.r3 * s.r3 / ds3) * d3;
    //return r2v + b.a3 * (s3v - r2v);
}
// Original behavior: all three passes at full strength (slow mode disabled).
static inline Vector3D composeSigma(const Vector3D& p, const SigmaParams& s) {
    return composeSigmaBlended(p, s, {1.0, 1.0, 1.0});
}

// Capture the pristine identity lattice ONCE, before any deformation.
static void captureIdentityLattice() {
    cube.fillVertexLattice(g_identityPositions);  // n^3*8*3 floats — same call as main17_gpu
    cout << "[anim] identity lattice captured: " << (g_identityPositions.size() / 3)
         << " verts (" << (g_identityPositions.size() * sizeof(float) / (1024.0 * 1024.0))
         << " MB)\n";
}

// Recompose the three animated inversions from the identity lattice and push the
// deformed vertices back into the Cube's subcells (same iteration order
// fillVertexLattice() uses: physical [0,n)^3, vertex 0..7). In slow mode each
// pass is blended by currentAlphas(g_invProgress); in original mode all three
// passes run at full strength (== the original composeSigma(p, sp)).
static void updateAnimatedGeometry() {
    SigmaParams sp = currentSigmaParams(g_animTime);
    const InvBlend blend = g_slowMode ? currentAlphas(g_invProgress)
                                      : InvBlend{1.0, 1.0, 1.0};
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
                    cube.updateSubCellVertex(lx, ly, lz, l, composeSigmaBlended(p, sp, blend));
                    idx += 3;
                }
            }
}

// ~30 Hz animation clock. Advances the inversion time and requests a redraw.
// The lattice deformation itself runs in display(); [Space] pauses/resumes.
// In slow mode also advances the phase clock g_invProgress by speed-scaled
// ticks (currentAlphas fmods it into [0,6), so no manual wrap is needed).
static void animTimer(int /*value*/) {
    if (!g_animPaused) {
        g_animTime += g_timeStep;
        if (g_slowMode)
            g_invProgress = fmod(g_invProgress + INV_RATE_PER_TICK * g_invSpeed, 6.0);  // wraps 6->0 so the cycle-completion detector (g_invProgress < g_prevProgress) fires
        g_animStepped = true;   // a new animation step occurred (drives per-frame capture)
    }
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
        captureIdentityLattice();
        updateAnimatedGeometry();
        cube.writeSTL_s_octagonal("/home/ubuntu/Downloads/mesh_output_truncated.stl", "MyTruncatedCube", "checkerboard", ii, kk, jj, g_hollow);

        // GPU path: build the octagonal IBO + size the VBO for the current
        // <ii,jj,kk,hollow> selection, then expand the t=0 deformed subcells into
        // the VBO and upload. Each frame display() re-expands + re-uploads.
        rebuildBuffers();
        updateOctagonalGeometry();
        g_selectionDirty = false;   // Setup already built the buffers; no rebuild on first display
    }
}

///////////////////     DRAW (GPU)     ///////////////////////
// The CPU immediate-mode Facet loop is replaced by drawGPU(): one indexed
// glDrawElements in two passes (lit fill + black outline) through the shader.
// drawFacetMainCStyle is no longer used for the render path (the STL snapshot
// path still builds Facets via getCheckerboardFacetsOctagonal on demand).
void Draw() {
    // Intentionally empty — the octagonal mesh is rasterized by drawGPU() below.
}

void ProcessingProto() {
    drawGPU();
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

//==============================================================================
// GPU geometry buffers (octagonal checkerboard lattice).
// rebuildBuffers(): build the static octagonal IBO + (re)allocate the VBO to the
//   current selection size. Independent of deformation — rebuild only when
//   <ii,jj,kk> or g_hollow changes (g_selectionDirty).
// updateOctagonalGeometry(): expand the (just-deformed) subcells into the flat
//   octagonal vertex buffer and upload it to the VBO. Called every animated frame.
// drawGPU(): one indexed glDrawElements, two passes (lit fill + black outline),
//   using the fixed-function PROJECTION/MODELVIEW matrices display() already set.
//==============================================================================
static void rebuildBuffers() {
    std::vector<unsigned int> indices;
    cube.fillCheckerboardIndicesOctagonal(ii, kk, jj, indices, g_hollow);
    g_indexCount = (GLsizei)indices.size();

    if (g_ibo == 0) glGenBuffers(1, &g_ibo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int),
                 indices.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    // VBO size = selectedCount * 126 verts * 3 floats. selectedCount derives from
    // the triangle count (152 solid / 104 hollow per cell) — identical to the
    // count fillOctagonalVertexLattice will emit (same predicate + active filter).
    const size_t trisPerCell = g_hollow ? 104 : 152;
    const size_t totalTris    = trisPerCell ? (indices.size() / 3) / trisPerCell * trisPerCell : 0;
    const size_t selectedCount = trisPerCell ? totalTris / trisPerCell : 0;
    const size_t floats = selectedCount * 126 * 3;
    g_deformedOctagonal.assign(floats, 0.0f);   // sized for glBufferSubData; refilled each frame

    if (g_vbo == 0) glGenBuffers(1, &g_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, g_vbo);
    glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)(floats * sizeof(float)),
                 g_deformedOctagonal.data(), GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    std::cout << "[trunc-gpu] IBO: " << (g_indexCount / 3) << " triangles ("
              << selectedCount << " cells); VBO: " << (floats / 3) << " verts ("
              << (floats * sizeof(float) / (1024.0 * 1024.0)) << " MB)\n";
}

static void updateOctagonalGeometry() {
    cube.fillOctagonalVertexLattice(ii, kk, jj, g_deformedOctagonal, g_hollow);
    if (g_deformedOctagonal.empty()) return;
    glBindBuffer(GL_ARRAY_BUFFER, g_vbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0,
                    (GLsizeiptr)(g_deformedOctagonal.size() * sizeof(float)),
                    g_deformedOctagonal.data());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

static void drawGPU() {
    if (!g_program || !g_vbo || !g_ibo || g_indexCount == 0) return;

    // uMVP = projection * modelview (read back from the fixed-function stack that
    // display() set up with the orbit camera + fisheye — same as main17_gpu).
    float proj[16], view[16], mvp[16];
    glGetFloatv(GL_PROJECTION_MATRIX, proj);
    glGetFloatv(GL_MODELVIEW_MATRIX, view);
    matMul(proj, view, mvp);

    // Camera world position and light direction. The light is fixed in eye space
    // (set with identity modelview in initGL); invR maps them back to world space.
    // display()'s modelview is Translate(panX,panY,-5*effZoom) * RotX * RotY, so the
    // world camera is invR * (-panX, -panY, 5*effZoom) (only affects the spec term).
    float invR[16]; buildInvRotation(view, invR);
    float effectiveZoom = g_surfaceZoomFactor > 0.0f ? (g_zoom / g_surfaceZoomFactor) : g_zoom;
    float camPos[3];   xformPoint(invR, -g_panX, -g_panY, 5.0f * effectiveZoom, camPos);
    float lightDir[3]; xformPoint(invR, 1.0f, 1.0f, 0.25f, lightDir);

    glUseProgram(g_program);
    glUniformMatrix4fv(g_locMVP, 1, GL_FALSE, mvp);
    glUniform3f(g_locLightDir, lightDir[0], lightDir[1], lightDir[2]);
    glUniform3f(g_locCamPos,   camPos[0],   camPos[1],   camPos[2]);

    glBindBuffer(GL_ARRAY_BUFFER, g_vbo);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_ibo);

    // Fill pass — gray lit triangles, pushed back by polygon offset so the
    // outline wins the z-fight (matches the original glPolygonOffset(1,1)).
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0f, 1.0f);
    glUniform1i(g_locMode, 0);
    glDrawElements(GL_TRIANGLES, g_indexCount, GL_UNSIGNED_INT, (void*)0);
    glDisable(GL_POLYGON_OFFSET_FILL);

    // Outline pass — black triangle edges (GL_LINE on GL_TRIANGLES draws each
    // triangle's 3 edges, matching the original GL_LINE_LOOP per facet).
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glLineWidth(1.0f);
    glUniform1i(g_locMode, 1);
    glDrawElements(GL_TRIANGLES, g_indexCount, GL_UNSIGNED_INT, (void*)0);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    // Restore fixed-function state for axes/HUD.
    glDisableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glUseProgram(0);
}

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
        else if (a == "--hollow")                             g_hollow = true;   // omit inner octagon (hollow)
        else { char* end = nullptr; double v = strtod(argv[i], &end);  // positional => speed
               if (end != argv[i] && *end == '\0') g_invSpeed = v; }
    }
    if (g_invSpeed < 1e-4) g_invSpeed = 1e-4;   // never stall
    if (g_invSpeed > 10.0) g_invSpeed = 10.0;

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
            g_captureDir = base + "/cycle_truncated_" + ts;
        }
        ensureDir(g_captureDir);
        g_capturing = true;                                  // start capturing on the first display
        long expect = (long)(6.0 / (INV_RATE_PER_TICK * g_invSpeed) + 0.5);
        printf("[capture] per-frame %s STL for one cycle -> %s\n",
               g_binarySTL ? "BINARY" : "ASCII", g_captureDir.c_str());
        printf("[capture] ~%ld frames expected (speed=%.3f); program exits when the cycle wraps.\n",
               expect, g_invSpeed);
    }

    printf("[slow-inv] mode=%s  speed=%.3f  (per-inversion ~%.1fs)  hollow=%s\n",
           g_slowMode ? "SLOW" : "ORIGINAL", g_invSpeed,
           g_slowMode ? 1.0 / (INV_RATE_PER_TICK * 30.0 * g_invSpeed) : 0.0,
           g_hollow ? "ON" : "OFF");
}

int main(int argc, char** argv)
{
    srand((unsigned)time(nullptr));
    glutInit(&argc, argv);
    parseArgs(argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(720, 720);
    glutCreateWindow(" JAZ 4D Enhanced Camera U.U  (Truncated Cube)");

    // Load post-GL-1.1 entry points + build the shader program (needs a context).
    loadGL();
    g_program = buildProgram();
    if (g_program) {
        g_locMVP      = glGetUniformLocation(g_program, "uMVP");
        g_locLightDir = glGetUniformLocation(g_program, "uLightDir");
        g_locCamPos   = glGetUniformLocation(g_program, "uCamPos");
        g_locMode     = glGetUniformLocation(g_program, "uMode");
    } else {
        std::cerr << "[trunc-gpu] WARNING: shader program failed to build; "
                     "no mesh will be drawn.\n";
    }

    // Enable smoothing & blending by default
    ProcessMenu(1);
    initGL();
    createUI();
    
    // Objects setup (captures the identity lattice, deforms to t=0)
    Setup();

    // Start the ~30 Hz inversion-animation clock (advances g_animTime)
    glutTimerFunc(33, animTimer, 0);

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
        "[C]: Save STL Snapshot (~/Downloads)",
        "[Space]: Pause/Resume Anim",
        "[M]: Toggle Slow/Original Inversion",
        "[s/S]: Inversion speed +/-",
        "[O]: Toggle Hollow (octagonal hole per face)",
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

    // Animation status: inversion clock + pause state
    {
        y -= 0.05f;
        char animStatus[120];
        sprintf(animStatus, "Anim: %s  t=%.3f  c1.z=%.4f",
                g_animPaused ? "PAUSED" : "RUNNING", g_animTime, 0.1 * sin(g_animTime));
        glColor3f(0.0f, 0.5f, 0.0f);  // green
        glRasterPos2f(0.02f, y);
        for(const char* c=animStatus; *c; ++c)
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);
    }

    // Slow-inversion status: mode, speed, phase, and the three per-pass blends
    {
        y -= 0.05f;
        char invStatus[160];
        const char* ph = g_slowMode ? phaseName(g_invProgress) : "orig";
        InvBlend b = g_slowMode ? currentAlphas(g_invProgress) : InvBlend{1.0, 1.0, 1.0};
        sprintf(invStatus, "Inv: %s  speed=%.2f  phase=%s  a=(%.2f,%.2f,%.2f)",
                g_slowMode ? "SLOW" : "ORIG", g_invSpeed, ph, b.a1, b.a2, b.a3);
        glColor3f(0.6f, 0.3f, 0.6f);  // purple
        glRasterPos2f(0.02f, y);
        for(const char* c=invStatus; *c; ++c)
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);
    }

    // Hollow-truncated-cube status (octagonal holes per face)
    {
        y -= 0.05f;
        char hollowStatus[100];
        sprintf(hollowStatus, "Hollow: %s", g_hollow ? "ON (octagonal holes)" : "OFF (solid)");
        glColor3f(0.2f, 0.6f, 0.6f);  // teal
        glRasterPos2f(0.02f, y);
        for(const char* c=hollowStatus; *c; ++c)
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
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    // Apply surface zoom to the base zoom
    float effectiveZoom = g_zoom / g_surfaceZoomFactor;
    
    // Move back and zoom (with surface zoom effect)
    glTranslatef(g_panX, g_panY, -5.0f * effectiveZoom);

    // Apply rotations
    glRotatef(g_angleX, 1, 0, 0);
    glRotatef(g_angleY, 0, 1, 0);

    // Drive the realtime inversion animation: the timer advances g_animTime;
    // here we deform the lattice to the current frame before rendering it
    // (mirrors main17_gpu, which deforms in display() and advances in the timer).
    if (!g_animPaused) updateAnimatedGeometry();

    // GPU buffers: rebuild the octagonal IBO + size the VBO when the selection
    // (ii/jj/kk or hollow) changed, then expand the deformed subcells into the
    // VBO and upload. Upload runs every animated frame, and also after a selection
    // change so a paused mesh still shows the correct (reselected) geometry.
    bool selectionChanged = g_selectionDirty;
    if (selectionChanged) { rebuildBuffers(); g_selectionDirty = false; }
    if (!g_animPaused || selectionChanged) updateOctagonalGeometry();

    // Per-frame STL capture for one cycle (if --capture-cycle). Done after the
    // mesh is deformed so the STL matches the rendered frame.
    if (g_capturing) maybeCaptureFrame();

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
        FacetBox sel = cube.getCheckerboardFacetsOctagonal(ii, kk, jj, g_hollow);  // matches Draw()/captureSTLSnapshot()
        writeFacetBox(sel, path, "MyTruncatedCube");
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
// cube.getCheckerboardFacets(ii, kk, jj) — predicate (i%ii==0 && k%jj==0) || j%kk==0,
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
    snprintf(fname, sizeof(fname), "mesh_snapshot_truncated_%s_%03d.stl", ts, g_snapshotCount++);

    std::string path = dir + "/" + fname;

    try {
        // EXACTLY the render's selection: Draw() does getCheckerboardFacets(ii, kk, jj).
        // getCheckerboardFacets() reads the current (animated) subcell vertices
        // directly, so this matches the last rendered frame — no refreshTriangulation
        // needed, and the STL matches what OpenGL shows triangle-for-triangle.
        FacetBox sel = cube.getCheckerboardFacetsOctagonal(ii, kk, jj, g_hollow);
        writeFacetBox(sel, path, "MyTruncatedCube");
        cout << "[snapshot] rendered selection saved -> " << path
             << "  (facets: " << sel.size()
             << ", ii/jj/kk=" << ii << "/" << jj << "/" << kk << ")"
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
            
        // Switch 'k' module
        case 'k':
            kk += 1;
            cout << "modules <ii, jj, kk> = <" << ii << ", " << jj << ", " << kk << ">\n";
            g_selectionDirty = true;
            glutPostRedisplay();
            break;
        case 'K':
            kk -= 1;
            cout << "modules <ii, jj, kk> = <" << ii << ", " << jj << ", " << kk << ">\n";
            g_selectionDirty = true;
            glutPostRedisplay();
            break;

        // Switch 'i' module
        case 'i':
            ii += 1;
            cout << "modules <ii, jj, kk> = <" << ii << ", " << jj << ", " << kk << ">\n";
            g_selectionDirty = true;
            glutPostRedisplay();
            break;
        case 'I':
            ii -= 1;
            cout << "modules <ii, jj, kk> = <" << ii << ", " << jj << ", " << kk << ">\n";
            g_selectionDirty = true;
            glutPostRedisplay();
            break;

        // Switch 'j' module
        case 'j':
            jj += 1;
            cout << "modules <ii, jj, kk> = <" << ii << ", " << jj << ", " << kk << ">\n";
            g_selectionDirty = true;
            glutPostRedisplay();
            break;
        case 'J':
            jj -= 1;
            cout << "modules <ii, jj, kk> = <" << ii << ", " << jj << ", " << kk << ">\n";
            g_selectionDirty = true;
            glutPostRedisplay();
            break;

        // Toggle slow-sequential vs original simultaneous inversion mode
        case 'm':
        case 'M':
            if (g_capturing) {
                printf("Inversion mode locked to SLOW during cycle capture.\n");
                break;
            }
            g_slowMode = !g_slowMode;
            printf("Inversion mode: %s\n", g_slowMode ? "SLOW sequential" : "ORIGINAL simultaneous");
            glutPostRedisplay();
            break;

        // Inversion speed: 's' faster, 'S' slower (lower=+ matches i/I, j/J, k/K)
        case 's':
            g_invSpeed = fmin(10.0, g_invSpeed + g_invSpeedStep);
            printf("Inversion speed: %.3f (per-inversion ~%.1fs)\n", g_invSpeed,
                   1.0 / (INV_RATE_PER_TICK * 30.0 * g_invSpeed));
            glutPostRedisplay();
            break;
        case 'S':
            g_invSpeed = fmax(1e-4, g_invSpeed - g_invSpeedStep);
            printf("Inversion speed: %.3f\n", g_invSpeed);
            glutPostRedisplay();
            break;

        // Toggle hollow truncated cubes: omit the 8 center-fan triangles per face,
        // leaving an octagonal hole in each face so the cube is hollow.
        case 'o':
        case 'O':
            g_hollow = !g_hollow;
            printf("Hollow truncated cubes: %s\n", g_hollow ? "ON (octagonal holes)" : "OFF (solid)");
            g_selectionDirty = true;
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
