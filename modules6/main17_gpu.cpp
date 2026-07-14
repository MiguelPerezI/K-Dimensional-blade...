#include <vector>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <math.h>
#include <GL/glut.h>
#include <GL/glx.h>     // glXGetProcAddressARB (manual modern-GL loader, no GLEW needed)
#include "Vector3D.cpp"
#include "Vector4D.cpp"
#include "Quaternion.cpp"
#include "Facet.cpp"
#include "FacetBox.hpp"
#include "Dodecahedron.hpp"
#include "Cube.hpp"

//////////////////////////////////////
//
//  main17_gpu.cpp  —  car camera that drives on the mesh surface.
//
//  Builds on main16_gpu.cpp (the realtime-animated inversion build). The mesh is a
//  3D checkerboard of small subcell-cubes warped each frame by two sphere-inversions
//  (composeSigma) into an inflating/contracting blob. The inversion is applied to a
//  STATIC identity lattice (g_identityPositions); the Cube object never moves.
//
//  New in main17: a "car camera". The car lives in IDENTITY space — it drives on the
//  flat top face of the identity cube (y = yMax, x,z within the cube extent). Each
//  frame, the SAME composeSigma that deforms the mesh maps the car's identity-space
//  position + orientation (via a finite-difference Jacobian) into world space, so the
//  camera rides the blob's outer skin exactly — sticking to the surface through every
//  inflation/contraction with no mesh-topology/adjacency work.
//
//    Controls:  W/S throttle | A/D steer | C toggles car<->orbit cam | Esc quit
//    View: cockpit/first-person on the skin; starts at rest in the centre of the face.
//
//  Inherited from main16_gpu (unchanged): the NVIDIA GLX backend force, the dynamic
//  VBO upload path, the index buffer rebuilt only on ii/jj/kk change, the #version 460
//  core shader, and the fixed-function matrix stack that the shader's uMVP is read back
//  from (so the car camera lands correctly with no shader changes).
//
//////////////////////////////////////

//////////////////////////////////////
//  (legacy main16 header below)
//////////////////////////////////////
//
//  Changes vs main15_gpu.cpp:
//    1. The two sphere-inversion (sigma) passes that main15 baked once in Setup()
//       (lines 306/307) are now driven per frame from the IDENTITY lattice:
//         - pass 1 center oscillates in z:  Vector3D(0, 0, sin(t)), radius 0.5
//         - pass 2 center/radius walk a curated list of lattice corners via
//           cube.getSubcellCenter(i,j,k) / cube.getSubcellRadius(i,j,k)
//       The cube itself stays pristine (identity) forever; each frame the
//       composed sigma is computed from a cached identity vertex buffer into the
//       VBO (now GL_DYNAMIC_DRAW). This avoids the ~608k-Facet cost of
//       refreshTriangulation()/subdivide() per frame and the frame-over-frame
//       compounding that re-applying sigma in place would cause.
//    2. A glutTimerFunc drives the animation clock; [Space] pauses/resumes.
//
//  Inherited from main15_gpu: the NVIDIA GLX backend force, the static-style VBO
//  upload path (now dynamic), the index buffer rebuilt only on ii/jj/kk change,
//  and the #version 460 core shader. The fixed-function matrix stack is kept
//  exactly as in main15 (reshape sets the projection, display sets the modelview).
//
//////////////////////////////////////

//////////////////////////////////////
//  (legacy header from main15_gpu below)
//////////////////////////////////////
//
//  GPU-accelerated build of main15.
//
//  Two changes vs main15.cpp:
//    1. Force the NVIDIA GLX backend (setenv before glutInit) so the Tesla T4
//       rasterizes instead of Mesa llvmpipe (software).
//    2. Replace the per-frame ~264K glBegin/glEnd batches + ~132K Facet
//       constructions with: a static VBO (uploaded once after Setup()), an
//       index buffer rebuilt ONLY when ii/jj/kk change, and a #version 460 core
//       shader that computes per-triangle face normals via dFdx/dFdy.
//
//  The fixed-function matrix stack is kept exactly as in main15 (reshape sets
//  the projection, display sets the modelview). The shader's uMVP is built by
//  reading those matrices back with glGetFloatv, so the GPU geometry lands at
//  the same screen positions as the original, and the axes/HUD (which stay in
//  fixed-function) are unchanged.
//
//////////////////////////////////////

//////////////////////////////////////
//  VARIABLES GLOBALES PARA EL TECLADO
//////////////////////////////////////

int ciclo = 0;
double count = 0.25 * M_PI;
double count2 = 0.25 * 3.14159265358979;
double count3 = 1e-2;
double rad = 5.0;

bool g_glassMode = false;

static bool g_fisheyeMode = false;
static float g_fisheyeStrength = 1.0f;

static float g_minDistanceToSurface = 0.1f;
static float g_surfaceZoomFactor = 1.0f;
static bool g_enableSurfaceZoom = true;

int g_currentOrientation = 0;
int g_currentLayer = 4;
int g_maxLayers = 8;

//==============================================================================
// Manual loader for post-GL-1.1 functions (no GLEW dependency).
// GL 1.x functions (glDrawElements, glClear, glPolygonMode, ...) link directly.
//==============================================================================
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

//==============================================================================
// GPU render state
//==============================================================================
static GLuint g_vbo = 0;       // vertex positions (n^3 * 8 verts, static)
static GLuint g_ibo = 0;       // triangle indices (rebuilt on selection change)
static GLuint g_program = 0;   // linked shader program
static GLint  g_locMVP = -1, g_locLightDir = -1, g_locCamPos = -1, g_locMode = -1;
static GLsizei g_indexCount = 0;
static bool g_selectionDirty = true;  // (re)build IBO on first frame & on key press

//==============================================================================
// Realtime inversion animation (main16)
//   pass 1 center = Vector3D(0, 0, sin(t)),  radius = 0.5  (oscillates in z)
//   pass 2 center/radius fixed at subcell (1,1,1) — as in main15_gpu Setup():
//            cube.getSubcellCenter(1,1,1) / cube.getSubcellRadius(1,1,1)
// The cube stays pristine (identity) forever; each frame the composed sigma is
// computed from g_identityPositions into g_deformedPositions and uploaded.
//==============================================================================
static std::vector<float> g_identityPositions;    // pristine identity lattice (filled once)
static std::vector<float> g_deformedPositions;    // recomputed each frame, uploaded to VBO
static double g_animTime    = 0.0;                 // animation clock (advanced by timer)
static bool   g_animPaused  = false;
static const double g_timeStep    = 0.02;          // animTime advance per timer tick (~33ms)

// Interactive camera state (declared early so drawGPU can read it).
static float g_angleX = 20.0f, g_angleY = -30.0f;
static float g_zoom = 1.0f;
static float g_panX = 0.0f, g_panY = 0.0f;
static int   g_lastX = 0, g_lastY = 0;
static bool  g_leftDown = false, g_middleDown = false, g_rightDown = false;
static bool  g_showHelp = true;

//==============================================================================
// Car camera (main17). The car lives in IDENTITY space — it drives on the flat
// top face of the identity cube. Each frame the same composeSigma that deforms the
// mesh maps the car into world space, so it sticks to the blob's outer skin through
// every inflation/contraction. See display() for the mapping + Jacobian.
//==============================================================================
static bool   g_carMode    = true;     // ON by default; 'C' toggles back to orbit cam
static double g_carX       = 0.0;      // position on identity top face (y = g_idYmax)
static double g_carZ       = 0.0;
static double g_carHeading = 0.0;      // radians, measured in the X-Z plane
static double g_carSpeed   = 0.0;      // identity-space units / second
// identity-space bounding box of the lattice (filled once in initAnimatedBuffers)
static float  g_idXmin = 0, g_idXmax = 0, g_idYmax = 0, g_idZmin = 0, g_idZmax = 0;
// drive keys (set in keyboard(), cleared in keyboardUp()) — consumed by animTimer()
static bool   g_keyW = false, g_keyS = false, g_keyA = false, g_keyD = false;
// world-space eye computed each frame in display(), read by drawGPU() for uCamPos
static float  g_carEye[3] = {0.0f, 0.0f, 0.0f};

static void drawAxes(float length);
void Setup();
void initGL();
void reshape(int w, int h);
void display();
void mouseButton(int button, int state, int x, int y);
void mouseMotion(int x, int y);
void keyboard(unsigned char key, int x, int y);
void keyboardUp(unsigned char key, int x, int y);
void drawHUD();
void createUI();
void ProcessMenu(int value);

//==============================================================================
// Shaders: #version 460 core.
//  - VS passes world-space position (model matrix is identity: the cube is
//    already in world space, and uMVP already includes the view+projection).
//  - FS computes a flat per-triangle face normal from screen-space derivatives
//    of the world position (dFdx/dFdy), then does Blinn-Phong matching the
//    fixed-function light/material from initGL(). uMode=1 => black outline.
//==============================================================================
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
    "uniform vec3  uLightDir;   // world space (not necessarily normalized)\n"
    "uniform vec3  uCamPos;      // world space\n"
    "out vec4 fragColor;\n"
    "void main() {\n"
    "    if (uMode == 1) { fragColor = vec4(0.0, 0.0, 0.0, 1.0); return; }\n"
    "    vec3 N = normalize(cross(dFdx(vWorldPos), dFdy(vWorldPos)));\n"
    "    vec3 V = normalize(uCamPos - vWorldPos);\n"
    "    if (dot(N, V) < 0.0) N = -N;        // two-sided lighting\n"
    "    vec3 L = normalize(uLightDir);\n"
    "    vec3 H = normalize(L + V);\n"
    "    vec3 base = vec3(200.0 / 255.0);   // gray material (matches glColor3ub(200,200,200))\n"
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
        std::cerr << "[main16_gpu] " << (type == GL_VERTEX_SHADER ? "VERTEX" : "FRAGMENT")
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
    // Bind attribute 0 to aPos (matches glEnableVertexAttribArray(0) / layout).
    glBindAttribLocation(prog, 0, "aPos");
    glLinkProgram(prog);
    glDeleteShader(vs);
    glDeleteShader(fs);
    GLint ok = 0; glGetProgramiv(prog, GL_LINK_STATUS, &ok);
    if (!ok) {
        char log[4096] = {0}; glGetProgramInfoLog(prog, sizeof(log), nullptr, log);
        std::cerr << "[main16_gpu] program link failed:\n" << log << std::endl;
        glDeleteProgram(prog);
        return 0;
    }
    g_locMVP      = glGetUniformLocation(prog, "uMVP");
    g_locLightDir = glGetUniformLocation(prog, "uLightDir");
    g_locCamPos   = glGetUniformLocation(prog, "uCamPos");
    g_locMode     = glGetUniformLocation(prog, "uMode");
    return prog;
}

// ---- small column-major mat4 helpers (OpenGL convention) ----
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
// transform a point (w=1) by column-major M
static void xformPoint(const float* M, float x, float y, float z, float out[3]) {
    out[0] = M[0]*x + M[4]*y + M[8]*z  + M[12];
    out[1] = M[1]*x + M[5]*y + M[9]*z  + M[13];
    out[2] = M[2]*x + M[6]*y + M[10]*z + M[14];
}
// build a pure-rotation inverse-transpose-3x3 of the upper-left 3x3 of `view`,
// padded to a 4x4 (translation = 0). Used to map eye-space directions/points
// (the light direction and camera origin) back into world space.
static void buildInvRotation(const float* view, float invR[16]) {
    memset(invR, 0, sizeof(float) * 16);
    for (int c = 0; c < 3; ++c)
        for (int r = 0; r < 3; ++r)
            invR[c * 4 + r] = view[r * 4 + c];  // transpose of 3x3
    invR[15] = 1.0f;
}

//==============================================================================
// Fisheye projection (unchanged from main15)
//==============================================================================
void applyFisheyeProjection(int w, int h) {
    if (!g_fisheyeMode) {
        gluPerspective(50.0, (double)w / h, 0.001, 1000.0);
        return;
    }
    glLoadIdentity();
    float aspect = (float)w / h;
    float fisheyeFov = fminf(50.0f * (1.0f + g_fisheyeStrength * 2.5f), 170.0f);
    gluPerspective(fisheyeFov, aspect, 0.001, 1000.0);
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

//==============================================================================
// Global setup parameters
//==============================================================================
int cube_dim = 2.0;
int N = 37;
Cube cube(cube_dim, Vector3D{0, 0.0001, 0}, N);

double scal = 0.0;
double radi = 0.0;

// ii, jj, kk : module-action integers for the subcube lattice.
int ii = 4;
int jj = 5;
int kk = 9;

void applySigmaTransformationToCube(Cube& cube, const Vector3D& center, double radius) {
    if (!cube.hasSubdivision()) return;
    int n = cube.getSubdivisionLevels();
    int centerIdx = n / 2;
    for (int i = -centerIdx; i < centerIdx + 1; i++)
        for (int j = -centerIdx; j < centerIdx + 1; j++)
            for (int k = -centerIdx; k < centerIdx + 1; k++)
                for (int l = 0; l < 8; l++) {
                    Vector3D v_p = cube.getSubCell(i, j, k).vertices[l];
                    Vector3D new_pos = sigma(v_p, center, radius);
                    cube.updateSubCellVertex(i, j, k, l, new_pos);
                }
    cube.refreshTriangulation();
    std::cout << "  Triangulation refreshed for cube\n";
}

void Setup() {
    if (ciclo != 0) return;

    std::cout << "\n———————————————————————————————————————————————————————————————————————\n";
    std::cout <<   "|- CUBEs  (GPU build — NVIDIA T4, realtime inversion) —————————————————\n";
    std::cout <<   "———————————————————————————————————————————————————————————————————————\n\n";
    std::cout << "modules <ii, jj, kk> = <" << ii << ", " << jj << ", " << kk << ">\n";

    // One-time deformed STL export at t=0 (matches the old main15_gpu behavior), then
    // reset the cube to the pristine identity lattice so initAnimatedBuffers() captures
    // identity. The realtime loop never mutates the cube again.
    //   t=0  => pass-1 center = (0,0,sin 0) = (0,0,0), radius 0.5
    //   pass-2 = subcell (1,1,1) — same as main15_gpu Setup()
    applySigmaTransformationToCube(cube, Vector3D{0, 0, 0}, 0.5);
    applySigmaTransformationToCube(cube, cube.getSubcellCenter(1, 1, 1),
                                         cube.getSubcellRadius(1, 1, 1));
    cube.writeSTL_s("/home/mike666/Downloads/mesh_output_gpu.stl", "MyCube", "checkerboard", ii, kk, jj);
    cube.subdivide(N);   // reset to pristine identity lattice
    std::cout << "[main16_gpu] cube reset to identity lattice for realtime loop\n";
}

//==============================================================================
// GPU geometry buffers
//==============================================================================
static void initAnimatedBuffers() {
    // Capture the pristine identity lattice ONCE. The cube is never deformed during
    // the realtime loop, so this stays valid; each frame we sigma-transform a copy of
    // it into the VBO.
    cube.fillVertexLattice(g_identityPositions);  // n^3 * 8 * 3 floats
    g_deformedPositions.resize(g_identityPositions.size());
    glGenBuffers(1, &g_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, g_vbo);
    glBufferData(GL_ARRAY_BUFFER, g_identityPositions.size() * sizeof(float),
                 g_identityPositions.data(), GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    std::cout << "[main16_gpu] VBO (dynamic): " << (g_identityPositions.size() / 3)
              << " identity verts captured ("
              << (g_identityPositions.size() * sizeof(float) / (1024.0 * 1024.0))
              << " MB)\n";

    // Identity-space bounding box of the lattice. The car drives on the top face
    // (y = g_idYmax), clamped to the x/z extents, so the bounds must match the
    // actual cube placement (cube_dim / center) rather than being hard-coded.
    g_idXmin = g_idXmax = g_identityPositions[0];
    g_idYmax = g_identityPositions[1];
    g_idZmin = g_idZmax = g_identityPositions[2];
    float yMin = g_idYmax;
    for (size_t i = 0; i + 2 < g_identityPositions.size(); i += 3) {
        float x = g_identityPositions[i], y = g_identityPositions[i+1], z = g_identityPositions[i+2];
        if (x < g_idXmin) g_idXmin = x; if (x > g_idXmax) g_idXmax = x;
        if (y < yMin)     yMin     = y; if (y > g_idYmax) g_idYmax = y;
        if (z < g_idZmin) g_idZmin = z; if (z > g_idZmax) g_idZmax = z;
    }
    std::cout << "[main17_gpu] identity bounds X[" << g_idXmin << "," << g_idXmax
              << "] Y[" << yMin << "," << g_idYmax << "] Z[" << g_idZmin << "," << g_idZmax
              << "]  car drives on top face y=" << g_idYmax << "\n";
}

static void rebuildIndexBuffer() {
    std::vector<unsigned int> indices;
    cube.fillCheckerboardIndices(ii, kk, jj, indices);  // same call order as getCheckerboardFacets(ii,kk,jj)
    g_indexCount = (GLsizei)indices.size();
    if (g_ibo == 0) glGenBuffers(1, &g_ibo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int),
                 indices.data(), GL_DYNAMIC_DRAW);
    g_selectionDirty = false;
    std::cout << "[main16_gpu] IBO rebuilt: " << (g_indexCount / 3)
              << " triangles (ii,jj,kk=" << ii << "," << jj << "," << kk << ")\n";
}

//==============================================================================
// Composed sphere-inversion (the two Setup passes, now animated).
// Inlined with a singular-center guard: the free sigma() in Vector3D.cpp throws
// when a vertex coincides with the sphere center, which would crash the render
// loop. Here we just leave such a vertex untouched.
//==============================================================================
static inline Vector3D composeSigma(const Vector3D& p,
                                    const Vector3D& c1, double r1,
                                    const Vector3D& c2, double r2) {
    Vector3D d1 = p - c1;  double ds1 = d1 * d1;
    Vector3D q  = (ds1 < 1e-12) ? p : c1 + (r1 * r1 / ds1) * d1;   // pass 1 (was line 306)
    Vector3D d2 = q - c2;  double ds2 = d2 * d2;
    return (ds2 < 1e-12) ? q : c2 + (r2 * r2 / ds2) * d2;         // pass 2 (was line 307)
}

//==============================================================================
// Shared sigma parameters for the animated frame. Both the mesh update
// (updateAnimatedGeometry) and the car-camera mapping (display) call this so they
// use IDENTICAL inversion centres/radii each frame — that is what makes the car
// stick to the mesh exactly as it inflates/contracts.
//==============================================================================
struct SigmaParams { Vector3D c1; double r1; Vector3D c2; double r2; };
static SigmaParams currentSigmaParams(double t) {
    SigmaParams s;
    s.c1 = Vector3D(0.0, 0.0, 0.1 * sin(t));
    s.r1 = 0.5;
    s.c2 = cube.getSubcellCenter(1, 1, 1);
    s.r2 = cube.getSubcellRadius(1, 1, 1);
    return s;
}
static inline Vector3D composeSigmaAt(const Vector3D& p, const SigmaParams& s) {
    return composeSigma(p, s.c1, s.r1, s.c2, s.r2);
}

//==============================================================================
// Per-frame lattice computation + VBO upload.
// Reads the pristine identity lattice, applies the two inversions, and streams
// the result into the VBO. The cube object is never mutated, so
// getSubcellCenter/getSubcellRadius always return stable identity values.
//==============================================================================
static void updateAnimatedGeometry() {
    // Same two animated inversions as main16, now via the shared helper so the car
    // camera (display()) sees the identical deformation this frame.
    SigmaParams sp = currentSigmaParams(g_animTime);

    const size_t n = g_identityPositions.size();
    for (size_t i = 0; i + 2 < n; i += 3) {
        Vector3D p(g_identityPositions[i], g_identityPositions[i+1], g_identityPositions[i+2]);
        Vector3D q = composeSigmaAt(p, sp);
        g_deformedPositions[i]   = (float)q.x();
        g_deformedPositions[i+1] = (float)q.y();
        g_deformedPositions[i+2] = (float)q.z();
    }
    glBindBuffer(GL_ARRAY_BUFFER, g_vbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, (GLsizeiptr)(n * sizeof(float)),
                    g_deformedPositions.data());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

//==============================================================================
// DRAW — replaces ProcessingProto()/Draw()/drawFacetMainCStyle() loop.
//==============================================================================
static void drawGPU() {
    if (!g_program || !g_vbo) return;
    if (g_selectionDirty) rebuildIndexBuffer();
    if (g_indexCount == 0) return;

    // uMVP = projection * modelview (read back from the fixed-function stack,
    // which display() has already set exactly as main15 did).
    float proj[16], view[16], mvp[16];
    glGetFloatv(GL_PROJECTION_MATRIX, proj);
    glGetFloatv(GL_MODELVIEW_MATRIX, view);
    matMul(proj, view, mvp);

    // Camera world position and light direction (light is fixed in eye space,
    // set with identity modelview in initGL — see comment in main).
    float invR[16]; buildInvRotation(view, invR);
    float camPos[3];
    if (g_carMode) {
        // display() already computed the car eye in world space.
        camPos[0] = g_carEye[0]; camPos[1] = g_carEye[1]; camPos[2] = g_carEye[2];
    } else {
        // Orbit camera: recover the world-space eye from the orbit modelview.
        xformPoint(invR, -g_panX, -g_panY, 5.0f * g_zoom, camPos);
    }
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

    // Outline pass — black triangle edges (glPolygonMode LINE on GL_TRIANGLES
    // draws each triangle's 3 edges, matching the original GL_LINE_LOOP).
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

void drawAxes(float length) {
    glLineWidth(2.0f);
    glBegin(GL_LINES);
      glColor3f(1, 0, 0); glVertex3f(-length, 0, 0); glVertex3f(length, 0, 0);
      glColor3f(0, 1, 0); glVertex3f(0, -length, 0); glVertex3f(0, length, 0);
      glColor3f(0, 0, 1); glVertex3f(0, 0, -length); glVertex3f(0, 0, length);
    glEnd();
}

//////////////////////////////////////
//  OPENGL AS BACKGROUND
//////////////////////////////////////

//-----------------------------------------------------------------------------
// Animation timer — advances the realtime inversion clock and requests a redraw.
// Re-registers itself for a ~30 Hz tick. While paused, t holds so the geometry
// freezes (display() also skips the lattice update when paused).
//-----------------------------------------------------------------------------
static void animTimer(int /*value*/) {
    if (!g_animPaused) g_animTime += g_timeStep;

    // Car physics, integrated in IDENTITY space (the flat top face). Time-based via
    // g_timeStep so behaviour is independent of the frame rate. The world-space
    // mapping happens later in display() via composeSigma, so this stays simple.
    if (g_carMode && !g_animPaused) {
        const double accel = 0.6, drag = 0.9, maxSpeed = 0.6, turnRate = 1.6;
        double throttle = (g_keyW ? 1.0 : 0.0) - (g_keyS ? 1.0 : 0.0);
        double steer    = (g_keyA ? 1.0 : 0.0) - (g_keyD ? 1.0 : 0.0);
        g_carSpeed += throttle * accel * g_timeStep;
        g_carSpeed *= drag;                                   // friction / natural stop
        g_carSpeed = fmax(-maxSpeed, fmin(maxSpeed, g_carSpeed));
        if (fabs(g_carSpeed) > 1e-4)                          // only steer while moving
            g_carHeading += steer * turnRate * g_timeStep * (g_carSpeed > 0 ? 1.0 : -1.0);
        g_carX += g_carSpeed * cos(g_carHeading) * g_timeStep;
        g_carZ += g_carSpeed * sin(g_carHeading) * g_timeStep;
        // clamp to the top face; zero speed on hit so the car "bumps" the edge
        if (g_carX < g_idXmin) { g_carX = g_idXmin; g_carSpeed = 0; }
        if (g_carX > g_idXmax) { g_carX = g_idXmax; g_carSpeed = 0; }
        if (g_carZ < g_idZmin) { g_carZ = g_idZmin; g_carSpeed = 0; }
        if (g_carZ > g_idZmax) { g_carZ = g_idZmax; g_carSpeed = 0; }
    }

    glutPostRedisplay();
    glutTimerFunc(33, animTimer, 0);
}

//-----------------------------------------------------------------------------
// Main
//-----------------------------------------------------------------------------
int main(int argc, char** argv) {
    // --- Force the NVIDIA GLX backend BEFORE creating the GL context ---
    // Without this the app falls back to Mesa llvmpipe (CPU software raster).
    setenv("__GLX_VENDOR_LIBRARY_NAME", "nvidia", 1);

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(720, 720);
    glutCreateWindow(" JAZ 4D GPU (Tesla T4) — realtime inversion ");

    loadGL();

    // Confirm which renderer we actually got.
    std::cout << "[main16_gpu] GL VENDOR   = " << glGetString(GL_VENDOR)   << "\n";
    std::cout << "[main16_gpu] GL RENDERER = " << glGetString(GL_RENDERER) << "\n";
    std::cout << "[main16_gpu] GL VERSION  = " << glGetString(GL_VERSION)  << "\n";
    std::string rendererStr = (const char*)glGetString(GL_RENDERER);
    if (rendererStr.find("llvmpipe") != std::string::npos) {
        std::cerr << "[main16_gpu] WARNING: still on llvmpipe (software). "
                     "Run via ./run_gpu.sh or set __GLX_VENDOR_LIBRARY_NAME=nvidia.\n";
    }

    ProcessMenu(1);   // smoothing/blending (same as main15)
    initGL();
    createUI();
    Setup();          // one-time: t=0 sigma + STL, then reset cube to identity

    g_program = buildProgram();
    if (!g_program) {
        std::cerr << "[main16_gpu] shader program build failed; aborting.\n";
        return 1;
    }
    initAnimatedBuffers();   // capture identity lattice once; create dynamic VBO

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouseButton);
    glutMotionFunc(mouseMotion);
    glutKeyboardFunc(keyboard);
    glutKeyboardUpFunc(keyboardUp);
    glutTimerFunc(33, animTimer, 0);   // ~30 Hz animation clock

    glutMainLoop();
    return 0;
}

//-----------------------------------------------------------------------------
// Setup OpenGL lights, materials, and default projection
//-----------------------------------------------------------------------------
void initGL() {
    GLfloat ambient[]  = {0.3f, 0.3f, 0.3f, 1.0f};
    GLfloat diffuse[]  = {0.7f, 0.7f, 0.7f, 1.0f};
    GLfloat specular[] = {1.0f, 1.0f, 1.0f, 1.0f};
    GLfloat position[] = {1.0f, 1.0f, 0.25f, 0.0f};  // directional, set with identity modelview => eye space

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glLightfv(GL_LIGHT0, GL_AMBIENT,  ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE,  diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
    glLightfv(GL_LIGHT0, GL_POSITION, position);

    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
    GLfloat specref[] = {1.0f, 1.0f, 1.0f, 1.0f};
    GLfloat shine[]   = {128.0f};
    glMaterialfv(GL_FRONT, GL_SPECULAR, specref);
    glMaterialfv(GL_FRONT, GL_SHININESS, shine);

    glEnable(GL_DEPTH_TEST);
    glFrontFace(GL_CCW);
    glEnable(GL_NORMALIZE);
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
}

//-----------------------------------------------------------------------------
// Handle window size changes with fisheye support (same as main15)
//-----------------------------------------------------------------------------
void reshape(int w, int h) {
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    applyFisheyeProjection(w, h);
    glMatrixMode(GL_MODELVIEW);
}

//-----------------------------------------------------------------------------
// Draw help overlay
//-----------------------------------------------------------------------------
void drawHUD() {
    if (!g_showHelp) return;
    const char* lines[] = {
        "[C]: Toggle Car / Orbit camera",
        "Car: W/S throttle | A/D steer (steer while moving)",
        "L-drag: Rotate (orbit cam)",
        "M-drag: Pan (orbit cam)",
        "R-drag/Wheel: Zoom (orbit cam)",
        "[H]: Toggle Help",
        "[G]: Toggle Glass Mode",
        "[F]: Toggle Fisheye Mode",
        "[]: Fisheye Strength -/+",
        "[i/I j/J k/K]: Change <ii,jj,kk> (GPU rebuild)",
        "[Space]: Pause/Resume inversion animation",
        "Right-click: UI Menu",
        "[Esc]: Quit"
    };

    glMatrixMode(GL_PROJECTION); glPushMatrix();
    glLoadIdentity(); glOrtho(0, 1, 0, 1, -1, 1);
    glMatrixMode(GL_MODELVIEW); glPushMatrix(); glLoadIdentity();

    glColor3f(0, 0, 0);
    float y = 0.95f;
    for (auto& ln : lines) {
        glRasterPos2f(0.02f, y);
        for (const char* c = ln; *c; ++c)
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);
        y -= 0.05f;
    }

    // status line
    y -= 0.05f;
    char status[160];
    sprintf(status, "<ii,jj,kk> = <%d,%d,%d>  tris=%d",
            ii, jj, kk, (int)(g_indexCount / 3));
    glColor3f(0.8f, 0.2f, 0.2f);
    glRasterPos2f(0.02f, y);
    for (const char* c = status; *c; ++c)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);

    // animation status line — pass-1 oscillation z + run state
    y -= 0.05f;
    char anim[160];
    sprintf(anim, "inv: center=(0,0,%.2f)  pass2=subcell(1,1,1)  %s",
            0.1 * sin(g_animTime), g_animPaused ? "PAUSED" : "RUNNING");
    glColor3f(0.2f, 0.6f, 0.2f);
    glRasterPos2f(0.02f, y);
    for (const char* c = anim; *c; ++c)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);

    if (g_fisheyeMode) {
        y -= 0.05f;
        char fs[100];
        sprintf(fs, "Fisheye: ON (Strength: %.1f)", g_fisheyeStrength);
        glColor3f(0.2f, 0.2f, 0.8f);
        glRasterPos2f(0.02f, y);
        for (const char* c = fs; *c; ++c)
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);
    }

    // car status line — speed / heading / identity-space position
    if (g_carMode) {
        y -= 0.05f;
        char car[180];
        sprintf(car, "CAR  speed=%.2f  heading=%.0f deg  pos=(%.2f, %.2f) on skin y=%.2f",
                g_carSpeed, g_carHeading * 180.0 / M_PI, g_carX, g_carZ, g_idYmax);
        glColor3f(0.8f, 0.4f, 0.0f);
        glRasterPos2f(0.02f, y);
        for (const char* c = car; *c; ++c)
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);
    }

    glMatrixMode(GL_MODELVIEW); glPopMatrix();
    glMatrixMode(GL_PROJECTION); glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
}

//-----------------------------------------------------------------------------
// Main display
//-----------------------------------------------------------------------------
void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Recompute the animated inversion lattice and stream it into the VBO before
    // drawing. When paused, the last computed frame is kept (we skip the update so
    // the frozen geometry stays exactly as-is, and the timer stops advancing t).
    if (!g_animPaused) updateAnimatedGeometry();

    // Fixed-function modelview. This is read back to build the shader's uMVP, so
    // whatever camera we set here (car or orbit) lands on screen with no shader
    // changes. The HUD is drawn in its own ortho projection, unaffected.
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    if (g_carMode) {
        // --- Car camera (cockpit/first-person on the blob's outer skin) ---
        // The car lives in identity space on the top face. Map it to world space
        // with the SAME composeSigma that deformed the mesh this frame, and lift
        // its orientation through the deformation's Jacobian (finite differences).
        SigmaParams sp = currentSigmaParams(g_animTime);
        Vector3D pId(g_carX, (double)g_idYmax, g_carZ);
        Vector3D pW = composeSigmaAt(pId, sp);                 // world point on the skin

        const double e = 1e-4;
        Vector3D jx = (composeSigmaAt(pId + Vector3D(e, 0, 0), sp) - pW) / e;
        Vector3D jy = (composeSigmaAt(pId + Vector3D(0, e, 0), sp) - pW) / e;
        Vector3D jz = (composeSigmaAt(pId + Vector3D(0, 0, e), sp) - pW) / e;

        // identity up = +Y (top-face normal) -> world skin normal via the Jacobian.
        Vector3D worldUp = unit(jy);
        // identity forward = heading in the X-Z plane -> world, then re-orthogonalise
        // against the skin normal so it stays tangent to the surface.
        Vector3D fwdId(cos(g_carHeading), 0.0, sin(g_carHeading));
        Vector3D worldFwd = fwdId.x() * jx + fwdId.z() * jz;   // J * fwdId (y term ~ 0)
        worldFwd = worldFwd - (worldFwd * worldUp) * worldUp;   // (* = dot; scalar*vec)
        worldFwd = unit(worldFwd);

        Vector3D eye = pW + 0.01 * worldUp;                    // cockpit: on the skin
        g_carEye[0] = (float)eye.x(); g_carEye[1] = (float)eye.y(); g_carEye[2] = (float)eye.z();
        Vector3D aim = eye + worldFwd;
        gluLookAt(eye.x(), eye.y(), eye.z(),
                  aim.x(), aim.y(), aim.z(),
                  worldUp.x(), worldUp.y(), worldUp.z());
    } else {
        // --- Orbit camera (original main15/main16 behaviour) ---
        glTranslatef(g_panX, g_panY, -5.0f * (g_zoom / g_surfaceZoomFactor));
        glRotatef(g_angleX, 1, 0, 0);
        glRotatef(g_angleY, 0, 1, 0);
    }

    //drawAxes(10.0f);
    drawGPU();   // shader + VBO/IBO replaces the ~264K glBegin/glEnd loop
    drawHUD();

    glutSwapBuffers();
}

//-----------------------------------------------------------------------------
// Mouse (same as main15)
//-----------------------------------------------------------------------------
void mouseButton(int button, int state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON)        g_leftDown   = (state == GLUT_DOWN);
    else if (button == GLUT_MIDDLE_BUTTON) g_middleDown = (state == GLUT_DOWN);
    else if (button == GLUT_RIGHT_BUTTON)  g_rightDown  = (state == GLUT_DOWN);
    else if (button == 3) { g_zoom *= 1.05f; glutPostRedisplay(); }
    else if (button == 4) { g_zoom /= 1.05f; glutPostRedisplay(); }
    g_lastX = x; g_lastY = y;
}

void mouseMotion(int x, int y) {
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

//-----------------------------------------------------------------------------
// Keyboard — same as main15, but i/I/j/J/k/K now mark the selection dirty so the
// IBO is rebuilt once on the next frame instead of every frame.
//-----------------------------------------------------------------------------
void keyboard(unsigned char key, int x, int y) {
    switch (key) {
        case ' ':
            g_animPaused = !g_animPaused;
            printf("Animation: %s\n", g_animPaused ? "PAUSED" : "RUNNING");
            glutPostRedisplay();
            break;
        case 'g': case 'G':
            g_glassMode = !g_glassMode;
            printf("Glass mode: %s\n", g_glassMode ? "ON" : "OFF");
            glutPostRedisplay();
            break;
        case 'h': case 'H':
            g_showHelp = !g_showHelp;
            glutPostRedisplay();
            break;
        case 'f': case 'F':
            g_fisheyeMode = !g_fisheyeMode;
            printf("Fisheye mode: %s\n", g_fisheyeMode ? "ON" : "OFF");
            reshape(glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
            glutPostRedisplay();
            break;
        case '[':
            g_fisheyeStrength = fmaxf(0.0f, g_fisheyeStrength - 0.1f);
            printf("Fisheye strength: %.1f\n", g_fisheyeStrength);
            if (g_fisheyeMode) reshape(glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
            glutPostRedisplay();
            break;
        case ']':
            g_fisheyeStrength = fminf(2.0f, g_fisheyeStrength + 0.1f);
            printf("Fisheye strength: %.1f\n", g_fisheyeStrength);
            if (g_fisheyeMode) reshape(glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
            glutPostRedisplay();
            break;

        case 'k': kk += 1; g_selectionDirty = true; std::cout << "modules <ii, jj, kk> = <" << ii << ", " << jj << ", " << kk << ">\n"; glutPostRedisplay(); break;
        case 'K': kk -= 1; g_selectionDirty = true; std::cout << "modules <ii, jj, kk> = <" << ii << ", " << jj << ", " << kk << ">\n"; glutPostRedisplay(); break;
        case 'i': ii += 1; g_selectionDirty = true; std::cout << "modules <ii, jj, kk> = <" << ii << ", " << jj << ", " << kk << ">\n"; glutPostRedisplay(); break;
        case 'I': ii -= 1; g_selectionDirty = true; std::cout << "modules <ii, jj, kk> = <" << ii << ", " << jj << ", " << kk << ">\n"; glutPostRedisplay(); break;
        case 'j': jj += 1; g_selectionDirty = true; std::cout << "modules <ii, jj, kk> = <" << ii << ", " << jj << ", " << kk << ">\n"; glutPostRedisplay(); break;
        case 'J': jj -= 1; g_selectionDirty = true; std::cout << "modules <ii, jj, kk> = <" << ii << ", " << jj << ", " << kk << ">\n"; glutPostRedisplay(); break;

        // --- Car camera (main17) ---
        case 'c': case 'C':
            g_carMode = !g_carMode;
            if (g_carMode) {  // entering car mode: re-centre on the top face, at rest
                g_carX = 0.5 * (g_idXmin + g_idXmax);
                g_carZ = 0.5 * (g_idZmin + g_idZmax);
                g_carHeading = 0.0;
                g_carSpeed   = 0.0;
            }
            printf("Car camera: %s\n", g_carMode ? "ON (cockpit on skin)" : "OFF (orbit)");
            glutPostRedisplay();
            break;
        case 'w': case 'W': g_keyW = true; break;
        case 's': case 'S': g_keyS = true; break;
        case 'a': case 'A': g_keyA = true; break;
        case 'd': case 'D': g_keyD = true; break;

        case 27: exit(0); break;
    }
}

//-----------------------------------------------------------------------------
// Key release — clears the held drive keys so throttle/steer stop cleanly.
//-----------------------------------------------------------------------------
void keyboardUp(unsigned char key, int x, int y) {
    (void)x; (void)y;
    switch (key) {
        case 'w': case 'W': g_keyW = false; break;
        case 's': case 'S': g_keyS = false; break;
        case 'a': case 'A': g_keyA = false; break;
        case 'd': case 'D': g_keyD = false; break;
    }
}

//-----------------------------------------------------------------------------
// UI menu (same as main15)
//-----------------------------------------------------------------------------
void MenuHandler(int choice) {
    switch (choice) {
        case 1: g_showHelp = !g_showHelp; break;
        case 2: g_angleX = 20; g_angleY = -30; g_zoom = 1; g_panX = g_panY = 0;
                g_fisheyeMode = false; g_fisheyeStrength = 1.0f;
                g_enableSurfaceZoom = true; break;
        case 3: exit(0); break;
    }
    glutPostRedisplay();
}

void createUI() {
    int menu = glutCreateMenu(MenuHandler);
    glutAddMenuEntry("Toggle Help", 1);
    glutAddMenuEntry("Reset View",  2);
    glutAddMenuEntry("Quit",        3);
    glutAttachMenu(GLUT_RIGHT_BUTTON);
}

void ProcessMenu(int value) {
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
