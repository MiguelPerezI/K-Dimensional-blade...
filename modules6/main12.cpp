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

Vector3D d_ui(1.0, 0.0, 0.0);
Vector3D e_ui(0.0, 1.0, 0.0);
Vector3D f_ui(0.0, 0.0, 1.0);
Facet f_1(d_ui, e_ui, f_ui);

// ==================== CUBE INSTANCES ====================
// Basic cube: radius=1.0 (side length=2.0), positioned right of center
// This creates a simple cube with 12 triangular faces (2 triangles per face × 6 faces)
Cube cube1(1.0, Vector3D{2.5,0,0});

// Subdivided cube: radius=0.8, positioned left of center, with 2 subdivision levels
// This creates 2³=8 subcubes, each triangulated with regular 12 triangle mesh
// Total faces: 8 subcubes × 12 triangular faces = 96 triangles
Cube cube2(1.6, Vector3D{0,0,0}, 8);

// Demonstration cube for advanced subdivision features
Cube cube3(1.6, Vector3D{0,0,0}, 8);  // 3³=27 subcells for detailed examples
Cube cube4(1.6, Vector3D{0,0,0}, 8);
Cube cube5(1.6, Vector3D{0,0,0}, 8);
Cube cube6(1.6, Vector3D{0,0,0}, 8);
Cube cube7(1.6, Vector3D{0,0,0}, 8);

// ==================== FACET STORAGE ====================
// Storage for cube triangle collections
FacetBox cube_facets_1;  // Basic cube triangulation
FacetBox cube_facets_2;  // Subdivided cube triangulation
FacetBox cube_facets_3;  // Advanced subdivision demos

// Advanced subdivision control demos
FacetBox single_subcell;     // Single subcell rendering demo
FacetBox plane_subcells;     // Plane rendering demo
FacetBox plane_subcells_2;
FacetBox selected_subcells;  // Selective rendering demo

double ra = 0.0;
Vector3D ce = Vector3D{0,0,0};

//==============================================================================
// OPTION A: Real-Time Multi-Reflection System
//==============================================================================
class ReflectionManager {
public:
    struct ReflectionPoint {
        Vector3D center;
        double radius;
        bool active = true;

        ReflectionPoint(const Vector3D& c, double r) : center(c), radius(r) {}
    };

    std::vector<ReflectionPoint> reflection_points;

    void addReflectionPoint(const Vector3D& center, double radius) {
        reflection_points.emplace_back(center, radius);
    }

    void clearReflectionPoints() {
        reflection_points.clear();
    }

    // Generate reflected cube facets on-the-fly for rendering
    FacetBox getReflectedCubeFacets(const Cube& original, int reflection_index) {
        if (reflection_index < 0 || reflection_index >= reflection_points.size()) {
            return FacetBox(); // Return empty
        }

        const auto& rp = reflection_points[reflection_index];
        if (!rp.active) {
            return FacetBox(); // Return empty if disabled
        }

        FacetBox reflected_facets;
        const FacetBox& original_facets = original.getFacets();

        try {
            for (size_t i = 0; i < original_facets.size(); ++i) {
                const Facet& face = original_facets[i];

                // Reflect each vertex of the triangle
                Vector3D v0 = sigma(face[0], rp.center, rp.radius);
                Vector3D v1 = sigma(face[1], rp.center, rp.radius);
                Vector3D v2 = sigma(face[2], rp.center, rp.radius);

                reflected_facets.push(v0, v1, v2);
            }
        } catch (const std::runtime_error& e) {
            // Skip reflections that fail (e.g., point at sphere center)
            std::cerr << "Reflection warning: " << e.what() << std::endl;
        }

        return reflected_facets;
    }

    void setReflectionActive(int index, bool active) {
        if (index >= 0 && index < reflection_points.size()) {
            reflection_points[index].active = active;
        }
    }
};

//==============================================================================
// Simple Cube Reflection Function
//==============================================================================

/**
 * @brief Reflect a subdivided cube using checkboard pattern (i+j+k)%2
 * @param source_cube The cube to reflect
 * @param reflection_center Center of the reflection sphere
 * @param reflection_radius Radius of the reflection sphere
 * @returns FacetBox containing reflected triangles from checkboard subcells
 */
void reflectCubeWithCheckboard(Cube& source_cube, const Vector3D& reflection_center, double reflection_radius) {
    FacetBox reflected_facets;

    if (!source_cube.hasSubdivision()) {

        int n =source_cube.getSubdivisionLevels();
        for (int i = 0; i < n; i++) 
        for (int j = 0; j < n; j++) 
        for (int k = 0; k < n; k++) {
        
            for (int l = 0; l < 8; l++) {
                Vector3D v_p = source_cube.getSubCell(i, j, k).vertices[l];
        
                // Move it upward by 0.2 units
                Vector3D new_pos = sigma(v_p, reflection_center, reflection_radius);
                source_cube.updateSubCellVertex(i, j, k, l, new_pos);
            }
        
            //cout << "  Moved vertex  of subcell (0,0,0) from " << original_pos << " to " << new_pos << "\n";
        }
        // Refresh triangulation
        source_cube.refreshTriangulation();
        //cube_facets_3 = cube3.getFacets();  // Update facets
        //plane_subcells = cube3.getPlaneFacets(2, 4);
        //cout << "  Triangulation refreshed\n";
    }
}

// Global instances
ReflectionManager g_reflectionManager;

// Global reflection parameters for simple approach
Vector3D g_reflection_center{0, 0, 0};
double g_reflection_radius = 1.0;
bool g_show_reflection = true;

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
        cout << "                                     ~[Cube Class Example]\n\n";

        // ==================== CUBE SETUP ====================
        // Extract triangular facets from cube instances
        // cube1: Basic cube with 12 triangles
        cube_facets_1 = cube1.getFacets();

        // cube2: Subdivided cube with 96 triangles
        cube_facets_2 = cube2.getFacets();

        // cube3: Advanced subdivision demos (3³ = 27 subcells)
        cube_facets_3 = cube3.getFacets();

        cout << "Cube 1 (basic): " << cube_facets_1.size() << " triangular faces\n";
        cout << "Cube 2 (subdivided): " << cube_facets_2.size() << " triangular faces\n";
        cout << "Cube 3 (advanced): " << cube_facets_3.size() << " triangular faces\n";

        // ==================== ADVANCED SUBDIVISION DEMOS ====================
        cout << "\n--- Advanced Subdivision Control Demos ---\n";

        // 1. Individual subcell access demo
        cout << "Demo 1: Individual Subcell Access\n";
        if (cube2.hasSubdivision()) {
            const Cube::SubCell& center_cell = cube2.getSubCell(0, 0, 0);  // Center subcell of 8x8x8
            cout << "  Center subcell (1,1,1) position: " << center_cell.center << "\n";
            cout << "  Active status: " << (center_cell.active ? "enabled" : "disabled") << "\n";

            // Extract facets for this specific subcell
            single_subcell = cube2.getSubCellFacets(0, 0, 0);
            cout << "  Subcell triangles: " << single_subcell.size() << "\n";
        }

        // 2. Plane access demo
        cout << "Demo 2: Plane Access\n";
        if (cube2.hasSubdivision()) {
            auto xy_plane = cube2.getPlane(2, 4);  // XY plane at Z=1
            cout << "  XY plane at Z=1 contains " << xy_plane.size() << " subcells\n";

            // Extract all facets for this plane
            plane_subcells_2 = cube2.getPlaneFacets(2, 4);
            cout << "  Plane triangles: " << plane_subcells_2.size() << "\n";
        }

        // 3. Selective rendering demo - create "hollow frame"
        cout << "Demo 3: Selective Rendering (Hollow Frame)\n";
        if (cube3.hasSubdivision()) {
            int n = cube3.getSubdivisionLevels();
            int disabled_count = 0;

            // Disable interior subcells to create hollow effect
            int center = n / 2;
            for (int i = -center; i < center; i++) {
                for (int j = -center; j < center; j++) {
                    for (int k = -center; k < center; k++) {

                        if ((i + j + k) % 2 == 0) {
                            cube3.setSubCellActive(i, j, k, false);
                            cube2.setSubCellActive(i, j, k, false);
                            disabled_count++;
                       }
                    }
                }
            }

            cout << "  Disabled " << disabled_count << " interior subcells\n";
            cout << "  Creating hollow frame effect\n";

            // Refresh triangulation after modifications
            cube3.refreshTriangulation();
            cube2.refreshTriangulation();
            cube_facets_3 = cube3.getFacets();  // Update facets
            cube_facets_2 = cube2.getFacets();  // Update facets
            cout << "  Frame triangles: " << cube_facets_3.size() << "\n";
        }

        // 4. Real-time vertex manipulation demo
        cout << "Demo 4: Vertex Manipulation\n";
        if (cube3.hasSubdivision()) {
            int n = cube3.getSubdivisionLevels();
            int center = n / 2;
            double ra = cube2.getSubcellRadius(0, 0, 0);
            Vector3D ce = cube2.getSubcellCenter(0, 0, 0);
            // Get original position of vertex (l) in subcell (i,j,k)
            for (int i = -center; i < center; i++) 
            for (int j = -center; j < center; j++) 
            for (int k = -center; k < center; k++) {

                for (int l = 0; l < 8; l++) {
                    Vector3D v_p = cube3.getSubCell(i, j, k).vertices[l];

                    // Move it upward by 0.2 units
                    Vector3D new_pos = sigma(v_p, ce, ra);
                    cube3.updateSubCellVertex(i, j, k, l, new_pos);
                }
            
                //cout << "  Moved vertex  of subcell (0,0,0) from " << original_pos << " to " << new_pos << "\n";
            }
            // Refresh triangulation
            cube3.refreshTriangulation();                                                                                                                                               
            //cube_facets_3 = cube3.getFacets();  // Update facets
            //plane_subcells = cube3.getPlaneFacets(2, 4);
            cout << "  Triangulation refreshed\n";
        }




        if (cube4.hasSubdivision()) {                                                                                                                                                   
            int n = cube4.getSubdivisionLevels();
            int center = n / 2;
            double ra = cube2.getSubcellRadius(0, 0, 1);
            Vector3D ce = cube2.getSubcellCenter(0, 0, 1);
            // Get original position of vertex (l) in subcell (i,j,k)
            for (int i = -center; i < center; i++) 
            for (int j = -center; j < center; j++) 
            for (int k = -center; k < center; k++) {

                for (int l = 0; l < 8; l++) {
                    Vector3D v_p = cube4.getSubCell(i, j, k).vertices[l];

                    // Move it upward by 0.2 units
                    Vector3D new_pos = sigma(v_p, ce, ra);
                    cube4.updateSubCellVertex(i, j, k, l, new_pos);
                }
            
                //cout << "  Moved vertex  of subcell (0,0,0) from " << original_pos << " to " << new_pos << "\n";
            }
            // Refresh triangulation
            cube4.refreshTriangulation();                                                                                                                                               
            cout << "  Triangulation refreshed\n";
        }

        


        if (cube5.hasSubdivision()) {
            int n = cube5.getSubdivisionLevels();
            int center = n / 2;
            double ra = cube2.getSubcellRadius(0, 1, 0);
            Vector3D ce = cube2.getSubcellCenter(0, 1, 0);
            // Get original position of vertex (l) in subcell (i,j,k)
            for (int i = -center; i < center; i++)
            for (int j = -center; j < center; j++)
            for (int k = -center; k < center; k++) {

                for (int l = 0; l < 8; l++) {
                    Vector3D v_p = cube5.getSubCell(i, j, k).vertices[l];

                    // Move it upward by 0.2 units
                    Vector3D new_pos = sigma(v_p, ce, ra);
                    cube5.updateSubCellVertex(i, j, k, l, new_pos);
                }

                //cout << "  Moved vertex  of subcell (0,0,0) from " << original_pos << " to " << new_pos << "\n";
            }
            // Refresh triangulation
            cube5.refreshTriangulation();
            cout << "  Triangulation refreshed\n";
        }



        if (cube6.hasSubdivision()) {                                                                                                                                                   
            int n = cube6.getSubdivisionLevels();
            int center = n / 2;
            double ra = cube2.getSubcellRadius(1, 0, 0);
            Vector3D ce = cube2.getSubcellCenter(1, 0, 0);
            // Get original position of vertex (l) in subcell (i,j,k)
            for (int i = -center; i < center; i++) 
            for (int j = -center; j < center; j++) 
            for (int k = -center; k < center; k++) {

                for (int l = 0; l < 8; l++) {
                    Vector3D v_p = cube6.getSubCell(i, j, k).vertices[l];

                    // Move it upward by 0.2 units
                    Vector3D new_pos = sigma(v_p, ce, ra);
                    cube6.updateSubCellVertex(i, j, k, l, new_pos);
                }
            
                //cout << "  Moved vertex  of subcell (0,0,0) from " << original_pos << " to " << new_pos << "\n";
            }
            // Refresh triangulation
            cube6.refreshTriangulation();                                                                                                                                               
            cout << "  Triangulation refreshed\n";
        }



        if (cube7.hasSubdivision()) {
            int n = cube7.getSubdivisionLevels();
            int center = n / 2;
            double ra = cube2.getSubcellRadius(1, 1, 1);
            Vector3D ce = cube2.getSubcellCenter(1, 1, 1);
            // Get original position of vertex (l) in subcell (i,j,k)
            for (int i = -center; i < center; i++)
            for (int j = -center; j < center; j++)
            for (int k = -center; k < center; k++) {

                for (int l = 0; l < 8; l++) {
                    Vector3D v_p = cube7.getSubCell(i, j, k).vertices[l];

                    // Move it upward by 0.2 units
                    Vector3D new_pos = sigma(v_p, ce, ra);
                    cube7.updateSubCellVertex(i, j, k, l, new_pos);
                }

                //cout << "  Moved vertex  of subcell (0,0,0) from " << original_pos << " to " << new_pos << "\n";
            }
            // Refresh triangulation
            cube7.refreshTriangulation();
            cout << "  Triangulation refreshed\n";
        }



        // ==================== NEW STRING-BASED getPlaneFacets DEMO ====================
        cout << "Demo 5: New String-based getPlaneFacets Method\n";
        if (cube2.hasSubdivision()) {
            try {
                // Test the three examples mentioned in the request

                // getPlaneFacets("x", "y", "0") - center XY plane (perpendicular to Z axis)
                FacetBox xy_center = cube2.getPlaneFacets("x", "y", "0");
                cout << "  getPlaneFacets(\"x\", \"y\", \"0\") -> " << xy_center.size() << " triangles (center XY plane)\n";

                // getPlaneFacets("x", "0", "z") - center XZ plane (perpendicular to Y axis)
                FacetBox xz_center = cube2.getPlaneFacets("x", "0", "z");
                cout << "  getPlaneFacets(\"x\", \"0\", \"z\") -> " << xz_center.size() << " triangles (center XZ plane)\n";

                // getPlaneFacets("0", "y", "z") - center YZ plane (perpendicular to X axis)
                FacetBox yz_center = cube2.getPlaneFacets("0", "y", "z");
                cout << "  getPlaneFacets(\"0\", \"y\", \"z\") -> " << yz_center.size() << " triangles (center YZ plane)\n";

                // Test with different layers
                FacetBox xy_top = cube2.getPlaneFacets("x", "y", "2");
                cout << "  getPlaneFacets(\"x\", \"y\", \"2\") -> " << xy_top.size() << " triangles (top XY plane)\n";

                FacetBox xy_bottom = cube2.getPlaneFacets("x", "y", "-2");
                cout << "  getPlaneFacets(\"x\", \"y\", \"-2\") -> " << xy_bottom.size() << " triangles (bottom XY plane)\n";

                // Store one for rendering (replace the old plane_subcells assignment)
                plane_subcells = xy_center;
                cout << "  String-based method working correctly!\n";

            } catch (const std::exception& e) {
                cout << "  Error testing string-based getPlaneFacets: " << e.what() << "\n";
            }
        }

    }

}

///////////////////     DRAW       ///////////////////////
void Draw() {
    
    extern void Setup(); // assume you define this elsewhere
    static int ciclo = 1;  // or however you manage visibility
	if (ciclo > 0) {
        /*Draw here with OpenGL*/	
        // ==================== CUBE RENDERING ====================

        //// Draw first cube (basic cube, positioned right side) - Main.c color scheme
        //// This demonstrates the basic cube triangulation with 12 faces
        //size_t cube1_total = cube_facets_1.size();
        //for(size_t i = 0; i < cube1_total; ++i) {
        //    drawFacetMainCStyle(cube_facets_1[i], i);
        //}

        //// Draw second cube (subdivided cube with vertex manipulation, positioned left side)
        //// This demonstrates cube subdivision with real-time vertex modification - Main.c colors
        //size_t cube2_total = cube_facets_2.size();
        //for(size_t i = 0; i < cube2_total; ++i) {
        //    drawFacetMainCStyle(cube_facets_2[i], i + 100);  // Offset index for color variation
        //}

        //// Draw third cube (hollow frame effect, positioned above) - Main.c colors
        //// This demonstrates selective subcell rendering
        //size_t cube3_total = cube_facets_3.size();
        //for(size_t i = 0; i < cube3_total; ++i) {
        //    drawFacetMainCStyle(cube_facets_3[i], i + 200);  // Different offset for variety
        //}

        //// ==================== ADVANCED SUBDIVISION DEMOS RENDERING ====================

        //// Demo 1: Highlight single subcell using Main.c color scheme
        //for(size_t i = 0; i < single_subcell.size(); ++i) {
        //    drawFacetMainCStyle(single_subcell[i], i + 300);  // Different offset
        //}

        // Draw original cube2 plane subcells based on current UI selection
        //plane_subcells_2 = cube2.getPlaneFacets(g_currentOrientation, g_currentLayer);
        //for(size_t i = 0; i < plane_subcells_2.size(); ++i) {
        //    drawFacetMainCStyle(plane_subcells_2[i], i + 100);  // Offset index for color variation
        //}

        plane_subcells = cube3.getCheckerboardFacets();//getPlaneFacets(g_currentOrientation, g_currentLayer);
        // Draw reflected cube using checkboard pattern
        for (size_t i = 0; i < plane_subcells.size(); ++i) {
            drawFacetMainCStyle(plane_subcells[i], (int)(i + 500));
        }

        plane_subcells = cube4.getCheckerboardFacets();//getPlaneFacets(g_currentOrientation, g_currentLayer);
        // Draw reflected cube using checkboard pattern
        for (size_t i = 0; i < plane_subcells.size(); ++i) {
            drawFacetMainCStyle(plane_subcells[i], (int)(i + 500));
        }

        plane_subcells = cube5.getCheckerboardFacets();//getPlaneFacets(g_currentOrientation, g_currentLayer);
        // Draw reflected cube using checkboard pattern
        for (size_t i = 0; i < plane_subcells.size(); ++i) {
            drawFacetMainCStyle(plane_subcells[i], (int)(i + 500));
        }

        plane_subcells = cube6.getCheckerboardFacets();//getPlaneFacets(g_currentOrientation, g_currentLayer);
        // Draw reflected cube using checkboard pattern
        for (size_t i = 0; i < plane_subcells.size(); ++i) {
            drawFacetMainCStyle(plane_subcells[i], (int)(i + 500));
        }

        plane_subcells = cube7.getCheckerboardFacets();//getPlaneFacets(g_currentOrientation, g_currentLayer);
        // Draw reflected cube using checkboard pattern
        for (size_t i = 0; i < plane_subcells.size(); ++i) {
            drawFacetMainCStyle(plane_subcells[i], (int)(i + 500));
        }

        // ==================== CHECKERBOARD SUBCELLS DEMO ====================
        // Example: Draw only subcells where (i+j+k)%2 == 0 creating a 3D checkerboard pattern
        // This demonstrates the new getCheckerboardFacets() method in Cube.hpp
        /*
        if (cube2.hasSubdivision()) {
            // Get all checkerboard pattern subcells (where (i+j+k)%2 == 0)
            FacetBox checkerboard_facets = cube2.getCheckerboardFacets();

            // Draw checkerboard pattern with distinctive color
            for (size_t i = 0; i < checkerboard_facets.size(); ++i) {
                // Use red color to highlight the checkerboard pattern
                glColor3ub(255, 100, 100);  // Light red
                drawFacetMainCStyle(checkerboard_facets[i], (int)(i + 800));
            }

            // Alternative: Get subcells directly for more control
            auto checkerboard_subcells = cube2.getCheckerboardSubcells();
            cout << "Checkerboard pattern contains " << checkerboard_subcells.size()
                 << " subcells out of " << (cube2.getSubdivisionLevels() * cube2.getSubdivisionLevels() * cube2.getSubdivisionLevels()) << " total\n";

            // You can also process individual checkerboard subcells:
            for (const auto& subcell_ref : checkerboard_subcells) {
                const Cube::SubCell& cell = subcell_ref.get();
                Vector3D center = cell.center;
                // Do something with each checkerboard subcell center, e.g., draw a sphere
                // drawSphere(center, 0.05f);
            }
        }
        */

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
        "[G]: Toggle Glass Mode",
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

    // Display reflection status
    y -= 0.05f;
    char reflectionStatus[150];
    sprintf(reflectionStatus, "Reflection: %s (Checkboard)", g_show_reflection ? "ON" : "OFF");
    glColor3f(0.2f, 0.8f, 0.2f);  // Green color for reflections
    glRasterPos2f(0.02f, y);
    for(const char* c=reflectionStatus; *c; ++c)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);

    // Show reflection parameters
    y -= 0.03f;
    char refParams[100];
    sprintf(refParams, "  Center: (%.2f, %.2f, %.2f)", g_reflection_center.x(), g_reflection_center.y(), g_reflection_center.z());
    glColor3f(0.0f, 0.6f, 0.0f);
    glRasterPos2f(0.02f, y);
    for(const char* c=refParams; *c; ++c)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, *c);

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

    // Draw a test sphere to demonstrate glass effect
    drawSphere(ce, ra, 16, 16);

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
// Keyboard handler for glass mode toggle
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

        // Plane orientation controls
        case 'o':
        case 'O':
            g_currentOrientation = (g_currentOrientation + 1) % 3;
            printf("Orientation: %s\n", g_currentOrientation == 0 ? "XY" :
                                       g_currentOrientation == 1 ? "YZ" : "XZ");
            glutPostRedisplay();
            break;

        // Layer controls
        case '+':
        case '=':
            g_currentLayer = (g_currentLayer + 1) % g_maxLayers;
            printf("Layer: %d/%d\n", g_currentLayer, g_maxLayers-1);
            glutPostRedisplay();
            break;
        case '-':
        case '_':
            g_currentLayer = (g_currentLayer - 1 + g_maxLayers) % g_maxLayers;
            printf("Layer: %d/%d\n", g_currentLayer, g_maxLayers-1);
            glutPostRedisplay();
            break;

        // Reflection controls
        case 'r':
        case 'R':
            g_show_reflection = !g_show_reflection;
            printf("Reflection: %s\n", g_show_reflection ? "ON" : "OFF");
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

