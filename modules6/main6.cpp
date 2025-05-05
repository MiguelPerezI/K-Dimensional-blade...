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

Vector3D d_ui(1.0, 0.0, 0.0);
Vector3D e_ui(0.0, 1.0, 0.0);
Vector3D f_ui(0.0, 0.0, 1.0);
Facet f_1(d_ui, e_ui, f_ui);


void Setup() {

	if (ciclo == 0) {
        
        cout << "\n\n";
        cout << "|———————————————————————————————————|\n";
        cout << "|          _                        |\n"; 
        cout << "|          \\`*-.                    |\n"; 
        cout << "|           )  _`-.                 |\n"; 
        cout << "|          .  : `. .                |\n"; 
        cout << "|          : _   '  \\               |\n"; 
        cout << "|          ; *` _.   `*-._          |\n"; 
        cout << "|          `-.-'          `-.       |\n"; 
        cout << "|            ;       `       `.     |\n"; 
        cout << "|            :.       .        \\    |\n"; 
        cout << "|            . \\  .   :   .-'   .   |\n"; 
        cout << "|            '  `+.;  ;  '      :   |\n"; 
        cout << "|            :  '  |    ;       ;-. |\n"; 
        cout << "|            ; '   : :`-:     _.`* ;|\n";
        cout << "|[K-Blade] .*' /  .*' ; .*`- +'  `*'|\n";
        cout << "|         `*-*   `*-*  `*-*'        |\n";
        cout << "|———————————————————————————————————|\n";
        cout << "\n\n";
        cout << "———————————————————————————————————————————————————————————————————————\n";
        cout << "|- Testing Vector3D Class——————————————————————————————————————————————\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n";
        
        // --- Construction ---------------------------------------------------
        cout << "- List of constructors ————————————————————————————————————————————————\n";
        Vector3D p;                       // Default constructor
        cout << "Default constructor:       Vector3D p;              → p = " << p << "\n";
        
        Vector3D q(1.0, 2.0, 3.0);        // Constructor with components
        cout << "Component constructor:     Vector3D q(1.0, 2.0, 3.0) → q = " << q << "\n";
        
        Vector3D r = q;                   // Copy constructor
        cout << "Copy constructor:          Vector3D r = q;           → r = " << r << "\n";
        
        Vector3D r1{0.1, 0.1, 0.12};      // Brace-initialized constructor
        cout << "Brace constructor:         Vector3D r1{0.1,0.1,0.12}          → r1 = " << r1 << "\n";
       
        Vector3D sup = Vector3D(r1);  // Copy constructor (explicit expression)
        cout << "Copy constructor (explicit):   Vector3D sup = Vector3D(r1)  → sup = " << sup << "\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n"; 

        // --- Access & Modification -----------------------------------------
        cout << "- Access & Modification Examples ——————————————————————————————————————\n";
        cout << "";
        q[0] = 4.0;                       // Write to q[0]
        double zy = q[2];                 // Read from q[2]
        cout << "After q[0] = 4.0:          q = " << q << "\n";
        cout << "Reading q[2]:              q[2] = " << zy << "\n";
        
        r += Vector3D(2, 2, 2);           // Compound addition
        cout << "After r += Vector3D(2,2,2): r = " << r << "\n";
        
        r -= Vector3D(0.5, 0.5, 0.5);     // Compound subtraction
        cout << "After r -= Vector3D(0.5,...): r = " << r << "\n";
        
        r1 /= 0.5;                        // Compound scalar division
        cout << "After r1 /= 0.5: r1 = " << r1 << "\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n";

        // --- Arithmetic -----------------------------------------------------
        cout << "- Arithmetic ——————————————————————————————————————————————————————————\n";
        Vector3D u = q + r;               // Addition
        cout << "Vector addition:           u = q + r                → u = " << u << "\n";
        
        Vector3D v = r - q;               // Subtraction
        cout << "Vector subtraction:        v = r - q                → v = " << v << "\n";
        
        Vector3D w = 2.5 * u;             // Scalar multiplication
        cout << "Scalar multiplication:     w = 2.5 * u              → w = " << w << "\n";
        
        Vector3D n = unit(w);             // Normalization
        cout << "Normalized vector:         n = unit(w)              → n = " << n << "\n";
        
        double d  = u * r;                // Dot product
        cout << "Dot product:               d = u · r                → d = " << d << "\n";
        
        Vector3D c = u % r;               // Cross product
        cout << "Cross product:             c = u × r                → c = " << c << "\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n";
        
        // --- Geometry helpers ----------------------------------------------
        cout << "- Geometry helpers ————————————————————————————————————————————————————\n";
        Vector3D mid = line(0.5, q, r);   // Midpoint
        cout << "Midpoint of q and r:       mid = line(0.5, q, r)    → mid = " << mid << "\n";
        
        // --- Summary Output ------------------------------------------------
        cout << "\n--- Summary ---\n";
        cout << "q = " << q << "\n";
        cout << "r = " << r << "\n";
        cout << "u = " << u << "\n";
        cout << "v = " << v << "\n";
        cout << "w = " << w << "\n";
        cout << "n = " << n << "\n";
        cout << "u·r = " << d << "\n";
        cout << "u×r = " << c << "\n";
        cout << "midpoint = " << mid << "\n";
    
        cout << "———————————————————————————————————————————————————————————————————————\n";
        cout << "|- Testing Quaternion Class————————————————————————————————————————————\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n";

        // --- Construction ---------------------------------------------------
        cout << "- List of constructors ————————————————————————————————————————————————\n";
        Quaternion q0;  // default constructor
        cout << "Default constructor:             Quaternion q0;              \t→ q0 = " << q0 << "\n";
        
        Vector3D v1{1.0, 2.0, 3.0};
        cout << "- Let v1 = " << v1 << "\n";
        Quaternion q1(5.0, v1);  // scalar-vector constructor
        cout << "Parameterized constructor:        Quaternion(5, v1)          \t→ q1 = " << q1 << "\n";
        
        Quaternion q2 = q1;  // copy constructor
        cout << "Copy constructor:                 Quaternion q2 = q1         \t→ q2 = " << q2 << "\n";
        
        Vector3D v2{0.1, 0.2, 0.3};
        cout << "- Let v2 = " << v2 << "\n";
        Quaternion q3(Vector3D{0.1, 0.2, 0.3});  // pure-vector constructor
        cout << "Pure-vector constructor:          Quaternion(v2)             \t→ q3 = " << q3 << "\n";
        
        Quaternion q4(3.1415169265358979);  // pure-scalar constructor
        cout << "Pure-scalar constructor:        Quaternion(3.1415...)        \t→ q4 = " << q4 << "\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n";
        // --- Element Access -------------------------------------------------
        cout << "———————————————————————————————————————————————————————————————————————\n";
        cout << "Element access q1[0]:             q1[0] (vector part)        \t→ " << q1[0] << "\n";
        //cout << "Element access q1[1]:             q1[1] (scalar part)         \t→ " << q1[1] << "\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n";
        
        // --- Inspectors -----------------------------------------------------
        cout << "- Inspectors ——————————————————————————————————————————————————————————\n";
        cout << "Scalar (r) part of q1:            q1.r()                      \t→ " << q1.r() << "\n";
        cout << "Vector (v) part of q1:            q1.V()                      \t→ " << q1.V() << "\n";
        cout << "Components: i j k                 q1.i() q1.j() q1.k()        \t→ <i="
                  << q1.i() << ", j=" << q1.j() << ", k=" << q1.k() << ">\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n"; 
        
        // --- Conversion -----------------------------------------------------
        cout << "- Conversion ——————————————————————————————————————————————————————————\n";
        cout << "Convert q1 to Vector4D:           q1.v4()                     \t→ " << q1.v4() << "\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n";

        // --- Arithmetic Operators -------------------------------------------
        cout << "- Arithmetic Operators ————————————————————————————————————————————————\n";
        Quaternion qAdd = q1 + q3;
        cout << "Addition:                         q1 + q3                     \t→ " << qAdd << "\n";
        
        Quaternion qSub = q1 - q3;
        cout << "Subtraction:                      q1 - q3                     \t→ " << qSub << "\n";
        
        Quaternion qScaled = 2.0 * q3;
        cout << "Scalar multiplication:           2.0 * q3                    \t→ " << qScaled << "\n";
        
        Quaternion qMul = q1 * q3;
        cout << "Quaternion multiplication:        q1 * q3                    \t→ " << qMul << "\n";
        
        Quaternion Q1(1.11, Vector3D(0.121, 4.123, 1.1));
        Quaternion Q2(1.11, Vector3D(0.121, 4.123, 1.1));
        if (q1 == q2) {
        cout << "Quaternion comparison:             Q1==Q2                    \t→  Are Equal\n";
        } else {
        cout << "Quaternion comparison:             Q1==Q2                    \t→  Are Distinct\n";
        }

        cout << "———————————————————————————————————————————————————————————————————————\n\n";

        // --- Compound assignment --------------------------------------------
        cout << "- Compound assignment —————————————————————————————————————————————————\n";
        q1 += q3;
        cout << "Compound += :                     q1 += q3                    \t→ q1 = " << q1 << "\n";
        
        q1 -= q3;
        cout << "Compound -= :                     q1 -= q3                    \t→ q1 = " << q1 << "\n";
        
        q1 /= 2.0;
        cout << "Compound /= :                     q1 /= 2.0                   \t→ q1 = " << q1 << "\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n";
        
        // --- Utility --------------------------------------------------------
        cout << "- Utility —————————————————————————————————————————————————————————————\n";
        Quaternion qc = q1.conjugate();
        cout << "Conjugate:                        q1.conjugate()              \t→ " << qc << "\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n";

        // --- Free Functions -------------------------------------------------
        cout << "- Free Funcionts ——————————————————————————————————————————————————————\n";
        Quaternion qRot = Qan(3.1415 / 2, Vector3D(0, 0, 1));
        cout << "Axis-angle (Qan):                 Qan(pi/2, z-axis)           \t→ " << qRot << "\n";
        
        Quaternion qCross = cross(q1, q3);
        cout << "Quaternion cross product:         cross(q1, q3)               \t→ " << qCross << "\n";
        
        Quaternion qRotated = rotate(q1, Vector3D(1,0,0), Vector3D(0,1,0), qRot);
        cout << "Rotate q1 around axis:            rotate(q1, a, b, qRot)      \t→ " << qRotated << "\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n";


        cout << "\n———————————————————————————————————————————————————————————————————————\n";
        cout <<   "|- Testing Facet Class ———————————————————————————————————————————\n";
        cout <<   "———————————————————————————————————————————————————————————————————————\n\n";
        cout << "               [A] _________[B]\n";
        cout << "                  \\^========/\n";
        cout << "                   \\^=[N]==/     \n";                                                                                                                                                                            
        cout << "                    \\^====/\n";
        cout << "                     \\^==/ [Facet Class]\n";
        cout << "                      \\^/\n";
        cout << "                      [C]\n";
        cout << "                \n";
        // --- Construction ---------------------------------------------------
        cout << "- List of constructors ————————————————————————————————————————————————\n";
        
        Facet f0;  // Default constructor
        cout << "Default constructor: Facet f0; → " << f0 << "\n";
        
        // Define three points for Vector3D-based constructor
        Vector3D a(0.0, 0.0, 0.0);
        Vector3D b(1.0, 0.0, 0.0);
        Vector3D c0(0.0, 1.0, 0.0);
        cout << "- Let a = " << a << "\n      b = " << b << "\n      c = " << c0 << "\n";
        
        Facet f1(a, b, c);  // Vector3D constructor
        cout << "Vector3D constructor:\n - Facet(a, b, c) → " << f1 << "\n";

        // Quaternion-based constructor
        Quaternion q_a(0.0, a), q_b(0.0, b), q_c(0.0, c);
        Facet f2(q_a, q_b, q_c);
        cout << "Quaternion constructor:\n - Facet(q_a, q_b, q_c) → " << f2 << "\n";

        // Copy constructor
        Facet f3(f2);
        cout << "Copy constructor:\n - Facet f3(f2) → " << f3 << "\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n";

        // --- Element Access -------------------------------------------------
        cout << "- Element access (vertex positions) ————————————————————————————————————\n";
        cout << "f1[0]= " << f1[0] << "\nf1[1]= " << f1[1] << "\nf1[2]= " << f1[2] << "\n";
        cout << "(Normal) f1[3]=" << f1[3] << "\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n";

        // --- Update Facet ---------------------------------------------------
        cout << "- Update Facet ————————————————————————————————————————————————————————\n";
        Vector3D d_u(1.0, 1.0, 0.0);
        Vector3D e_u(2.0, 1.0, 0.0);
        Vector3D f_u(1.0, 2.0, 0.0);
        cout << " - Before f1 = " << f1 << "\n";
        f1.updateFacet(d_u, e_u, f_u);
        cout << " - After f1.updateFacet(d_u,e_u,f_u): f1 = " << f1 << "\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n";

        // --- Translation ----------------------------------------------------
        cout << "- Geometric Translation  ——————————————————————————————————————————————\n";
        Vector3D offset(0.0, 0.0, 1.0);
        cout << "- offset: " << offset << "\n";
        cout << " - Before translate: f1 = " << f1 << "\n";
        f1.translate(offset);
        cout << " - After f1.translate(offset): f1 = " << f1 << "\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n";

        // --- Crunch (scale) --------------------------------------------------
        cout << " - Crunch (scale) —————————————————————————————————————————————————————\n";
        // Geometrically: scale the triangle's vertices by factor 0.5 about the pivot at the origin.
        // Each vertex is pulled halfway closer to (0,0,0), shrinking the facet uniformly.
        cout << " - Before crunch: f2 = " << f2 << "\n";
        double factor = 0.25;
        Vector3D origin{0,0,0}; 
        cout << " - factor = " << factor << "        origin = " << origin << "\n";
        f2.crunch(factor, origin);
        cout << " - After f2.crunch(factor, origin): f2 = " << f2 << "\n";
        // The facet's center and normal have been recomputed to match the scaled geometry.
        cout << "———————————————————————————————————————————————————————————————————————\n\n";

        /*
         * Parsing operator>>: We have an input stream containing three
         * 3D points in the form "(x0,y0,z0) (x1,y1,z1) (x2,y2,z2)".
         * Each triplet defines one vertex of the triangular facet.
         * The Facet extraction operator reads these three Vector3D values,
         * assigns them to A, B, C respectively, and recomputes the facet's
         * normal N as cross(B-A, C-A).
         *
         * Geometric interpretation:
         *  - "(0,0,0) (1,1,0) (1,0,1)" defines a triangle whose vertices
         *    are P0=(0,0,0), P1=(1,1,0), P2=(1,0,1).
         *  - After parsing, f4 represents that triangle in space.
         *  - The facet normal N is perpendicular to the plane containing
         *    P0,P1,P2, computed via N = (P1-P0) × (P2-P0).
         */
        cout << "- Parsing operator>> ——————————————————————————————————————————————————\n";
        {
            istringstream iss("(0,0,0) (1,1,0) (1,0,1)");
            Facet f4;
            cout << " - Constructed default facet f4 = " << f4 << "\n";
            iss >> f4;
            cout << " - istringstream iss('(0,0,0) (1,1,0) (1,0,1)');\n";
            cout << "- Parsed f4 from stream: f4 = " << f4 << "\n";
        }
        cout << "———————————————————————————————————————————————————————————————————————\n\n";
    

        cout << "\n———————————————————————————————————————————————————————————————————————\n";
        cout <<   "|- Testing FacetBox Class ———————————————————————————————————————————\n";
        cout <<   "———————————————————————————————————————————————————————————————————————\n\n";
cout << "*—————————————————————————————————————————————————————|*\n";
cout << "| [A] _________[B]  [A] _________[B]   [A] _________[B]|\n";      
cout << "|    \\^========/       \\^========/        \\^========/  |\n";    
cout << "|     \\^=[N]==/         \\^=[N]==/          \\^=[N]==/   |\n";    
cout << "|      \\^====/           \\^====/            \\^====/    |\n";        
cout << "|       \\^==/             \\^==/              \\^==/     | ~ [FacetBox Class]\n";
cout << "|        \\^/               \\^/                \\^/      |\n";        
cout << "|        [C]               [C]                [C]      |\n";            
cout << "*—————————————————————————————————————————————————————|*\n\n";

        // --- Construction ---------------------------------------------------
        cout << "- List of constructors ————————————————————————————————————————————————\n";
        // ————————————————————————————————————————————————
        cout << "1) Default constructor: starts empty\n";
        // ————————————————————————————————————————————————
        FacetBox boxDefault;
        cout << "FacetBox boxDefault (default):\n" << boxDefault << "\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n";


        // Prepare a single Facet to use with the next constructors
        Vector3D A(0,0,0), B(1,0,0), C(0,1,0);
        Facet singleFacet(A, B, C);
        cout << " - Let Facet singleFacet = " << singleFacet << "\n\n";


        cout << "———————————————————————————————————————————————————————————————————————\n";
        // ————————————————————————————————————————————————
        cout << "2) Explicit single-facet constructor";
        // ————————————————————————————————————————————————
        FacetBox boxSingle(singleFacet);
        std::cout << "FacetBox boxSingle(singleFacet):\n" << boxSingle << "\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n";


        cout << "———————————————————————————————————————————————————————————————————————\n";
        // ————————————————————————————————————————————————
        cout << "3) Copy constructor\n";
        // ————————————————————————————————————————————————
        FacetBox boxCopy(boxSingle);
        std::cout << "FacetBox boxCopy(boxSingle):\n" << boxCopy << "\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n";


        cout << "———————————————————————————————————————————————————————————————————————\n";
        // ————————————————————————————————————————————————
        cout << "4) Move constructor\n";
        // ————————————————————————————————————————————————
        FacetBox boxMoved(std::move(boxSingle));
        std::cout << "boxMoved (moved from boxSingle):\n" << boxMoved << "\n";
        std::cout << "boxSingle now (should be empty):\n" << boxSingle << "\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n";


        cout << "———————————————————————————————————————————————————————————————————————\n";
        // ————————————————————————————————————————————————
        cout << "5) push() to extend a box";
        // ————————————————————————————————————————————————
        boxDefault.push(A, B, Vector3D(0,0,1));
        std::cout << "boxDefault after push(A,B,C):\n" << boxDefault << "\n";
        boxMoved.push(A, B, Vector3D(0.2,0,1)); 
        std::cout << "boxMoved after push(A,B,C):\n" << boxMoved << "\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n";


        cout << "———————————————————————————————————————————————————————————————————————\n";
        // ————————————————————————————————————————————————
        // 6) copy-assignment
        // ————————————————————————————————————————————————
        boxDefault = boxCopy;
        std::cout << "boxDefault after copy-assignment = boxCopy:\n" << boxDefault << "\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n";


        cout << "———————————————————————————————————————————————————————————————————————\n";
        // ————————————————————————————————————————————————
        // 7) move-assignment
        // ————————————————————————————————————————————————
        boxDefault = std::move(boxMoved);
        std::cout << "boxDefault after move-assignment = boxMoved:\n" << boxDefault << "\n";
        std::cout << "boxMoved now (should be empty):\n" << boxMoved << "\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n";


        // 1) Create an empty FacetBox
        FacetBox box;
       

        cout << "———————————————————————————————————————————————————————————————————————\n"; 
        cout << " - push method\n";
        // 2) push(A,B,C): add a triangle (0,0,0)-(1,0,0)-(0,1,0)
        Vector3D A1(1,1,-1), B1(1,0,-1), C1(0,1,-1), D1(0,0,-1);
        Vector3D A2(1,1,1), B2(1,0,1), C2(0,1,1), D2(0,0,1);
        box.push(A1, B1, C1);
        box.push(C1, D1, A1);
        box.push(A2, B2, C2);
        box.push(C2, D2, A2);
        std::cout << "After push: box.size() = " << box.size() << "\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n";


        cout << "———————————————————————————————————————————————————————————————————————\n";
        // 3) center(): compute centroid of all vertices
        Vector3D ctr = box.center();
        std::cout << "Center of box = " << ctr << "\n";
        cout << "———————————————————————————————————————————————————————————————————————\n\n";
    }

}

///////////////////     DRAW       ///////////////////////
void Draw() {
    
    extern void Setup(); // assume you define this elsewhere
    static int ciclo = 1;  // or however you manage visibility
	if (ciclo > 0) {
        /*Draw here with OpenGL*/	
        drawFacet(f_1, 200, 10, 40, 0.75f);
        drawSphere(f_1[0], 0.1f, 6, 6);
        drawSphere(f_1[1], 0.1f, 6, 6);
        drawSphere(f_1[2], 0.1f, 6, 6);
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
                      float alpha = 0.75f)
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
    glLineWidth(1.5f);
    glBegin(GL_LINE_LOOP);
        glVertex3f(A.x(), A.y(), A.z());
        glVertex3f(B.x(), B.y(), B.z());
        glVertex3f(C.x(), C.y(), C.z());
    glEnd();

    // 6) Restore state
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

