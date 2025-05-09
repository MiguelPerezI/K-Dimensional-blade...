# Table of Contents

- [K-Dimensional Blade](#k-dimensional-blade)
- [Vector3D Usage Guide](#vector3d-usage-guide)
  - [What each line demonstrates](#what-each-line-demonstrates)
  - [Floating-Point Comparison: `operator==`](#floating-point-comparison-operator)
  - [Building the Sample](#building-the-sample)
  - [Next Steps](#next-steps)
- [Quaternion Class Usage Guide](#quaternion-class-usage-guide)
  - [Overview](#overview)
  - [Constructors](#constructors)
  - [Element Access](#element-access)
  - [Inspectors](#inspectors)
  - [Conversion](#conversion)
  - [Arithmetic Operators](#arithmetic-operators)
  - [Compound Assignment](#compound-assignment)
  - [Utility Methods](#utility-methods)
  - [Free Functions](#free-functions)
   * [rotate() - explanation](#rotate-explanation)
  - [Output Streaming](#output-streaming)
  - [Notes](#notes)
  - [Header Dependencies](#header-dependencies)
- [Facet Class Usage Guide](#facet-class-usage-guide)
  - [Overview](#overview-1)
  - [Constructors](#constructors-1)
  - [Element Access](#element-access-1)
  - [Method: `updateFacet()`](#method-updatefacet)
  - [Method: `translate()`](#method-translate)
  - [Method: `crunch()`](#method-crunch)
  - [Stream Operators](#stream-operators)
  - [Geometry Helpers](#geometry-helpers)
- [FacetBox Usage Guide](#facetbox-class-usage-guide)
  - [Including](#including)
  - [Methods](#methods)
   * [push](#void-pushconst-vector3d-a-const-vector3d-b-const-vector3d-c)
   * [replace](#void-replace-size_t-idx-const-vector3d-a-const-vector3d-b-const-vector3d-c)
   * [clear](#void-clear-noexcept)
   * [center](#vector3d-center-const)
   * [translate](#void-translateconst-vector3d-offset)
   * [crunch](#void-crunchdouble-t-const-vector3d-pivot)
  - [Full Example](#full-example)
  - [Building](#building)

---

```cpp

/*                          |———————————————————————————————————|
                            |          _                        | 
                            |          \`*-.                    | 
                            |           )  _`-.                 | 
                            |          .  : `. .                | 
                            |          : _   '  \               | 
                            |          ; *` _.   `*-._          | 
                            |          `-.-'          `-.       | 
                            |            ;       `       `.     | 
                            |            :.       .        \    | 
                            |            . \  .   :   .-'   .   | 
                            |            '  `+.;  ;  '      :   | 
                            |            :  '  |    ;       ;-. | 
                            |            ; '   : :`-:     _.`* ;|
                            |[K-Blade] .*' /  .*' ; .*`- +'  `*'|
                            |         `*-*   `*-*  `*-*'        |
                            |———————————————————————————————————|
*/
```

---

## Introduction

Yes, I know you can download Blender and quickly design 3D objects. And yes, I’m aware of powerful libraries in Python or C++ like NumPy, SciPy, Eigen, CGAL, and OpenMesh.

But this isn’t about that.

This is my personal geometric and topological sandbox — a "coding diary" where I experiment with abstract math ideas for fun and exploration. It’s a space for recreational mathematics: What happens when a dodecahedron moves through 3D or even 4D hyperbolic space? How can we build intuitive tools for understanding the geometry of complex surfaces?

Like a chemist tinkering in their lab, I use this toolkit to model, test, and visualize higher-dimensional constructs. The goal isn’t production code — it’s clarity, curiosity, and a deepening of geometric intuition.


---


# K-Dimensional Blade

A C++ geometry toolkit for exploring **3‑manifolds**, surfaces, and vectors, with ongoing work to generalize the core to **quaternions**, 3‑D simplices (Facets), and higher‑order polytopes such as dodecahedra. Quaternions serve as the common backbone, enabling easy manipulation, rotation, and composition of every geometric object.

---

# Vector3D Usage Guide

Below is a quick reference demonstrating common ways to construct, manipulate, and inspect a `Vector3D` object. Feel free to copy‑paste the snippet into `main.cpp` and explore the output.

```cpp
// --- Construction ---------------------------------------------------
Vector3D p;                       // (0,0,0)
Vector3D q(1.0, 2.0, 3.0);        // explicit components
Vector3D r = q;                   // copy constructor
Vector3D w{0,0,0};                // brace initialization

// --- Access & modification -----------------------------------------
q[0] = 4.0;                       // write via operator[]
double zy = q[2];                 // read via operator[]
r += Vector3D(1,1,1);             // compound addition
r -= Vector3D(0.5,0.5,0.5);       // compound subtraction

// --- Arithmetic -----------------------------------------------------
Vector3D u = q + r;               // vector addition
Vector3D v = r - q;               // subtraction
Vector3D s = 2.5 * u;             // scalar multiplication
Vector3D n = unit(s);             // normalization

double d  = u * r;                // dot product
Vector3D c = u % r;               // cross product

// --- Geometry helpers ----------------------------------------------
Vector3D mid = line(0.5, q, r);   // midpoint of segment

// --- Output ---------------------------------------------------------
std::cout << "u = " << u << '\n'
          << "n (unit s) = " << n << '\n'
          << "u·r = " << d << '\n'
          << "u×r = " << c << '\n'
          << "midpoint = " << mid << '\n';
```

## What each line demonstrates

| Section                   | Purpose                                                | Key Functionality                                      |
| ------------------------- | ------------------------------------------------------ | ------------------------------------------------------ |
| **Construction**          | Default, parameterized, copy, and brace initialization | Constructors                                           |
| **Access & Modification** | Indexed assignment and compound operations             | `operator[]`, `operator+=`, `operator-=`               |
| **Arithmetic**            | Vector math (add, subtract, scale, normalize)          | `operator+`, `operator-`, scalar `operator*`, `unit()` |
| **Dot & Cross**           | Inner and outer products                               | `operator*` (dot), `operator%` (cross)                 |
| **Geometry Helpers**      | Linear interpolation                                   | `line()`                                               |
| **Output**                | Streaming to `std::ostream`                            | `operator<<`                                           |
| **Vector Length**         | Euclidean norm in $\mathbb{R}^3$                       | `abs()`                                                |
| **Inifity Norm**          | Euclidean inifinty norm in $\mathbb{R}^3$              | `infty()`                                              |
| **Normal**                | Mapping to $B^3 = \\{\,x \in \mathbb{R}^3 : \|x\| < 1\\}$| `unit()`                                             |


---

### Floating-Point Comparison: `operator==`

The equality comparison operator `==` for `Vector3D` objects uses a tolerance-based approach inspired by a result from **real analysis**:

> **Theorem (Bartle & Sherbert, _Introduction to Real Analysis_):**  
> $Let a\in\mathbb{R}\text{ and for each }\epsilon>0\text{ define }V_\epsilon(a)=\\{\,x\in\mathbb{R}:\;|x-a|<\epsilon\\}.\text{ If }x\in V_\epsilon(a)\quad\text{for every }\epsilon>0,\text{ then }x=a.$

This guides the implementation of `operator==` using a neighborhood epsilon approach. In `Vector3D`, two vectors `a` and `b` are considered equal if their Euclidean distance is less than a small epsilon:

```cpp
constexpr double epsilon = 1e-12;
```

This accounts for small numerical inaccuracies common in floating-point arithmetic.

---
## Notes

* Vector3D equality uses tolerance `1e-12`

---

### Building the Sample

Assuming all sources are in the current directory:

```bash
g++ -std=c++11 -I. main.cpp Vector3D.cpp -o vector_demo
./vector_demo
```

You should see something like:

```
u = (5.0, 5.0, 6.0)
n (unit s) = (0.57, 0.57, 0.68)
u·r = 43
u×r = (-3, 6, -5)
midpoint = (2.5, 3.5, 3.5)
```

*(Exact numbers may vary depending on your edits.)*

---

### Next Steps

- Explore additional helper functions like `centerM4`, `cPenta`, `centerM8`, and `cDodeca`.
- Integrate `Vector3D` into your geometry pipeline or simulation projects.
- Consider migrating storage to `std::array<double,3>` for built-in bounds checking or extended operator overloads.


# Quaternion Class Usage Guide

This guide documents the `Quaternion` class provided in the C++ geometry module. It demonstrates the class's core functionality and showcases usage with clear examples.

---

## Overview

A quaternion is represented as:

```
q = u + v
```

Where:

* `u` is the **scalar** part (a `double`)
* `v` is the **vector** part (a `Vector3D`)

---

## Constructors

```cpp
Quaternion q0;                      // Default constructor → (0, (0, 0, 0))
// Define a 3D vector v1
Vector3D v1(1, 2, 3);
Quaternion q1(5, v1);               // Parameterized constructor
Quaternion q2 = q1;                 // Copy constructor
Vector3D v2(0.1, 0.2, 0.3);
// Define a 3D vector v2
Quaternion q3(v2);                  // Pure-vector constructor (scalar=0)
Quaternion q4(3.1415);              // Pure-scalar constructor (vector=(0,0,0))
```

---

## Element Access

```cpp
q1[0];

// Example usage
Vector3D vp = Vector3D(q1[0]);       // Returns vector part (1, 2, 3)
```

---

## Inspectors

```cpp
double r = q1.r();      // Scalar part → 5
Vector3D w = q1.V();      // Vector part → (1, 2, 3)
double i = q1.i();      // x-component → 1
double j = q1.j();      // y-component → 2
double k = q1.k();      // z-component → 3
```

---

## Conversion

```cpp
q1.v4();     // Converts to Vector4D → (5, 1, 2, 3)

// Example usage:
Vector4D vQ = Vector4D(q1.v4());
```

---

## Arithmetic Operators

```cpp
// Suppose q1, q2, and q3 are Quaternions
Quaternion qAdd = q1 + q3;     // → (5, (1.1, 2.2, 3.3))
Quaternion qSub = q1 - q3;     // → (5, (0.9, 1.8, 2.7))
Quaternion qScaled = 2.0 * q3;    // → (0, (0.2, 0.4, 0.6))
Quaternion qMul = q1 * q3;     // Quaternion multiplication → (-1.4, (0.5, 1, 1.5))
q1 == q2;    // Quaternion comparison → true if close

if (q1 == q2) {
    cout << "Quaternion comparison:             Q1==Q2                    \t→  Are Equal\n";
} else {
    cout << "Quaternion comparison:             Q1==Q2                    \t→  Are Distinct\n";
}
```

---

## Compound Assignment

```cpp
q1 += q3;    // q1 modified in-place
q1 -= q3;
q1 /= 2.0;
```

---

## Utility Methods

```cpp
Quaternion qc = q1.conjugate();  // Returns conjugate → (scalar, -vector)
```

---

## Free Functions

```cpp
Quaternion qRot = Qan(M_PI/2, Vector3D(0, 0, 1));  // Axis-angle quaternion
Quaternion qCross = cross(q1, q3);                 // Quaternion-like cross product
Quaternion qRotated = rotate(q1, a, b, qRot);      // Rotate q1 from point a to b using qRot axis
```

---

### Qan() - explanation

**Purpose**: Construct a unit‐quaternion representing a rotation of $\theta$ radians about the given 3D `axis`.
1. Any rotation of a $\mathbf{p}$ in $\mathbb{R}^3$ by angle θ around a unit vector **n** = $(nx, ny, nz)$ is encoded as:
    
    $Q = \Bigl(\cos\tfrac\theta2,\sin\tfrac\theta2\,\mathbf n\Bigr) = \bigl(w\,\mathbf v\bigr)$

where $w = \cos\bigl(\tfrac{\theta}{2}\bigr)$ and $\mathbf v = \sin(\tfrac\theta2)\,\mathbf n$.

Acting on a pure-vector quaternion $p=\bigl(0, \mathbf{p} \bigr)$, the rotated vector is $p´= QpQ_{-1}$.

> **Small-Angle Shortcut**: If θ is extremely small (near zero), we skip the trigonometric calls and return the “pure‐vector” quaternion $(0,\mathbf n)$, which represents an infinitesimal rotation in the axis direction.
> For $|\theta| << 1 $, we get $\cos \tfrac{\theta}{2} \approx 1$ and $\sin \tfrac{\theta}{2} \approx \tfrac{\theta}{2}$, so $Q \approx \bigl( 1, \tfrac{\theta}{2}\mathbf{n}\bigr)$. Dropping the scalar part altogether and returning $\bigl(0, \mathbf{n} \bigr)$ captures the axis direction without catastrophic cancellation.

2. **Output**: A unit quaternion $Q=(w, \mathbf{v})$ such that rotating any pure-vector quaternion $\mathbf{p}$ via $p´= Q\mathbf{p}Q_{-q}$ applies exactly the desired rotation in $\mathbb{R}^3$.

---
### rotate() - explanation

Rotate a point‐quaternion `p` from position `a` to position `b` around a given axis.
We treat `p` as a “point” in 3D space encoded as a pure‐vector quaternion,
then:

1. Determine the rotation axis (τ) and angle (φ) needed to turn the ray from **b → a** 
   onto the global “up” direction $(0, 0, 1)$.

2. **Colinearity check**: does $a-b$ already align with *z-axis*?

3. Compute the rotation axis $\tau = normalize(normal × difference)$:
    - normal.V(): axis we want to rotate around (as Vector3D)
    - difference.V(): direction to align with $z$.

4. Compute the rotation angle $\phi$ between the normal and difference:
    - $\phi=arccos(normal-difference)$

5. Build the unit quaternion representing rotation of φ about τ:
    - $Q_{\mathrm{an}} = (cos(φ/2),  sin(φ/2)\tau)$
    - $\tau$ = eigenvector for space

3. Build the corresponding rotation quaternion  
   
   $Q = Q_{\mathrm{an}}(\varphi,\tau) = \left(\cos\frac{\varphi}{2},\, \sin\frac{\varphi}{2}\,\tau\right)$

4. Apply that rotation to `p`:

    $p_1 = Q\,p\,Q^{-1}$

5. Finally translate `p₁` so it’s moved back to “around point `b`”:

    $p_{\mathrm{out}} = p_1 + (0,\mathbf{b})$

---

## Hyperbolic $n-space$ - $\mathbb{H}^n$

The *Lorentzian inner product* in $\mathbb{R}^{n+1}$ is defined as
$\langle x, y \rangle = -x_{n+1}y_{n+1} + \sum_{i=1}^{n} x_{i}y_{i}$.

Denote $\mathbb{R}^n$ as $\mathbb{R}^n \bigcup \\{0\\}$. Let $x \in D^{n}=\\{x \in \mathbb{R}^{n} : \|x\| < 1\\}$
and $e_{n+1} \in \mathbb{R}^{n+1}$, thus $x+e_{n+1} \in \mathbb{R}^{n+1}$.

Since $x=(x_{1}, x_{2}, \cdots, x_{n}. 0)$, lets define $\bar{x} = (x_{1}, x_{2}, \cdots, x_{n})$. With this notation
the *Lorentzian norm* in $\mathbb{R}^{n+1}$ is defined as $| \|x\| |= \langle x, y \rangle^{\tfrac12}= \bigl(-x_{n+1} + \sum_{i=1}^{n} x_i^2 \bigr)^{\tfrac12}= \bigl(-x_{n+1} + \|\bar x\|^{2}\bigr)^{\tfrac12}$.

## Output Streaming

```cpp
std::cout << q1;  // Outputs: (scalar, (x, y, z))
```

---

## Notes

* Scalar equality uses tolerance `1e-12`
* Vector equality uses `Vector3D::operator==`
* Designed to follow C++ Rule of Zero with noexcept-safe operations

---

## Header Dependencies

* Requires: `Vector3D.hpp`, `Vector4D.hpp`
* No external dependencies

---

This documentation provides a solid foundation for anyone looking to integrate or extend quaternion support in a C++ geometry pipeline.

## Facet Class Overview

The `Facet` class models a triangle in 3D space using three `Quaternion` vertices: A, B, and C. It also maintains a normal vector `N` computed as the cross product of vectors (B - A) and (C - A). This structure is ideal for geometric applications, mesh processing, and graphics.

```
     [A] _________ [B]
        \^========/
         \^=[N]==/   [Facet]
          \^====/
           \^==/
            \^/
            [C]
```

### Constructors

```cpp
Facet f0;
```

Default constructor initializes all vertices and normal to zero.

```cpp
Vector3D a(0,0,0), b(1,0,0), c(0,1,0);
Facet f1(a, b, c);
```

Constructs the facet using 3 `Vector3D` points.

```cpp
Quaternion qa(0, a), qb(0, b), qc(0, c);
Facet f2(qa, qb, qc);
```

Constructs the facet from 3 `Quaternion` points.

```cpp
Facet f3(f2);
```

* **Copy constructor**: Creates a new `Facet` by duplicating the vertices and normal of an existing facet. This results in an exact geometric copy.

---

### Element Access

```cpp
f1[0] // returns vertex A as a Vector3D
f1[1] // returns vertex B as a Vector3D
f1[2] // returns vertex C as a Vector3D
f1[3] // returns normal N as a Vector3D
```

Access vertices and the precomputed normal using `operator[]`.

---

### Update Facet

```cpp
f1.updateFacet(d, e, f);
```

Updates the facet's vertices and recomputes the normal.

---

### Translate

```cpp
Vector3D offset(0, 0, 1);
f1.translate(offset);
```

Shifts all three vertices by a given Vector3D offset.

---

### Method: `crunch()`

```cpp
double factor = 0.25;
Vector3D origin(0, 0, 0);
f2.crunch(factor, origin);
```

**Purpose**: Scales the triangle toward or away from a pivot point.

#### Explanation:

* Each vertex is repositioned using the formula:

  ```cpp
  newP = pivot + t * (P - pivot);
  ```
* **t = 1.0** → No change
* **0 < t < 1** → Vertex moves closer to `pivot`
* **t > 1** → Vertex moves farther away from `pivot`
* This uniformly scales the triangle about the `pivot`. After scaling, the triangle’s **centroid** and **normal vector** are updated to reflect the new geometry.

This is especially useful in animations or mesh deformation workflows, where you need to dynamically "shrink" or "expand" triangles around a point.

---

### Stream Operators

```cpp
std::istringstream iss("(0,0,0) (1,1,0) (1,0,1)");
Facet f4;
iss >> f4;
```

**`operator>>`** reads three `(x,y,z)` coordinates to construct a facet.

```cpp
std::cout << f4;
```

**`operator<<`** outputs all three vertices.

---

Prints the facet in a readable format.

---

### Geometry Helpers

* `getNormal()` → Returns the current normal vector.
* `getCenter()` → Computes the centroid:

  ```cpp
  center = (A + B + C) / 3;
  ```

---


This `Facet` class provides a lightweight and expressive interface for 3D triangle modeling. Combined with the `Vector3D` and `Quaternion` classes, it forms a robust foundation for geometric computation and graphics programming.

# FacetBox Usage Guide

`FacetBox` is a lightweight container for managing a dynamic collection of `Facet` objects (triangles). It uses `std::vector` internally and provides convenience methods for common geometry operations.

---

## Including

```cpp
#include "FacetBox.hpp"
```

Make sure you also have:

```cpp
#include "Vector3D.hpp"
#include "Quaternion.hpp"
#include "Facet.hpp"
```

---

## Methods

### `void push(const Vector3D& A, const Vector3D& B, const Vector3D& C)`

Appends a new `Facet` constructed from points **A**, **B**, **C**.

```cpp
FacetBox box;
Vector3D A(0,0,0), B(1,0,0), C(0,1,0);
box.push(A, B, C);
// box.size() == 1, contains facet with vertices A, B, C
```

---

### `void replace(size_t idx, const Vector3D& A, const Vector3D& B, const Vector3D& C)`

Replaces the facet at index `idx` by a new `Facet(A,B,C)`. Throws `std::out_of_range` if `idx >= size()`.

```cpp
Vector3D D(1,1,0), E(2,1,0), F(1,2,0);
box.replace(0, D, E, F);
// box[0] now has vertices D, E, F
```

---

### `void clear() noexcept`

Removes all facets.

```cpp
box.clear();
// box.size() == 0
```

---

### `Vector3D center() const`

Computes the centroid of *all* vertices across *all* stored facets:

$$
\text{center} = \frac{1}{3N} \sum_{i=0}^{N-1} \bigl(A_i + B_i + C_i\bigr)
$$

```cpp
Vector3D ctr = box.center();
std::cout << "Center = " << ctr << "\n";
```

---

### `void translate(const Vector3D& offset)`

Translates every facet by adding `offset` to each vertex.

```cpp
Vector3D offset(0,0,1);
box.translate(offset);
// every facet’s vertices have z–coordinate increased by 1
```

---

### `void crunch(double t, const Vector3D& pivot)`

Scales (“crunches”) each facet toward (or away from) `pivot` by factor `t`:

$$
P' = \text{pivot} + t\,(P - \text{pivot})
$$

* `t == 1.0`: no change
* `0 < t < 1.0`: moves toward pivot
* `t > 1.0`: moves away (expansion)

```cpp
box.crunch(0.5, Vector3D(0,0,0));
// each vertex is now halfway between its old position and the origin
```

---

## Full Example

```cpp
#include <iostream>
#include "FacetBox.hpp"
#include "Vector3D.hpp"

int main() {
    // 1) Create an empty FacetBox
    FacetBox box;

    // 2) push(A,B,C): add a triangle (0,0,0)-(1,0,0)-(0,1,0)
    Vector3D A(0,0,0), B(1,0,0), C(0,1,0);
    box.push(A, B, C);
    std::cout << "After push: box.size() = " << box.size() << "\n";

    // 3) center(): compute centroid of all vertices
    Vector3D ctr = box.center();
    std::cout << "Center of box = " << ctr << "\n";

    // 4) translate(offset): move every facet by (0,0,1)
    Vector3D offset(0,0,1);
    box.translate(offset);
    std::cout << "After translate by " << offset
              << ": box.center() = " << box.center() << "\n";

    // 5) crunch(t, pivot): scale each facet toward the pivot
    box.crunch(0.5, Vector3D(0,0,0));
    std::cout << "After crunch(0.5, origin): box.center() = "
              << box.center() << "\n";

    // 6) replace(idx, A',B',C'): replace the first triangle
    Vector3D D(1,1,0), E(2,1,0), F(1,2,0);
    box.replace(0, D, E, F);
    std::cout << "After replace(0): new center = "
              << box.center() << "\n";

    // 7) clear(): remove all facets
    box.clear();
    std::cout << "After clear: box.size() = " << box.size() << "\n";

    return 0;
}
```

---

## Building

```bash
g++ -std=c++11 -I. main.cpp \
    Vector3D.cpp Quaternion.cpp Facet.cpp FacetBox.cpp \
    -o facetbox_demo
./facetbox_demo
```

This guide covers every public method of **`FacetBox`**, demonstrating how to add, replace, clear, query, and transform your collection of facets in a modern, safe C++ style.

