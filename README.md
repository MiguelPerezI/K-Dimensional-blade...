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

---

### Floating-Point Comparison: `operator==`

The equality comparison operator `==` for `Vector3D` objects uses a tolerance-based approach inspired by a result from **real analysis**:

> **Theorem (Bartle & Sherbert, _Introduction to Real Analysis_):**  
> Let a in R (real number set). If x in V(a) (neighborhood of a) for every epsilon > 0, then x = a.

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

