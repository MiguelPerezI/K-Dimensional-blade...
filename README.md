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
Vector3D v1(1, 2, 3);
Quaternion q1(5, v1);               // Parameterized constructor
Quaternion q2 = q1;                 // Copy constructor
Vector3D v2(0.1, 0.2, 0.3);
Quaternion q3(v2);                  // Pure-vector constructor (scalar=0)
Quaternion q4(3.1415);              // Pure-scalar constructor (vector=(0,0,0))
```

---

## Element Access

```cpp
q1[0];       // Returns vector part (1, 2, 3)
q1[1];       // Returns scalar part as (5, 5, 5)
```

---

## Inspectors

```cpp
q1.r();      // Scalar part → 5
q1.V();      // Vector part → (1, 2, 3)
q1.i();      // x-component → 1
q1.j();      // y-component → 2
q1.k();      // z-component → 3
```

---

## Conversion

```cpp
q1.v4();     // Converts to Vector4D → (5, 1, 2, 3)
```

---

## Arithmetic Operators

```cpp
q1 + q3;     // → (5, (1.1, 2.2, 3.3))
q1 - q3;     // → (5, (0.9, 1.8, 2.7))
2.0 * q3;    // → (0, (0.2, 0.4, 0.6))
q1 * q3;     // Quaternion multiplication → (-1.4, (0.5, 1, 1.5))
q1 == q2;    // Quaternion comparison → true if close
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
q1.conjugate();  // Returns conjugate → (scalar, -vector)
```

---

## Free Functions

```cpp
Qan(M_PI/2, Vector3D(0, 0, 1));  // Axis-angle quaternion
cross(q1, q3);                   // Quaternion-like cross product
rotate(q1, a, b, qRot);          // Rotate q1 from point a to b using qRot axis
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


