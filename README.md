# K-Dimensional Blade: A Geometric Computation Framework

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

## Table of Contents

- [Introduction](#introduction)
- [Project Overview](#project-overview)
- [Mathematical Foundations](#mathematical-foundations)
- [Core Components](#core-components)
  - [Vector3D](#vector3d-class)
  - [Vector4D](#vector4d-class)
  - [Quaternion](#quaternion-class)
  - [Facet](#facet-class)
  - [FacetBox](#facetbox-class)
  - [Dodecahedron](#dodecahedron-class)
- [Advanced Features](#advanced-features)
  - [Hyperbolic Geometry](#hyperbolic-geometry)
  - [Mesh Refinement Algorithms](#mesh-refinement-algorithms)
- [Usage Examples](#usage-examples)
- [Building & Installation](#building--installation)
- [API Reference](#api-reference)
- [Contributing](#contributing)

---

## Introduction

Yes, I know you can download Blender and quickly design 3D objects. And yes, I'm aware of powerful libraries in Python or C++ like NumPy, SciPy, Eigen, CGAL, and OpenMesh.

But this isn't about that.

This is my personal geometric and topological sandbox — a "coding diary" where I experiment with abstract math ideas for fun and exploration. It's a space for recreational mathematics: What happens when a dodecahedron moves through 3D or even 4D hyperbolic space? How can we build intuitive tools for understanding the geometry of complex surfaces?

Like a chemist tinkering in their lab, I use this toolkit to model, test, and visualize higher-dimensional constructs. The goal isn't production code — it's clarity, curiosity, and a deepening of geometric intuition.

---

## Project Overview

**K-Dimensional Blade** is a C++ geometry toolkit for exploring **3-manifolds**, surfaces, and vectors, with ongoing work to generalize the core to **quaternions**, 3-D simplices (Facets), and higher-order polytopes such as dodecahedra. Quaternions serve as the common backbone, enabling easy manipulation, rotation, and composition of every geometric object.

### Key Features

- **Multi-dimensional Vector Algebra**: Complete implementations for 3D and 4D vector spaces
- **Quaternion Mathematics**: Full quaternion algebra for rotations and transformations
- **Mesh Generation**: Triangle-based mesh system with multiple subdivision algorithms
- **Hyperbolic Geometry**: Support for hyperbolic space transformations
- **Interactive Visualization**: Real-time OpenGL rendering with intuitive controls
- **Geometric Primitives**: Regular polyhedra with configurable tessellation

---

## Mathematical Foundations

### Vector Spaces

The framework implements vector spaces $\mathbb{R}^3$ and $\mathbb{R}^4$ with complete algebraic operations.

#### 3D Vector Operations

For vectors $\vec{v}, \vec{w} \in \mathbb{R}^3$:

- **Dot Product**: $\vec{v} \cdot \vec{w} = \sum_{i=1}^{3} v_i w_i$
- **Cross Product**: $\vec{v} \times \vec{w} = (v_2w_3 - v_3w_2, v_3w_1 - v_1w_3, v_1w_2 - v_2w_1)$
- **Euclidean Norm**: $\|\vec{v}\|_2 = \sqrt{\sum_{i=1}^{3} v_i^2}$
- **Infinity Norm**: $\|\vec{v}\|_\infty = \max_{i \in \{1,2,3\}} |v_i|$

#### 4D Vector Extension

The 4D cross product requires three vectors $\vec{a}, \vec{b}, \vec{c} \in \mathbb{R}^4$ to produce a vector orthogonal to all three:

$$\text{Cross}(\vec{a}, \vec{b}, \vec{c}) = \vec{d} \text{ where } \vec{d} \perp \vec{a}, \vec{b}, \vec{c}$$

### Quaternion Algebra

Quaternions extend complex numbers to four dimensions, represented as:

$$q = s + x\mathbf{i} + y\mathbf{j} + z\mathbf{k}$$

where $\mathbf{i}^2 = \mathbf{j}^2 = \mathbf{k}^2 = \mathbf{ijk} = -1$.

#### Quaternion Multiplication

For $q_1 = s_1 + \vec{v}_1$ and $q_2 = s_2 + \vec{v}_2$:

$$q_1 q_2 = (s_1 s_2 - \vec{v}_1 \cdot \vec{v}_2) + (s_1 \vec{v}_2 + s_2 \vec{v}_1 + \vec{v}_1 \times \vec{v}_2)$$

#### Rotation Representation

A rotation by angle $\theta$ around unit axis $\hat{n}$ is represented as:

$$q = \cos\left(\frac{\theta}{2}\right) + \sin\left(\frac{\theta}{2}\right)\hat{n}$$

To rotate a vector $\vec{v}$, we compute: $\vec{v}' = q\vec{v}q^{-1}$

---

## Core Components

### Vector3D Class

The `Vector3D` class provides a complete implementation of 3-dimensional vector algebra.

#### Key Features

- **Constructors**: Default, parameterized, copy, and brace initialization
- **Operators**: Element access, arithmetic operations, compound assignments
- **Geometric Operations**: Dot product, cross product, normalization
- **Utilities**: Linear interpolation, collinearity testing, centroid calculations

#### Mathematical Operations

```cpp
Vector3D a(1, 2, 3), b(4, 5, 6);
double dot = a * b;              // Dot product: a·b
Vector3D cross = a % b;          // Cross product: a×b
Vector3D normalized = unit(a);   // Unit vector: a/|a|
double length = abs(a);          // Euclidean norm: |a|
```

#### Floating-Point Comparison

The equality operator uses tolerance-based comparison inspired by real analysis:

> **Theorem**: Let $a \in \mathbb{R}$ and for each $\epsilon > 0$ define $V_\epsilon(a) = \{x \in \mathbb{R} : |x-a| < \epsilon\}$. If $x \in V_\epsilon(a)$ for every $\epsilon > 0$, then $x = a$.

Implementation uses $\epsilon = 10^{-12}$ for numerical stability.

### Vector4D Class

Extends 3D operations to 4-dimensional space with additional functionality for 4D cross products and transformations.

```cpp
Vector4D v(1, 2, 3, 4);
Vector4D a(1,0,0,0), b(0,1,0,0), c(0,0,1,0);
Vector4D orthogonal = Cross(a, b, c);  // 4D cross product
```

### Quaternion Class

Implements Hamilton's quaternions for 3D rotations and hyperbolic transformations.

#### Representation

A quaternion $q = u + \vec{v}$ where:
- $u$ is the scalar part (real component)
- $\vec{v}$ is the vector part (imaginary components)

#### Key Operations

```cpp
// Rotation quaternion: 90° around z-axis
Quaternion q = Qan(M_PI/2, Vector3D(0, 0, 1));

// Quaternion multiplication
Quaternion q1(1, Vector3D(1,2,3));
Quaternion q2(2, Vector3D(4,5,6));
Quaternion product = q1 * q2;

// Conjugate: q* = u - v
Quaternion conj = q1.conjugate();
```

#### Rotation Function

The `rotate()` function implements rotation of a point around an arbitrary axis:

1. Translate point to origin relative to axis
2. Apply quaternion rotation: $p' = qpq^{-1}$
3. Translate back to original position

### Facet Class

Represents a triangular face in 3D space with vertices and normal vector.

#### Structure

```
     [A] _________ [B]
        \^========/
         \^=[N]==/   
          \^====/
           \^==/
            \^/
            [C]
```

#### Mathematical Definition

For vertices $A, B, C \in \mathbb{R}^3$, the normal vector is:

$$\vec{N} = \frac{(B - A) \times (C - A)}{|(B - A) \times (C - A)|}$$

#### Transformations

- **Translation**: $T(\vec{v}): P \mapsto P + \vec{v}$
- **Scaling**: $S(t, \vec{p}): P \mapsto \vec{p} + t(P - \vec{p})$
- **Hyperbolic Projection**: Maps to hyperboloid model

```cpp
Facet triangle(Vector3D(0,0,0), Vector3D(1,0,0), Vector3D(0,1,0));
triangle.translate(Vector3D(0, 0, 1));     // Move up
triangle.crunch(0.5, Vector3D(0,0,0));     // Scale by 50%
triangle.applyHyperboloid();               // Project to hyperbolic space
```

### FacetBox Class

Container for managing collections of triangular facets with advanced subdivision algorithms.

#### Subdivision Modes

1. **Centroid3**: Subdivides using centroid
   $$C = \frac{A + B + C}{3}$$

2. **Midpoint4**: 4-way subdivision using edge midpoints
   $$M_{AB} = \frac{A + B}{2}, \quad M_{BC} = \frac{B + C}{2}, \quad M_{CA} = \frac{C + A}{2}$$

3. **Midpoint6**: 6-way subdivision combining centroid and midpoints

4. **Sierpinski**: Fractal subdivision pattern

#### Operations

```cpp
FacetBox mesh;
mesh.push(v1, v2, v3);                    // Add triangle
mesh += otherMesh;                        // Merge meshes
FacetBox refined = mesh.refine(3, FacetBox::SubdivisionMode::Midpoint4);
refined.applyHyperboloid();               // Transform to hyperbolic space
```

### Dodecahedron Class

Generates regular dodecahedra with configurable triangulation.

#### Mathematical Properties

- **Vertices**: 20 vertices positioned using golden ratio $\phi = \frac{1 + \sqrt{5}}{2}$
- **Faces**: 12 pentagonal faces (triangulated into 36 triangles)
- **Edges**: 30 edges

#### Vertex Coordinates

Based on permutations of:
- $(0, \pm\frac{1}{\phi}, \pm\phi)$
- $(\pm 1, \pm 1, \pm 1)$
- $(\pm\frac{1}{\phi}, \pm\phi, 0)$

```cpp
Dodecahedron dodec(2.0, Vector3D(0,0,0), FaceMode::Pentagons);
FacetBox faces = dodec.getFacets();
```

---

## Advanced Features

### Hyperbolic Geometry

The framework implements the hyperboloid model of hyperbolic $n$-space $\mathbb{H}^n$.

#### Lorentzian Geometry

The Lorentzian inner product in $\mathbb{R}^{n+1}$ is defined as:

$$\langle x, y \rangle = -x_{n+1}y_{n+1} + \sum_{i=1}^{n} x_i y_i$$

#### Projections

1. **Gnomonic Projection** $\mu: D^n \to \mathbb{H}^n$:
   $$\mu(x) = \frac{x + e_{n+1}}{\sqrt{1 - \|x\|^2}}$$

2. **Stereographic Projection** $\lambda: \mathbb{H}^n \to B^n$:
   $$\lambda(x) = \left(\frac{2x_1}{1-|x|^2}, \ldots, \frac{2x_n}{1-|x|^2}, \frac{1+|x|^2}{1-|x|^2}\right)$$

#### Implementation

The `toHyperboloid()` method implements the composition $\lambda \circ \mu$:

```cpp
Quaternion Quaternion::toHyperboloid() const {
    double sq = v * v;
    if (sq >= 1.0)
        throw std::domain_error("toHyperboloid: vector norm must be < 1");
    
    double s = 1.0 / std::sqrt(1.0 - sq);
    Quaternion ret{s, s * v};
    return Quaternion{0.0, (1.0 / (1.0 + ret.r())) * ret.V()};
}
```

### Mesh Refinement Algorithms

The framework provides multiple subdivision schemes that preserve geometric properties while increasing resolution:

- **Centroid subdivision**: Preserves area ratios
- **Midpoint subdivision**: Maintains edge length ratios
- **Sierpinski subdivision**: Creates self-similar fractal patterns

Each algorithm can be applied iteratively:

```cpp
FacetBox mesh = dodecahedron.getFacets();
for (int i = 0; i < levels; ++i) {
    mesh = mesh.refine(1, FacetBox::SubdivisionMode::Sierpinski);
}
```

---

## Usage Examples

### Example 1: Creating and Manipulating Vectors

```cpp
// Basic vector operations
Vector3D v1(1, 2, 3), v2(4, 5, 6);
Vector3D sum = v1 + v2;              // (5, 7, 9)
double dot = v1 * v2;                // 32
Vector3D cross = v1 % v2;            // (-3, 6, -3)
Vector3D unit_v1 = unit(v1);        // Normalized vector

// Linear interpolation
Vector3D midpoint = line(0.5, v1, v2);
```

### Example 2: Quaternion Rotations

```cpp
// Create rotation quaternion: 45° around y-axis
Quaternion rot = Qan(M_PI/4, Vector3D(0, 1, 0));

// Rotate a point
Vector3D point(1, 0, 0);
Quaternion p_quat(point);
Quaternion rotated = rotate(p_quat, Vector3D(0,0,0), Vector3D(0,1,0), rot);
Vector3D result = rotated.V();
```

### Example 3: Building and Refining Meshes

```cpp
// Create dodecahedron
Dodecahedron dodec(1.0, Vector3D(0,0,0));
FacetBox mesh = dodec.getFacets();

// Apply iterative refinement
mesh = mesh.refine(4, FacetBox::SubdivisionMode::Midpoint4);

// Transform to hyperbolic space
mesh.applyHyperboloid();
```

### Example 4: Facet Manipulation

```cpp
// Create and transform a triangle
Facet triangle(
    Vector3D(0, 0, 0),
    Vector3D(1, 0, 0),
    Vector3D(0, 1, 0)
);

// Apply transformations
triangle.translate(Vector3D(0, 0, 1));  // Move up
triangle.crunch(0.5, Vector3D(0,0,0));  // Scale by 50%

// Get properties
Vector3D normal = triangle.getNormal();
Vector3D center = triangle.getCenter();
```

---

## Building & Installation

### Requirements

- C++ compiler with C++11 support
- OpenGL
- GLUT (FreeGLUT recommended)
- GLU

### Basic Build

```bash
# Compile a basic example
g++ -std=c++11 main.cpp Vector3D.cpp Vector4D.cpp Quaternion.cpp \
    Facet.cpp FacetBox.cpp Dodecahedron.cpp \
    -o k_blade -lGL -lGLU -lglut

# Run
./k_blade
```

### Building Specific Examples

The project includes multiple example programs (`main.cpp` through `main9.cpp`):

```bash
# Build the comprehensive class demonstration (main7.cpp)
g++ -std=c++11 main7.cpp Vector3D.cpp Vector4D.cpp Quaternion.cpp \
    Facet.cpp FacetBox.cpp Dodecahedron.cpp \
    -o demo -lGL -lGLU -lglut

# Build the hyperbolic dodecahedron demo (main8.cpp)
g++ -std=c++11 main8.cpp Vector3D.cpp Vector4D.cpp Quaternion.cpp \
    Facet.cpp FacetBox.cpp Dodecahedron.cpp \
    -o hyperbolic_demo -lGL -lGLU -lglut
```

### Interactive Controls

- **Left Mouse Drag**: Rotate view
- **Middle Mouse Drag**: Pan camera
- **Right Mouse Drag / Scroll**: Zoom
- **H**: Toggle help display
- **Right Click**: UI menu

---

## API Reference

### Vector3D

#### Constructors & Access
| Method | Description | Example |
|--------|-------------|---------|
| `Vector3D()` | Default constructor | `Vector3D v;` → (0,0,0) |
| `Vector3D(x,y,z)` | Component constructor | `Vector3D v(1,2,3);` |
| `Vector3D{x,y,z}` | Brace initialization | `Vector3D v{1,2,3};` |
| `operator[]` | Element access | `v[0]` returns x |
| `x()`, `y()`, `z()` | Component accessors | `v.x()` returns x component |

#### Operators
| Method | Description | Mathematical Operation |
|--------|-------------|------------------------|
| `operator+` | Vector addition | $\vec{a} + \vec{b}$ |
| `operator-` | Vector subtraction/negation | $\vec{a} - \vec{b}$, $-\vec{a}$ |
| `operator*` | Dot product / Scalar multiplication | $\vec{a} \cdot \vec{b}$, $c\vec{a}$ |
| `operator/` | Scalar division | $\vec{a}/c$ |
| `operator%` | Cross product | $\vec{a} \times \vec{b}$ |
| `operator+=` | Compound addition | $\vec{a} = \vec{a} + \vec{b}$ |
| `operator-=` | Compound subtraction | $\vec{a} = \vec{a} - \vec{b}$ |
| `operator/=` | Compound division | $\vec{a} = \vec{a}/c$ |
| `operator==` | Equality test (with tolerance) | $\|\vec{a} - \vec{b}\| < \epsilon$ |
| `operator<<`, `operator>>` | Stream I/O | Input/output formatting |

#### Functions
| Method | Description | Mathematical Operation |
|--------|-------------|------------------------|
| `abs()` | Euclidean norm | $\|\vec{v}\|_2 = \sqrt{x^2 + y^2 + z^2}$ |
| `infty()` | Infinity norm | $\|\vec{v}\|_\infty = \max(\|x\|,\|y\|,\|z\|)$ |
| `unit()` | Normalization | $\vec{v}/\|\vec{v}\|$ |
| `line()` | Linear interpolation | $(1-t)\vec{a} + t\vec{b}$ |
| `cruz()` | Cross product (function form) | $\vec{a} \times \vec{b}$ |
| `areColinear()` | Test collinearity of 3 points | Returns true if collinear |
| `centerM4()` | Center of 4 points | $(\vec{a}+\vec{b}+\vec{c}+\vec{d})/4$ |
| `cPenta()` | Center of 5 points | Average of 5 vectors |
| `cDodeca()` | Center of 20 points | Average of dodecahedron vertices |
| `centerM8()` | Center of 8 points | Average of 8 vectors |

### Vector4D

#### Constructors & Access
| Method | Description | Example |
|--------|-------------|---------|
| `Vector4D()` | Default constructor | `Vector4D v;` → (0,0,0,0) |
| `Vector4D(x,y,z,t)` | Component constructor | `Vector4D v(1,2,3,4);` |
| `operator[]` | Element access | `v[3]` returns t |
| `x()`, `y()`, `z()`, `t()` | Component accessors | `v.t()` returns t component |

#### Operators & Functions
| Method | Description | Mathematical Operation |
|--------|-------------|------------------------|
| `operator+`, `-`, `*`, `/` | Arithmetic operations | Standard vector operations |
| `operator+=`, `/=` | Compound operations | In-place modifications |
| `operator==` | Equality test | Component-wise comparison |
| `operator<` | Lexicographic ordering | For use in sorted containers |
| `Cross()` | 4D cross product (3 vectors) | Produces orthogonal vector |
| `abs()` | Euclidean norm | $\|\vec{v}\|_2$ in $\mathbb{R}^4$ |
| `infty()` | Infinity norm | $\max(\|x\|,\|y\|,\|z\|,\|t\|)$ |
| `unit()` | Normalization | $\vec{v}/\|\vec{v}\|$ |
| `zero()` | Zero vector | Returns (0,0,0,0) |
| `rect()` | Linear interpolation | $(1-t)\vec{a} + t\vec{b}$ |

### Quaternion

#### Constructors & Access
| Method | Description | Example |
|--------|-------------|---------|
| `Quaternion()` | Default constructor | `Quaternion q;` → (0, (0,0,0)) |
| `Quaternion(s, v)` | Scalar-vector constructor | `Quaternion q(1, Vector3D(0,0,1));` |
| `Quaternion(v)` | Pure vector constructor | `Quaternion q(Vector3D(1,2,3));` |
| `Quaternion(s)` | Pure scalar constructor | `Quaternion q(3.14);` |
| `operator[]` | Access vector part | `q[0]` returns vector part |
| `r()` | Get scalar part | Returns real component |
| `V()` | Get vector part | Returns Vector3D |
| `i()`, `j()`, `k()` | Component accessors | Individual vector components |

#### Operations
| Method | Description | Mathematical Operation |
|--------|-------------|------------------------|
| `operator+`, `-` | Addition/subtraction | Component-wise |
| `operator*` | Quaternion/scalar multiplication | $q_1 q_2$ or $cq$ |
| `operator/` | Scalar division | $q/c$ |
| `operator+=`, `-=`, `/=` | Compound operations | In-place modifications |
| `operator==` | Equality test | With tolerance |
| `conjugate()` | Quaternion conjugate | $q^* = s - \vec{v}$ |
| `v4()` | Convert to Vector4D | $(s, x, y, z)$ |
| `toHyperboloid()` | Hyperbolic projection | $\lambda \circ \mu$ |

#### Free Functions
| Function | Description | Mathematical Operation |
|----------|-------------|------------------------|
| `Qan()` | Axis-angle to quaternion | $\cos(\theta/2) + \sin(\theta/2)\hat{n}$ |
| `cross()` | Quaternion cross product | Special quaternion product |
| `rotate()` | Apply rotation | $qpq^{-1}$ with translation |
| `lerp()` | Linear interpolation | $(1-t)q_1 + tq_2$ |

### Facet

#### Constructors & Access
| Method | Description |
|--------|-------------|
| `Facet()` | Default constructor |
| `Facet(A,B,C)` | From three Vector3D points |
| `Facet(qA,qB,qC)` | From three Quaternions |
| `operator[]` | Access vertices (0=A, 1=B, 2=C, 3=N) |

#### Operations
| Method | Description |
|--------|-------------|
| `getNormal()` | Get face normal vector |
| `getCenter()` | Compute triangle centroid |
| `updateFacet()` | Update all vertices |
| `translate()` | Apply translation |
| `crunch()` | Scale relative to pivot |
| `applyHyperboloid()` | In-place hyperbolic projection |
| `hyperboloid()` | Return hyperbolic-transformed copy |
| `operator>>`, `<<` | Stream I/O |

### FacetBox

#### Constructors & Factory Methods
| Method | Description |
|--------|-------------|
| `FacetBox()` | Default constructor (empty) |
| `FacetBox(facet)` | Single facet constructor |
| `fromFacet()` | Static: subdivide facet into 3 around centroid |
| `subdivide4()` | Static: 4-way midpoint subdivision |
| `subdivide6()` | Static: 6-way subdivision |
| `sierpinski()` | Static: Sierpinski pattern |

#### Operations
| Method | Description |
|--------|-------------|
| `operator[]` | Access facet by index |
| `push()` | Add new facet from 3 points |
| `replace()` | Replace facet at index |
| `clear()` | Remove all facets |
| `size()` | Number of facets |
| `center()` | Compute mesh centroid |
| `translate()` | Translate all facets |
| `crunch()` | Scale all facets |
| `merge()` | Append another FacetBox |
| `operator+=` | Merge operator |
| `operator+` | Union (creates new FacetBox) |
| `refine()` | Apply subdivision (with mode) |
| `applyHyperboloid()` | In-place hyperbolic transform |
| `hyperboloid()` | Return hyperbolic copy |
| `operator<<` | Stream output |

### Dodecahedron

#### Constructors
| Method | Description |
|--------|-------------|
| `Dodecahedron()` | Default constructor |
| `Dodecahedron(radius, center)` | Basic constructor |
| `Dodecahedron(radius, center, mode)` | With face mode (Triangles/Pentagons) |

#### Operations
| Method | Description |
|--------|-------------|
| `faceCount()` | Number of triangular faces (36) |
| `operator[]` | Access face by index |
| `center()` | Compute geometric center |
| `translate()` | Move by offset vector |
| `scale()` | Scale relative to pivot |
| `getFacets()` | Get underlying FacetBox |
| `operator<<` | Stream output |

---

## Contributing

This is a personal exploration project, but contributions that enhance mathematical clarity or add interesting geometric algorithms are welcome. Please ensure:

1. Mathematical operations are clearly documented with LaTeX notation
2. Code follows existing style conventions
3. New features include example usage
4. Complex algorithms include references to mathematical literature

---

## License

[Specify your license here]

---

## Acknowledgments

This project is inspired by the beauty of geometric algebra and the elegance of quaternion mathematics. Special thanks to the mathematical community for centuries of geometric insights.