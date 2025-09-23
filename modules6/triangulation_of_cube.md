### **A Standard Triangulation of the Unit Cube \([0,1]^d\)**

This document formalizes the construction of a simplicial complex that triangulates the \(d\)-dimensional unit cube, generalizing the standard decompositions of the square and cube.

#### **1. Preliminaries**

Let \(d, n \in \mathbb{N}\) with \(d \geq 1\). The goal is to triangulate the \(d\)-dimensional unit cube, \([0,1]^d\).

*   **Grid Parameters:** The parameter \(n\) defines the granularity of the subdivision. The cube is divided into \(n^d\) smaller congruent cubes.
*   **Index Set:** We define the index set \(I_n = \{0, 1, \dots, n-1\}\).
*   **Vertex Coordinates:** A vertex is identified by a multi-index \(\mathbf{i} = (i_1, i_2, \dots, i_d) \in I_n^d\). Its geometric position in \([0,1]^d\) is given by the scaling map:
    \[
    \varphi(\mathbf{i}) = \left( \frac{i_1}{n}, \frac{i_2}{n}, \dots, \frac{i_d}{n} \right).
    \]

#### **2. Triangulation of the Unit Square (\(d=2\))**

The 2-dimensional case serves as the foundational example.

1.  **Subdivision:** The unit square \([0,1]^2\) is subdivided into an \(n \times n\) grid of smaller squares. Each square is indexed by \((i, j) \in I_n^2\) and is defined as:
    \[
    Q_{i,j} = \left[\frac{i}{n}, \frac{i+1}{n}\right] \times \left[\frac{j}{n}, \frac{j+1}{n}\right].
    \]

2.  **Cell Triangulation:** Each square \(Q_{i,j}\) is subdivided into two 2-simplices (triangles) by adding a diagonal. A standard choice is the diagonal from the vertex \((i, j)\) to \((i+1, j+1)\). This yields the following two triangles:
    \[
    \begin{aligned}
    \sigma^{(1)}_{i,j} &= \left( (i,j),\ (i+1,j),\ (i,j+1) \right), \\
    \sigma^{(2)}_{i,j} &= \left( (i+1,j),\ (i+1,j+1),\ (i,j+1) \right).
    \end{aligned}
    \]
    Here, each triple denotes an oriented 2-simplex.

3.  **Global Complex:** The full simplicial complex \(K_n^{(2)}\) triangulating \([0,1]^2\) is the union of all such simplices:
    \[
    K_n^{(2)} = \left\{ \sigma^{(1)}_{i,j},\ \sigma^{(2)}_{i,j} \mid (i,j) \in I_n^2 \right\}.
    \]
    The geometric realization is obtained by applying the vertex map \(\varphi\).

#### **3. A Specific Triangulation of the Unit Cube (\(d=3\))**

The construction generalizes to three dimensions. We now describe a specific triangulation of the unit cube into 12 tetrahedra.

1.  **Subdivision:** The unit cube \([0,1]^3\) is subdivided into an \(n \times n \times n\) grid of smaller cubes. Each cube is indexed by \((i, j, k) \in I_n^3\):
    \[
    C_{i,j,k} = \left[\frac{i}{n}, \frac{i+1}{n}\right] \times \left[\frac{j}{n}, \frac{j+1}{n}\right] \times \left[\frac{k}{n}, \frac{k+1}{n}\right].
    \]

2.  **Cell Triangulation:** Each cube \(C_{i,j,k}\) is subdivided into twelve 3-simplices (tetrahedra). The following is a valid triangulation, with simplices defined by their vertices (indices are relative to the cube's lower-left-front corner \((i, j, k)\)):
    \[
    \begin{aligned}
    \sigma_1 &= \left( (i,j,k),\ (i+1,j,k),\ (i,j+1,k),\ (i,j,k+1) \right) \\
    \sigma_2 &= \left( (i,j,k+1),\ (i+1,j,k+1),\ (i+1,j+1,k),\ (i,j,k) \right) \\
    \sigma_3 &= \left( (i,j+1,k),\ (i,j,k),\ (i,j+1,k+1),\ (i+1,j+1,k) \right) \\
    \sigma_4 &= \left( (i,j,k+1),\ (i,j+1,k+1),\ (i,j+1,k),\ (i,j,k) \right) \\
    \sigma_5 &= \left( (i,j+1,k),\ (i+1,j+1,k+1),\ (i+1,j,k+1),\ (i+1,j,k) \right) \\
    \sigma_6 &= \left( (i+1,j+1,k+1),\ (i+1,j+1,k),\ (i,j+1,k+1),\ (i,j+1,k) \right) \\
    \sigma_7 &= \left( (i,j,k+1),\ (i+1,j+1,k+1),\ (i+1,j,k+1),\ (i,j,k) \right) \\
    \sigma_8 &= \left( (i+1,j+1,k+1),\ (i+1,j,k+1),\ (i+1,j,k),\ (i,j,k) \right) \\
    \sigma_9 &= \left( (i+1,j+1,k+1),\ (i+1,j+1,k),\ (i+1,j,k),\ (i,j,k) \right) \\
    \sigma_{10} &= \left( (i,j,k+1),\ (i+1,j,k+1),\ (i+1,j,k),\ (i,j,k) \right) \\
    \sigma_{11} &= \left( (i,j,k),\ (i,j+1,k),\ (i,j,k+1),\ (i+1,j+1,k+1) \right) \\
    \sigma_{12} &= \left( (i,j,k),\ (i+1,j,k),\ (i,j,k+1),\ (i+1,j+1,k+1) \right)
    \end{aligned}
    \]
    *Note: This list corrects apparent typos from the original (e.g., in \(\sigma_8\) and \(\sigma_9\)) to ensure all lists contain four distinct vertices that form a non-degenerate tetrahedron, consistent with the triangulation's logic.*

3.  **Global Complex:** The full simplicial complex \(K_n^{(3)}\) is the union of all such tetrahedra over all indices:
    \[
    K_n^{(3)} = \left\{ \sigma_{m}(i,j,k) \mid (i,j,k) \in I_n^3,\ m = 1, \dots, 12 \right\}.
    \]
    The geometric realization is given by the vertex map \(\varphi\).

