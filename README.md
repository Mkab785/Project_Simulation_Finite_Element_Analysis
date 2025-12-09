# 🏗️ 2D Truss Analysis using the Finite Element Method (FEM)

## Project Overview

This repository contains a full, **from-scratch implementation** of the **Finite Element Method (FEM)** algorithm designed for static analysis of **2D Truss Structures**.

The primary goal of this solver is to establish and solve the global equilibrium equation, $[K]\{q\} = \{F\}$, to determine the unknown nodal **displacements** ($\{q\}$) under external loads ($\{F\}$) and imposed boundary conditions.

---

## ⚠️ Core Engineering Assumption (Truss Element)

This solver is based on the fundamental principles of truss analysis, which relies on a critical assumption:

* **Pure Axial Loading:** Each bar element is modeled exclusively in **tension or compression**. The element's stiffness matrix ($[K_e]$) is derived assuming that loads are applied only along the element's axis.
* **Neglect of Bending:** The model **does not account for bending moments, shear forces, or torsion** within the elements. This is suitable for structures where members are connected by idealized pins (hinges).

This specific focus validates the use of simplified, yet powerful, 2D bar elements in the finite element framework. 

---

## 🛠️ Technical Implementation Steps

The solution demonstrates expertise in numerical methods and scientific computing by implementing all key stages of the FEM:

1.  **Input Data Handling:** Interactive collection and validation of geometric data (Nodal Coordinates), material properties (Young's Modulus, $E$, and Cross-sectional Area, $A$), and member connectivity.
2.  **Element Geometry:** Calculation of the **length ($L$)** and **orientation angle ($\theta$)** for every bar.
3.  **Stiffness Matrix Assembly:**
    * Calculation of the elemental stiffness matrix ($[K_e]$) in the global coordinate system.
    * Direct assembly of the full **Global Stiffness Matrix ($[K]$)**.
4.  **Boundary Conditions (BCs):** Application of displacement constraints (fixed or roller supports) and known external forces ($\{F\}$).
5.  **System Reduction:** Identification and elimination of constrained Degrees of Freedom (DOFs) to formulate the **Reduced System** ($[K_{red}]\{q_{red}\} = \{F_{red}\}$).
6.  **Solution:** Solving the linear system for the unknown nodal displacements ($\{q_{red}\}$) using optimized linear algebra routines (e.g., the backslash operator (`\`) in MATLAB).
7.  **Post-Processing:** Visualization of the structure's behavior, comparing the initial geometry with the **deformed geometry** (magnified for clarity).

---

## 💻 Languages & Availability

The original implementation was developed in MATLAB, and was subsequently translated to Python to demonstrate language flexibility and adaptability to engineering environments. The Python version also allows for easier integration into modern data science workflows and simplified visualization capabilities.

| Language | File | Key Numerical Libraries |
| :--- | :--- | :--- |
| **Python** | `MiniProject_FEA_LEGROS_KABIR.ipynb` | `NumPy`, `Matplotlib` |
| **MATLAB** | `truss_fem_solver.m` | Native Matrix Operations (`\`, `zeros`) |

### How to Run:

1.  Clone this repository.
2.  Open the file (`.ipynb` for Python or `.m` for MATLAB) in the respective environment.
3.  Run the script and follow the interactive prompts to define your truss structure.
