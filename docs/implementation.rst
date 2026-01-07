Implementation Details
======================

This section describes the software architecture and numerical implementation of the BEM solver in Panelyze.

Core Architecture
-----------------
The solver is divided into four main components:

1. **Material Definition** (`panelyze.analysis.material.OrthotropicMaterial`): Handles the compliance and stiffness matrices and solves for the characteristic roots :math:`\mu_k`.
2. **Kernels** (`panelyze.analysis.kernels.BEMKernels`): Implements the fundamental solutions from NASA CR-1934.
3. **Geometry** (`panelyze.analysis.geometry`): Simple boundary discretization into constant linear elements.
4. **Solver** (`panelyze.analysis.solver.BEMSolver`): Assembles the global matrices :math:`H` and :math:`G` and solves the system.

Matrix Assembly
---------------
The boundary is discretized into :math:`M` elements. For each element :math:`i`, we collocate the boundary integral equation at its center. This leads to a system of equations:

.. math::
   [H] \{u\} = [G] \{t\}

where:
- :math:`\{u\}` is the vector of boundary displacements.
- :math:`\{t\}` is the vector of boundary tractions.

Singular Integrals
~~~~~~~~~~~~~~~~~~
When the collocation point :math:`i` lies on the element :math:`j` being integrated, the kernels become singular (:math:`1/r` for :math:`T` and :math:`\ln(r)` for :math:`U`).
- The singular :math:`G` matrix components (integrals of :math:`\ln z`) are handled analytically using the identity :math:`\int \ln z dz = z(\ln z - 1)`.
- The singular :math:`H` matrix components involve a "jump term" (typically 0.5) and the Cauchy Principal Value (CPV). Panelyze uses the **Rigid-Body Motion trick** to determine the diagonal terms of :math:`H`:

  .. math::
     H_{ii} = - \sum_{j \neq i} H_{ij}

This ensures that a rigid body translation results in zero tractions, which is a physical requirement of the system.

Numerical Integration
~~~~~~~~~~~~~~~~~~~~~
For non-singular elements, the kernels are integrated analytically using the endpoint complex coordinates :math:`z_1` and :math:`z_2`. This approach is more robust than Gauss quadrature for elements very close to each other (e.g., near sharp cutouts).

Internal Point Evaluation
-------------------------
Displacements and stresses at interior points are computed using the Somigliana identity after the boundary values are fully solved.

.. math::
   u_i(\xi) = \sum_{j=1}^M \left( \int_{\Gamma_j} U_{ki} d\Gamma \right) t_k(x) - \sum_{j=1}^M \left( \int_{\Gamma_j} T_{ki} d\Gamma \right) u_k(x)

Note that in the implementation, **Reciprocity** is accounted for. The fundamental solution :math:`U_{ki}(\xi, x)` represents the displacement in direction $i$ at $x$ due to a load at $\xi$. To get the displacement at $\xi$, the kernel results are transposed.

Stress Calculation
~~~~~~~~~~~~~~~~~~
Stresses at interior points are calculated by taking the derivative of the displacement identity. Direct differentiation of the kernels is used to avoid numerical differentiation errors.

.. math::
   \sigma_{ab}(\xi) = C_{abcd} \frac{\partial u_c(\xi)}{\partial x_d}
