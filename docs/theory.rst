Theory of the Boundary Element Method for Orthotropic Panels
============================================================

This section describes the mathematical foundation of the Boundary Element Method (BEM) as applied to orthotropic elasticity in the Panelyze project.

Boundary Element Method (BEM)
-----------------------------
The BEM is a numerical technique based on the transformation of the governing partial differential equations into a surface (or boundary) integral equation. For 2D linear elasticity, this is achieved through **Somigliana's Identity**, which relates the displacement at any point :math:`\xi` (inside the domain or on the boundary) to the boundary displacements :math:`u` and tractions :math:`t`.

The displacement at a point :math:`\xi` is given by:

.. math::
   c_{ji}(\xi) u_i(\xi) + \oint_{\Gamma} T_{ji}(\xi, x) u_i(x) d\Gamma(x) = \oint_{\Gamma} U_{ji}(\xi, x) t_i(x) d\Gamma(x)

where:
- :math:`u_i(x)` and :math:`t_i(x)` are the displacements and tractions on the boundary :math:`\Gamma`.
- :math:`U_{ji}` and :math:`T_{ji}` are the fundamental displacement and traction kernels.
- :math:`c_{ji}` is a coefficient that depends on the local geometry at :math:`\xi`. For points inside the domain, :math:`c_{ji} = \delta_{ji}`. For points on a smooth boundary, :math:`c_{ji} = 0.5 \delta_{ji}`.

Fundamental Solutions for Anisotropic Elasticity
------------------------------------------------
The kernels :math:`U_{ji}` and :math:`T_{ji}` used in this solver are derived from the complex variable formulation for anisotropic elasticity. Unlike isotropic materials, where the fundamental solutions are relatively simple, anisotropic fundamental solutions depend on the roots of a characteristic equation related to the material's elastic constants.

Characteristic Roots (:math:`\mu_k`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
For a plane stress problem, the characteristic equation is:

.. math::
   \beta_{11} \mu^4 - 2 \beta_{16} \mu^3 + (2\beta_{12} + \beta_{66}) \mu^2 - 2 \beta_{26} \mu + \beta_{22} = 0

The roots :math:`\mu_k` (where :math:`k=1,2`) are always complex and occur in conjugate pairs. Panelyze uses the roots with positive imaginary parts.

Displacement Kernel (:math:`U_{ji}`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Following **NASA CR-1934** (Equation 23), the displacement kernel is:

.. math::
   U_{ji} = 2 \text{Re} \{ p_{j1} A_{i1} \ln(z_1) + p_{j2} A_{i2} \ln(z_2) \}

where:
- :math:`z_k = x_1 + \mu_k x_2`.
- :math:`p_{jk}` and :math:`A_{ik}` are complex constants derived from the material's compliance matrix and the characteristic roots.

Traction Kernel (:math:`T_{ji}`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The traction kernel is obtained by differentiating the displacement kernel to find the stress state and then projecting it onto the boundary normal :math:`n`:

.. math::
   T_{ji} = \sigma_{j1} n_1 + \sigma_{j2} n_2

Theoretical References
----------------------
The implementation in Panelyze is primarily based on the following technical reports:

1. **NASA CR-1934 (1971)**: *Boundary-Integral Equation Method for Elasticity and Thermal Elasticity*. This report provides the core formulation for the fundamental solutions used in the `BEMKernels` class.
2. **NASA-CR-125596**: *Advanced Fracture Mechanics Analysis*. Cruse and Besuner detail the application of BEM to complex stress concentration and fracture problems, which is the primary use case for Panelyze (analyzing cutouts in panels).
3. **Cruse, T. A. (1988)**: *Boundary Element Analysis in Computational Fracture Mechanics*. This text provides the theoretical background for the internal point evaluation and the treatment of singular integrals.

Verification (NASA 5.3.1)
-------------------------
The validity of these kernels and their implementation is verified in the `nasa_531.py` validation script, which reproduces the results for a circular hole in an infinite orthotropic plate. The theoretical stress concentration factor (SCF) for an infinite plate is compared against the BEM results to ensure convergence.
