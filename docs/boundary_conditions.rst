Boundary Conditions
===================

This section outlines the default boundary conditions used in the ``panelyze`` BEM solver and their impact on simulation results.

Overview
--------

To simplify usage and ensure models are properly constrained against rigid-body motion (RBM), the ``BEMSolver.solve`` method implements a set of default boundary conditions when explicit values are not provided. These defaults are tailored for rectangular panel analysis under standard loading conditions (tension, shear).

Default Constraints (RBM Suppression)
-------------------------------------

BEM simulations require sufficient displacement constraints to prevent the model from "flying away" due to numerical or physical imbalances (Rigid-Body Motion). By default, the following pinning logic is applied to the outer boundary:

- **Lower-Left Corner**: Fixed in both global X and Y directions ($u_x=0, u_y=0$).
- **Lower-Right Corner**: Fixed in the global Y direction ($u_y=0$).

These constraints allow the panel to expand or contract freely in the X-direction while preventing translation and rotation.

.. note::
   The solver identifies "corners" by searching for the elements closest to the geometric Extents (min/max X and Y) of the elements tagged as ``"outer"``.

Default Load Convention
-----------------------

When the ``qx``, ``qy``, or ``qxy`` parameters are used in ``solve()``, the following running load (force per unit length) convention is adopted, matching standard in-plane element formulations (e.g., Nastran CQUAD4):

+-----------------+---------------------------+---------------------------------------------+
| Load Parameter  | Description               | Boundary Application                        |
+=================+===========================+=============================================+
| ``qx``          | Uniaxial Tension in X     | $+F_x$ on Max-X edge, $-F_x$ on Min-X edge  |
+-----------------+---------------------------+---------------------------------------------+
| ``qy``          | Uniaxial Tension in Y     | $+F_y$ on Max-Y edge, $-F_y$ on Min-Y edge  |
+-----------------+---------------------------+---------------------------------------------+
| ``qxy``         | In-plane Shear            | See Shear Convention below                  |
+-----------------+---------------------------+---------------------------------------------+

Shear Convention (``qxy``)
~~~~~~~~~~~~~~~~~~~~~~~~~~

The shear load ``qxy`` is applied as a balanced set of tractions on all four outer edges:
- **Right Edge (Max-X)**: $+t_y$
- **Left Edge (Min-X)**: $-t_y$
- **Top Edge (Max-Y)**: $+t_x$
- **Bottom Edge (Min-Y)**: $-t_x$

Impact of Assumptions
---------------------

- **Geometry Alignment**: The default BC logic assumes the panel edges are aligned with the global X and Y axes. For rotated panels, the coordinate-based search for "horizontal" and "vertical" edges may yield unexpected results.
- **Homogeneity**: The RBM suppression logic fixes single elements. In a perfectly balanced theoretical model, this has negligible impact on the internal stress field (Saint-Venant's Principle). However, close to the pinned corners, local stress artifacts may be present.
- **Applicability**: These defaults are highly applicable to standard coupon-level analysis (e.g., hole-in-plate problems). For complex built-up structures where the panel interfaces with other components, users should define explicit boundary conditions.

Future Refactor Path
--------------------

The current implementation uses coordinate-based element identification. As the codebase evolves, we plan to move towards a more robust nodal or geometric set-based boundary condition system. Code comments in ``solver.py`` highlight these areas for ease of future refactoring.
