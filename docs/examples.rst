Getting Started: Isotropic Plate Example
========================================

This example walks through setting up a simple stress analysis for a square isotropic plate with a circular cutout in the center, subjected to uniaxial tension.

Problem Description
-------------------

We will analyze a 100mm x 100mm square plate with a 10mm diameter hole at the center. The material is aluminum (isotropic), and we apply a 100 MPa tensile stress to the left and right edges.

Free Body Diagram (FBD)
-----------------------

.. mermaid::

   graph LR
       subgraph Plate
       A[Left Edge: -100 MPa] --> B[Square Panel]
       B --> C[Right Edge: +100 MPa]
       B -- contains -- D((Circular Hole r=5mm))
       end

Step-by-Step Implementation
---------------------------

1. Define the Material
~~~~~~~~~~~~~~~~~~~~~~

For an isotropic material, we set :math:`E_1 = E_2` and calculate :math:`G_{12}` from :math:`E` and :math:`\nu`.

.. code-block:: python

   import numpy as np
   from panelyze.analysis.material import OrthotropicMaterial

   E = 70000.0  # MPa (Aluminum)
   nu = 0.33
   G = E / (2 * (1 + nu))

   # Isotropic material represented as orthotropic
   mat = OrthotropicMaterial(e1=E, e2=E, nu12=nu, g12=G)

2. Create Geometry and Mesh
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Define the panel dimensions and add a circular cutout.

.. code-block:: python

   from panelyze.analysis.geometry import PanelGeometry, CircularCutout

   W, H = 100.0, 100.0
   radius = 5.0
   geom = PanelGeometry(W, H)
   geom.add_cutout(CircularCutout(x_center=W/2, y_center=H/2, radius=radius))

   # Discretize the boundary (coarse mesh for speed)
   elements = geom.discretize(num_elements_per_side=10, num_elements_cutout=40)

3. Assemble and Solve
~~~~~~~~~~~~~~~~~~~~~

Using the `BEMKernels` and `BEMSolver` to find the unknown boundary values.

.. code-block:: python

   from panelyze.analysis.kernels import BEMKernels
   from panelyze.analysis.solver import BEMSolver

   kernels = BEMKernels(mat)
   solver = BEMSolver(kernels, elements)
   solver.assemble()

   # Define Boundary Conditions (BCs)
   # indices: 0:2*N_panel: BCs for outer boundary, then for cutouts
   # Let's apply tension on x-edges (sides)
   num_dofs = 2 * len(elements)
   bc_type = np.zeros(num_dofs, dtype=int)  # 0 = Traction given
   bc_value = np.zeros(num_dofs)

   # Find elements on left and right edges
   for i, el in enumerate(elements):
       if np.isclose(el.center[0], 0.0): # Left Edge
           bc_value[2*i] = -100.0 # Tension out
       elif np.isclose(el.center[0], W): # Right Edge
           bc_value[2*i] = 100.0

   # Constrain a point to prevent rigid body motion
   # (In BEM, sometimes single point constraints are needed if H is singular)
   # For this example, we assume the solver's internal rigid body trick handles it

   u, t = solver.solve(bc_type, bc_value)

4. Extract Results
~~~~~~~~~~~~~~~~~~

Evaluate the stress at the stress concentration point (tip of the hole at :math:`\theta = 90^\circ`).

.. code-block:: python

   # Point at theta=90 deg relative to hole center
   eval_pt = np.array([[W/2, H/2 + 5.05]])
   stresses = solver.compute_stress(eval_pt, u, t)

   print(f"Sigma_xx at hole pole: {stresses[0, 0]:.2f} MPa")
   # Expected SCF for isotropic infinite plate is 3.0.
   # For finite plate, it will be slightly higher.

Verification with xdoctest
--------------------------

The following block is a testable example that can be run with `nox -s docs`.

.. code-block:: python

   >>> import numpy as np
   >>> from panelyze.analysis.material import OrthotropicMaterial
   >>> from panelyze.analysis.geometry import PanelGeometry, CircularCutout
   >>> from panelyze.analysis.kernels import BEMKernels
   >>> from panelyze.analysis.solver import BEMSolver
   >>> E, nu = 70000.0, 0.33
   >>> G = E / (2 * (1 + nu))
   >>> mat = OrthotropicMaterial(e1=E, e2=E, nu12=nu, g12=G)
   >>> geom = PanelGeometry(100.0, 100.0)
   >>> geom.add_cutout(CircularCutout(50.0, 50.0, 5.0))
   >>> elements = geom.discretize(num_elements_per_side=20, num_elements_cutout=80)
   >>> solver = BEMSolver(BEMKernels(mat), elements)
   >>> solver.assemble()
   >>> bc_type = np.zeros(2 * len(elements), dtype=int)
   >>> bc_value = np.zeros(2 * len(elements))
   >>> for i, el in enumerate(elements):
   ...     if np.isclose(el.center[0], 0.0): bc_value[2*i] = -100.0
   ...     if np.isclose(el.center[0], 100.0): bc_value[2*i] = 100.0
   >>> u, t = solver.solve(bc_type, bc_value)
   >>> eval_pts = np.array([[50.0, 55.05]])
   >>> stress = solver.compute_stress(eval_pts, u, t)
   >>> # Check if SCF is around 3.0 (approximate due to finite size)
   >>> scf = stress[0, 0] / 100.0
   >>> print(f"{scf:.1f}")
   3.1
