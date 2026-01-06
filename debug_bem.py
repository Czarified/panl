import numpy as np

from panelyze.analysis.geometry import PanelGeometry
from panelyze.analysis.kernels import BEMKernels
from panelyze.analysis.material import OrthotropicMaterial
from panelyze.analysis.solver import BEMSolver


def debug_solver():
    # Orthotropic properties
    E1, E2, nu12, G12 = 10000.0, 5000.0, 0.3, 3000.0
    mat = OrthotropicMaterial(e1=E1, e2=E2, nu12=nu12, g12=G12)
    kernels = BEMKernels(mat)

    print(f"Material Roots (mu): {kernels.mu1}, {kernels.mu2}")
    print(f"Coefficients A (x-load):\n{kernels.A[0]}")
    print(f"Coefficients A (y-load):\n{kernels.A[1]}")

    # 10x10 panel to keep it small
    W, H = 10.0, 10.0
    geom = PanelGeometry(W, H)
    # No cutout for now, just simple tension

    elements = geom.discretize(num_elements_per_side=1)
    solver = BEMSolver(kernels, elements)
    solver.assemble()

    print(f"H matrix:\n{solver.H}")
    print(f"G matrix:\n{solver.G}")

    # 4 elements total
    bc_type = np.zeros(2 * len(elements), dtype=int)
    bc_value = np.zeros(2 * len(elements))

    # Right edge (element 1): tx = 100
    bc_value[2] = 100.0
    # Left edge (element 3): tx = -100
    bc_value[6] = -100.0

    # Pin element 0 (Bottom)
    bc_type[0:2] = 1
    bc_value[0:2] = 0.0

    u, t = solver.solve(bc_type, bc_value)

    print(f"Boundary Displacements (u):\n{u}")
    print(f"Boundary Tractions (t):\n{t}")

    # Check interior stress at center (5,5)
    stress = solver.compute_stress(np.array([[5.0, 5.0]]), u, t)[0]
    print(f"Stress at center: {stress}")
    print("Expected: [100.0, 0.0, 0.0]")


if __name__ == "__main__":
    debug_solver()
