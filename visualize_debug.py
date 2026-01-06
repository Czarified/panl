import numpy as np

from panelyze.analysis.geometry import PanelGeometry
from panelyze.analysis.kernels import BEMKernels
from panelyze.analysis.material import OrthotropicMaterial
from panelyze.analysis.solver import BEMSolver


def visualize_stress():
    # NASA Example 5.3.1
    E1, E2, nu12, G12 = 10000.0, 5000.0, 0.3, 3000.0
    mat = OrthotropicMaterial(e1=E1, e2=E2, nu12=nu12, g12=G12)
    kernels = BEMKernels(mat)

    W, H = 100.0, 100.0
    geom = PanelGeometry(W, H)

    N = 40  # per side
    elements = geom.discretize(num_elements_per_side=N)
    solver = BEMSolver(kernels, elements)
    solver.assemble()

    bc_type = np.zeros(2 * len(elements), dtype=int)
    bc_value = np.zeros(2 * len(elements))

    # Tension X: Right (idx 80 to 159), Left (idx 240 to 319)
    # N=40. Right: ele 40-79. idx 80-159.
    # Left: ele 120-159. idx 240-319.
    bc_value[80:160:2] = 100.0
    bc_value[240:320:2] = -100.0

    # Pinning:
    # Center of Bottom edge (ele 20): idx 40, 41 (u,v=0)
    bc_type[40:42] = 1
    bc_value[40:42] = 0.0
    # Center of Top edge (ele 100): idx 200 (u=0)
    bc_type[200] = 1
    bc_value[200] = 0.0

    print(f"Material Roots (mu):\n{kernels.mu1}, {kernels.mu2}")

    # Gauss Law Check
    p_int = np.array([50.0, 50.0])
    sum_t11 = 0.0
    for el in elements:
        dx, dy = el.center[0] - p_int[0], el.center[1] - p_int[1]
        T = kernels.traction_kernel(dx, dy, el.nx, el.ny)
        sum_t11 += T[0, 0] * el.length
    print(f"Gauss Law Integral (sum T11 * L): {sum_t11:.6f}")

    u, t = solver.solve(bc_type, bc_value)

    # Right: ele N to 2N-1. idx 2N to 4N-1.
    # Left: ele 3N to 4N-1. idx 6N to 8N-1.
    ux_right = u[2 * N : 4 * N : 2]
    ux_left = u[6 * N : 8 * N : 2]
    uy_top = u[4 * N + 1 : 6 * N : 2]
    print(f"Mean ux Right: {np.mean(ux_right):.6f}")
    print(f"Mean ux Left: {np.mean(ux_left):.6f}")
    print(f"Delta ux: {np.mean(ux_right) - np.mean(ux_left):.6f}")
    print(f"Mean uy Top: {np.mean(uy_top):.6f}")

    # Internal profile
    pts = np.array(
        [[25.0, 50.0], [40.0, 50.0], [50.0, 50.0], [60.0, 50.0], [75.0, 50.0]]
    )
    stresses = solver.compute_stress(pts, u, t)
    print("Sigma_xx profile along Y=50:")
    print(stresses[:, 0])
    print("Expected: Uniform 100.0")


if __name__ == "__main__":
    visualize_stress()
