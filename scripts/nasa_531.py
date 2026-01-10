import numpy as np

from panelyze.analysis.geometry import CircularCutout, PanelGeometry
from panelyze.analysis.kernels import BEMKernels
from panelyze.analysis.material import OrthotropicMaterial
from panelyze.analysis.solver import BEMSolver


def nasa_5_3_1_validation():
    # Properties (NASA 5.3.1)
    E1, E2, nu12, G12 = 10000.0, 5000.0, 0.3, 3000.0
    mat = OrthotropicMaterial(e1=E1, e2=E2, nu12=nu12, g12=G12)
    kernels = BEMKernels(mat)

    W, H = 200.0, 200.0
    radius = 5.0
    geom = PanelGeometry(W, H)
    geom.add_cutout(CircularCutout(W / 2, H / 2, radius))

    # Mesh (Ultra Dense)
    elements = geom.discretize(num_elements_per_side=20, num_elements_cutout=400)
    solver = BEMSolver(kernels, elements)
    solver.assemble()

    # BCs
    bc_type = np.zeros(2 * len(elements), dtype=int)
    bc_value = np.zeros(2 * len(elements))
    bc_value[2 * 20 : 2 * 40 : 2] = 100.0
    bc_value[2 * 60 : 2 * 80 : 2] = -100.0
    bc_type[20:22] = 1
    bc_value[20:22] = 0.0
    bc_type[100] = 1
    bc_value[100] = 0.0

    u, t = solver.solve(bc_type, bc_value)

    # Peak Profiling at r=5.05 (1.01*R)
    angles = np.array([0, np.pi / 4, np.pi / 2])
    pts_circ = np.array(
        [[W / 2 + 5.05 * np.cos(a), H / 2 + 5.05 * np.sin(a)] for a in angles]
    )

    st_circ = solver.compute_stress(pts_circ, u, t)
    print("Sigma_xx profile at r=5.05:")
    for i, a in enumerate(angles):
        print(f"  theta={np.degrees(a):.0f}: {st_circ[i, 0]:.4f}")

    print(f"Calculated Peak SCF: {st_circ[2, 0]/100.0:.4f}")

    expected_scf = 1 + np.sqrt(2 * (np.sqrt(E1 / E2) - nu12) + E1 / G12)
    print(f"Theoretical Infinite SCF: {expected_scf:.4f}")


if __name__ == "__main__":
    nasa_5_3_1_validation()
