import numpy as np

from panelyze.analysis.geometry import CircularCutout, PanelGeometry
from panelyze.analysis.kernels import BEMKernels
from panelyze.analysis.material import OrthotropicMaterial
from panelyze.analysis.postproc import calculate_scf
from panelyze.analysis.solver import BEMSolver


def test_isotropic_scf():
    # Isotropic material properties
    E = 10000.0
    nu = 0.3
    G = E / (2 * (1 + nu))

    # Define as quasi-orthotropic (slightly perturbed to avoid singularity)
    mat = OrthotropicMaterial(e1=10000.0, e2=10005.0, nu12=nu, g12=G)
    kernels = BEMKernels(mat)

    # Geometry: 100x100 panel with R=5 hole at center
    W, H = 100.0, 100.0
    geom = PanelGeometry(W, H)
    geom.add_cutout(CircularCutout(50.0, 50.0, 5.0))

    # Discretize
    elements = geom.discretize(num_elements_per_side=40, num_elements_cutout=80)
    solver = BEMSolver(kernels, geom)
    solver.assemble()

    # Boundary Conditions: Uniform tension in X (sigma = 100)
    bc_type = np.zeros(2 * len(elements), dtype=int)
    bc_value = np.zeros(2 * len(elements))

    for i, el in enumerate(elements):
        if i < 40:  # Bottom
            bc_type[2 * i] = 0
            bc_type[2 * i + 1] = 0
        elif i < 80:  # Right (applied tension sigma_x = 100)
            bc_type[2 * i] = 0
            bc_value[2 * i] = 100.0
            bc_type[2 * i + 1] = 0
        elif i < 120:  # Top
            bc_type[2 * i] = 0
            bc_type[2 * i + 1] = 0
        elif i < 160:  # Left (tx = -100)
            bc_type[2 * i] = 0
            bc_value[2 * i] = -100.0
            bc_type[2 * i + 1] = 0
        else:  # Cutout boundary (Free)
            bc_type[2 * i] = 0
            bc_type[2 * i + 1] = 0

    # Pin the model
    bc_type[0] = 1
    bc_value[0] = 0.0
    bc_type[1] = 1
    bc_value[1] = 0.0

    u, t = solver.solve(bc_type, bc_value)

    scf = calculate_scf(solver, u, t, nominal_stress=100.0)
    print(f"Calculated SCF: {scf:.4f}")
    print("Expected SCF (Isotropic): ~3.0")


if __name__ == "__main__":
    test_isotropic_scf()
