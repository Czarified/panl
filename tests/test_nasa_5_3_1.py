import numpy as np
import pytest

from panelyze.analysis.geometry import CircularCutout, PanelGeometry
from panelyze.analysis.kernels import BEMKernels
from panelyze.analysis.material import OrthotropicMaterial
from panelyze.analysis.postproc import calculate_scf_at_points
from panelyze.analysis.solver import BEMSolver


def test_orthotropic_scf_5_3_1():
    """
    NASA Example 5.3.1: Orthotropic plate with a circular hole.
    Using properties: E1=10000, E2=5000, nu12=0.3, G12=3000.
    """
    E1, E2, nu12, G12 = 10000.0, 5000.0, 0.3, 3000.0
    mat = OrthotropicMaterial(e1=E1, e2=E2, nu12=nu12, g12=G12)
    kernels = BEMKernels(mat)

    # Geometry: 100x100 panel with R=5 hole at center
    W, H = 100.0, 100.0
    geom = PanelGeometry(W, H)
    geom.add_cutout(CircularCutout(50.0, 50.0, 5.0))

    # Discretize: Use dense mesh to ensure accuracy
    elements = geom.discretize(num_elements_per_side=40, num_elements_cutout=80)
    solver = BEMSolver(kernels, geom)
    solver.assemble()

    # Boundary Conditions: Uniform tension in X (sigma = 100)
    bc_type = np.zeros(2 * len(elements), dtype=int)
    bc_value = np.zeros(2 * len(elements))

    # Outer bounds: 0-159 elements
    for i in range(160):
        if i < 40:  # Bottom (y=0)
            bc_type[2 * i] = 0
            bc_type[2 * i + 1] = 0
        elif i < 80:  # Right (x=100)
            bc_type[2 * i] = 0
            bc_value[2 * i] = 100.0
            bc_type[2 * i + 1] = 0
        elif i < 120:  # Top (y=100)
            bc_type[2 * i] = 0
            bc_type[2 * i + 1] = 0
        elif i < 160:  # Left (x=0)
            bc_type[2 * i] = 0
            bc_value[2 * i] = -100.0
            bc_type[2 * i + 1] = 0

    # Fix rigid body motion: Pin bottom-left
    bc_type[0] = 1
    bc_value[0] = 0.0
    bc_type[1] = 1
    bc_value[1] = 0.0

    u, t = solver.solve(bc_type, bc_value)

    # Theoretical SCF for tension in fiber direction (1)
    expected_scf = 1 + np.sqrt(E1 / G12 - 2 * nu12 + 2 * np.sqrt(E1 / E2))

    # Evaluate stresses around the hole boundary
    # Shift slightly into the domain
    theta = np.linspace(0, 2 * np.pi, 100)
    eval_pts = []
    r_eval = 5.0 + 0.1  # R + delta
    for t_val in theta:
        eval_pts.append([50.0 + r_eval * np.cos(t_val), 50.0 + r_eval * np.sin(t_val)])

    eval_pts = np.array(eval_pts)
    scf = calculate_scf_at_points(solver, u, t, eval_pts, nominal_stress=100.0)

    print("\nNASA Example 5.3.1 - Orthotropic Circular Hole")
    print(f"Expected SCF: {expected_scf:.4f}")
    print(f"Calculated SCF: {scf:.4f}")

    # Assert within 10% tolerance
    assert scf == pytest.approx(expected_scf, rel=0.1)


if __name__ == "__main__":
    test_orthotropic_scf_5_3_1()
