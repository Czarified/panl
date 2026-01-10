import numpy as np
import pytest

from panelyze.analysis.geometry import PanelGeometry
from panelyze.analysis.kernels import BEMKernels
from panelyze.analysis.material import OrthotropicMaterial
from panelyze.analysis.solver import BEMSolver


@pytest.fixture
def sample_setup():
    # Setup for a simple unit tension case
    E, nu = 10.0e6, 0.3
    thickness = 0.1
    mat = OrthotropicMaterial(
        e1=E, e2=E * 1.001, nu12=nu, g12=E / (2 * (1 + nu)), thickness=thickness
    )
    kernels = BEMKernels(mat)

    W, H = 2.0, 2.0
    geom = PanelGeometry(W, H)
    # Simple mesh (4 elements)
    elements = geom.discretize(num_elements_per_side=1)
    solver = BEMSolver(kernels, elements)
    solver.assemble()

    return solver, W, H, thickness


def test_solve_running_load_vs_stress(sample_setup):
    solver, W, H, h = sample_setup
    elements = solver.elements

    # 1. Solve with running loads (default)
    bc_type = np.zeros(2 * len(elements), dtype=int)
    bc_value = np.zeros(2 * len(elements))
    q_applied = 500.0  # lbf/in

    for i, el in enumerate(elements):
        if np.isclose(el.center[0], 0.0):
            bc_value[2 * i] = -q_applied
        if np.isclose(el.center[0], W):
            bc_value[2 * i] = q_applied

    # Simple displacement constraints to avoid singular A
    bc_type[0:2] = 1  # u,v fixed
    bc_value[0:2] = 0.0
    bc_type[3] = 1  # v fixed at (W,0)
    bc_value[3] = 0.0

    u1, t1 = solver.solve(bc_type, bc_value, is_stress=False)

    # 2. Solve with stress
    bc_stress_value = np.copy(bc_value)
    # Convert running load indices to stress
    for i in range(len(bc_value)):
        if bc_type[i] == 0:
            bc_stress_value[i] = bc_value[i] / h

    # Note: bc_value in solver.solve(is_stress=True) is treated as traction directly
    u2, t2 = solver.solve(bc_type, bc_stress_value, is_stress=True)

    # Results should be identical
    np.testing.assert_allclose(u1, u2)
    np.testing.assert_allclose(t1, t2)


def test_compute_resultants(sample_setup):
    solver, W, H, h = sample_setup
    elements = solver.elements

    bc_type = np.zeros(2 * len(elements), dtype=int)
    bc_value = np.zeros(2 * len(elements))
    q = 100.0
    for i, el in enumerate(elements):
        if np.isclose(el.center[0], W):
            bc_value[2 * i] = q

    # Constraints
    bc_type[0:2] = 1
    bc_type[3] = 1

    u, t = solver.solve(bc_type, bc_value)

    pts = np.array([[W / 2, H / 2]])
    stresses = solver.compute_stress(pts, u, t)
    resultants = solver.compute_resultants(pts, u, t)

    np.testing.assert_allclose(resultants, stresses * h)


def test_compute_displacement(sample_setup):
    solver, W, H, h = sample_setup
    elements = solver.elements

    bc_type = np.zeros(2 * len(elements), dtype=int)
    bc_value = np.zeros(2 * len(elements))
    bc_type[0:2] = 1  # fix corner

    u, t = solver.solve(bc_type, bc_value)

    pts = np.array([[W / 2, H / 2]])
    disp = solver.compute_displacement(pts, u, t)
    assert disp.shape == (1, 2)


def test_kernels_direct(sample_setup):
    solver, _, _, _ = sample_setup
    kernels = solver.kernels

    # Target kernels.py coverage
    U = kernels.displacement_kernel(0.5, 0.5)
    assert U.shape == (2, 2)

    T = kernels.traction_kernel(0.5, 0.5, 1.0, 0.0)
    assert T.shape == (2, 2)


def test_elliptical_cutout():
    from panelyze.analysis.geometry import EllipticalCutout

    ell = EllipticalCutout(0, 0, 2.0, 1.0, 45.0)
    els = ell.discretize(10)
    assert len(els) == 10
    assert np.isclose(els[0].p1[0], 2.0 * np.cos(np.radians(45.0)))
