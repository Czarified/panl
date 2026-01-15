import matplotlib.pyplot as plt
import numpy as np
import pytest

from panelyze.analysis.geometry import CircularCutout, PanelGeometry
from panelyze.analysis.kernels import BEMKernels
from panelyze.analysis.material import OrthotropicMaterial
from panelyze.analysis.plotting import plot_deformed_shape, plot_results
from panelyze.analysis.solver import BEMSolver


@pytest.fixture
def solved_system():
    # Setup a simple system: Square panel with a circular cutout under tension
    E, nu = 10.0e6, 0.3
    thickness = 0.1
    mat = OrthotropicMaterial(
        e1=E, e2=E * 1.01, nu12=nu, g12=E / (2 * (1 + nu)), thickness=thickness
    )
    kernels = BEMKernels(mat)

    W, H = 10.0, 10.0
    geom = PanelGeometry(W, H)
    geom.add_cutout(CircularCutout(W / 2, H / 2, 1.0))

    elements = geom.discretize(num_elements_per_side=4, num_elements_cutout=16)
    solver = BEMSolver(kernels, geom)
    solver.assemble()

    # Boundary Conditions: Tension in X
    bc_type = np.zeros(2 * len(elements), dtype=int)
    bc_value = np.zeros(2 * len(elements))
    q = 1000.0

    for i, el in enumerate(elements):
        if np.isclose(el.center[0], W):
            bc_value[2 * i] = q
        elif np.isclose(el.center[0], 0.0):
            bc_value[2 * i] = -q

    # Fix a point to prevent rigid body motion
    bc_type[0:2] = 1  # fix (0,0)
    bc_value[0:2] = 0.0

    u, t = solver.solve(bc_type, bc_value)
    return solver, u, t


def test_plot_results_basic(solved_system):
    solver, u, t = solved_system
    fig = plot_results(solver, u, t, show_mesh=True, show_labels=True)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_deformed_shape(solved_system):
    solver, u, t = solved_system
    fig = plot_deformed_shape(solver, u, scale=100.0)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_results_save(solved_system, tmp_path):
    solver, u, t = solved_system
    save_file = tmp_path / "test_plot.png"
    fig = plot_results(solver, u, t, save_path=str(save_file))
    assert save_file.exists()
    plt.close(fig)


def test_plot_results_no_labels(solved_system):
    solver, u, t = solved_system
    fig = plot_results(solver, u, t, show_labels=False)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


@pytest.mark.parametrize("s_type", ["vm", "xx", "yy", "xy", "principal"])
def test_plot_results_stress_types(solved_system, s_type):
    solver, u, t = solved_system
    fig = plot_results(solver, u, t, stress_type=s_type)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_results_invalid_type(solved_system):
    solver, u, t = solved_system
    with pytest.raises(ValueError, match="Unknown stress_type"):
        plot_results(solver, u, t, stress_type="invalid")
