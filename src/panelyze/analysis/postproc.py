from typing import Tuple

import numpy as np

from .solver import BEMSolver


def calculate_scf(
    solver: BEMSolver,
    u_boundary: np.ndarray,
    t_boundary: np.ndarray,
    nominal_stress: float,
) -> float:
    """
    Calculates the Stress Concentration Factor around the cutout.

    Args:
        solver: The BEM solver instance.
        u_boundary: Solved boundary displacements.
        t_boundary: Solved boundary tractions.
        nominal_stress: The far-field stress to normalize by.

    Returns:
        float: The calculated Stress Concentration Factor (SCF).
    """
    # Assuming the cutout elements are at the end of the elements list
    # (This depends on how geometry.discretize is implemented)
    # Better: Identify cutout elements.

    # Evaluate stresses at points slightly shifted from the boundary into the domain
    # to avoid boundary singularities in the displacement gradient.
    # For midpoint elements, we can shift along the inward normal.

    max_vm_stress = 0.0

    for el in solver.elements:
        # Check if element is part of a cutout (simple heuristic for now)
        # In a real tool, we'd tag elements by boundary type.
        # For now, let's just check all elements and find max stress.

        # Shift point slightly inward (normal points OUT, so inward is -n)
        P = el.center - 0.1 * el.length * np.array([el.nx, el.ny])

        stress = solver.compute_stress(np.array([P]), u_boundary, t_boundary)[0]
        sxx, syy, txy = stress

        # Von Mises stress (Plane Stress)
        vm = np.sqrt(sxx**2 - sxx * syy + syy**2 + 3 * txy**2)
        if vm > max_vm_stress:
            max_vm_stress = vm

    return max_vm_stress / nominal_stress


def get_stress_field(
    solver: BEMSolver,
    u_boundary: np.ndarray,
    t_boundary: np.ndarray,
    x_range: Tuple[float, float],
    y_range: Tuple[float, float],
    grid_res: int = 50,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Computes stress field components on a grid.

    Args:
        solver: The BEM solver instance.
        u_boundary: Solved boundary displacements.
        t_boundary: Solved boundary tractions.
        x_range: Tuple of (min, max) x-coordinates.
        y_range: Tuple of (min, max) y-coordinates.
        grid_res: Resolution of the grid.

    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray]: X grid, Y grid,
            and stress tensor field.
    """
    xv = np.linspace(x_range[0], x_range[1], grid_res)
    yv = np.linspace(y_range[0], y_range[1], grid_res)
    X, Y = np.meshgrid(xv, yv)

    points = np.vstack([X.ravel(), Y.ravel()]).T
    stresses = solver.compute_stress(points, u_boundary, t_boundary)

    return X, Y, stresses.reshape((grid_res, grid_res, 3))


def calculate_scf_at_points(
    solver: BEMSolver,
    u_boundary: np.ndarray,
    t_boundary: np.ndarray,
    points: np.ndarray,
    nominal_stress: float,
) -> float:
    """
    Calculates max Von Mises stress at specific points.

    Args:
        solver: The BEM solver instance.
        u_boundary: Solved boundary displacements.
        t_boundary: Solved boundary tractions.
        points: Array of points to evaluate.
        nominal_stress: The far-field stress to normalize by.

    Returns:
        float: The calculated Peak Stress Concentration Factor (SCF).
    """
    stresses = solver.compute_stress(points, u_boundary, t_boundary)
    vms = np.sqrt(
        stresses[:, 0] ** 2
        - stresses[:, 0] * stresses[:, 1]
        + stresses[:, 1] ** 2
        + 3 * stresses[:, 2] ** 2
    )
    return np.max(vms) / nominal_stress
