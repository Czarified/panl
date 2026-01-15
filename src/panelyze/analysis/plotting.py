from typing import List, Optional

import matplotlib.pyplot as plt
import numpy as np

from .geometry import BoundaryElement
from .solver import BEMSolver


def plot_results(
    solver: BEMSolver,
    u_boundary: np.ndarray,
    t_boundary: np.ndarray,
    deform_scale: float = 0.0,
    stress_type: str = "vm",
    show_mesh: bool = True,
    show_labels: bool = True,
    title: Optional[str] = None,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """
    Plots the panel geometry, mesh, and results.

    Args:
        solver: The BEM solver instance.
        u_boundary: Solved boundary displacements.
        t_boundary: Solved boundary tractions.
        deform_scale: Scaling factor for displacements. If 0, no deformation is shown.
        stress_type: Type of stress to evaluate (vm, xx, yy, xy, principal).
        show_mesh: Whether to show the boundary element discretization.
        show_labels: Whether to highlight and label max stress points.
        title: Optional title for the plot.
        save_path: Optional path to save the figure to.

    Returns:
        plt.Figure: The generated figure.
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    # 1. Draw Geometry and Mesh
    _draw_boundary(ax, solver.elements, u_boundary, deform_scale, show_mesh)

    # 2. Highlight Max Stress Points
    if show_labels:
        _highlight_max_stresses(ax, solver, u_boundary, t_boundary, stress_type)

    # 3. Aesthetics
    ax.set_aspect("equal")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    if title:
        ax.set_title(title)
    else:
        ax.set_title("Panel Analysis Results")

    ax.grid(True, linestyle="--", alpha=0.6)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig


def _draw_boundary(
    ax: plt.Axes,
    elements: List[BoundaryElement],
    u_boundary: np.ndarray,
    deform_scale: float,
    show_mesh: bool,
):
    """Draws the Undeformed and Deformed (if scale > 0) boundaries.

    Args:
        ax: Matplotlib axes object.
        elements: List of boundary elements.
        u_boundary: Solved boundary displacements.
        deform_scale: Scaling factor for displacements.
        show_mesh: Whether to show the boundary element discretization.
    """

    # Organize elements into contiguous boundary loops
    # (For now, we assume elements are already ordered by PanelGeometry.discretize)
    # We'll just draw each element separately for simplicity, but could be optimized.

    for i, el in enumerate(elements):
        # Undeformed
        p1 = el.p1
        p2 = el.p2

        # Deformed
        if deform_scale > 0:
            # Constant element: use midpoint displacement for both ends?
            # Or use displacement values at nodes if we had them.
            # Since it's a constant element BEM, the solved u is at the center.
            # To show a mesh-like deformation, we'd need nodal values.
            # For constant elements, we can shift the entire segment.
            u_j = u_boundary[2 * i : 2 * i + 2]
            p1_def = p1 + u_j * deform_scale
            p2_def = p2 + u_j * deform_scale

            # Draw deformed (Dashed red)
            ax.plot(
                [p1_def[0], p2_def[0]], [p1_def[1], p2_def[1]], "r-", lw=1.5, alpha=0.8
            )

        # Draw undeformed (Solid black)
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], "k-", lw=1.0)

        if show_mesh:
            # Draw nodes/midpoints
            ax.plot(el.center[0], el.center[1], "k.", markersize=4)
            if deform_scale > 0:
                u_j = u_boundary[2 * i : 2 * i + 2]
                c_def = el.center + u_j * deform_scale
                ax.plot(c_def[0], c_def[1], "r.", markersize=4)


def _highlight_max_stresses(
    ax: plt.Axes,
    solver: BEMSolver,
    u_boundary: np.ndarray,
    t_boundary: np.ndarray,
    stress_type: str,
):
    """Identifies and labels max stress points on cutout boundaries.

    Args:
        ax: Matplotlib axes object.
        solver: The BEM solver instance.
        u_boundary: Solved boundary displacements.
        t_boundary: Solved boundary tractions.
        stress_type: Type of stress to evaluate (vm, xx, yy, xy, principal).

    Raises:
        ValueError: If stress_type is not recognized.
    """
    # Find elements that likely belong to cutouts.
    # Better: use tagging if available. For now, check all.

    stress_values = []
    points = []

    boundary_stresses = solver.compute_boundary_stress(u_boundary, t_boundary)

    for i, el in enumerate(solver.elements):
        if el.tag == "cutout":
            sxx, syy, txy = boundary_stresses[i]

            if stress_type.lower() == "vm":
                val = np.sqrt(sxx**2 - sxx * syy + syy**2 + 3 * txy**2)
                label = "Max VM"
            elif stress_type.lower() == "xx":
                val = sxx
                label = "Max σxx"
            elif stress_type.lower() == "yy":
                val = syy
                label = "Max σyy"
            elif stress_type.lower() == "xy":
                val = abs(txy)
                label = "Max |τxy|"
            elif stress_type.lower() == "principal":
                val = (sxx + syy) / 2 + np.sqrt(((sxx - syy) / 2) ** 2 + txy**2)
                label = "Max P1"
            else:
                raise ValueError(f"Unknown stress_type: {stress_type}")

            stress_values.append(val)
            points.append(el.center)

    stress_values = np.array(stress_values)
    max_idx = np.argmax(stress_values)
    max_val = stress_values[max_idx]
    max_pt = points[max_idx]

    # Label the max stress point
    ax.annotate(
        f"{label}: {max_val:.2f}",
        xy=(max_pt[0], max_pt[1]),
        xytext=(
            max_pt[0] + 0.1 * solver.kernels.mat.thickness,
            max_pt[1] + 0.1 * solver.kernels.mat.thickness,
        ),
        arrowprops=dict(facecolor="black", shrink=0.05, width=1, headwidth=5),
        bbox=dict(boxstyle="round,pad=0.3", fc="yellow", ec="k", lw=1, alpha=0.7),
        fontsize=9,
        fontweight="bold",
    )

    # Highlight the point
    ax.plot(
        max_pt[0],
        max_pt[1],
        "go",
        markersize=8,
        markeredgecolor="k",
        label="Max Stress",
    )


def plot_deformed_shape(
    solver: BEMSolver, u_boundary: np.ndarray, scale: float = 100.0, **kwargs
) -> plt.Figure:
    """Convenience wrapper for deformed plotting.

    Args:
        solver: The BEM solver instance.
        u_boundary: Solved boundary displacements.
        scale: Scaling factor for displacements.
        **kwargs: Additional arguments passed to plot_results.

    Returns:
        plt.Figure: The generated figure.
    """
    return plot_results(
        solver, u_boundary, np.zeros_like(u_boundary), deform_scale=scale, **kwargs
    )
