"""Results Validation Against Peterson's Textbook - Example 4.24."""

import matplotlib.pyplot as plt
import numpy as np
from common import get_material, get_output_dir, peterson_4_24
from rich.console import Console
from rich.table import Table

from panl.analysis import plot_results
from panl.analysis.geometry import CircularCutout, PanelGeometry
from panl.analysis.kernels import BEMKernels
from panl.analysis.solver import BEMSolver


def run():
    """Run Peterson's Example 4.24."""
    # Setup
    mat = get_material()
    output_dir = get_output_dir("peterson_validation")

    # Example 4.24, Basic Panel Dimensions
    d = 2.0
    height = d * 10
    width = d * 15

    # Constant and Sweep parameters
    q_x = 500
    q_y = 500
    n_side = 40
    aspect_ratios = [1.5, 2.0, 2.5, 3.0, 4.0, 5.0]
    n_cutout_elements = [44, 220]
    shapes = ["o", "x"]

    # Plot the results
    results = []
    fig, ax = plt.subplots(figsize=(10, 6))

    # Peterson's SCFs
    k_tgs = []
    k_tns = []
    for ar in aspect_ratios:
        _tg, _tn = peterson_4_24(ar)
        k_tgs.append(_tg)
        k_tns.append(_tn)
    ax.plot(aspect_ratios, k_tgs, linestyle="-", color="k")
    ax.plot(aspect_ratios, k_tns, linestyle="--", color="k")

    for ar in aspect_ratios:
        for n_cutout, shape in zip(n_cutout_elements, shapes):
            # Create the panel geometry
            geom = PanelGeometry(width, height)
            a_pos = width / 2 - d * ar / 2
            b_pos = a_pos + d * ar
            geom.add_cutout(CircularCutout(a_pos, height / 2, d / 2))
            geom.add_cutout(CircularCutout(b_pos, height / 2, d / 2))

            # Discretize the panel
            geom.discretize(num_elements_per_side=n_side, num_elements_cutout=n_cutout)
            elem_B = geom.find_element_by_angle("B", 180)

            # Solve the system
            solver = BEMSolver(BEMKernels(mat), geom)
            solver.assemble()
            u, t = solver.solve(qx=q_x, qy=q_y)

            # Get the stress data
            stress_data = solver.cutout_stress_table(u, t)

            # Find the maximum stress in xx direction
            max_stress_idx = np.where(stress_data == elem_B)[0][0]
            local_xx = stress_data[max_stress_idx, 3]
            local_yy = stress_data[max_stress_idx, 4]
            local_xy = stress_data[max_stress_idx, 5]

            local_max_principal = (local_xx + local_yy) / 2 + np.sqrt(
                (local_xx - local_yy) ** 2 / 4 + local_xy**2
            )
            local_von_mises = np.sqrt(
                local_xx**2 - local_xx * local_yy + local_yy**2 + 3 * local_xy**2
            )

            # Calculate the stress concentration
            q_sigma = q_x / mat.thickness
            scf = local_max_principal / q_sigma

            ax.plot(
                ar,
                scf,
                marker=shape,
                linestyle="None",
                color="r",
            )

            # Store the results
            results.append(
                {
                    "aspect_ratio": ar,
                    "n_cutout": n_cutout,
                    "f_i": q_sigma,
                    "f_xx": local_xx,
                    "f_1": local_max_principal,
                    "f_vm": local_von_mises,
                    "scf_panl": scf,
                    "K_tg": peterson_4_24(ar)[0],
                    "K_tn": peterson_4_24(ar)[1],
                }
            )

            fig_res = plot_results(
                solver,
                u,
                t,
                stress_type="xx",
                # save_path=output_dir / f"4_24_{ar}_{n_cutout}.png",
            )
            ax_1 = fig_res.axes[0]

            _ = ax_1.annotate(
                f"{scf:.2f}x",
                xy=(stress_data[max_stress_idx, 1], stress_data[max_stress_idx, 2]),
                xytext=(
                    0.25 * width,
                    0.75 * height,
                ),
                arrowprops=dict(facecolor="black", width=1, headwidth=5),
                fontsize=9,
            )
            fig_res.savefig(output_dir / f"4_24_{ar}_{n_cutout}_annotated.png")
            plt.close(fig_res)

    # Labels and Title
    ax.set_xlabel("Aspect Ratio (l/d)")
    ax.set_ylabel("Stress Concentration Factor (K_t)")
    ax.set_title("Stress Concentration Factor vs. Aspect Ratio")
    ax.legend(["K_tg", "K_tn", "Panl"])
    ax.grid(True)

    # Print the results
    console = Console()
    table = Table(title="Peterson's 4.24 Stress Concentrations")
    table.add_column("Aspect Ratio", style="cyan", justify="right")
    table.add_column("N_Cutout", style="green", justify="right")
    table.add_column("f_i", style="green", justify="right")
    table.add_column("f_xx", style="green", justify="right")
    table.add_column("f_1", style="green", justify="right")
    table.add_column("f_vm", style="green", justify="right")
    table.add_column("SCF_Panl", style="green", justify="right")
    table.add_column("K_tg", style="green", justify="right")
    table.add_column("K_tn", style="green", justify="right")

    for r in results:
        table.add_row(
            f"{r['aspect_ratio']:.2f}",
            f"{r['n_cutout']:.0f}",
            f"{r['f_i']:.0f}",
            f"{r['f_xx']:.0f}",
            f"{r['f_1']:.0f}",
            f"{r['f_vm']:.0f}",
            f"{r['scf_panl']:.2f}",
            f"{r['K_tg']:.2f}",
            f"{r['K_tn']:.2f}",
        )

    console.print(table)


if __name__ == "__main__":
    run()
