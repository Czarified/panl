"""Results Validation Against Peterson's Textbook."""

from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
from rich.console import Console
from rich.table import Table

from panl.analysis import plot_results
from panl.analysis.geometry import CircularCutout, PanelGeometry
from panl.analysis.kernels import BEMKernels
from panl.analysis.material import OrthotropicMaterial
from panl.analysis.solver import BEMSolver


def peterson_4_1(aspect_ratio: int) -> Tuple[float, float]:
    """Returns the Peterson's stress concentrations from chart 4.1.

    Chart 4.1 is a table of stress concentrations for a single hole in a
    plate of infinite extent, but finite height. The K_tg value is defined
    as the ratio of the maximum stress to the applied stress, based on the
    gross area of the cross section. The K_tn value is defined as the ratio
    of the maximum stress to the applied stress, based on the net area of
    the cross section.

    Pilkey advises that if the stress gradient is of concern, as in certain
    fatigue problems, the proper factor to use is K_tn.

    Args:
        aspect_ratio (int): Hole diameter over panel height.

    Returns:
        Tuple[float, float]: The Peterson's K_tg and K_tn values.
    """
    K_tg = (
        0.284
        + 2 / (1 - aspect_ratio)
        - 0.6 * (1 - aspect_ratio)
        + 1.32 * (1 - aspect_ratio) ** 2
    )
    K_tn = (
        2
        + 0.284 * (1 - aspect_ratio)
        - 0.6 * (1 - aspect_ratio) ** 2
        + 1.32 * (1 - aspect_ratio) ** 3
    )
    return K_tg, K_tn


def peterson_4_3(aspect_ratio: float, eccentricity: float) -> Tuple[float, float]:
    """Returns the Peterson's stress concentrations from chart 4.3.

    Chart 4.3 is a table of stress concentrations for a single eccentrically
    located hole in a plate of infinite extent, but finite height. The K_tg
    value is defined as the ratio of the maximum stress to the applied stress,
    based on the gross area of the cross section. The K_tn value is defined
    as the ratio of the maximum stress to the applied stress, based on the
    net area of the cross section.

    Pilkey advises that if the stress gradient is of concern, as in certain
    fatigue problems, the proper factor to use is K_tn.

    Args:
        aspect_ratio (float): Hole radius over hole offset from bottom edge.
        eccentricity (float): Hole offset from bottom edge over offset from top edge.

    Returns:
        Tuple[float, float]: The Peterson's K_tg and K_tn values.
    """
    # K_tg constants
    C_1 = 2.9969 - 0.0090 * (eccentricity) + 0.01338 * (eccentricity) ** 2
    C_2 = 0.1217 + 0.5180 * (eccentricity) - 0.5297 * (eccentricity) ** 2
    C_3 = 0.5565 + 0.7215 * (eccentricity) + 0.6153 * (eccentricity) ** 2
    C_4 = 4.082 + 6.0146 * (eccentricity) - 3.9815 * (eccentricity) ** 2
    K_tg = (
        C_1
        + C_2 * (1 / aspect_ratio)
        + C_3 * (1 / aspect_ratio) ** 2
        + C_4 * (1 / aspect_ratio) ** 3
    )

    # K_tn constants
    D_1 = 2.989 - 0.0064 * eccentricity
    D_2 = -2.872 + 0.095 * eccentricity
    D_3 = 2.348 + 0.196 * eccentricity
    K_tn = D_1 + D_2 * aspect_ratio + D_3 * aspect_ratio**2
    return K_tg, K_tn


# # #   G L O B A L S   # # #

# Setup a basic material for the parametric study.
# Note that the thickness is defined here,
# but not used in the calculation of the stress concentrations.
E, nu = 10.5e6, 0.33
G = E / (2 * (1 + nu))
thickness = 0.080
mat = OrthotropicMaterial(e1=E, e2=E * 1.001, nu12=nu, g12=G, thickness=thickness)


# # #   E X A M P L E   4 . 1   # # #

# Example 4.1, Basic Panel Dimensions
# To approximate an infinite panel, the width will be set to 3 time the height.
height = 10.0
width = 3 * height

# Constant and Sweep parameters
# Use the same Uniaxial Tension load for each case.
q_x = 500
n_side = 40
aspect_ratios = [0.1, 0.2, 0.3, 0.4, 0.5]
n_cutout_elements = [44, 88, 160, 220]

# Storage for results
results = []

# Loop through each aspect ratio
for n_cutout in n_cutout_elements:
    for ar in aspect_ratios:
        # Calculate the cutout diameter
        d = ar * height

        # Create the panel geometry
        geom = PanelGeometry(width, height)
        geom.add_cutout(CircularCutout(width / 2, height / 2, d / 2))

        # Discretize the panel
        elements = geom.discretize(
            num_elements_per_side=n_side, num_elements_cutout=n_cutout
        )

        # Solve the system
        solver = BEMSolver(BEMKernels(mat), geom)
        solver.assemble()
        u, t = solver.solve(qx=q_x)

        # Get the stress data
        stress_data = solver.cutout_stress_table(u, t)

        # Find the maximum stress in xx direction
        max_stress_idx = np.argmax(stress_data[:, 3])
        local_xx = stress_data[max_stress_idx, 3]
        local_yy = stress_data[max_stress_idx, 4]
        local_xy = stress_data[max_stress_idx, 5]

        local_max_principal = np.sqrt((local_xx - local_yy) ** 2 / 4 + local_xy**2)
        local_von_mises = np.sqrt(
            local_xx**2 - local_xx * local_yy + local_yy**2 + 3 * local_xy**2
        )

        # Calculate the stress concentration
        q_sigma = q_x / thickness
        scf = local_xx / q_sigma

        # Store the results
        results.append(
            {
                "aspect_ratio": ar,
                "n_cutout": n_cutout,
                "diameter": d,
                "f_i": q_sigma,
                "f_xx": local_xx,
                "f_1": local_max_principal,
                "f_vm": local_von_mises,
                "scf_panl": scf,
                "K_tg": peterson_4_1(ar)[0],
                "K_tn": peterson_4_1(ar)[1],
            }
        )

        # Save the plot
        fig = plot_results(
            solver,
            u,
            t,
            stress_type="xx",
            save_path=f"scripts/fig/peterson_validation_{ar}_{n_cutout}.png",
        )
        plt.close(fig)

# Print the results
console = Console()

# Create a table
table = Table(title="Peterson's 4.1 Stress Concentrations")
table.add_column("Aspect Ratio", style="cyan", justify="right")
table.add_column("N_Cutout", style="green", justify="right")
table.add_column("Diameter", style="green", justify="right")
table.add_column("f_i", style="green", justify="right")
table.add_column("f_xx", style="green", justify="right")
table.add_column("f_1", style="green", justify="right")
table.add_column("f_vm", style="green", justify="right")
table.add_column("SCF_Panl", style="green", justify="right")
table.add_column("K_tg", style="green", justify="right")
table.add_column("K_tn", style="green", justify="right")

# Add the results to the table
for r in results:
    table.add_row(
        f"{r['aspect_ratio']:.2f}",
        f"{r['n_cutout']:.0f}",
        f"{r['diameter']:.2f}",
        f"{r['f_i']:.2f}",
        f"{r['f_xx']:.2f}",
        f"{r['f_1']:.2f}",
        f"{r['f_vm']:.2f}",
        f"{r['scf_panl']:.2f}",
        f"{r['K_tg']:.2f}",
        f"{r['K_tn']:.2f}",
    )

console.print(table)

fig, ax = plt.subplots(figsize=(10, 6))

# Peterson's SCFs
k_tgs = []
k_tns = []
for ar in aspect_ratios:
    _tg, _tn = peterson_4_1(ar)
    k_tgs.append(_tg)
    k_tns.append(_tn)
ax.plot(aspect_ratios, k_tgs, linestyle="-", color="k")  # Peterson's K_tg
ax.plot(aspect_ratios, k_tns, linestyle="--", color="k")  # Peterson's K_tn

# Panl's SCF
ax.plot(
    [r["aspect_ratio"] for r in results],
    [r["scf_panl"] for r in results],
    marker="o",
    linestyle="None",
    color="r",
)

# Labels and Title
ax.set_xlabel("Aspect Ratio (d/h)")
ax.set_ylabel("Stress Concentration Factor (K_t)")
ax.set_title("Stress Concentration Factor vs. Aspect Ratio")
ax.legend(["K_tg", "K_tn", "Panl"])
ax.grid(True)


# # #   E X A M P L E   4 . 3   # # #

# Example 4.3, Basic Panel Dimensions
# To approximate an infinite panel, the width will be set to 3 time the height.
height = 10.0
width = 3 * height

# Constant and Sweep parameters
# Use the same Uniaxial Tension load for each case.
q_x = 500
n_side = 40
aspect_ratios = [0.1, 0.2, 0.3, 0.4, 0.5]
# Define the sweep for eccentricity as a percentaget of the height.
# We'll calculate the ratio in the loop.
eccentricities = [0.33, 0.25]
n_cutout_elements = [44, 220]

# Storage for results
results = []

for ar in aspect_ratios:
    for ec in eccentricities:
        c = ec * height
        e = height - c
        e_ratio = c / e
        for n_cutout in n_cutout_elements:
            # Calculate the cutout radius
            a = ar * c

            # Create the panel geometry
            geom = PanelGeometry(width, height)
            geom.add_cutout(CircularCutout(width / 2, c, a))

            # Discretize the panel
            elements = geom.discretize(
                num_elements_per_side=n_side, num_elements_cutout=n_cutout
            )

            # Solve the system
            solver = BEMSolver(BEMKernels(mat), geom)
            solver.assemble()
            u, t = solver.solve(qx=q_x)

            # Get the stress data
            stress_data = solver.cutout_stress_table(u, t)

            # Find the maximum stress in xx direction
            max_stress_idx = np.argmax(stress_data[:, 3])
            local_xx = stress_data[max_stress_idx, 3]
            local_yy = stress_data[max_stress_idx, 4]
            local_xy = stress_data[max_stress_idx, 5]

            local_max_principal = np.sqrt((local_xx - local_yy) ** 2 / 4 + local_xy**2)
            local_von_mises = np.sqrt(
                local_xx**2 - local_xx * local_yy + local_yy**2 + 3 * local_xy**2
            )

            # Calculate the stress concentration
            q_sigma = q_x / thickness
            scf = local_xx / q_sigma

            # Store the results
            results.append(
                {
                    "aspect_ratio": ar,
                    "n_cutout": n_cutout,
                    "diameter": d,
                    "eccentricity": ec,
                    "f_i": q_sigma,
                    "f_xx": local_xx,
                    "f_1": local_max_principal,
                    "f_vm": local_von_mises,
                    "scf_panl": scf,
                    "K_tg": peterson_4_3(ar, e_ratio)[0],
                    "K_tn": peterson_4_3(ar, e_ratio)[1],
                }
            )

            # Save the plot
            # fig = plot_results(
            #     solver,
            #     u,
            #     t,
            #     stress_type="xx",
            #     save_path=f"scripts/fig/peterson_validation_{ar}_{n_cutout}.png",
            # )
            # plt.close(fig)

# Print the results
console = Console()

# Create a table
table = Table(title="Peterson's 4.3 Stress Concentrations")
table.add_column("Aspect Ratio", style="cyan", justify="right")
table.add_column("N_Cutout", style="green", justify="right")
table.add_column("Diameter", style="green", justify="right")
table.add_column("Eccentricity", style="green", justify="right")
table.add_column("f_i", style="green", justify="right")
table.add_column("f_xx", style="green", justify="right")
table.add_column("f_1", style="green", justify="right")
table.add_column("f_vm", style="green", justify="right")
table.add_column("SCF_Panl", style="green", justify="right")
table.add_column("K_tg", style="green", justify="right")
table.add_column("K_tn", style="green", justify="right")

# Add the results to the table
for r in results:
    table.add_row(
        f"{r['aspect_ratio']:.2f}",
        f"{r['n_cutout']:.0f}",
        f"{r['diameter']:.2f}",
        f"{r['eccentricity']:.2f}",
        f"{r['f_i']:.2f}",
        f"{r['f_xx']:.2f}",
        f"{r['f_1']:.2f}",
        f"{r['f_vm']:.2f}",
        f"{r['scf_panl']:.2f}",
        f"{r['K_tg']:.2f}",
        f"{r['K_tn']:.2f}",
    )

console.print(table)

# Plot the results
fig, ax = plt.subplots(figsize=(10, 6))

# Peterson's SCFs
k_tgs = []
k_tns = []
ec = 0.25
for ar in aspect_ratios:
    _tg, _tn = peterson_4_3(ar, ec)
    k_tgs.append(_tg)
    k_tns.append(_tn)
ax.plot(aspect_ratios, k_tgs, linestyle="-", color="k")  # Peterson's K_tg
ax.plot(aspect_ratios, k_tns, linestyle="--", color="k")  # Peterson's K_tn

# Panl's SCF
ax.plot(
    [r["aspect_ratio"] for r in results if r["eccentricity"] == ec],
    [r["scf_panl"] for r in results if r["eccentricity"] == ec],
    marker="o",
    linestyle="None",
    color="r",
)

# Labels and Title
ax.set_xlabel("Aspect Ratio (d/h)")
ax.set_ylabel("Stress Concentration Factor (K_t)")
ax.set_title("Stress Concentration Factor vs. Aspect Ratio (Eccentricity = 0.25)")
ax.legend(["K_tg", "K_tn", "Panl"])
ax.grid(True)


plt.show()
