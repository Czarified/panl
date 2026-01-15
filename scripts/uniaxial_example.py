import matplotlib.pyplot as plt
import numpy as np

from panelyze.analysis import plot_results
from panelyze.analysis.geometry import CircularCutout, PanelGeometry
from panelyze.analysis.kernels import BEMKernels
from panelyze.analysis.material import OrthotropicMaterial
from panelyze.analysis.solver import BEMSolver

E, nu = 10.5e6, 0.33
G = E / (2 * (1 + nu))
thickness = 0.080
mat = OrthotropicMaterial(e1=E, e2=E * 1.001, nu12=nu, g12=G, thickness=thickness)

W, H = 30.0, 15.0
geom = PanelGeometry(W, H)
geom.add_cutout(CircularCutout(W / 2, H / 2, 1.5))

n_side = 40
elements = geom.discretize(num_elements_per_side=n_side, num_elements_cutout=88)

solver = BEMSolver(BEMKernels(mat), elements)
solver.assemble()

bc_type = np.zeros(2 * len(elements), dtype=int)
bc_value = np.zeros(2 * len(elements))

q_applied = 500
for i, el in enumerate(elements):
    if np.isclose(el.center[0], 0.0):
        bc_value[2 * i] = -q_applied
    if np.isclose(el.center[0], W):
        bc_value[2 * i] = q_applied

bc_type[0:2] = 1
bc_value[0:2] = 0.0
bc_type[2 * (n_side - 1) + 1] = 1
bc_value[2 * (n_side - 1) + 1] = 0.0

u, t = solver.solve(bc_type, bc_value)

stress_data = solver.print_cutout_stress_table(u, t)
# The array contains [id, x, y, sxx, syy, txy]
# Find the maximum sigma_xx (column index 3)
max_stress_idx = np.argmax(stress_data[:, 3])
local_stress = stress_data[max_stress_idx, 3]

# Does displacement make sense?
idx = n_side + n_side // 2
u_x_end = u[2 * idx]
u_y_end = u[2 * idx + 1]
print(f"Displacement at end: {u_x_end:.3f}, {u_y_end:.4f}")

# Expect SCF ~3.14 based on external sources
q_sigma = q_applied / thickness
scf = local_stress / q_sigma
print(f"Far-Field Stress: {q_sigma:.0f} [psi]")
print(f"Max Local Stress: {local_stress:.0f} [psi]")
print(f"K_t: {scf:.2f}")

# Peterson's Factors
# Use Chart 4.1
d = 3
h = 15
x = d / h
K_tg = 0.284 + 2 / (1 - x) - 0.6 * (1 - x) + 1.32 * (1 - x) ** 2
print(f"Peterson's K_tg: {K_tg:.2f}")

# After solving the system
fig = plot_results(
    solver,
    u,
    t,
    deform_scale=150,
    title="Circular Cutout under X-Tension",
    stress_type="xx",
)
ax = fig.get_axes()[0]
_ = ax.annotate(
    f"Peterson's K_tg: {K_tg:.2f}, f_xx = {K_tg*q_sigma:.0f}[psi]",
    xy=(W / 2, H / 2 + 1.51),
    xytext=(0.6 * W, 0.75 * H),
    arrowprops=dict(facecolor="black", shrink=0.05, width=1, headwidth=5),
    fontsize=9,
)
plt.show()
