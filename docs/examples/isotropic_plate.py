import numpy as np

from panelyze.analysis.geometry import CircularCutout, PanelGeometry
from panelyze.analysis.kernels import BEMKernels
from panelyze.analysis.material import OrthotropicMaterial
from panelyze.analysis.solver import BEMSolver


def verify_isotropic_plate_running_loads():
    # 1. Setup Material (Imperial: psi)
    # Al 6061-T6 approx: E = 10 Msi, nu = 0.33
    E, nu = 10.0e6, 0.33
    G = E / (2 * (1 + nu))
    # Using a professional thickness for aerospace panels
    thickness = 0.25
    mat = OrthotropicMaterial(e1=E, e2=E * 1.001, nu12=nu, g12=G, thickness=thickness)

    # 2. Setup Geometry (Imperial: inches)
    W, H = 10.0, 10.0
    radius = 0.5
    geom = PanelGeometry(W, H)
    geom.add_cutout(CircularCutout(W / 2, H / 2, radius))

    # 3. Discretize
    n_side = 20
    elements = geom.discretize(num_elements_per_side=n_side, num_elements_cutout=80)

    # 4. Assemble
    solver = BEMSolver(BEMKernels(mat), geom)
    solver.assemble()

    # 5. Boundary Conditions (Running Loads lbf/in)
    bc_type = np.zeros(2 * len(elements), dtype=int)
    bc_value = np.zeros(2 * len(elements))

    # Tension q = 250 lbf/in (corresponds to sigma = 1000 psi for t=0.25)
    q_applied = 250.0
    for i, el in enumerate(elements):
        if np.isclose(el.center[0], 0.0):  # Left
            bc_value[2 * i] = -q_applied
        if np.isclose(el.center[0], W):  # Right
            bc_value[2 * i] = q_applied

    # Corner Constraints
    bc_type[0:2] = 1  # u=0, v=0 at (0,0)
    bc_value[0:2] = 0.0

    k_br = n_side - 1  # (W,0)
    bc_type[2 * k_br + 1] = 1  # v=0
    bc_value[2 * k_br + 1] = 0.0

    # Solve using DEFAULT (running loads)
    u, t = solver.solve(bc_type, bc_value, is_stress=False)

    # 6. Evaluate stress and resultants at hole pole (r=0.51 in)
    eval_pts = np.array([[W / 2, H / 2 + 0.51]])
    stresses = solver.compute_stress(eval_pts, u, t)
    resultants = solver.compute_resultants(eval_pts, u, t)

    sigma_applied = q_applied / thickness
    scf_stress = stresses[0, 0] / sigma_applied
    scf_resultant = resultants[0, 0] / q_applied

    print(f"Applied Sigma: {sigma_applied:.1f} psi")
    print(f"Applied Load q: {q_applied:.1f} lbf/in")
    print(f"Stress sigma_xx: {stresses[0, 0]:.2f} psi")
    print(f"Resultant Nx: {resultants[0, 0]:.2f} lbf/in")
    print(f"SCF (from stress): {scf_stress:.3f}")
    print(f"SCF (from resultant): {scf_resultant:.3f}")

    return scf_stress


if __name__ == "__main__":
    verify_isotropic_plate_running_loads()
