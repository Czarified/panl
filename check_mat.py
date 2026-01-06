import numpy as np

from panelyze.analysis.material import OrthotropicMaterial


def check_material():
    E1, E2, nu12, G12 = 10000.0, 5000.0, 0.3, 3000.0
    mat = OrthotropicMaterial(e1=E1, e2=E2, nu12=nu12, g12=G12)

    print(f"Compliance Matrix (beta):\n{mat.beta}")
    print(f"Stiffness Matrix (C):\n{mat.C}")

    # Check Poisson relation
    nu21 = mat.nu21
    print(f"nu21: {nu21}")

    # Expected beta[0, 1] = -nu12 / E1
    print(f"Expected beta[0, 1]: {-nu12/E1}")

    # Expected beta[1, 1] = 1/E2
    print(f"Expected beta[1, 1]: {1/E2}")

    # Check strain-stress
    # sigma = [100, 0, 0]
    sigma = np.array([100.0, 0.0, 0.0])
    # epsilon = beta * sigma
    epsilon = mat.beta @ sigma
    print(f"Strains for sigma_xx=100:\n{epsilon}")

    # Check back-stress
    sigma_recon = mat.C @ epsilon
    print(f"Reconstructed Stress:\n{sigma_recon}")


if __name__ == "__main__":
    check_material()
