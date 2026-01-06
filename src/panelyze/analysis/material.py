from typing import Tuple

import numpy as np


class OrthotropicMaterial:
    """
    Represents an orthotropic material and handles the calculation of compliance
    matrices and characteristic equation roots for BEM.
    """

    def __init__(
        self,
        e1: float,
        e2: float,
        nu12: float,
        g12: float,
        theta_deg: float = 0.0,
        thickness: float = 1.0,
    ):
        """
        Initialize orthotropic material properties.

        Args:
            e1: Young's modulus in the fiber direction (1)
            e2: Young's modulus in the transverse direction (2)
            nu12: Poisson's ratio
            g12: Shear modulus
            theta_deg: Angle of orientation in degrees relative to x-axis
            thickness: Panel thickness
        """
        self.e1 = e1
        self.e2 = e2
        self.nu12 = nu12
        self.g12 = g12
        self.theta = np.radians(theta_deg)
        self.thickness = thickness

        # nu21 = nu12 * E2 / E1
        self.nu21 = nu12 * e2 / e1

        # Calculate compliance in principal directions (1, 2)
        # S11 S12 0
        # S12 S22 0
        # 0   0   S66
        self.s11 = 1.0 / e1
        self.s22 = 1.0 / e2
        self.s12 = -nu12 / e1
        self.s66 = 1.0 / g12

        # Compliance matrix in principal directions
        self.S_principal = np.array(
            [[self.s11, self.s12, 0.0], [self.s12, self.s22, 0.0], [0.0, 0.0, self.s66]]
        )

        # Compliance matrix in global (x, y) directions
        self.beta = self._transform_compliance(self.theta)

        # Roots of the characteristic equation
        self.mu1, self.mu2 = self._solve_characteristic_roots()

        # Stiffness matrix [C] = [S]^-1
        self.C = np.linalg.inv(self.beta)

    def _transform_compliance(self, theta: float) -> np.ndarray:
        """
        Transform the compliance matrix from principal directions (1, 2)
        to global coordinates (x, y) using rotation angle theta.

        Args:
            theta: Rotation angle in radians.

        Returns:
            np.ndarray: Transformed compliance matrix beta.
        """
        c = np.cos(theta)
        s = np.sin(theta)

        # beta coefficients for plane stress (anisotropic)
        # Formulae for transformed compliance coefficients beta_ij
        b11 = (
            self.s11 * c**4 + (2 * self.s12 + self.s66) * s**2 * c**2 + self.s22 * s**4
        )
        b22 = (
            self.s11 * s**4 + (2 * self.s12 + self.s66) * s**2 * c**2 + self.s22 * c**4
        )
        b12 = (self.s11 + self.s22 - self.s66) * s**2 * c**2 + self.s12 * (s**4 + c**4)
        b66 = (
            4 * (self.s11 + self.s22 - 2 * self.s12) * s**2 * c**2
            + self.s66 * (c**2 - s**2) ** 2
        )
        b16 = (2 * self.s11 - 2 * self.s12 - self.s66) * s * c**3 - (
            2 * self.s22 - 2 * self.s12 - self.s66
        ) * s**3 * c
        b26 = (2 * self.s11 - 2 * self.s12 - self.s66) * s**3 * c - (
            2 * self.s22 - 2 * self.s12 - self.s66
        ) * s * c**3

        return np.array([[b11, b12, b16], [b12, b22, b26], [b16, b26, b66]])

    def _solve_characteristic_roots(self) -> Tuple[complex, complex]:
        """
        Solve the characteristic equation:
        b11*mu^4 - 2*b16*mu^3 + (2*b12 + b66)*mu^2 - 2*b26*mu + b22 = 0

        Returns:
            Tuple[complex, complex]: Two roots with positive imaginary parts.

        Raises:
            ValueError: If the equation does not yield valid complex roots.
        """
        b11 = self.beta[0, 0]
        b12 = self.beta[0, 1]
        b16 = self.beta[0, 2]
        b22 = self.beta[1, 1]
        b26 = self.beta[1, 2]
        b66 = self.beta[2, 2]

        coeffs = [b11, -2 * b16, (2 * b12 + b66), -2 * b26, b22]
        roots = np.roots(coeffs)

        # We need two roots mu_1, mu_2 with positive imaginary parts.
        # Typically they appear as pairs mu1, bar(mu1) and mu2, bar(mu2).
        pos_im_roots = [r for r in roots if np.imag(r) > 1e-12]

        if len(pos_im_roots) != 2:
            # Handle cases with multiple roots or purely real roots?
            # For orthotropic materials in physical ranges, they should be complex.
            # If they are double roots, we might need to handle it.
            # Sort by imaginary part or just take the first two?
            pos_im_roots = sorted(pos_im_roots, key=lambda r: np.imag(r), reverse=True)[
                :2
            ]

        if len(pos_im_roots) < 2:
            raise ValueError(
                "Characteristic equation did not yield two complex roots "
                "with positive imaginary parts. Check material properties."
            )

        return pos_im_roots[0], pos_im_roots[1]

    def get_stiffness_matrix(self) -> np.ndarray:
        """
        Returns the constitutive matrix [C] = [S]^-1.

        Returns:
            np.ndarray: Stiffness matrix.
        """
        return np.linalg.inv(self.beta)
