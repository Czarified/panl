import numpy as np

from .material import OrthotropicMaterial


class BEMKernels:
    """
    Implements the fundamental solution kernels for anisotropic elasticity
    as formulated in NASA CR-1934.
    """

    def __init__(self, material: OrthotropicMaterial):
        self.mat = material
        self.mu1 = material.mu1
        self.mu2 = material.mu2

        # Calculate p_k, q_k (Eq 16)
        b = material.beta
        b11, b12, b16 = b[0, 0], b[0, 1], b[0, 2]
        b22, b26 = b[1, 1], b[1, 2]

        self.p1 = b11 * self.mu1**2 + b12 - b16 * self.mu1
        self.p2 = b11 * self.mu2**2 + b12 - b16 * self.mu2
        self.q1 = b12 * self.mu1 + b22 / self.mu1 - b26
        self.q2 = b12 * self.mu2 + b22 / self.mu2 - b26

        # Solve for coefficients A_ik (Eq 25, 26)
        # Unit load in direction i (1=x, 2=y)
        # 2 * Re [ sum A_ik ] = -Y_i/pi
        # 2 * Re [ sum mu_k A_ik ] = X_i/pi
        # 2 * Re [ sum p_k A_ik ] = 0
        # 2 * Re [ sum q_k A_ik ] = 0

        self.A = self._solve_A_coeffs()

    def _solve_A_coeffs(self) -> np.ndarray:
        """
        Solves the 4x4 real system for the complex coefficients A_i1, A_i2.

        Returns:
            np.ndarray: A (2, 2) complex array where A[i, k-1] is A_{i+1, k}.
        """
        # We have 4 complex unknowns: A_11, A_12, A_21, A_22.
        # Actually 2 separate 4x4 real systems (one for i=1, one for i=2).
        # Let's write the system matrix M such that M * {Re(A_k), Im(A_k)} = RHS
        # Row 1: 2 * Re(A_1 + A_2) = 2*Re(A_1) + 2*Re(A_2)
        # Row 2: 2 * Re(mu_1*A_1 + mu_2*A_2) = 2*(mu1_R*Re(A_1) - mu1_I*Im(A_1) + ...)

        m1r, m1i = self.mu1.real, self.mu1.imag
        m2r, m2i = self.mu2.real, self.mu2.imag
        p1r, p1i = self.p1.real, self.p1.imag
        p2r, p2i = self.p2.real, self.p2.imag
        q1r, q1i = self.q1.real, self.q1.imag
        q2r, q2i = self.q2.real, self.q2.imag

        # System matrix M for columns [Re(A1), Im(A1), Re(A2), Im(A2)]
        # Row 1: 2*Im(A1 + A2) = 2*b1 + 2*b2
        # Row 2: 2*Im(mu1*A1 + mu2*A2) = 2*(mu1_I*a1 + mu1_R*b1 + ...)
        # Row 3: 2*Im(p1*A1 + p2*A2) = 2*(p1_I*a1 + p1_R*b1 + ...)
        # Row 4: 2*Im(q1*A1 + q2*A2) = 2*(q1_I*a1 + q1_R*b1 + ...)
        M = np.array(
            [
                [0, 2, 0, 2],
                [2 * m1i, 2 * m1r, 2 * m2i, 2 * m2r],
                [2 * p1i, 2 * p1r, 2 * p2i, 2 * p2r],
                [2 * q1i, 2 * q1r, 2 * q2i, 2 * q2r],
            ],
            dtype=float,
        )

        # RHS for unit x-load (i=1 => X=1, Y=0)
        # Row 2: 2*Im(sum mu_k A_k) = 1/(2*pi)
        rhs_x = np.array([0, 1.0 / (2.0 * np.pi), 0, 0])
        # RHS for unit y-load (i=2 => X=0, Y=1)
        # Row 1: 2*Im(sum A_k) = -1/(2*pi)
        rhs_y = np.array([-1.0 / (2.0 * np.pi), 0, 0, 0])

        sol_x = np.linalg.solve(M, rhs_x)
        sol_y = np.linalg.solve(M, rhs_y)

        A = np.zeros((2, 2), dtype=complex)
        A[0, 0] = sol_x[0] + 1j * sol_x[1]
        A[0, 1] = sol_x[2] + 1j * sol_x[3]
        A[1, 0] = sol_y[0] + 1j * sol_y[1]
        A[1, 1] = sol_y[2] + 1j * sol_y[3]

        return A

    def displacement_kernel(self, dx: float, dy: float) -> np.ndarray:
        """
        Calculates the displacement tensor U_ji(dx, dy).

        Args:
            dx: x-distance between source and field point.
            dy: y-distance between source and field point.

        Returns:
            np.ndarray: A (2, 2) array where U[j, i] is displacement in j due
                to load in i.
        """
        # z_k = dx + mu_k * dy
        z1 = dx + self.mu1 * dy
        z2 = dx + self.mu2 * dy

        # log(z) is complex. NumPy's log handles complex.
        # Note: BEM often requires careful branch cut handling,
        # but for simple panels usually standard log is fine.
        ln_z1 = np.log(z1)
        ln_z2 = np.log(z2)

        U = np.zeros((2, 2))
        # U_ji = 2 * Re { P_j1 A_i1 log(z1) + P_j2 A_i2 log(z2) }
        # where P = [[p1, p2], [q1, q2]]
        P = np.array([[self.p1, self.p2], [self.q1, self.q2]])

        for i in range(2):  # load direction (x=0, y=1)
            for j in range(2):  # displacement direction (x=0, y=1)
                val = P[j, 0] * self.A[i, 0] * ln_z1 + P[j, 1] * self.A[i, 1] * ln_z2
                U[j, i] = 2.0 * val.real

        return U

    def traction_kernel(self, dx: float, dy: float, nx: float, ny: float) -> np.ndarray:
        """
        Calculates the traction tensor T_ji(dx, dy) with normal (nx, ny).

        Args:
            dx: x-distance between source and field point.
            dy: y-distance between source and field point.
            nx: x-component of the normal vector.
            ny: y-component of the normal vector.

        Returns:
            np.ndarray: A (2, 2) array where T[j, i] is traction in j due
                to load in i.
        """
        z1 = dx + self.mu1 * dy
        z2 = dx + self.mu2 * dy

        # Terms for T_ji: (mu_k * ny - nx) * A_ik / z_k

        T = np.zeros((2, 2))
        for i in range(2):  # load direction
            for j in range(2):  # traction direction
                if j == 0:  # t1 = sigma_x * nx + tau_xy * ny
                    q1 = self.mu1**2 * nx - self.mu1 * ny
                    q2 = self.mu2**2 * nx - self.mu2 * ny
                else:  # t2 = tau_xy * nx + sigma_y * ny
                    q1 = -self.mu1 * nx + ny
                    q2 = -self.mu2 * nx + ny

                term1 = q1 * self.A[i, 0] / z1
                term2 = q2 * self.A[i, 1] / z2
                T[j, i] = 2.0 * (term1 + term2).real

        return T
