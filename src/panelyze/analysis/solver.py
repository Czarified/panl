from typing import List, Tuple

import numpy as np

from .geometry import BoundaryElement
from .kernels import BEMKernels


class BEMSolver:
    """
    Boundary Element Method solver for anisotropic elasticity.
    """

    def __init__(self, kernels: BEMKernels, elements: List[BoundaryElement]):
        self.kernels = kernels
        self.elements = elements
        self.num_elements = len(elements)
        self.M = self.num_elements

        # System matrices
        self.H = np.zeros((2 * self.M, 2 * self.M))
        self.G = np.zeros((2 * self.M, 2 * self.M))

        # Solve for coefficients p, q, A, mu from kernels
        self.mu1 = kernels.mu1
        self.mu2 = kernels.mu2
        self.p1, self.p2 = kernels.p1, kernels.p2
        self.q1, self.q2 = kernels.q1, kernels.q2
        self.A = kernels.A  # A[i, k] where i load dir, k root index

    def assemble(self):
        """Assembles G and H matrices using constant elements."""
        for i in range(self.M):  # Source element (collocation point at center)
            Pi = self.elements[i].center

            for j in range(self.M):  # Field element
                el_j = self.elements[j]

                # Compute integrals of U and T over el_j with source at Pi
                du, dt = self._integrate_kernels(Pi, el_j, i == j)

                # Fill matrices
                # row indices: 2*i, 2*i+1 (x and y directions)
                # col indices: 2*j, 2*j+1
                # Fill matrices (Transposed due to reciprocity)
                # Reciprocity: U_ji(P, Q) = U_ij(Q, P).
                # My integrate(P, Q) gives response at Q due to load at P.
                # So we transpose to get response at P due to condition at Q.
                self.G[2 * i : 2 * i + 2, 2 * j : 2 * j + 2] = du.T
                self.H[2 * i : 2 * i + 2, 2 * j : 2 * j + 2] = dt.T

        # Rigid-body sum trick for diagonal of H
        # This replaces the 0.5 jump term and ensures H * 1 = 0
        for i in range(self.M):
            h_diag = np.zeros((2, 2))
            for k in range(self.M):
                if i != k:
                    h_diag += self.H[2 * i : 2 * i + 2, 2 * k : 2 * k + 2]
            self.H[2 * i : 2 * i + 2, 2 * i : 2 * i + 2] = -h_diag

    def _integrate_kernels(
        self, source_pt: np.ndarray, el: BoundaryElement, is_singular: bool
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Integrates U_ji and T_ji kernels over a straight line element.

        Args:
            source_pt: Point where load is applied.
            el: Boundary element to integrate over.
            is_singular: True if source_pt is the center of el.

        Returns:
            Tuple[np.ndarray, np.ndarray]: Integrated G and H matrices.
        """
        # Element properties
        p1 = el.p1 - source_pt
        p2 = el.p2 - source_pt

        nx, ny = el.nx, el.ny
        # Normal angle alpha
        alpha = np.arctan2(ny, nx)

        # Complex coords z_k at endpoints
        z1_1 = p1[0] + self.mu1 * p1[1]
        z1_2 = p2[0] + self.mu1 * p2[1]
        z2_1 = p1[0] + self.mu2 * p1[1]
        z2_2 = p2[0] + self.mu2 * p2[1]

        # Denominator term (mu_k * cos(alpha) - sin(alpha))
        # This matches the NASA doc where alpha is normal angle.
        den1 = self.mu1 * np.cos(alpha) - np.sin(alpha)
        den2 = self.mu2 * np.cos(alpha) - np.sin(alpha)

        # Robust integration using log-difference identity
        # integral(log z) = z(log z - 1)
        # diff = z2(log z2 - 1) - z1(log z1 - 1) = (z2-z1)(log z1 - 1) + z2 log(z2/z1)
        # where log(z2/z1) = log|z2/z1| + i*angle_change.

        # Continuous change in log
        dln1 = np.log(z1_2 / z1_1)
        dln2 = np.log(z2_2 / z2_1)

        # Integral of 1/z (Traction kernel derivative)
        di21 = dln1 / den1
        di22 = dln2 / den2

        # Integral of log z (Displacement kernel)
        di11 = ((z1_2 - z1_1) * (np.log(z1_1) - 1.0) + z1_2 * dln1) / den1
        di12 = ((z2_2 - z2_1) * (np.log(z2_1) - 1.0) + z2_2 * dln2) / den2

        if is_singular:
            # Singular T is 0 in CPV.
            di21 = 0j
            di22 = 0j
            # Singular U is fine with the general formula if we are careful,
            # but using z1_2 = -z1_1 works.
            # My general formula above handles it correctly.

        # Eq 55: DU_ji = 2 Re { P_j1 A_i1 Di11 + P_j2 A_i2 Di12 }
        # Note: j is disp dir, i is load dir. My kernel code used U[j, i].
        du = np.zeros((2, 2))
        dt = np.zeros((2, 2))

        # Displacement kernel assembly
        P = np.array([[self.p1, self.p2], [self.q1, self.q2]])
        # Traction kernel coefficients Q_jk * (mu_k * nx - ny)
        Q = np.array([[self.mu1, self.mu2], [-1.0, -1.0]], dtype=complex)

        for i in range(2):  # load dir (x,y)
            for j in range(2):  # component dir (x,y)
                # Displacement
                val_u = P[j, 0] * self.A[i, 0] * di11 + P[j, 1] * self.A[i, 1] * di12
                du[j, i] = 2.0 * val_u.real

                # Traction
                # Eq 56: DT_ji = 2 Re { Q_j1 (mu1 nx - ny) A_i1 Di21 + ... }
                # Wait, NASA 56 says Q_i1 (mu1 n1 - n2) A_j1?
                # Let's re-check the indices. j component, i load.
                term1 = Q[j, 0] * (self.mu1 * nx - ny) * self.A[i, 0] * di21
                term2 = Q[j, 1] * (self.mu2 * nx - ny) * self.A[i, 1] * di22
                dt[j, i] = 2.0 * (term1 + term2).real

        return du, dt

    def solve(
        self, bc_type: np.ndarray, bc_value: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Solves the BEM system H*u = G*t for unknown boundary values.

        Args:
            bc_type: Array of 0 (traction known) or 1 (displacement known).
            bc_value: Value of the prescribed boundary condition.

        Returns:
            Tuple[np.ndarray, np.ndarray]: Complete u and t boundary arrays.
        """
        size = 2 * self.M
        A = np.zeros((size, size))
        b = np.zeros(size)

        # H u = G t
        # If u is known: move H*u to RHS. Matrix column becomes -G.
        # If t is known: move G*t to RHS. Matrix column remains H.

        for j in range(size):  # column index
            if (
                bc_type[j] == 1
            ):  # Displacement unknown? No, bc_type defines what is GIVEN.
                pass  # logic below is better

        # Let x be the vector of unknowns (u if t given, t if u given).
        # H u - G t = 0
        # For each DOF k:
        # If u[k] given: t[k] unknown. col k of A = -G[:, k]. RHS -= H[:, k] * u[k]
        # If t[k] given: u[k] unknown. col k of A = H[:, k]. RHS += G[:, k] * t[k]

        for k in range(size):
            if bc_type[k] == 1:  # Displacement u[k] is GIVEN
                A[:, k] = -self.G[:, k]
                b -= self.H[:, k] * bc_value[k]
            else:  # Traction t[k] is GIVEN
                A[:, k] = self.H[:, k]
                b += self.G[:, k] * bc_value[k]

        x = np.linalg.solve(A, b)

        # Recover u and t
        u = np.zeros(size)
        t = np.zeros(size)
        for k in range(size):
            if bc_type[k] == 1:
                u[k] = bc_value[k]
                t[k] = x[k]
            else:
                u[k] = x[k]
                t[k] = bc_value[k]

        return u, t

    def compute_displacement(
        self, points: np.ndarray, u_boundary: np.ndarray, t_boundary: np.ndarray
    ) -> np.ndarray:
        """
        Computes displacements at interior points using Somigliana Identity.

        Args:
            points: (N, 2) array of interior points.
            u_boundary: Solved boundary displacements.
            t_boundary: Solved boundary tractions.

        Returns:
            np.ndarray: Computed displacements at each point.
        """
        N = points.shape[0]
        results = np.zeros((N, 2))

        for i in range(N):
            P = points[i]
            u_pt = np.zeros(2)

            for j in range(self.M):
                el = self.elements[j]
                # Integrate kernels over element j with source at P
                du, dt = self._integrate_kernels(P, el, False)

                # t_j = [t_x, t_y], u_j = [u_x, u_y]
                tj = t_boundary[2 * j : 2 * j + 2]
                uj = u_boundary[2 * j : 2 * j + 2]

                # Reciprocity: Transpose kernels
                # du[res_Q, load_P] -> du_T[res_P, load_Q]
                u_pt += -du.T @ tj + dt.T @ uj

            results[i] = u_pt

        return results

    def compute_stress(
        self, points: np.ndarray, u_boundary: np.ndarray, t_boundary: np.ndarray
    ) -> np.ndarray:
        """
        Calculates stress components (sigma_xx, sigma_yy, tau_xy) at interior points.

        Args:
            points: Array of points (N, 2) where stress is computed.
            u_boundary: Solved boundary displacements.
            t_boundary: Solved boundary tractions.

        Returns:
            np.ndarray: Array of stresses (N, 3).
        """
        N = points.shape[0]
        stresses = np.zeros((N, 3))  # sigma_xx, sigma_yy, tau_xy

        for i in range(N):
            P = points[i]
            # Grad u_i,k = sum G_ijk * t_k - sum H_ijk * u_k
            # where G_ijk = integral of U_ij,k
            grad_u = self._compute_u_gradient(P, u_boundary, t_boundary)

            # Strain e_ij = 0.5 * (u_i,j + u_j,i)
            exx = grad_u[0, 0]
            eyy = grad_u[1, 1]
            gxy = grad_u[0, 1] + grad_u[1, 0]

            # Stress from Hooke's Law (Plane Stress)
            # sigma = [E] * epsilon
            # We use the stiffness matrix (inverse of beta)
            C = self.kernels.mat.C
            stresses[i, 0] = C[0, 0] * exx + C[0, 1] * eyy + C[0, 2] * gxy
            stresses[i, 1] = C[1, 0] * exx + C[1, 1] * eyy + C[1, 2] * gxy
            stresses[i, 2] = C[2, 0] * exx + C[2, 1] * eyy + C[2, 2] * gxy

        return stresses

    def _compute_u_gradient(
        self, P: np.ndarray, u_boundary: np.ndarray, t_boundary: np.ndarray
    ) -> np.ndarray:
        """
        Computes the displacement gradient matrix du_i/dx_j at point P.

        Args:
            P: The interior point (x, y) where the gradient is computed.
            u_boundary: Solved boundary displacements.
            t_boundary: Solved boundary tractions.

        Returns:
            np.ndarray: The 2x2 displacement gradient matrix.
        """
        grad = np.zeros((2, 2))  # grad[i, j] = du_i/dx_j

        for m in range(self.M):
            el = self.elements[m]
            nx, ny = el.nx, el.ny
            alpha = np.arctan2(ny, nx)

            p1 = el.p1 - P
            p2 = el.p2 - P

            z1_1 = p1[0] + self.mu1 * p1[1]
            z1_2 = p2[0] + self.mu1 * p2[1]
            z2_1 = p1[0] + self.mu2 * p1[1]
            z2_2 = p2[0] + self.mu2 * p2[1]

            den1 = self.mu1 * np.cos(alpha) - np.sin(alpha)
            den2 = self.mu2 * np.cos(alpha) - np.sin(alpha)

            # Integrals for gradients:
            dln1 = np.log(z1_2 / z1_1)
            dln2 = np.log(z2_2 / z2_1)

            di21 = dln1 / den1
            di22 = dln2 / den2

            di31 = (1.0 / z1_2 - 1.0 / z1_1) / den1
            di32 = (1.0 / z2_2 - 1.0 / z2_1) / den2

            tm = t_boundary[2 * m : 2 * m + 2]
            um = u_boundary[2 * m : 2 * m + 2]

            P_mat = np.array([[self.p1, self.p2], [self.q1, self.q2]])
            Q_mat = np.array([[self.mu1, self.mu2], [-1.0, -1.0]], dtype=complex)

            for j in range(2):  # displacement component
                for k in range(2):  # derivative wrt x_k
                    grad_z1 = 1.0 if k == 0 else self.mu1
                    grad_z2 = 1.0 if k == 0 else self.mu2

                    for i in range(2):  # load dir
                        # dG_ji,k = integral( d/dP_k U ) = - integral( d/dQ_k U )
                        # Transpose for reciprocity: load-i at Q, response-j at P
                        val_g = (
                            P_mat[i, 0] * self.A[j, 0] * di21 * grad_z1
                            + P_mat[i, 1] * self.A[j, 1] * di22 * grad_z2
                        )
                        grad[j, k] += (-2.0 * val_g.real) * tm[i]

                        # dH_ji,k = integral( d/dP_k T ) = - integral( d/dQ_k T )
                        # val_h = integral( d/dQ_k T )
                        # Transpose for reciprocity
                        val_h = (
                            Q_mat[i, 0]
                            * (self.mu1 * nx - ny)
                            * self.A[j, 0]
                            * di31
                            * grad_z1
                            + Q_mat[i, 1]
                            * (self.mu2 * nx - ny)
                            * self.A[j, 1]
                            * di32
                            * grad_z2
                        )
                        grad[j, k] += (2.0 * val_h.real) * um[i]

        return grad
