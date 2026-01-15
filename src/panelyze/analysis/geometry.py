from typing import List

import numpy as np


class BoundaryElement:
    """
    Represents a straight-line boundary element with constant
    traction/displacement assumption.
    """

    def __init__(self, p1: np.ndarray, p2: np.ndarray, tag: str = "outer"):
        self.p1 = p1  # Start point (x, y)
        self.p2 = p2  # End point (x, y)
        self.tag = tag  # 'outer' or 'cutout'
        self.center = (p1 + p2) / 2.0
        self.dx = p2[0] - p1[0]
        self.dy = p2[1] - p1[1]
        self.length = np.sqrt(self.dx**2 + self.dy**2)

        # Outward unit normal (assuming CCW orientation for outer boundary)
        # Normal to (dx, dy) is (dy, -dx) or (-dy, dx).
        # For CCW outer boundary, normal is (dy/L, -dx/L).
        self.nx = self.dy / self.length
        self.ny = -self.dx / self.length


class PanelGeometry:
    """
    Defines the geometry of the panel and any cutouts. handles discretization.
    """

    def __init__(self, width: float, height: float):
        self.width = width
        self.height = height
        self.cutouts: List["Cutout"] = []

    def add_cutout(self, cutout: "Cutout"):
        self.cutouts.append(cutout)

    def discretize(
        self, num_elements_per_side: int, num_elements_cutout: int = 20
    ) -> List[BoundaryElement]:
        """
        Discretizes the panel outer boundary and all cutouts.

        Args:
            num_elements_per_side: Number of elements for each side of the panel.
            num_elements_cutout: Number of elements for each cutout.

        Returns:
            List[BoundaryElement]: List of all boundary elements.
        """
        elements = []

        # Outer boundary (Rectangular)
        # 1. Bottom side: (0,0) to (W, 0)
        elements.extend(
            self._discretize_line(
                np.array([0, 0]), np.array([self.width, 0]), num_elements_per_side
            )
        )
        # 2. Right side: (W,0) to (W, H)
        elements.extend(
            self._discretize_line(
                np.array([self.width, 0]),
                np.array([self.width, self.height]),
                num_elements_per_side,
            )
        )
        # 3. Top side: (W, H) to (0, H)
        elements.extend(
            self._discretize_line(
                np.array([self.width, self.height]),
                np.array([0, self.height]),
                num_elements_per_side,
            )
        )
        # 4. Left side: (0, H) to (0, 0)
        elements.extend(
            self._discretize_line(
                np.array([0, self.height]),
                np.array([0, 0]),
                num_elements_per_side,
                tag="outer",
            )
        )

        # Cutouts
        for cutout in self.cutouts:
            elements.extend(cutout.discretize(num_elements_cutout))

        self.elements = elements
        return elements

    def _discretize_line(
        self, p1: np.ndarray, p2: np.ndarray, num_els: int, tag: str = "outer"
    ) -> List[BoundaryElement]:
        """
        Discretizes a straight line into elements.

        Args:
            p1: Start point.
            p2: End point.
            num_els: Number of elements.
            tag: Tag for the elements.

        Returns:
            List[BoundaryElement]: List of line elements.
        """
        els = []
        pts = np.linspace(p1, p2, num_els + 1)
        for i in range(num_els):
            els.append(BoundaryElement(pts[i], pts[i + 1], tag=tag))
        return els


class Cutout:
    """Base class for cutouts"""

    def discretize(self, num_elements: int) -> List[BoundaryElement]:
        raise NotImplementedError


class CircularCutout(Cutout):
    """Circular cutout definition."""

    def __init__(self, x_center: float, y_center: float, radius: float):
        self.xc = x_center
        self.yc = y_center
        self.r = radius

    def discretize(self, num_elements: int) -> List[BoundaryElement]:
        """
        Discretizes the circular cutout.

        Args:
            num_elements: Number of elements.

        Returns:
            List[BoundaryElement]: List of elements.
        """
        # Discretize CW for an internal boundary
        angles = np.linspace(0, -2 * np.pi, num_elements + 1)
        pts = []
        for a in angles:
            pts.append(
                np.array([self.xc + self.r * np.cos(a), self.yc + self.r * np.sin(a)])
            )

        els = []
        for i in range(num_elements):
            els.append(BoundaryElement(pts[i], pts[i + 1], tag="cutout"))
        return els


class EllipticalCutout(Cutout):
    """Elliptical cutout definition."""

    def __init__(
        self,
        x_center: float,
        y_center: float,
        a: float,
        b: float,
        theta_deg: float = 0.0,
    ):
        self.xc = x_center
        self.yc = y_center
        self.a = a
        self.b = b
        self.theta = np.radians(theta_deg)

    def discretize(self, num_elements: int) -> List[BoundaryElement]:
        """
        Discretizes the elliptical cutout.

        Args:
            num_elements: Number of elements.

        Returns:
            List[BoundaryElement]: List of elements.
        """
        # Discretize CW
        angles = np.linspace(0, -2 * np.pi, num_elements + 1)
        c, s = np.cos(self.theta), np.sin(self.theta)

        pts = []
        for alpha in angles:
            # Point in ellipse coords
            xi = self.a * np.cos(alpha)
            eta = self.b * np.sin(alpha)
            # Rotate and translate
            x = self.xc + xi * c - eta * s
            y = self.yc + xi * s + eta * c
            pts.append(np.array([x, y]))

        els = []
        for i in range(num_elements):
            els.append(BoundaryElement(pts[i], pts[i + 1], tag="cutout"))
        return els
