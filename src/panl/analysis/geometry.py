from typing import List

import numpy as np


class BoundaryElement:
    """
    Represents a straight-line boundary element with constant
    traction/displacement assumption.
    """

    def __init__(
        self, p1: np.ndarray, p2: np.ndarray, tag: str = "outer", label: str = None
    ):
        self.p1 = p1  # Start point (x, y)
        self.p2 = p2  # End point (x, y)
        self.tag = tag  # 'outer' or 'cutout'
        self.label = label
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
        if cutout.label is None:
            count = len(self.cutouts)
            cutout.label = chr(ord("A") + count)
        self.cutouts.append(cutout)

    def get_element_ids_by_label(self, label: str) -> List[int]:
        """
        Returns the indices of elements that match the given label.

        Args:
            label: The label to search for (e.g., 'A', 'B').

        Returns:
            List[int]: List of element indices.

        Raises:
            ValueError: If the geometry has not been discretized yet.
        """
        if not hasattr(self, "elements"):
            raise ValueError(
                "Geometry must be discretized before retrieving element IDs."
            )
        return [i for i, el in enumerate(self.elements) if el.label == label]

    def find_element_by_angle(self, label: str, angle_deg: float) -> int:
        """
        Finds the ID of the boundary element at a specific angle relative
        to the cutout center.

        Args:
            label: The label of the cutout to search within.
            angle_deg: Angle in degrees, measured CCW from positive X-axis.

        Returns:
            int: The index of the boundary element closest to the specified angle.

        Raises:
            ValueError: If the geometry is not discretized or the label is not found.
        """
        if not hasattr(self, "elements"):
            raise ValueError(
                "Geometry must be discretized before searching for elements."
            )

        # Find the cutout with matching label
        cutout = next((c for c in self.cutouts if c.label == label), None)
        if cutout is None:
            raise ValueError(f"No cutout found with label '{label}'.")

        # Get elements for this cutout
        indices = self.get_element_ids_by_label(label)
        if not indices:
            raise ValueError(f"No elements found for cutout label '{label}'.")

        target_angle_rad = np.radians(angle_deg % 360)

        best_idx = -1
        min_angle_diff = float("inf")

        center = cutout.center
        for idx in indices:
            el = self.elements[idx]
            # Vector from center to element center
            vec = el.center - center
            angle = np.arctan2(vec[1], vec[0]) % (2 * np.pi)

            # Shortest angular distance
            diff = abs(angle - target_angle_rad)
            diff = min(diff, 2 * np.pi - diff)

            if diff < min_angle_diff:
                min_angle_diff = diff
                best_idx = idx

        return best_idx

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

    def __init__(self, label: str = None):
        self.label = label

    def discretize(self, num_elements: int) -> List[BoundaryElement]:
        raise NotImplementedError

    @property
    def center(self) -> np.ndarray:
        """
        Returns the center of the cutout as (x, y).

        Raises:
            NotImplementedError: If not implemented in subclass.
        """
        raise NotImplementedError


class CircularCutout(Cutout):
    """Circular cutout definition."""

    def __init__(
        self, x_center: float, y_center: float, radius: float, label: str = None
    ):
        super().__init__(label=label)
        self.xc = x_center
        self.yc = y_center
        self.r = radius

    @property
    def center(self) -> np.ndarray:
        return np.array([self.xc, self.yc])

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
            els.append(
                BoundaryElement(pts[i], pts[i + 1], tag="cutout", label=self.label)
            )
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
        label: str = None,
    ):
        super().__init__(label=label)
        self.xc = x_center
        self.yc = y_center
        self.a = a
        self.b = b
        self.theta = np.radians(theta_deg)

    @property
    def center(self) -> np.ndarray:
        return np.array([self.xc, self.yc])

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
            els.append(
                BoundaryElement(pts[i], pts[i + 1], tag="cutout", label=self.label)
            )
        return els
