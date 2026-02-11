import pytest

from panl.analysis.geometry import CircularCutout, EllipticalCutout, PanelGeometry


def test_cutout_auto_labeling():
    geom = PanelGeometry(10, 10)
    c1 = CircularCutout(2, 2, 1)
    c2 = EllipticalCutout(5, 5, 2, 1)

    geom.add_cutout(c1)
    geom.add_cutout(c2)

    assert c1.label == "A"
    assert c2.label == "B"


def test_cutout_manual_labeling():
    geom = PanelGeometry(10, 10)
    c1 = CircularCutout(2, 2, 1, label="Custom")
    c2 = EllipticalCutout(5, 5, 2, 1)

    geom.add_cutout(c1)
    geom.add_cutout(c2)

    assert c1.label == "Custom"
    # Even if we manual label the first,
    # the next one should be "B" because len(cutouts) was 1
    assert c2.label == "B"


def test_get_element_ids_by_label():
    geom = PanelGeometry(10, 10)
    c1 = CircularCutout(2, 2, 1)
    c2 = EllipticalCutout(5, 5, 2, 1)

    geom.add_cutout(c1)
    geom.add_cutout(c2)

    # 4 elements per side = 16 (actually 4 per side = 4*4 = 16)
    # 5 elements per cutout
    num_per_side = 2
    num_per_cutout = 5

    geom.discretize(
        num_elements_per_side=num_per_side, num_elements_cutout=num_per_cutout
    )

    ids_a = geom.get_element_ids_by_label("A")
    ids_b = geom.get_element_ids_by_label("B")

    # Outer elements: 4 * 2 = 8
    # Cutout A: 5 elements
    # Cutout B: 5 elements
    # Total: 18 elements

    assert len(ids_a) == num_per_cutout
    assert len(ids_b) == num_per_cutout

    # Verify they are the correct indices
    # Outer is 0-7, A is 8-12, B is 13-17
    assert ids_a == [8, 9, 10, 11, 12]
    assert ids_b == [13, 14, 15, 16, 17]


def test_get_element_ids_before_discretize():
    geom = PanelGeometry(10, 10)
    geom.add_cutout(CircularCutout(2, 2, 1))

    with pytest.raises(ValueError, match="Geometry must be discretized"):
        geom.get_element_ids_by_label("A")


def test_find_element_by_angle():
    geom = PanelGeometry(10, 10)
    # Put a circular cutout at (5, 5) with radius 1
    # Note: CircularCutout discretizes CW starting from 0 (angle=0)
    # pts: np.array([self.xc + self.r * np.cos(a), self.yc + self.r * np.sin(a)])
    # angles = np.linspace(0, -2 * np.pi, num_elements + 1)
    # For 4 elements, angles are 0, -pi/2, -pi, -3pi/2, -2pi
    # centers are:
    # 0: angle = -pi/4 (315 deg)
    # 1: angle = -3pi/4 (225 deg)
    # 2: angle = -5pi/4 (135 deg)
    # 3: angle = -7pi/4 (45 deg)
    geom.add_cutout(CircularCutout(5, 5, 1, label="A"))
    geom.discretize(num_elements_per_side=1, num_elements_cutout=4)

    # Top dead center is 90 degrees.
    # Closest should be element 2 (135 deg) or 3 (45 deg).
    # Wait, 135-90 = 45. 90-45 = 45.
    # If we ask for 90, it might pick either. Let's ask for 135 precisely.
    idx = geom.find_element_by_angle("A", 135)
    # Outer is 4 elements (0-3). Cutout A is index 4-7.
    # Element 2 of cutout A is index 4+2 = 6.
    assert idx == 6

    # Bottom dead center is 270 degrees.
    # Element 1 is 225 deg, element 0 is 315 deg.
    # 270-225 = 45. 315-270 = 45.
    idx = geom.find_element_by_angle("A", 225)
    assert idx == 5

    # Test 0 degrees (315 and 45 are neighbors)
    idx = geom.find_element_by_angle("A", 0)
    assert idx == 4 or idx == 7
