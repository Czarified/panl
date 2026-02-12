"""Common logic for Peterson's validation studies."""

from pathlib import Path
from typing import Tuple

import numpy as np

from panl.analysis.material import OrthotropicMaterial

THICKNESS = 0.080

RAW_DATA_424 = np.array(
    [
        (1.2005958291956307, 3.9845059162136227),
        (1.2090367428003972, 3.7743347415548776),
        (1.2428003972194637, 3.5302602124114246),
        (1.2596822244289971, 3.3641551512295282),
        (1.3018867924528303, 3.208214532173093),
        (1.3440913604766633, 3.1065112012522507),
        (1.3862959285004965, 3.014977361856832),
        (1.445382323733863, 2.9098808341608735),
        (1.496027805362463, 2.8183453116321338),
        (1.5466732869910629, 2.7505386026627168),
        (1.5973187686196624, 2.689511554710248),
        (1.6564051638530286, 2.628482823624459),
        (1.7070506454816285, 2.5844049282143633),
        (1.7576961271102285, 2.54541177856698),
        (1.7999006951340613, 2.5131999730698666),
        (1.8336643495531282, 2.482684765960311),
        (1.9011916583912611, 2.448772995808998),
        (1.9602780536246278, 2.4148629087910045),
        (2.0362462760675273, 2.377559624997896),
        (2.1122144985104274, 2.3487309174759727),
        (2.183962264150943, 2.3199030515207104),
        (2.2810327706057594, 2.2910701361654855),
        (2.3907646474677255, 2.252065204584855),
        (2.5469215491559076, 2.218135761533671),
        (2.694637537239325, 2.184208001615808),
        (2.9056603773584904, 2.1587421944692236),
        (3.099801390268123, 2.136669584097756),
        (3.344587884806355, 2.1162817901806004),
        (3.555610724925522, 2.1026803898136768),
        (3.83416087388282, 2.0890655243801857),
        (4.138033763654419, 2.072055779038257),
        (4.585402184707051, 2.0584072509383464),
        (5.074975173783516, 2.0447503071718307),
        (5.497020854021847, 2.0344966589803577),
        (6.1385302879841115, 2.024199249322539),
        (6.982621648460775, 2.007081783448067),
        (7.978649453823238, 2.0001035126992406),
    ]
)


def get_material() -> OrthotropicMaterial:
    """Setup a basic material for the parametric study.

    Args:
        None

    Returns:
        OrthotropicMaterial: The material for the parametric study.
    """
    E, nu = 10.5e6, 0.33
    G = E / (2 * (1 + nu))
    thickness = THICKNESS
    return OrthotropicMaterial(e1=E, e2=E * 1.001, nu12=nu, g12=G, thickness=thickness)


def get_output_dir(example_name: str) -> Path:
    """Create and return the output directory for an example.

    Args:
        example_name (str): The name of the example.

    Returns:
        Path: The output directory for the example.
    """
    output_dir = Path(f"scripts/validation_study/fig/{example_name}")
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def peterson_4_1(aspect_ratio: float) -> Tuple[float, float]:
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
        aspect_ratio (float): Hole diameter over panel height.

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
    K_tg = C_1 + C_2 * aspect_ratio + C_3 * aspect_ratio**2 + C_4 * aspect_ratio**3

    # K_tn constants
    D_1 = 2.989 - 0.0064 * eccentricity
    D_2 = -2.872 + 0.095 * eccentricity
    D_3 = 2.348 + 0.196 * eccentricity
    K_tn = D_1 + D_2 * aspect_ratio + D_3 * aspect_ratio**2
    return K_tg, K_tn


def peterson_4_24(aspect_ratio: float) -> Tuple[float, float]:
    """Returns the Peterson's stress concentrations from chart 4.24.

    Chart 4.24 is a table of stress concentrations for two holes in a
    plate of infinite extents. The K_tg value is defined
    as the ratio of the maximum stress to the applied stress, based on the
    gross area of the cross section. The K_tn value is defined as the ratio
    of the maximum stress to the applied stress, based on the net area of
    the cross section.

    Pilkey advises that if the stress gradient is of concern, as in certain
    fatigue problems, the proper factor to use is K_tn.

    Args:
        aspect_ratio (float): Hole spacing over diameter.

    Returns:
        Tuple[float, float]: The Peterson's K_tg and K_tn values.
    """
    K_tn = (
        2.000
        - 2.119 * aspect_ratio**-1
        + 2.493 * aspect_ratio**-2
        - 1.372 * aspect_ratio**-3
    )
    K_tg = K_tn * np.sqrt(1 - aspect_ratio**-2) / (1 - aspect_ratio**-1)
    return K_tg, K_tn
