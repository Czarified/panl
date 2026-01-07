Validation: NASA 5.3.1
=======================

The Panelyze BEM implementation is validated against Example 5.3.1 from **NASA CR-125596** (also referenced in subsequent BEM documentation).

Problem Description
-------------------
A circular hole of radius :math:`R=5.0` is located at the center of a large orthotropic plate. The plate is subjected to a remote uniaxial stress :math:`\sigma_{\infty} = 100` in the x-direction.

Material Properties (Orthotropic)
---------------------------------
The material properties used for verification are:

- :math:`E_1 = 10,000`
- :math:`E_2 = 5,000`
- :math:`\nu_{12} = 0.3`
- :math:`G_{12} = 3,000`

Theoretical Solution
--------------------
For an infinite plate, the theoretical stress concentration factor (SCF) at the pole of the hole (:math:`\theta = 90^\circ` or :math:`270^\circ`) is given by:

.. math::
   K = 1 + \sqrt{2 \left( \sqrt{\frac{E_1}{E_2}} - \nu_{12} \right) + \frac{E_1}{G_{12}}}

For the given properties:

.. math::
   K = 1 + \sqrt{2 \left( \sqrt{\frac{10000}{5000}} - 0.3 \right) + \frac{10000}{3000}} \approx 3.3583

BEM Results
-----------
The validation script `nasa_531.py` discretizes the hole into 400 linear elements and uses a 200x200 panel (width = :math:`40R`) to approximate an infinite plate.

At a distance :math:`r = 1.01R` from the center, the calculated stresses are:

+------------------+-----------------------------+-------------------+-------+
| Angle (:math:`\theta`) | Theoretical :math:`\sigma_{xx}` | BEM :math:`\sigma_{xx}` | Error |
+==================+=============================+===================+=======+
| :math:`90^\circ` | 335.83                      | 329.14            | 2.0%  |
+------------------+-----------------------------+-------------------+-------+
| :math:`45^\circ` | -                           | 84.14             | -     |
+------------------+-----------------------------+-------------------+-------+
| :math:`0^\circ`  | 0.0                         | -6.70             | -     |
+------------------+-----------------------------+-------------------+-------+

Discussion
----------
The 2% error is attributable to two primary factors:

1. **Domain Size**: The BEM models a finite panel (:math:`W=40R`), whereas the theoretical solution is for an infinite sheet.
2. **Evaluation Point**: The stress is evaluated at :math:`r=1.01R` to ensure numerical stability when using interior point kernels.

Despite these approximations, the BEM solver correctly captures the high-gradient stress concentration and the directionality of the orthotropic response.
