"""Rigorous Coupled Wave Analysis example."""

import numpy as np
import matplotlib.pyplot as plt

import EMpy
from EMpy.materials import IsotropicMaterial, RefractiveIndex


alpha = EMpy.utils.deg2rad(30.)
delta = EMpy.utils.deg2rad(45.)
psi = EMpy.utils.deg2rad(0.)  # TE
phi = EMpy.utils.deg2rad(90.)

wls = np.linspace(1.5495e-6, 1.550e-6, 101)

LAMBDA = 1.e-6  # grating periodicity
n = 3  # orders of diffraction

Top = IsotropicMaterial("Top", n0=RefractiveIndex(n0_const=1.))
Bottom = IsotropicMaterial("Bottom", n0=RefractiveIndex(n0_const=3.47))

multilayer = EMpy.utils.Multilayer(
    [
        EMpy.utils.Layer(Top, np.inf),
        EMpy.utils.BinaryGrating(Top, Bottom, .4, LAMBDA, .01),
        EMpy.utils.Layer(Bottom, np.inf),
    ]
)

solution = EMpy.RCWA.IsotropicRCWA(multilayer, alpha, delta, psi, phi, n).solve(wls)

plt.plot(
    wls,
    solution.DE1[n, :],
    "ko-",
    wls,
    solution.DE3[n, :],
    "ro-",
    wls,
    solution.DE1[n - 1, :],
    "kx-",
    wls,
    solution.DE3[n - 1, :],
    "rx-",
    wls,
    solution.DE1[n + 1, :],
    "k.-",
    wls,
    solution.DE3[n + 1, :],
    "r.-",
)
plt.xlabel("wavelength /m")
plt.ylabel("diffraction efficiency")
plt.legend(("DE1:0", "DE3:0", "DE1:-1", "DE3:-1", "DE1:+1", "DE3:+1"))
plt.axis("tight")
plt.ylim([0, 1])
plt.show()
