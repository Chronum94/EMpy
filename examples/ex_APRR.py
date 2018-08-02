"""All-pass ring resonator example."""

import EMpy
import numpy as np
import matplotlib.pyplot as plt

wls = np.linspace(1.5e-6, 1.6e-6, 1000)
K = EMpy.devices.Coupler(wls, np.sqrt(0.08), 1.)
l = 2 * np.pi * 5e-6
SWG = EMpy.devices.SWG(400, 220, 125).solve(wls)
APRR = EMpy.devices.APRR(K, SWG.neff, l).solve()

plt.plot(wls, np.unwrap(np.angle(APRR.THRU)), "r.-")
plt.axis("tight")
plt.xlabel("wavelength /m")
plt.ylabel("phase")
plt.show()
