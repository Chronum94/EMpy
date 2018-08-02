"""Single ring resonator example."""

import EMpy
import numpy as np
import matplotlib.pyplot as plt

wls = np.linspace(1.53e-6, 1.56e-6, 1000)
K1 = EMpy.devices.Coupler(wls, np.sqrt(0.08), 1.)
K2 = EMpy.devices.Coupler(wls, np.sqrt(0.08), 1.)
l1 = np.pi * 5e-6
l2 = np.pi * 5e-6
SWG = EMpy.devices.SWG(488, 220, 25).solve(wls)
SRR = EMpy.devices.SRR(K1, K2, SWG.neff, l1, l2).solve()

plt.plot(wls, np.absolute(SRR.THRU), "r.-", wls, np.absolute(SRR.DROP), "g.-")
plt.axis("tight")
plt.ylim([0, 1])
plt.xlabel("wavelength /m")
plt.ylabel("power")
plt.legend(("THRU", "DROP"))
plt.show()
