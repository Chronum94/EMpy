"""N-Ring resonators example."""

import EMpy
import numpy as np
import matplotlib.pyplot as plt

wls = np.linspace(1.53e-6, 1.57e-6, 1000)

Ks = [
    EMpy.devices.Coupler(wls, np.sqrt(0.08), 1.),
    EMpy.devices.Coupler(wls, np.sqrt(0.008), 1.),
    EMpy.devices.Coupler(wls, np.sqrt(0.006), 1.),
    EMpy.devices.Coupler(wls, np.sqrt(0.09), 1.),
]

R = 5e-6
l1s = [np.pi * R, np.pi * R, np.pi * R]
l2s = [np.pi * R, np.pi * R, np.pi * R]

SWG = EMpy.devices.SWG(400, 220, 125).solve(wls)
neffs = [SWG.neff, SWG.neff, SWG.neff]

NRR = EMpy.devices.NRR(Ks, neffs, l1s, l2s).solve()

plt.plot(
    wls,
    20 * np.log10(np.absolute(NRR.THRU)),
    "r.-",
    wls,
    20 * np.log10(np.absolute(NRR.DROP)),
    "g.-",
)
plt.axis("tight")
plt.ylim([-30, 0])
plt.xlabel("wavelength /m")
plt.ylabel("power /dB")
plt.legend(("THRU", "DROP"))
plt.show()
