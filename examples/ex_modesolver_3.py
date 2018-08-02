"""Fully vectorial finite-difference mode solver example."""

import numpy as np
import EMpy
import matplotlib.pyplot as plt

plt.rcParams["image.cmap"] = "coolwarm"


def epsfunc(x_, y_):
    """Similar to ex_modesolver.py, but using anisotropic eps."""
    eps = np.zeros((len(x_), len(y_), 5))
    for ix, xx in enumerate(x_):
        for iy, yy in enumerate(y_):
            if abs(xx - 1.24e-6) <= .24e-6 and abs(yy - 1.11e-6) <= .11e-6:
                a = 3.4757 ** 2
                b = 1  # some xy value
                # eps_xx, xy, yx, yy, zz
                eps[ix, iy, :] = [a, b, b, a, a]
            else:
                a = 1.446 ** 2
                # isotropic
                eps[ix, iy, :] = [a, 0, 0, a, a]
    return eps


wl = 1.55e-6
x = np.linspace(0, 2.48e-6, 125)
y = np.linspace(0, 2.22e-6, 112)

neigs = 2
tol = 1e-8
boundary = "0000"

solver = EMpy.modesolvers.FD.VFDModeSolver(wl, x, y, epsfunc, boundary).solve(
    neigs, tol
)

fig = plt.figure()
fig.add_subplot(1, 3, 1)
plt.contourf(abs(solver.modes[0].Ex), 50)
plt.title("Ex")
fig.add_subplot(1, 3, 2)
plt.contourf(abs(solver.modes[0].Ey), 50)
plt.title("Ey")
fig.add_subplot(1, 3, 3)
plt.contourf(abs(solver.modes[0].Ez), 50)
plt.title("Ez")
plt.show()
