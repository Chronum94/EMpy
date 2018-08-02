"""Semi-vectorial finite-difference mode solver example."""

import numpy as np
import EMpy
import matplotlib.pyplot as plt


def epsfunc(x_, y_):
    """Return a matrix describing a 2d material.

    :param x_: x values
    :param y_: y values
    :return: 2d-matrix
    """
    xx, yy = np.meshgrid(x_, y_)
    return np.where(
        (np.abs(xx.T - 1.24e-6) <= .24e-6) * (np.abs(yy.T - 1.11e-6) <= .11e-6),
        3.4757 ** 2,
        1.446 ** 2,
    )


wl = 1.55e-6
x = np.linspace(0, 2.48e-6, 125)
y = np.linspace(0, 2.22e-6, 112)

neigs = 2
tol = 1e-8
boundary = "0000"

solver = EMpy.modesolvers.FD.SVFDModeSolver(wl, x, y, epsfunc, boundary).solve(
    neigs, tol
)

fig = plt.figure()
fig.add_subplot(1, 2, 1)
plt.contourf(abs(solver.Ex[0]), 50)
plt.title("Ex first mode")
fig.add_subplot(1, 2, 2)
plt.contourf(abs(solver.Ex[1]), 50)
plt.title("Ex second mode")
plt.show()
