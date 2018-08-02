import numpy as np
import EMpy
import matplotlib.pyplot as plt

# define the multilayer
epsilon = [
    1.0 ** 2 * EMpy.constants.eps0 * np.eye(3),
    EMpy.constants.eps0 * np.diag([2.1, 2.0, 1.9]),
    2.3 ** 2 * EMpy.constants.eps0 * np.eye(3),
    4.3 ** 2 * EMpy.constants.eps0 * np.eye(3),
    3.0 ** 2 * EMpy.constants.eps0 * np.eye(3),
]

d = np.array([np.inf, 1e-6, 2.3e-6, 0.1e-6, np.inf])

aniso_layers = EMpy.utils.Multilayer()
for i in xrange(len(epsilon)):
    eps = EMpy.materials.EpsilonTensor(epsilon[i] * np.eye(3))
    mat = EMpy.materials.AnisotropicMaterial("layer_%d" % i, eps)
    layer = EMpy.utils.Layer(mat, d[i])
    aniso_layers.append(layer)

# define the planewave
theta_inc_x = EMpy.utils.deg2rad(0.)
theta_inc_y = 0.
wls = np.linspace(1.4e-6, 1.7e-6, 100)

# solve
tm = EMpy.transfer_matrix.AnisotropicTransferMatrix(
    aniso_layers, theta_inc_x, theta_inc_y
)
solution_aniso = tm.solve(wls)

# plot
plt.figure()
plt.plot(
    wls,
    solution_aniso.R[0, 0, :],
    wls,
    solution_aniso.R[1, 0, :],
    wls,
    solution_aniso.R[0, 1, :],
    wls,
    solution_aniso.R[1, 1, :],
    wls,
    solution_aniso.T[0, 0, :],
    wls,
    solution_aniso.T[1, 0, :],
    wls,
    solution_aniso.T[0, 1, :],
    wls,
    solution_aniso.T[1, 1, :],
)
plt.legend(("Rss", "Rps", "Rsp", "Rpp", "Tss", "Tps", "Tsp", "Tpp"))
plt.title("Anisotropic Multilayer")
plt.xlabel("wavelength /m")
plt.ylabel("Power /dB")
plt.xlim(wls.min(), wls.max())
plt.show()
