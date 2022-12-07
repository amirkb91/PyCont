import numpy as np
from scipy.integrate import odeint


def massspring(Z, t):
    # ODE of the mass spring system. Z\dot(t) = g(Z(t))
    dZdt = np.stack((Z[2], Z[3], -2 * Z[0] + Z[1] - 0.5 * Z[0] ** 3, -2 * Z[1] + Z[0]))
    return dZdt


def sensitivity(dZdZ0, t, x1):
    # dZdZ0\dot(t) = dg(t)dZ * dZdZ0
    # problem has been flattenned so reshape to recover, return flatten
    d_dZdZ0_dt = np.array(
        [[0, 0, 1, 0], [0, 0, 0, 1], [-2 - 1.5 * x1**2, 1, 0, 0], [1, -2, 0, 0]]
    ) @ dZdZ0.reshape(4, 4)
    return d_dZdZ0_dt.flatten()


def zerofunction(T, Z, par):
    # unpack run parameters
    ns = par["shooting"]["npts"]

    # run ODE
    t = np.linspace(0, T[0], ns)
    Z_all = odeint(massspring, Z, t)

    # periodicity condition
    H = Z_all[-1, :] - Z_all[0, :]
    H = H.reshape(-1, 1)

    # Sensitvity analysis **CHECK WITH FINITE DIFFERENCE**
    # IC = eye. numpy ODE only works on 1D arrays so flatten problem then reshape
    dZdZ0 = odeint(sensitivity, np.eye(4).flatten(), t, args=(Z_all[-1, 0],))
    dZdZ0 = dZdZ0[-1, :].reshape(4, 4)
    Mm0 = dZdZ0 - np.eye(4)
    dHdt = massspring(Z_all[-1, :], ())

    # Energy
    M = np.eye(2)
    fnl = np.array([[0.5 * Z_all[-1, 0] ** 3], [0]])
    energy = 0.5 * (Z_all[-1, 2:].T @ M @ Z_all[-1, 2:] + Z_all[-1, :2].T @ fnl)
    outputs = {"energy": np.array([energy])}

    cvg = True

    return H, Mm0, dHdt, outputs, cvg


def initialguess_eig():
    # Continuation variables initial guess from eigenvalues
    M = np.eye(2)
    K = np.array([[2, -1], [-1, 2]])
    frq, eig = np.linalg.eig(np.linalg.inv(M) @ K)
    idx = np.argsort(frq)
    frq = frq[idx]
    eig = eig[:, idx]
    frq = np.sqrt(frq) / (2 * np.pi)

    nnm = 2
    scale = 1e-10
    x0 = eig[:, nnm - 1] * scale
    v0 = np.zeros(len(x0))
    Z0 = np.concatenate([x0, v0])
    T0 = np.array([1 / frq[nnm - 1]])

    return Z0, T0
