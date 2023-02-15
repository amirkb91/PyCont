import numpy as np
from scipy.integrate import odeint
import scipy.linalg as spl


class Duffing:
    # mass and linear stiffness matrices
    M = np.eye(2)
    K = np.array([[2, -1], [-1, 2]])
    mu = 0.5

    @classmethod
    def statespace(cls, Z, t):
        # ODE of the mass spring system. Z\dot(t) = g(Z(t))
        X = Z[:2]
        Xdot = Z[2:]
        Minv = spl.inv(cls.M)
        KX = cls.K @ X
        fnl = np.array([cls.mu * X[0] ** 3, 0])
        Zdot = np.concatenate((Xdot, -Minv @ (KX + fnl)))
        return Zdot

    @classmethod
    def sensitivity(cls, dZdZ0, t, x1):
        # dZdZ0\dot(t) = dg(t)dZ * dZdZ0
        # problem has been flattenned so reshape to recover, return flatten
        d_dZdZ0_dt = np.array(
            [[0, 0, 1, 0], [0, 0, 0, 1], [-2 - 1.5 * x1 ** 2, 1, 0, 0], [1, -2, 0, 0]]
        ) @ dZdZ0.reshape(4, 4)
        return d_dZdZ0_dt.flatten()

    @classmethod
    def time_solve(cls, T, Z0, cont_params):
        nperiod = cont_params["shooting"]["single"]["nperiod"]
        nsteps = cont_params["shooting"]["single"]["nsteps_per_period"]

        # run ODE
        t = np.linspace(0, T*nperiod, nsteps*nperiod)
        Z = odeint(cls.statespace, Z0, t)

        # periodicity condition
        H = Z[-1, :] - Z[0, :]
        H = H.reshape(-1, 1)

        # Energy
        M = np.eye(2)
        fnl = np.array([[0.5 * Z[-1, 0] ** 3], [0]])
        energy = 0.5 * (Z[-1, 2:].T @ M @ Z[-1, 2:] + Z[-1, :2].T @ fnl)
        outputs = {"energy": np.array([energy])}

        # Sensitvity analysis **CHECK WITH FINITE DIFFERENCE**
        # IC = eye. numpy ODE only works on 1D arrays so flatten problem then reshape
        dZdZ0 = odeint(cls.sensitivity, np.eye(4).flatten(), t, args=(Z[-1, 0],))
        dZdZ0 = dZdZ0[-1, :].reshape(4, 4)
        Mm0 = dZdZ0 - np.eye(4)
        dHdt = cls.statespace(Z[-1, :], ())

        cvg = True

        return H, Mm0, dHdt, outputs, cvg

    @classmethod
    def eigen_solve(cls, cont_params):
        # Continuation variables initial guess from eigenvalues

        frq, eig = spl.eigh(cls.K, cls.M)
        frq = np.sqrt(frq) / (2 * np.pi)

        # frq, eig = np.linalg.eig(np.linalg.inv(cls.M) @ cls.K)
        # idx = np.argsort(frq)
        # frq = frq[idx]
        # eig = eig[:, idx]
        # frq = np.sqrt(frq) / (2 * np.pi)

        nnm = cont_params["first_point"]["eig_start"]["NNM"]
        scale = cont_params["first_point"]["eig_start"]["scale"]
        X0 = eig[:, nnm - 1] * scale
        V0 = np.zeros_like(X0)
        Z0 = np.concatenate([X0, V0])
        T0 = 1 / frq[nnm - 1]
        # initial position, fixed at 0 since displacement is stored in variable X0
        pose0 = np.zeros_like(X0)

        return Z0, T0, pose0
