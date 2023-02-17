import numpy as np
from scipy.integrate import odeint
import scipy.linalg as spl


class Duffing:
    # parameters of nonlinear system EoM    MX'' + KX + fnl = 0
    M = np.eye(2)
    K = np.array([[2, -1], [-1, 2]])
    Knl = 0.5
    Minv = spl.inv(M)

    # finite element data, 2 dof system
    free_dof = np.array([0, 1])
    ndof_all = 2
    ndof_fix = 0
    ndof_free = 2

    @classmethod
    def statespace(cls, Z, t):
        # ODE of the mass spring system. Z\dot(t) = g(Z(t))
        X = Z[:2]
        Xdot = Z[2:]
        KX = cls.K @ X
        fnl = np.array([cls.Knl * X[0] ** 3, 0])
        Zdot = np.concatenate((Xdot, -cls.Minv @ (KX + fnl)))
        return Zdot

    @classmethod
    def monodromy(cls, dZdZ0, t, ZT1):
        M = cls.M
        K = cls.K
        knl = cls.Knl

        dgdz = np.array([[0, 0, 1, 0],
                         [0, 0, 0, 1],
                         [-1 / M[0, 0] * K[0, 0] + knl * 3 * ZT1 ** 2, -1 / M[0, 0] * K[0, 1], 0, 0],
                         [-1 / M[1, 1] * K[1, 0], -1 / M[1, 1] * K[1, 1], 0, 0]])

        dZdZ0dot = dgdz @ dZdZ0.reshape(4, 4)
        return dZdZ0dot.flatten()

    @classmethod
    def sensitivity(cls, dZdZ0, t, x1):
        # dZdZ0\dot(t) = dg(t)dZ * dZdZ0
        # problem has been flattenned so reshape to recover, return flatten
        d_dZdZ0_dt = np.array(
            [[0, 0, 1, 0], [0, 0, 0, 1], [-2 - 1.5 * x1 ** 2, 1, 0, 0], [1, -2, 0, 0]]
        ) @ dZdZ0.reshape(4, 4)
        return d_dZdZ0_dt.flatten()

    @classmethod
    def time_solve(cls, T, Z0, pose_base, cont_params):
        nperiod = cont_params["shooting"]["single"]["nperiod"]
        nsteps = cont_params["shooting"]["single"]["nsteps_per_period"]
        rel_tol = cont_params["shooting"]["rel_tol"]

        # time integration, store position and velocity
        t = np.linspace(0, T*nperiod, nsteps*nperiod)
        Z = np.array(odeint(cls.statespace, Z0, t, rtol=rel_tol))
        pose_time = Z[:, :2].T
        vel_time = Z[:, 2:].T

        # periodicity condition
        H = Z[-1, :] - Z[0, :]
        H = H.reshape(-1, 1)

        # Energy, conservative system so take initial Z values
        X = Z[0, :2]
        Xdot = Z[0, 2:]
        fnl = np.array([cls.Knl * X[0] ** 3, 0])
        energy = 0.5 * (Xdot.T @ cls.M @ Xdot + X.T @ cls.K @ X + X @ fnl)

        # Sensitivity and Monodromy
        dHdt = cls.statespace(Z[-1, :], None)
        dZdZ0 = np.array(odeint(cls.monodromy, np.eye(4).flatten(), t, args=(Z[-1, 0],)))
        dZdZ0 = dZdZ0[-1, :].reshape(4, 4)

        # eps = 1e-6
        # for i in range(4):
        #
        # zp = Z0 +
        # Z = np.array(odeint(cls.statespace, Z0, t, rtol=rel_tol)


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

        nnm = cont_params["first_point"]["eig_start"]["NNM"]
        scale = cont_params["first_point"]["eig_start"]["scale"]
        X0 = eig[:, nnm - 1] * scale
        V0 = np.zeros_like(X0)
        Z0 = np.concatenate([X0, V0])
        T0 = 1 / frq[nnm - 1]

        # all position variables are required for the Lie group code only, can be set to None
        pose0 = None

        return Z0, T0, pose0

    @classmethod
    def get_fe_data(cls):
        return {"free_dof": cls.free_dof, "ndof_all": cls.ndof_all, "ndof_fix": cls.ndof_fix,
                "ndof_free": cls.ndof_free}
