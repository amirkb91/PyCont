import numpy as np
from scipy.integrate import odeint
import scipy.linalg as spl
import scipy.interpolate as spi


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
    def system_ode(cls, t, Z):
        # ODE of the mass spring system. Zdot(t) = g(Z(t))
        X = Z[:2]
        Xdot = Z[2:]
        KX = cls.K @ X
        fnl = np.array([cls.Knl * X[0] ** 3, 0])
        Zdot = np.concatenate((Xdot, -cls.Minv @ (KX + fnl)))
        return Zdot

    @classmethod
    def monodromy_ode(cls, t, dZdZ0, x1_interp):
        # ODE to solve for Monodromy matrix. dZdZ0dot = dg(Z)dZ . dZdZ0
        # estimate the value of x1(t) from interpolated curve
        M = cls.M
        K = cls.K
        knl = cls.Knl
        x1_t = spi.splev(t, x1_interp)
        dgdz = np.array([[0, 0, 1, 0],
                         [0, 0, 0, 1],
                         [-1 / M[0, 0] * K[0, 0] + knl * 3 * x1_t ** 2, -1 / M[0, 0] * K[0, 1], 0, 0],
                         [-1 / M[1, 1] * K[1, 0], -1 / M[1, 1] * K[1, 1], 0, 0]])
        dZdZ0dot = dgdz @ dZdZ0.reshape(4, 4)
        return dZdZ0dot.flatten()

    @classmethod
    def monodromy_centdiff(cls, t, Z0):
        # central difference
        eps = 1e-8
        M = np.zeros([4, 4])
        for i in range(4):
            Z0plus = Z0.copy()
            Z0plus[i] += eps
            ZTplus = np.array(odeint(cls.system_ode, Z0plus, t, tfirst=True))[-1, :]
            Z0mins = Z0.copy()
            Z0mins[i] -= eps
            ZTmins = np.array(odeint(cls.system_ode, Z0mins, t, tfirst=True))[-1, :]
            m = (ZTplus - ZTmins) / (2 * eps)
            M[:, i] = m
        return M

    @classmethod
    def time_solve(cls, T, Z0, pose_base, cont_params):
        nperiod = cont_params["shooting"]["single"]["nperiod"]
        nsteps = cont_params["shooting"]["single"]["nsteps_per_period"]
        rel_tol = cont_params["shooting"]["rel_tol"]

        # time integration, position and velocity
        t = np.linspace(0, T*nperiod, nsteps*nperiod)
        Z = np.array(odeint(cls.system_ode, Z0, t, rtol=rel_tol, tfirst=True))
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

        # Monodromy and augmented Jacobian
        # interpolate x1 = Z[0] as needed in monodromy time integration
        # odeint selects time points automatically so we need to have x1 at any t during integration
        x1_interp = spi.splrep(t, Z[:, 0])
        M = np.array(odeint(cls.monodromy_ode, np.eye(4).flatten(), t, args=(x1_interp,), tfirst=True))
        M = M[-1, :].reshape(4, 4)
        # M_centdiff = cls.monodromy_centdiff(t, Z0)  # central difference to check values
        M -= np.eye(len(M))
        dHdt = cls.system_ode(None, Z[-1, :])
        J = np.concatenate((M, dHdt.reshape(-1, 1)), axis=1)


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
