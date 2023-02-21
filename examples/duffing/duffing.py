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
    def system_ode(cls, t, X):
        # ODE of the mass spring system. Xdot(t) = g(X(t))
        x = X[:2]
        xdot = X[2:]
        KX = cls.K @ x
        fnl = np.array([cls.Knl * x[0] ** 3, 0])
        Xdot = np.concatenate((xdot, -cls.Minv @ (KX + fnl)))
        return Xdot

    @classmethod
    def monodromy_ode(cls, t, dXdX0, x1_interp):
        # ODE to solve for Monodromy matrix. dXdX0dot = dg(X)dX . dXdX0
        # estimate the value of x1(t) from interpolated curve
        M = cls.M
        K = cls.K
        knl = cls.Knl
        x1_t = spi.splev(t, x1_interp)
        dgdz = np.array([[0, 0, 1, 0],
                         [0, 0, 0, 1],
                         [-1 / M[0, 0] * K[0, 0] + knl * 3 * x1_t ** 2, -1 / M[0, 0] * K[0, 1], 0, 0],
                         [-1 / M[1, 1] * K[1, 0], -1 / M[1, 1] * K[1, 1], 0, 0]])
        dXdX0dot = dgdz @ dXdX0.reshape(4, 4)
        return dXdX0dot.flatten()

    @classmethod
    def monodromy_centdiff(cls, t, X0):
        # central difference
        eps = 1e-8
        M = np.zeros([4, 4])
        for i in range(4):
            X0plus = X0.copy()
            X0plus[i] += eps
            XTplus = np.array(odeint(cls.system_ode, X0plus, t, tfirst=True))[-1, :]
            X0mins = X0.copy()
            X0mins[i] -= eps
            XTmins = np.array(odeint(cls.system_ode, X0mins, t, tfirst=True))[-1, :]
            m = (XTplus - XTmins) / (2 * eps)
            M[:, i] = m
        return M

    @classmethod
    def time_solve(cls, T, X0, pose_base, cont_params):
        nperiod = cont_params["shooting"]["single"]["nperiod"]
        nsteps = cont_params["shooting"]["single"]["nsteps_per_period"]
        rel_tol = cont_params["shooting"]["rel_tol"]

        # get total displacements from x and pose_base. positions stored in pose_base, increments in X0
        X0_total = X0.copy()
        X0_total[0:2] += pose_base.copy()

        # time integration, position and velocity
        t = np.linspace(0, T * nperiod, nsteps * nperiod)
        X = np.array(odeint(cls.system_ode, X0_total, t, rtol=rel_tol, tfirst=True))
        pose_time = X[:, :2].T
        vel_time = X[:, 2:].T

        # periodicity condition
        H = X[-1, :] - X[0, :]
        H = H.reshape(-1, 1)

        # Energy, conservative system so take initial X values
        x = X[0, :2]
        xdot = X[0, 2:]
        fnl = np.array([cls.Knl * x[0] ** 3, 0])
        energy = 0.5 * (xdot.T @ cls.M @ xdot + x.T @ cls.K @ x + x @ fnl)

        # Monodromy and augmented Jacobian
        # interpolate x1 = X[0] as needed in monodromy time integration
        # odeint selects time points automatically so we need to have x1 at any t during integration
        x1_interp = spi.splrep(t, X[:, 0])
        # M = np.array(odeint(cls.monodromy_ode, np.eye(4).flatten(), t, args=(x1_interp,), tfirst=True))
        # M = M[-1, :].reshape(4, 4)
        M = cls.monodromy_centdiff(t, X0_total)  # central difference to check values
        M -= np.eye(len(M))
        dHdt = cls.system_ode(None, X[-1, :])
        J = np.concatenate((M, dHdt.reshape(-1, 1)), axis=1)

        cvg = True
        # all position variables are required for the Lie group code only, can be set to None
        pose_base_plus_inc = X0_total[0:2].copy()

        return H, J, pose_time, vel_time, pose_base_plus_inc, energy, cvg

    @classmethod
    def eigen_solve(cls, cont_params):
        # Continuation variables initial guess from eigenvalues
        frq, eig = spl.eigh(cls.K, cls.M)
        frq = np.sqrt(frq) / (2 * np.pi)

        nnm = cont_params["first_point"]["eig_start"]["NNM"]
        scale = cont_params["first_point"]["eig_start"]["scale"]
        x0 = eig[:, nnm - 1] * scale
        v0 = np.zeros_like(x0)
        X0 = np.concatenate([x0, v0])
        T0 = 1 / frq[nnm - 1]

        # pose base taken as zero for this test case. Would be different for Lie Group FE solver.
        pose_base0 = np.zeros(2)

        return X0, T0, pose_base0

    @classmethod
    def get_fe_data(cls):
        return {"free_dof": cls.free_dof, "ndof_all": cls.ndof_all, "ndof_fix": cls.ndof_fix,
                "ndof_free": cls.ndof_free}
