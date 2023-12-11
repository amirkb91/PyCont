import numpy as np
from scipy.integrate import odeint


class Duffing:
    # parameters of the model EoM
    # xddot + delta*xdot + alpha*x + beta*x^3 = 0
    delta = 0.0
    alpha = 1.0
    beta = 1.0

    # finite element data, 1 dof model
    free_dof = np.array([0])
    ndof_all = 1
    ndof_fix = 0
    ndof_free = 1

    @classmethod
    def system_ode(cls, t, X):
        # ODE of model: Xdot(t) = g(X(t))
        x = X[0]
        xdot = X[1]
        f = cls.delta * xdot + cls.alpha * x + cls.beta * x**3
        Xdot = np.array([xdot, -f])
        return Xdot

    @classmethod
    def augsystem_ode(cls, t, X_aug):
        # Augemented ODE of model + Monodromy, to be solved together
        # System: Xdot(t) = g(X(t))
        # Monodromy: dXdX0dot = dg(X)dX . dXdX0
        X, dXdX0 = X_aug[:2], X_aug[2:]
        x = X[0]
        xdot = X[1]
        f = cls.delta * xdot + cls.alpha * x + cls.beta * x**3
        Xdot = np.array([xdot, -f])
        dgdX = np.array([[0, 1], [-cls.alpha - 3 * cls.beta * x**2, -cls.delta]])
        dXdX0dot = dgdX @ dXdX0.reshape(2, 2)

        return np.concatenate([Xdot, dXdX0dot.flatten()])

    @classmethod
    def eigen_solve(cls):
        frq = np.array([[cls.alpha]])  # natural frequency
        frq = np.sqrt(frq) / (2 * np.pi)
        eig = np.array([[1.0]])

        # initial position taken as zero
        pose0 = 0.0

        return eig, frq, pose0

    @classmethod
    def time_solve(
        cls, omega, tau, X, pose_base, cont_params, return_time=False, sensitivity=True
    ):
        nperiod = cont_params["shooting"]["single"]["nperiod"]
        nsteps = cont_params["shooting"]["single"]["nsteps_per_period"]
        rel_tol = cont_params["shooting"]["rel_tol"]
        T = tau / omega

        # Add position to increment and do time sim to get solution and Monodromy M
        X_total = X.copy()
        X_total[0] += pose_base
        t = np.linspace(0, T * nperiod, nsteps * nperiod + 1)
        initial_cond_aug = np.concatenate((X_total, np.eye(2).flatten()))
        Xsol_aug = np.array(
            odeint(cls.augsystem_ode, initial_cond_aug, t, rtol=rel_tol, tfirst=True)
        )
        Xsol, M = Xsol_aug[:, :2], Xsol_aug[-1, 2:].reshape(2, 2)

        # periodicity condition
        H = Xsol[-1, :] - Xsol[0, :]
        H = H.reshape(-1, 1)

        # Augmented Jacobian (dHdX0 and dHdt)
        dHdX0 = M - np.eye(2)
        dHdt = cls.system_ode(None, Xsol[-1, :]) * nperiod
        J = np.concatenate((dHdX0, dHdt.reshape(-1, 1)), axis=1)

        # solution pose and vel taken from time 0
        pose = Xsol[0, 0]
        vel = Xsol[0, 1]

        # Energy, conservative model so take mean of all time
        E = np.zeros(nsteps * nperiod + 1)
        for i in range(nsteps * nperiod + 1):
            x = Xsol[i, 0]
            xdot = Xsol[i, 1]
            Fnl = 0.25 * cls.beta * x**4
            E[i] = 0.5 * (xdot**2 + cls.alpha * x**2) + Fnl
        energy = np.mean(E)

        cvg = True
        return H, J, M, pose, vel, energy, cvg

    # @classmethod
    # def time_solve_multiple(cls, omega, tau, X, pose_base, cont_params):
    #     npartition = cont_params["shooting"]["multiple"]["npartition"]
    #     nsteps = cont_params["shooting"]["multiple"]["nsteps_per_partition"]
    #     rel_tol = cont_params["shooting"]["rel_tol"]
    #     N = cls.ndof_free
    #     twoN = 2 * N
    #     delta_S = 1 / npartition
    #     T = tau / omega

    #     # initialise
    #     J = np.zeros((npartition * twoN, npartition * twoN + 1))
    #     pose_time = np.zeros((cls.ndof_all, (nsteps + 1) * npartition))
    #     vel_time = np.zeros((cls.ndof_all, (nsteps + 1) * npartition))
    #     E = np.zeros([nsteps + 1, npartition])

    #     for ipart in range(npartition):
    #         # index values required for looping the partitions
    #         i = ipart * twoN
    #         i1 = (ipart + 1) * twoN
    #         j = (ipart + 1) % npartition * twoN
    #         j1 = ((ipart + 1) % npartition + 1) * twoN
    #         p = ipart * (nsteps + 1)
    #         p1 = (ipart + 1) * (nsteps + 1)

    #         # get total displacements from x and pose_base and do time integration
    #         X_total = X[i:i1].copy()
    #         X_total[:N] += pose_base[:, ipart].flatten().copy()
    #         t = np.linspace(0, T * delta_S, nsteps + 1)
    #         initial_cond_aug = np.concatenate((X_total, np.eye(4).flatten()))
    #         Xsol_aug = np.array(
    #             odeint(cls.augsystem_ode, initial_cond_aug, t, rtol=rel_tol, tfirst=True)
    #         )
    #         Xsol, M = Xsol_aug[:, :twoN], Xsol_aug[-1, twoN:].reshape(twoN, twoN)
    #         pose_time[:, p:p1] = Xsol[:, :N].T
    #         vel_time[:, p:p1] = Xsol[:, N:].T

    #         # Monodromy and augmented Jacobian
    #         dHdt = cls.system_ode(None, Xsol[-1, :]) * delta_S
    #         J[i:i1, i:i1] = M
    #         J[i:i1, j:j1] -= np.eye(twoN)
    #         J[i:i1, -1] = dHdt

    #         # Energy
    #         for k in range(nsteps + 1):
    #             x = Xsol[k, :N]
    #             xdot = Xsol[k, N:]
    #             Fnl = 0.25 * cls.Knl * x[0]**4
    #             E[k, ipart] = 0.5 * (xdot.T @ cls.M @ xdot + x.T @ cls.K @ x) + Fnl

    #     # time solution indicies which enclose each partition & order of the partitions for periodicity
    #     timesol_partition_index_start = nsteps * np.arange(npartition) + np.arange(npartition)
    #     timesol_partition_index_end = timesol_partition_index_start - 1
    #     block_order = (np.arange(npartition) + 1) % npartition

    #     # Periodicity condition for all partitions
    #     H1 = (
    #         pose_time[cls.free_dof][:, timesol_partition_index_end[block_order]] -
    #         pose_time[cls.free_dof][:, timesol_partition_index_start[block_order]]
    #     )
    #     H2 = (
    #         vel_time[cls.free_dof][:, timesol_partition_index_end[block_order]] -
    #         vel_time[cls.free_dof][:, timesol_partition_index_start[block_order]]
    #     )
    #     H = np.reshape(np.concatenate([H1, H2]), (-1, 1), order="F")

    #     # solution pose and vel taken from time 0 for each partition
    #     pose = pose_time[:, timesol_partition_index_start]
    #     vel = vel_time[:, timesol_partition_index_start]

    #     # Energy, conservative model so take mean of all time
    #     energy = np.mean(E)

    #     cvg = True
    #     return H, J, pose, vel, energy, cvg

    # @classmethod
    # def partition_singleshooting_solution(cls, omega, tau, X, pose_base, cont_params):
    #     npartition = cont_params["shooting"]["multiple"]["npartition"]
    #     nsteps = cont_params["shooting"]["multiple"]["nsteps_per_partition"]
    #     rel_tol = cont_params["shooting"]["rel_tol"]
    #     N = cls.ndof_free
    #     slicing_index = nsteps * np.arange(npartition)
    #     T = tau / omega

    #     # do time integration along whole orbit before slicing
    #     X_total = X.copy()
    #     X_total[:N] += pose_base.flatten().copy()
    #     t = np.linspace(0, T, nsteps * npartition + 1)
    #     Xsol = np.array(odeint(cls.system_ode, X_total, t, rtol=rel_tol, tfirst=True))
    #     pose_time = Xsol[:, :N].T
    #     vel_time = Xsol[:, N:].T
    #     pose = pose_time[:, slicing_index]
    #     V = vel_time[:, slicing_index]
    #     # set inc to zero as solution stored in pose, keep velocity
    #     Xsol = np.concatenate((np.zeros((N, npartition)), V))
    #     Xsol = np.reshape(Xsol, (-1), order="F")
    #     return Xsol, pose

    @classmethod
    def get_fe_data(cls):
        return {
            "free_dof": cls.free_dof,
            "ndof_all": cls.ndof_all,
            "ndof_fix": cls.ndof_fix,
            "ndof_free": cls.ndof_free,
        }

    # @classmethod
    # def monodromy_centdiff(cls, t, X0):
    #     # central difference calculation of the monodromy matrix
    #     # can be used to check values from ode
    #     eps = 1e-8
    #     M = np.zeros([4, 4])
    #     for i in range(4):
    #         X0plus = X0.copy()
    #         X0plus[i] += eps
    #         XTplus = np.array(odeint(cls.system_ode, X0plus, t, tfirst=True))[-1, :]
    #         X0mins = X0.copy()
    #         X0mins[i] -= eps
    #         XTmins = np.array(odeint(cls.system_ode, X0mins, t, tfirst=True))[-1, :]
    #         m = (XTplus - XTmins) / (2 * eps)
    #         M[:, i] = m
    #     return M
