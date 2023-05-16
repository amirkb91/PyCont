import numpy as np
from scipy.integrate import odeint
import scipy.linalg as spl
import scipy.interpolate as spi


class Cubic_Spring:
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
        x = X[:cls.ndof_free]
        xdot = X[cls.ndof_free:]
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
                         [-1 / M[0, 0] * (K[0, 0] + knl * 3 * x1_t ** 2), -1 / M[0, 0] * K[0, 1], 0, 0],
                         [-1 / M[1, 1] * K[1, 0], -1 / M[1, 1] * K[1, 1], 0, 0]])
        dXdX0dot = dgdz @ dXdX0.reshape(4, 4)
        return dXdX0dot.flatten()

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

        # initial position taken as zero for both dofs.
        pose0 = np.zeros(cls.ndof_free)

        return X0, T0, pose0    
    
    @classmethod
    def time_solve(cls, T, X, pose_base, cont_params):
        nperiod = cont_params["shooting"]["single"]["nperiod"]
        nsteps = cont_params["shooting"]["single"]["nsteps_per_period"]
        rel_tol = cont_params["shooting"]["rel_tol"]
        N = cls.ndof_free
        twoN = 2 * N

        # Add increment onto pose and do time sim
        X_total = X.copy()
        X_total[:N] += pose_base.flatten().copy()
        t = np.linspace(0, T * nperiod, nsteps * nperiod + 1)
        Xsol = np.array(odeint(cls.system_ode, X_total, t, rtol=rel_tol, tfirst=True))
        # solution pose and vel taken from time 0
        pose = Xsol[0, :N].T
        vel = Xsol[0, N:].T

        # periodicity condition
        H = Xsol[-1, :] - Xsol[0, :]
        H = H.reshape(-1, 1)
        # Energy, conservative system so take initial Xsol values
        x = Xsol[0, :N]
        xdot = Xsol[0, N:]
        fnl = np.array([cls.Knl * x[0] ** 3, 0])
        energy = 0.5 * (xdot.T @ cls.M @ xdot + x.T @ cls.K @ x + x @ fnl)

        # Monodromy and augmented Jacobian
        # interpolate x1 = Xsol[0] as needed in monodromy time integration
        # odeint selects time points automatically so we need to have x1 at any t during integration
        x1_interp = spi.splrep(t, Xsol[:, 0])
        M = np.array(odeint(cls.monodromy_ode, np.eye(4).flatten(), t, args=(x1_interp,), tfirst=True))
        M = M[-1, :].reshape(twoN, twoN)
        M -= np.eye(twoN)
        dHdt = cls.system_ode(None, Xsol[-1, :]) * nperiod
        J = np.concatenate((M, dHdt.reshape(-1, 1)), axis=1)

        cvg = True
        return H, J, pose, vel, energy, cvg

    @classmethod
    def time_solve_multiple(cls, T, X0, pose_base, cont_params):
        npartition = cont_params["shooting"]["multiple"]["npartition"]
        nsteps = cont_params["shooting"]["multiple"]["nsteps_per_partition"]
        rel_tol = cont_params["shooting"]["rel_tol"]
        N = cls.ndof_free
        twoN = 2 * N
        delta_S = 1 / npartition

        # initialise
        J = np.zeros((npartition * twoN, npartition * twoN + 1))
        pose_time = np.zeros((np.shape(pose_base)[0], (nsteps + 1) * npartition))
        vel_time = np.zeros((cls.ndof_all, (nsteps + 1) * npartition))

        for ipart in range(npartition):
            # index values required for looping the partitions
            i = ipart * twoN
            i1 = (ipart + 1) * twoN
            j = (ipart + 1) % npartition * twoN
            j1 = ((ipart + 1) % npartition + 1) * twoN
            p = ipart * (nsteps + 1)
            p1 = (ipart + 1) * (nsteps + 1)

            # get total displacements from x and pose_base and do time integration
            X0_total = X0[i:i1].copy()
            X0_total[:cls.ndof_free] += pose_base[:, ipart].flatten().copy()
            t = np.linspace(0, T * delta_S, nsteps + 1)
            X = np.array(odeint(cls.system_ode, X0_total, t, rtol=rel_tol, tfirst=True))
            pose_time[:, p:p1] = X[:, :N].T
            vel_time[:, p:p1] = X[:, N:].T

            # Monodromy and augmented Jacobian
            # interpolate x1 = X[0] as needed in monodromy time integration
            # odeint selects time points automatically so we need to have x1 at any t during integration
            x1_interp = spi.splrep(t, X[:, 0])
            M = np.array(odeint(cls.monodromy_ode, np.eye(4).flatten(), t, args=(x1_interp,), tfirst=True))
            M = M[-1, :].reshape(twoN, twoN)
            dHdt = cls.system_ode(None, X[-1, :]) * delta_S
            J[i:i1, i:i1] = M
            J[i:i1, j:j1] -= np.eye(twoN)
            J[i:i1, -1] = dHdt

        # time solution indicies which enclose each partition & order of the partitions for periodicity
        timesol_partition_index_start = nsteps * np.arange(npartition) + np.arange(npartition)
        timesol_partition_index_end = timesol_partition_index_start - 1
        block_order = (np.arange(npartition) + 1) % npartition

        # Periodicity condition for all partitions
        H1 = pose_time[cls.free_dof][:, timesol_partition_index_end[block_order]] - \
             pose_time[cls.free_dof][:, timesol_partition_index_start[block_order]]
        H2 = vel_time[cls.free_dof][:, timesol_partition_index_end[block_order]] - \
             vel_time[cls.free_dof][:, timesol_partition_index_start[block_order]]
        H = np.reshape(np.concatenate([H1, H2]), (-1, 1), order='F')

        # update pose_base with the increments included
        pose_base_plus_inc = pose_time[:, timesol_partition_index_start]

        # Energy, conservative system so take initial X values from final partition
        x = X[0, :N]
        xdot = X[0, N:]
        fnl = np.array([cls.Knl * x[0] ** 3, 0])
        energy = 0.5 * (xdot.T @ cls.M @ xdot + x.T @ cls.K @ x + x @ fnl)

        cvg = True
        return H, J, pose_time, vel_time, pose_base_plus_inc, energy, cvg



    @classmethod
    def partition_singleshooting_solution(cls, T, X0, pose_base, cont_params):
        npartition = cont_params["shooting"]["multiple"]["npartition"]
        nsteps = cont_params["shooting"]["multiple"]["nsteps_per_partition"]
        rel_tol = cont_params["shooting"]["rel_tol"]
        slicing_index = nsteps * np.arange(npartition)

        # nsteps has to equal to total steps for multiple shooting so solution can be partitioned correctly
        # get total displacements from x and pose_base and do time integration
        X0_total = X0.copy()
        X0_total[:cls.ndof_free] += pose_base.flatten().copy()
        t = np.linspace(0, T, nsteps * npartition + 1)
        X = np.array(odeint(cls.system_ode, X0, t, rtol=rel_tol, tfirst=True))
        pose_time = X[:, :cls.ndof_free].T
        vel_time = X[:, cls.ndof_free:].T
        V = vel_time[cls.free_dof][:, slicing_index]
        # update pose_base and set inc to zero
        pose_base = pose_time[:, slicing_index]
        X = np.concatenate((np.zeros((cls.ndof_free, npartition)), V))
        X = np.reshape(X, (-1), order='F')

        return X, pose_base

    @classmethod
    def get_fe_data(cls):
        return {"free_dof": cls.free_dof, "ndof_all": cls.ndof_all, "ndof_fix": cls.ndof_fix,
                "ndof_free": cls.ndof_free}

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
