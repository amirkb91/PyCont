import numpy as np
import scipy.linalg as spl


def first_point(self):
    print("Shooting first point.")
    print("Iter \t Residual")
    restart = self.prob.cont_params["first_point"]["restart"]["file_name"]
    fixF = self.prob.cont_params["first_point"]["restart"]["fixF"]

    if not restart or (restart and fixF):
        iter_firstpoint = 0
        while True:
            if iter_firstpoint > self.prob.cont_params["first_point"]["itermax"]:
                raise Exception("Max number of iterations reached without convergence.")

            [H, J, self.pose_time0, self.vel_time0, self.energy0, cvg_zerof] = \
                self.prob.zerofunction(self.T0, self.X0, self.pose_base0, self.prob.cont_params)

            if not cvg_zerof:
                raise Exception("Zero function failed.")

            residual = spl.norm(H)
            print(f"{iter_firstpoint} \t {residual:.5e}")

            if residual < self.prob.cont_params["continuation"]["tol"]:
                print("First point converged.")
                print("\n^-_-^-_-^-_-^-_-^-_-^-_-^-_-^-_-^-_-^\n")
                break

            iter_firstpoint += 1
            if not restart:
                # correct X0 and T0
                # Jacobian with orthogonality to linear solution
                J = np.vstack((J,
                               np.concatenate([self.h, np.zeros((self.nphase, 1))], axis=1),
                               np.concatenate([self.X0, np.zeros(1)])))
                hx = np.matmul(self.h, self.X0)
                H = np.vstack([H, hx.reshape(-1, 1), np.zeros(1)])
                dxt = spl.lstsq(J, -H, cond=None, check_finite=False, lapack_driver="gelsy")[0]
                self.T0 += dxt[-1, 0]
                dx = dxt[:-1, 0]
                self.X0 += dx

            elif restart and fixF:
                # correct X0 - ortho to restart solution?
                ortho = False
                if not ortho:
                    J = np.concatenate((M, self.h), axis=0)
                    hx = np.matmul(self.h, self.X0)
                    H = np.vstack([H, hx.reshape(-1, 1)])
                else:
                    J = np.concatenate((M, self.h, self.X0.reshape(-1, 1).T), axis=0)
                    hx = np.matmul(self.h, self.X0)
                    H = np.vstack([H, hx.reshape(-1, 1), np.zeros(1)])
                dx = spl.lstsq(J, -H, cond=None, check_finite=False, lapack_driver="gelsy")[0]
                self.X0 += dx[:, 0]

    elif restart and not fixF:
        self.prob.updatefunction(self.pose_base0)
        [H, M, dHdt, self.pose_time0, self.vel_time0, energy0, cvg_zerof] = \
            self.prob.zerofunction(self.T0, self.X0, self.prob.cont_params)
        residual = spl.norm(H)
        print(f"{0} \t {residual:.5e}")
        print("RESTARTING from previous run.")
        print("\n^-_-^-_-^-_-^-_-^-_-^-_-^-_-^-_-^-_-^\n")

    # Compute Tangent: Partition with Poincare sections for multiple shooting
    dofdata = self.prob.doffunction()
    N = dofdata["ndof_free"]
    twoN = 2 * dofdata["ndof_free"]
    npartition = self.prob.cont_params["shooting"]["npartition_multipleshooting"]
    nsteps = self.prob.cont_params["shooting"]["nsteps_per_period"]
    nsteps_per_partition = nsteps // npartition
    delta_S = 1 / npartition
    timesol0_partition_index = int(nsteps * delta_S) * np.arange(npartition)
    # overwrite pose_base so INC=0 as pose_base contains initial conditions
    self.pose_base0 = self.pose_time0[:, timesol0_partition_index]
    V = self.vel_time0[dofdata["free_dof"]][:, timesol0_partition_index]
    X = np.concatenate((np.zeros((N, npartition)), V))
    self.X0 = np.reshape(X, (-1), order='F')

    # tangent matrix: set T component to 1 and solve overdetermined system
    J = np.zeros((npartition * twoN + self.nphase + 1, npartition * twoN + 1))
    self.pose_time0 = np.zeros((np.shape(self.pose_base0)[0], (nsteps_per_partition + 1) * npartition))
    self.vel_time0 = np.zeros((dofdata["ndof_all"], (nsteps_per_partition + 1) * npartition))
    for ipart in range(npartition):
        i = ipart * twoN
        i1 = (ipart + 1) * twoN
        j = (ipart + 1) % npartition * twoN
        j1 = ((ipart + 1) % npartition + 1) * twoN
        p = ipart * (nsteps_per_partition + 1)
        p1 = (ipart + 1) * (nsteps_per_partition + 1)

        self.prob.updatefunction(self.pose_base0[:, ipart])
        [_, M, dHdt, self.pose_time0[:, p:p1], self.vel_time0[:, p:p1], _, _] = \
            self.prob.zerofunction(self.T0 * delta_S, X[:, ipart], self.prob.cont_params, mult=True)
        J[i:i1, i:i1] = M
        J[i:i1, j:j1] -= np.eye(twoN)
        J[i:i1, -1] = dHdt * delta_S
    J[npartition * twoN:npartition * twoN + self.nphase, :twoN] = self.h
    J[-1, -1] = 1
    Z = np.zeros((np.shape(J)[0], 1))
    Z[-1] = 1
    self.tgt0 = spl.lstsq(J, Z, cond=None, check_finite=False, lapack_driver="gelsy")[0][:, 0]
    self.tgt0 /= spl.norm(self.tgt0)

    # store solution in logger
    self.log.store(sol_X=self.X0.copy(), sol_T=self.T0.copy(), sol_tgt=self.tgt0.copy(),
                   sol_pose_time=self.pose_time0.copy(), sol_vel_time=self.vel_time0.copy(),
                   sol_pose_base=self.pose_base0.copy().reshape(-1, order='F'), sol_energy=self.energy0.copy())
