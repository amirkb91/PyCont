import numpy as np
import scipy.linalg as spl
from ._phase_condition import phase_condition


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

            # residual and Jacobian with orthogonality to linear solution
            [H, J, self.pose_time0, self.vel_time0, self.energy0, cvg_zerof] = \
                self.prob.zerofunction_firstpoint(self.T0, self.X0, self.pose_base0, self.prob.cont_params)
            J = np.block([
                [J],
                [self.h, np.zeros((self.nphase, 1))],
                [self.X0, np.zeros(1)]])

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
                hx = np.matmul(self.h, self.X0)
                H_all = np.vstack([H, hx.reshape(-1, 1), np.zeros(1)])
                dxt = spl.lstsq(J, -H_all, cond=None, check_finite=False, lapack_driver="gelsy")[0]
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

    # Compute Tangent
    if self.prob.cont_params["shooting"]["method"] == "single":
        pass
    elif self.prob.cont_params["shooting"]["method"] == "multiple":
        # partition solution
        X_partition, pose_base0_partition = \
         self.prob.partitionfunction(self.T0, self.X0, self.pose_base0, self.prob.cont_params)
        self.X0 = X_partition
        self.pose_base0 = pose_base0_partition
        [H, J, self.pose_time0, self.vel_time0, _, _] = \
            self.prob.zerofunction(self.T0, self.X0, self.pose_base0, self.prob.cont_params)
        # size of X0 has changed so reconfigure phase condition matrix
        phase_condition(self)

    J = np.block([
        [J],
        [self.h, np.zeros((self.nphase, 1))],
        [np.zeros(len(self.X0)), np.ones(1)]])
    Z = np.zeros((np.shape(J)[0], 1))
    Z[-1] = 1
    self.tgt0 = spl.lstsq(J, Z, cond=None, check_finite=False, lapack_driver="gelsy")[0][:, 0]
    self.tgt0 /= spl.norm(self.tgt0)

    # store solution in logger
    self.log.store(sol_X=self.X0.copy(), sol_T=self.T0.copy(), sol_tgt=self.tgt0.copy(),
                   sol_pose_time=self.pose_time0.copy(), sol_vel_time=self.vel_time0.copy(),
                   sol_pose_base=self.pose_base0.copy().reshape(-1, order='F'), sol_energy=self.energy0.copy())
