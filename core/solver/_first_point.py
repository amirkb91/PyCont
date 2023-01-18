import numpy as np
import scipy.linalg as spl
from ._phase_condition import phase_condition


def first_point(self):
    print("Shooting first point.")
    print("Iter \t Residual")
    restart = self.prob.cont_params["first_point"]["restart"]["file_name"]
    dofdata = self.prob.doffunction()
    N = dofdata["ndof_free"]

    if not restart:
        iter_firstpoint = 0
        while True:
            if iter_firstpoint > self.prob.cont_params["first_point"]["itermax"]:
                raise Exception("Max number of iterations reached without convergence.")

            # residual and Jacobian with orthogonality to linear solution
            [H, J, self.pose_time0, self.vel_time0, pose_base_plus_inc, self.energy0, cvg_zerof] = \
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

            # correct X0 and T0
            iter_firstpoint += 1
            hx = np.matmul(self.h, self.X0)
            Z = np.vstack([H, hx.reshape(-1, 1), np.zeros(1)])
            dxt = spl.lstsq(J, -Z, cond=None, check_finite=False, lapack_driver="gelsy")[0]
            self.T0 += dxt[-1, 0]
            dx = dxt[:-1, 0]
            self.X0 += dx

    elif restart:
        pass

    # Compute Tangent
    if self.prob.cont_params["shooting"]["method"] == "single":
        # update pose_base and set inc to zero
        self.pose_base0 = pose_base_plus_inc
        self.X0[:N] = 0.0
        J[-1, :] = np.zeros(np.shape(J)[1])
    elif self.prob.cont_params["shooting"]["method"] == "multiple":
        # partition solution
        self.X0, self.pose_base0 = self.prob.partitionfunction(self.T0, self.X0, self.pose_base0, self.prob.cont_params)
        [_, J, self.pose_time0, self.vel_time0, _, self.energy0, _] = \
            self.prob.zerofunction(self.T0, self.X0, self.pose_base0, self.prob.cont_params)
        # size of X0 has changed so reconfigure phase condition matrix
        phase_condition(self)
        J = np.block([
            [J],
            [self.h, np.zeros((self.nphase, 1))],
            [np.zeros(np.shape(J)[1])]])
    J[-1, -1] = 1
    Z = np.zeros((np.shape(J)[0], 1))
    Z[-1] = 1
    self.tgt0 = spl.lstsq(J, Z, cond=None, check_finite=False, lapack_driver="gelsy")[0][:, 0]
    self.tgt0 /= spl.norm(self.tgt0)

    self.log.store(sol_X=self.X0, sol_T=self.T0, sol_tgt=self.tgt0, sol_pose_time=self.pose_time0,
                   sol_vel_time=self.vel_time0, sol_pose_base=self.pose_base0, sol_energy=self.energy0)
