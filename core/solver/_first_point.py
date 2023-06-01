import numpy as np
import scipy.linalg as spl
from ._phase_condition import phase_condition


def first_point(self):
    restart = self.prob.cont_params["first_point"]["restart"]["file_name"]
    recompute_tangent = self.prob.cont_params["first_point"]["restart"]["recompute_tangent"]
    method = self.prob.cont_params["shooting"]["method"]
    dofdata = self.prob.doffunction()
    N = dofdata["ndof_free"]

    if not restart:
        iter_firstpoint = 0        
        linearsol = self.X0.copy()  # velocities are zero so no scaling needed

        while True:
            if iter_firstpoint > self.prob.cont_params["first_point"]["itermax"]:
                raise Exception("Max number of iterations reached without convergence.")

            # residual and Jacobian with orthogonality to linear solution
            [H, J, self.pose, self.vel, self.energy0, cvg_zerof] = self.prob.zerofunction_firstpoint(
                self.omega, self.tau, self.X0, self.pose0, self.prob.cont_params)
            J = np.block([[J], [self.h, np.zeros((self.nphase, 1))], [linearsol, np.zeros(1)]])
            if not cvg_zerof:
                raise Exception("Zero function failed.")

            residual = spl.norm(H)
            self.log.screenout(iter=0, correct=iter_firstpoint, res=residual,
                               freq=self.omega/self.tau, energy=self.energy0)

            if residual < self.prob.cont_params["continuation"]["tol"]:
                break

            # correct X0 and tau
            iter_firstpoint += 1
            hx = np.matmul(self.h, self.X0)
            Z = np.vstack([H, hx.reshape(-1, 1), np.zeros(1)])
            dxt = spl.lstsq(J, -Z, cond=None, check_finite=False, lapack_driver="gelsd")[0]
            self.tau += dxt[-1, 0]
            dx = dxt[:-1, 0]
            self.X0 += dx

        # set inc to zero as solution stored in pose, keep velocity
        self.X0[:N] = 0.0
        # Compute Tangent
        if method == "single":
            J[-1, :] = np.zeros(np.shape(J)[1])
        elif method == "multiple":
            # partition solution
            self.X0, self.pose = self.prob.partitionfunction(
                self.T0, self.X0, self.pose, self.prob.cont_params)
            [_, J, self.pose, self.vel, self.energy0, _] = \
                self.prob.zerofunction(self.T0, self.X0, self.pose, self.prob.cont_params)
            # size of X0 has changed so reconfigure phase condition matrix
            phase_condition(self)
            J = np.block([[J], [self.h, np.zeros((self.nphase, 1))], [np.zeros(np.shape(J)[1])]])
        J[-1, -1] = 1
        Z = np.zeros((np.shape(J)[0], 1))
        Z[-1] = 1
        self.tgt0 = spl.lstsq(J, Z, cond=None, check_finite=False, lapack_driver="gelsd")[0][:, 0]
        self.tgt0 /= spl.norm(self.tgt0)

        self.log.store(
            sol_pose=self.pose, sol_vel=self.vel, sol_T=self.tau/self.omega, sol_tgt=self.tgt0,
            sol_energy=self.energy0, sol_itercorrect=iter_firstpoint, sol_step=0)

    elif restart:
        if method == "single":
            # residual and Jacobian and Compute Tangent
            [H, J, self.pose, self.vel, self.energy0, cvg_zerof] = self.prob.zerofunction_firstpoint(
                self.T0, self.X0, self.pose0, self.prob.cont_params)
            residual = spl.norm(H)
            if recompute_tangent:
                J = np.block(
                    [[J], [self.h, np.zeros((self.nphase, 1))], [np.zeros(np.shape(J)[1])]])
                J[-1, -1] = 1
                Z = np.zeros((np.shape(J)[0], 1))
                Z[-1] = 1
                self.tgt0 = spl.lstsq(
                    J, Z, cond=None, check_finite=False, lapack_driver="gelsd")[0][:, 0]
                self.tgt0 /= spl.norm(self.tgt0)

            self.log.screenout(
                iter=0, correct=0, res=residual, freq=1 / self.T0, energy=self.energy0)
            self.log.store(sol_pose=self.pose, sol_vel=self.vel, sol_T=self.T0, sol_tgt=self.tgt0,
                           sol_energy=self.energy0, sol_itercorrect=0, sol_step=0)

        elif method == "multiple":
            pass
