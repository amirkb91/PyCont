import numpy as np
import scipy.linalg as spl
from ._cont_step import cont_step


def psacont(self):
    frml = self.prob.cont_params["continuation"]["tangent"].lower()
    dofdata = self.prob.doffunction()
    N = dofdata["ndof_free"]
    twoN = 2 * N

    # first point solution
    T = self.T0.copy()
    X = self.X0.copy()
    pose_base = self.pose.copy()
    tgt = self.tgt0.copy()
    energy = self.energy0.copy()

    # continuation step and direction
    step = self.prob.cont_params["continuation"]["s0"]
    direction = self.prob.cont_params["continuation"]["dir"]
    stepsign = -1 * direction * np.sign(tgt[-1])    # corrections are always added

    # continuation loop
    itercont = 1
    while True:
        # prediction step along tangent
        T_pred = T + tgt[-1] * step * stepsign
        X_pred = X + tgt[:-1] * step * stepsign

        if 1 / T_pred > self.prob.cont_params["continuation"]["fmax"] or \
                1 / T_pred < self.prob.cont_params["continuation"]["fmin"]:
            print("Frequency outside of specified boundary.")
            break

        # correction step
        itercorrect = 0
        while True:
            # residual and block Jacobian
            [H, J, pose, vel, energy_next, cvg_zerof] = \
                self.prob.zerofunction(T_pred, X_pred, pose_base, self.prob.cont_params)
            J = np.block([[J], [self.h, np.zeros((self.nphase, 1))], [tgt]])

            if not cvg_zerof:
                cvg_cont = False
                print("Zero function failed to converge.")
                break

            residual = spl.norm(H)
            if (residual < self.prob.cont_params["continuation"]["tol"] and
                    itercorrect >= self.prob.cont_params["continuation"]["itermin"]):
                cvg_cont = True
                break
            elif itercorrect >= self.prob.cont_params["continuation"]["itermax"]:
                cvg_cont = False
                break

            # apply corrections orthogonal to tangent
            self.log.screenout(
                iter=itercont,
                correct=itercorrect,
                res=residual,
                freq=1 / T_pred,
                energy=energy_next,
                step=step)
            itercorrect += 1
            hx = np.matmul(self.h, X_pred)
            Z = np.vstack([H, hx.reshape(-1, 1), np.zeros(1)])
            dxt = spl.lstsq(J, -Z, cond=None, check_finite=False, lapack_driver="gelsd")[0]
            T_pred += dxt[-1, 0]
            dx = dxt[:-1, 0]
            X_pred += dx

        if cvg_cont:
            # find new tangent with converged solution
            # peeters Jacobian is different for tangent update
            if frml == "peeters":
                J[-1, :] = np.zeros(np.shape(J)[1])
                J[-1, -1] = 1
            Z = np.zeros((np.shape(J)[0], 1))
            Z[-1] = 1
            tgt_next = spl.lstsq(
                J, Z, cond=None, check_finite=False, lapack_driver="gelsd")[0][:, 0]
            tgt_next /= spl.norm(tgt_next)
            # tgt_next /= spl.norm(tgt_next, np.inf)

            # calculate beta and check against betamax if requested, fail convergence if check fails
            beta = np.degrees(
                np.arccos((tgt_next.T @ tgt) / (spl.norm(tgt) * spl.norm(tgt_next))))
            if (self.prob.cont_params["continuation"]["betacontrol"] and
                    beta > self.prob.cont_params["continuation"]["betamax"]):
                print("Beta exceeds maximum angle.")
                cvg_cont = False
            else:
                # passed check, finalise and update for next step
                if frml == "peeters":
                    stepsign = np.sign(stepsign * tgt_next.T @ tgt)

                self.log.store(
                    sol_pose=pose,
                    sol_vel=vel,
                    sol_T=T_pred,
                    sol_tgt=tgt_next,
                    sol_energy=energy_next,
                    sol_beta=beta,
                    sol_itercorrect=itercorrect,
                    sol_step=step)
                self.log.screenout(
                    iter=itercont,
                    correct=itercorrect,
                    res=residual,
                    freq=1 / T_pred,
                    energy=energy_next,
                    step=step,
                    beta=beta)

                itercont += 1
                T = T_pred
                X = X_pred[:]
                tgt = tgt_next[:]
                energy = energy_next
                # update pose_base and set inc to zero (slice 0:N on each partition)
                pose_base = pose[:]
                X[np.mod(np.arange(X.size), twoN) < N] = 0.0

        # adaptive step size for next point
        if itercont > self.prob.cont_params["continuation"]["nadapt"] or not cvg_cont:
            step = cont_step(self, step, itercorrect, cvg_cont)

        if itercont > self.prob.cont_params["continuation"]["npts"]:
            print("Maximum number of continuation points reached.")
            break
        if energy > self.prob.cont_params["continuation"]["Emax"]:
            print("Energy exceeds Emax.")
            break
