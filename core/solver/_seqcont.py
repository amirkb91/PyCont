import numpy as np
import scipy.linalg as spl
from ._cont_step import cont_step


def seqcont(self):
    cont_params = self.prob.cont_params
    cont_params_cont = cont_params["continuation"]
    forced = cont_params_cont["forced"]
    dofdata = self.prob.doffunction()
    N = dofdata["ndof_free"]
    twoN = 2 * N

    # first point solution
    X = self.X0
    pose_base = self.pose
    omega = 1.0
    tau = self.T0
    if cont_params["shooting"]["scaling"]:
        omega = 1 / self.T0
        tau = 1.0

    # continuation parameters
    step = cont_params_cont["s0"]
    direction = cont_params_cont["dir"]
    stepsign = -1 * direction  # corrections are always added

    # boolean mask to select inc from X (has no effect on single shooting)
    inc_mask = np.mod(np.arange(X.size), twoN) < N

    # --- MAIN CONTINUATION LOOP
    itercont = 1
    while True:
        # increment period
        tau_pred = tau + step * stepsign
        X_pred = X.copy()

        if (
            omega / tau_pred > cont_params_cont["fmax"]
            or omega / tau_pred < cont_params_cont["fmin"]
        ):
            print(f"Frequency {omega / tau_pred:.2e} Hz outside of specified boundary.")
            break

        # correction step
        itercorrect = 0
        while True:
            if itercorrect % cont_params_cont["iterjac"] == 0:
                sensitivity = True
            else:
                sensitivity = False

            [H, Jsim, pose, vel, energy, cvg_zerof] = self.prob.zerofunction(
                omega, tau_pred, X_pred, pose_base, cont_params, sensitivity=sensitivity
            )
            if not cvg_zerof:
                cvg_cont = False
                print(f"Zero function failed to converge with step = {step:.3e}.")
                break

            residual = spl.norm(H)

            if sensitivity:
                J = np.block([[Jsim[:, :-1]], [self.h]])

            if residual < cont_params_cont["tol"] and itercorrect >= cont_params_cont["itermin"]:
                cvg_cont = True
                break
            elif itercorrect > cont_params_cont["itermax"] or residual > 1e10:
                cvg_cont = False
                break

            self.log.screenout(
                iter=itercont,
                correct=itercorrect,
                res=residual,
                freq=omega / tau_pred,
                energy=energy,
                step=step,
            )

            # correction
            itercorrect += 1
            hx = self.h @ X_pred
            Z = np.vstack([H, hx.reshape(-1, 1)])
            if not forced:
                dx = spl.lstsq(J, -Z, cond=None, check_finite=False, lapack_driver="gelsd")[0]
            elif forced:
                dx = spl.solve(J, -Z, check_finite=False)
            X_pred[:] += dx[:, 0]

        if cvg_cont:
            self.log.screenout(
                iter=itercont,
                correct=itercorrect,
                res=residual,
                freq=omega / tau_pred,
                energy=energy,
                step=step,
                beta=0.0,
            )
            self.log.store(
                sol_pose=pose,
                sol_vel=vel,
                sol_T=tau_pred / omega,
                sol_energy=energy,
                sol_itercorrect=itercorrect,
                sol_step=step,
            )

            itercont += 1
            tau = tau_pred
            X = X_pred.copy()
            # update pose_base and set inc to zero, pose will have included inc from current sol
            pose_base = pose.copy()
            X[inc_mask] = 0.0

            # if cont_params["shooting"]["scaling"]:
            #     # reset tau to 1.0
            #     omega = omega / tau
            #     tau = 1.0

        # adaptive step size for next point
        if itercont > cont_params_cont["nadapt"] or not cvg_cont:
            step = cont_step(self, step, itercorrect, cvg_cont)

        if itercont > cont_params_cont["npts"]:
            print("Maximum number of continuation points reached.")
            break
        if cvg_cont and energy and energy > cont_params_cont["Emax"]:
            print(f"Energy {energy:.5e} exceeds Emax.")
            break
