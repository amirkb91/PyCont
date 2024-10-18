import numpy as np
import scipy.linalg as spl
from collections import namedtuple
from ._cont_step import cont_step
from examples.beam_cpp.Frame import Frame

import warnings

warnings.filterwarnings("ignore", category=spl.LinAlgWarning)


def psacont(self):
    cont_params = self.prob.cont_params
    cont_params_cont = cont_params["continuation"]
    frml = cont_params_cont["tangent"].lower()
    forced = cont_params_cont["forced"]
    dofdata = self.prob.doffunction()
    N = dofdata["ndof_free"]
    twoN = 2 * N
    cvg_sol = namedtuple("converged_solution", "X, T, H, J")

    # first point converged solution
    X = self.X0
    pose_base = self.pose
    pose_ref = self.pose_ref  # undeformed pose
    tgt = self.tgt0
    omega = 1.0
    tau = self.T0
    if cont_params["shooting"]["scaling"]:
        omega = 1 / self.T0
        tau = 1.0

    # continuation parameters
    step = cont_params_cont["s0"]
    direction = cont_params_cont["dir"]
    stepsign = -1 * direction * np.sign(tgt[-1])  # corrections are always added

    # boolean masks to select inc and vel from X (has no effect on single shooting)
    inc_mask = np.mod(np.arange(X.size), twoN) < N
    vel_mask = ~inc_mask

    # --- MAIN CONTINUATION LOOP
    itercont = 1
    while True:
        # prediction step along tangent
        tau_pred = tau + tgt[-1] * step * stepsign
        X_pred = X + tgt[:-1] * step * stepsign

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
            residual = normalise_residual(residual, pose_base, pose_ref, dofdata)

            if not sensitivity:
                # Broyden's Jacobian update
                deltaX = (
                    np.append(X_pred, tau_pred / omega) - np.append(soldata.X, soldata.T)
                ).reshape(-1, 1)
                deltaf = H - soldata.H
                Jsim = soldata.J + 1 / spl.norm(deltaX) * (deltaf - soldata.J @ deltaX) @ deltaX.T

            J = np.block([[Jsim], [self.h, np.zeros((self.nphase, 1))], [tgt]])
            soldata = cvg_sol(X_pred.copy(), tau_pred / omega, H.copy(), Jsim.copy())

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
                step=stepsign * step,
            )

            # apply corrections orthogonal to tangent
            itercorrect += 1
            Jcr = J.copy()
            # Jcr[-1, twoN:-1] = 0.0  # ortho only 1st partition and T (no effect single shooting)
            hx = self.h @ X_pred
            Z = np.vstack([H, hx.reshape(-1, 1), np.zeros(1)])
            if not forced:
                dxt = spl.lstsq(Jcr, -Z, cond=None, check_finite=False, lapack_driver="gelsd")[0]
            elif forced:
                dxt = spl.solve(Jcr, -Z, check_finite=False)
            tau_pred += dxt[-1, 0]
            dx = dxt[:-1, 0]
            X_pred += dx

        if cvg_cont:
            # find new tangent with converged solution
            if frml == "secant":
                # Find difference between solutions at current and previous step
                # For non-Lie group formulations, this is already equal to X_pred, since pose-pose_base = X_pred[:N]
                # For Lie group formulations, relative inc between pose and pose_base should be found with log map
                tgt_next = np.concatenate(
                    (X_pred[inc_mask], (X_pred - X)[vel_mask], [tau_pred - tau])
                )
            else:
                if frml == "peeters":
                    # remove tgt from Jacobian and fix period component to 1
                    J[-1, :] = 0.0
                    J[-1, -1] = 1
                Z = np.zeros((np.shape(J)[0], 1))
                Z[-1] = 1.0
                if not forced:
                    tgt_next = spl.lstsq(
                        J, Z, cond=None, check_finite=False, lapack_driver="gelsd"
                    )[0][:, 0]
                elif forced:
                    tgt_next = spl.solve(J, Z, check_finite=False)[:, 0]
            tgt_next /= spl.norm(tgt_next)
            tgt /= spl.norm(tgt)
            # tgt_next /= spl.norm(tgt_next[-1])

            # calculate beta and check against betamax if requested, fail convergence if check fails
            dot_product = np.dot(tgt_next, tgt)
            dot_product = np.clip(dot_product, -1.0, 1.0)
            beta = np.degrees(np.arccos(dot_product))
            # beta below found using tangent of first partition + T only
            # beta_first_partition = np.degrees(
            #     np.arccos(
            #         (tgt_next[np.r_[0:twoN, -1]].T @ tgt[np.r_[0:twoN, -1]]) /
            #         (spl.norm(tgt[np.r_[0:twoN, -1]]) * spl.norm(tgt_next[np.r_[0:twoN, -1]]))
            #     )
            # )

            if cont_params_cont["betacontrol"] and beta > cont_params_cont["betamax"]:
                print("Beta exceeds maximum angle.")
                cvg_cont = False
            else:
                # Adjust tangent vector and stepsign to maintain correct direction
                if dot_product < 0:
                    tgt_next = -tgt_next
                    stepsign = -stepsign

                self.log.screenout(
                    iter=itercont,
                    correct=itercorrect,
                    res=residual,
                    freq=omega / tau_pred,
                    energy=energy,
                    step=stepsign * step,
                    beta=beta,
                )
                self.log.store(
                    sol_pose=pose,
                    sol_vel=vel,
                    sol_T=tau_pred / omega,
                    sol_tgt=tgt_next,
                    sol_energy=energy,
                    sol_beta=beta,
                    sol_itercorrect=itercorrect,
                    sol_step=stepsign * step,
                )
                
                tau = tau_pred
                X = X_pred.copy()
                tgt = tgt_next.copy()
                # update pose_base and set inc to zero, pose will have included inc from current sol
                pose_base = pose.copy()
                X[inc_mask] = 0.0
                itercont += 1

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
        self.log.screenline("-")


def normalise_residual(residual, pose_base, pose_ref, dofdata):
    ndof_all = dofdata["ndof_all"]
    n_nodes = dofdata["nnodes_all"]
    config_per_node = dofdata["config_per_node"]
    dof_per_node = dofdata["dof_per_node"]
    n_dim = dofdata["n_dim"]
    SEbeam = dofdata["SEbeam"]
    inc_from_ref = np.zeros((ndof_all))

    if SEbeam:
        for k in range(n_nodes):
            f = Frame.relative_frame(
                n_dim,
                pose_ref[k * config_per_node : (k + 1) * config_per_node],
                pose_base[k * config_per_node : (k + 1) * config_per_node],
            )
            inc_from_ref[k * dof_per_node : (k + 1) * dof_per_node] = (
                Frame.get_parameters_from_frame(n_dim, f)
            )
    else:
        inc_from_ref = pose_base - pose_ref
    return residual / spl.norm(inc_from_ref)


# def qrlinearsolver(A, b):
#     Q, R, P = spl.qr(A, pivoting=True, mode="economic")
#     Qt_b = np.dot(Q.T, b)
#     x_temp = spl.solve_triangular(R, Qt_b[:R.shape[0]])
#     x = np.zeros_like(x_temp)
#     x[P] = x_temp
#     return x
