import numpy as np
import scipy.linalg as spl
from ._cont_step import cont_step


def psacont_mult(self):
    print("Pseudo-arc length continuation started.")
    frml = self.prob.cont_params["continuation"]["tangent"].lower()
    print(f"{frml.title()} tangent formulation.")
    if self.prob.cont_params["continuation"]["betacontrol"]:
        print("++ Beta control is active. ++")

    # dof data
    dofdata = self.prob.doffunction()
    N = dofdata["ndof_free"]
    twoN = 2 * dofdata["ndof_free"]

    # partition the starting orbit with Poincare sections
    npartition = self.prob.cont_params["shooting"]["npartition_multipleshooting"]
    nsteps = self.prob.cont_params["shooting"]["nsteps_per_period"]
    nsteps_per_partition = nsteps // npartition
    delta_S = 1 / npartition
    sol_index = (nsteps * delta_S * np.arange(npartition)).astype(int)
    pose_base = self.pose_time0.copy()[:, sol_index]
    V = self.vel_time0.copy()[dofdata["free_dof"]][:, sol_index]
    T = self.T0.copy()

    # find initial tangent matrix: set time component to 1 and solve overdetermined system
    J = np.zeros((npartition * (twoN + self.nphase) + 1, npartition * twoN + 1))
    for ipart in range(npartition):
        i = ipart * twoN
        i1 = (ipart + 1) * twoN
        j = (ipart + 1) % npartition * twoN
        j1 = ((ipart + 1) % npartition + 1) * twoN
        k = npartition * twoN + ipart * self.nphase
        k1 = npartition * twoN + (ipart + 1) * self.nphase

        self.prob.updatefunction(pose_base[:, ipart])
        X = np.concatenate((np.zeros(N), V[:, ipart]))
        [_, M, dHdt, _, _, _, _] = self.prob.zerofunction(T * delta_S, X, self.prob.cont_params, mult=True)
        J[i:i1, i:i1] = M
        J[i:i1, j:j1] -= np.eye(twoN)
        J[i:i1, -1] = dHdt * delta_S
        J[k:k1, i:i1] = self.h
    J[-1, -1] = 1
    Z = np.zeros((np.shape(J)[0], 1))
    Z[-1] = 1
    tgt = spl.lstsq(J, Z, cond=None, check_finite=False, lapack_driver="gelsy")[0][:, 0]
    tgt /= spl.norm(tgt)

    # continuation step and direction
    step = self.prob.cont_params["continuation"]["s0"]
    direction = self.prob.cont_params["continuation"]["dir"]
    stepsign = -1 * direction * np.sign(tgt[-1])  # corrections are always added

    # continuation loop
    energy = self.energy0.copy()
    itercont = 1
    while True:
        print("\n**************************************\n")
        print(f"Continuation point {itercont}")
        print(f"Freq = {1 / T:.2f} -- Energy = {energy:.2f}")
        print(f"Step = {stepsign * step:.3e}")
        print("Iter \t Residual")
        if itercont > self.prob.cont_params["continuation"]["npts"]:
            print("Maximum number of continuation points reached.")
            break

        # prediction step (INC is 0 as config is updated)
        T_pred = T + tgt[-1] * step * stepsign
        tgt_reshape = np.reshape(tgt[:-1], (-1, npartition), order='F')
        INC_pred = tgt_reshape[:N, :] * step * stepsign
        VEL_pred = V + tgt_reshape[N:, :] * step * stepsign
        X_pred = np.concatenate((INC_pred, VEL_pred))
        if 1 / T_pred > self.prob.cont_params["continuation"]["fmax"]:
            print("Maximum frequency reached.")
            break

        # correction step
        itercorrect = 0
        while True:
            # pose_time and vel_time initialised with size required to stitch from partitions
            J = np.zeros((npartition * (twoN + self.nphase) + 1, npartition * twoN + 1))
            cvg_zerof = [None] * npartition
            H = np.zeros((twoN, npartition))
            pose_time = np.zeros((np.shape(pose_base)[0], (nsteps_per_partition + 1) * npartition))
            vel_time = np.zeros((dofdata["ndof_all"], (nsteps_per_partition + 1) * npartition))
            energy_next = np.zeros(npartition)

            # residual and block Jacobian
            for ipart in range(npartition):
                i = ipart * twoN
                i1 = (ipart + 1) * twoN
                j = (ipart + 1) % npartition * twoN
                j1 = ((ipart + 1) % npartition + 1) * twoN
                k = npartition * twoN + ipart * self.nphase
                k1 = npartition * twoN + (ipart + 1) * self.nphase
                p = ipart * (nsteps_per_partition + 1)
                p1 = (ipart + 1) * (nsteps_per_partition + 1)

                self.prob.updatefunction(pose_base[:, ipart])
                # periodicity target: pose_base and Vel of next Poincare section along orbit
                target = np.concatenate((pose_base[dofdata["free_dof"]][:, (ipart + 1) % npartition],
                                         X_pred[N:, (ipart + 1) % npartition]))
                # calculate residual and block Jacobian
                [H[:, ipart], M, dHdt, pose_time[:, p:p1], vel_time[:, p:p1], energy_next[ipart], cvg_zerof[ipart]] = \
                    self.prob.zerofunction(T_pred * delta_S, X_pred[:, ipart], self.prob.cont_params, mult=True, target=target)
                J[i:i1, i:i1] = M
                J[i:i1, j:j1] -= np.eye(twoN)
                J[i:i1, -1] = dHdt * delta_S
                J[k:k1, i:i1] = self.h

            J[-1, :] = tgt
            H = np.reshape(H, (-1, 1), order='F')

            if not all(cvg_zerof):
                cvg_cont = False
                print("Zero function failed to converge.")
                break

            residual = spl.norm(H)
            print(f"{itercorrect} \t {residual:.5e}")

            if (residual < self.prob.cont_params["continuation"]["tol"]
                    and itercorrect >= self.prob.cont_params["continuation"]["itermin"]):
                cvg_cont = True
                print("Solution converged.")
                break
            elif itercorrect >= self.prob.cont_params["continuation"]["itermax"]:
                cvg_cont = False
                print("Max number of iterations reached without convergence.")
                break

            # apply corrections orthogonal to tangent
            itercorrect += 1
            hx = np.matmul(self.h, X_pred)
            H = np.vstack([H, hx.reshape(-1, 1, order='F'), np.zeros(1)])
            dxt = spl.lstsq(J, -H, cond=None, check_finite=False, lapack_driver="gelsy")[0]
            T_pred += dxt[-1, 0]
            dx = np.reshape(dxt[:-1], (-1, npartition), order='F')
            X_pred += dx

        if cvg_cont:
            # find new tangent with converged solution
            # peeters Jacobian is different for tangent update
            if frml == "peeters":
                J[-1, :] = np.concatenate([np.zeros(npartition * twoN), np.ones(1)])
            Z = np.zeros((np.shape(J)[0], 1))
            Z[-1] = 1
            tgt_next = spl.lstsq(J, Z, cond=None, check_finite=False, lapack_driver="gelsy")[0][:, 0]
            tgt_next /= spl.norm(tgt_next)

            # calculate beta and check against betamax if requested, fail convergence if check fails
            beta = np.array([np.degrees(np.arccos(tgt_next.T @ tgt))])
            print(f"Beta = {beta[0]:.2f} deg")
            if (self.prob.cont_params["continuation"]["betacontrol"]
                    and beta[0] > self.prob.cont_params["continuation"]["betamax"]):
                print("Beta exceeds maximum angle, roll back and reduce continuation step.")
                cvg_cont = False
            else:
                # passed check, finalise and update for next step
                itercont += 1
                if frml == "peeters":
                    stepsign = np.sign(stepsign * tgt_next.T @ tgt)
                T = T_pred
                V = X_pred[N:]
                tgt = tgt_next
                energy = energy_next

                # store solution in logger
                self.log.store(sol_X=X_pred.copy(), sol_T=T_pred.copy(), sol_tgt=tgt_next.copy(),
                               sol_pose_time=pose_time.copy(), sol_vel_time=vel_time.copy(),
                               sol_pose_base=pose_base.copy(), sol_energy=energy_next.copy(), sol_beta=beta.copy())

                # pose_base for next step
                # points of poincare sections are repeated when stitching, modify sol_index
                pose_base = pose_time[:, sol_index + np.arange(npartition)]

        # adaptive step size for next point
        if itercont > self.prob.cont_params["continuation"]["nadapt"] or not cvg_cont:
            step = cont_step(self, step, itercorrect, cvg_cont)
