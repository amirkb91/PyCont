import numpy as np
import scipy.linalg as spl
from ._cont_step import cont_step


def psacont(self):
    print("Pseudo-arc length continuation started.")
    frml = self.prob.cont_params["continuation"]["tangent"].lower()
    print(f"{frml.title()} tangent formulation.")
    if self.prob.cont_params["continuation"]["betacontrol"]:
        print("++ Beta control is active. ++")

    # first point solution
    len_V = len(self.X0) // 2
    T = self.T0.copy()
    V = self.X0[len_V:].copy()
    tgt = self.tgt0.copy()
    pose_base = self.pose_base0.copy()
    energy = self.energy0.copy()

    # continuation step and direction
    step = self.prob.cont_params["continuation"]["s0"]
    direction = self.prob.cont_params["continuation"]["dir"]
    stepsign = -1 * direction * np.sign(tgt[-1])  # corrections are always added

    # continuation loop
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

        # update config
        self.prob.updatefunction(pose_base)

        # prediction step along tangent (INC is 0 as pose_base is updated)
        INC_pred = tgt[:len_V] * step * stepsign
        VEL_pred = V + tgt[len_V:-1] * step * stepsign
        T_pred = T + tgt[-1] * step * stepsign
        X_pred = np.concatenate((INC_pred, VEL_pred))
        if 1 / T_pred > self.prob.cont_params["continuation"]["fmax"]:
            print("Maximum frequency reached.")
            break

        # correction step
        itercorrect = 0
        while True:
            # calculate residual and block Jacobian
            [H, M, dHdt, pose_time, vel_time, energy_next, cvg_zerof] = \
                self.prob.zerofunction(T_pred, X_pred, self.prob.cont_params)
            M = M - np.eye(len(M))
            J = np.block([
                [M, dHdt.reshape(-1, 1)],
                [self.h, np.zeros((self.nphase, 1))],
                [tgt]])

            if not cvg_zerof:
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
            H = np.vstack([H, hx.reshape(-1, 1), np.zeros(1)])
            dxt = spl.lstsq(J, -H, cond=None, check_finite=False, lapack_driver="gelsy")[0]
            T_pred += dxt[-1, 0]
            dx = dxt[:-1, 0]
            X_pred += dx

        if cvg_cont:
            # find new tangent with converged solution
            # peeters Jacobian is different for tangent update
            if frml == "peeters":
                J = np.block([
                    [M, dHdt.reshape(-1, 1)],
                    [self.h, np.zeros((self.nphase, 1))],
                    [np.zeros([1, len(X_pred)]), np.ones(1)]])
            Z = np.vstack([np.zeros((len(X_pred), 1)), np.zeros((self.nphase, 1)), np.ones(1)])
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
                V = X_pred[len_V:]
                tgt = tgt_next
                energy = energy_next

                # store solution in logger
                self.log.store(sol_X=X_pred.copy(), sol_T=T_pred.copy(), sol_tgt=tgt_next.copy(),
                               sol_pose_time=pose_time.copy(), sol_vel_time=vel_time.copy(),
                               sol_pose_base=pose_base.copy(), sol_energy=energy_next.copy(), sol_beta=beta.copy())

                # pose_base for next step
                pose_base = pose_time[:, 0]

        # adaptive step size for next point
        if itercont > self.prob.cont_params["continuation"]["nadapt"] or not cvg_cont:
            step = cont_step(self, step, itercorrect, cvg_cont)
