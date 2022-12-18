import numpy as np


class ConX:
    # rcond for least squares SVD solver: ratio of largest SV
    svd_rcond = None

    def __init__(self, prob, start, log):
        self.h = None
        self.nphase = None
        self.prob = prob
        self.X0 = start.X0
        self.T0 = start.T0
        self.pose_base0 = start.pose_base0
        self.energy0 = start.energy0
        self.tgt0 = start.tgt0
        self.log = log

    def solve(self):
        # calculate phase condition matrix h
        self.phase_condition()

        # compute first point of the branch
        self.first_point()

        if self.prob.cont_params["continuation"]["method"].lower() == "seq":
            # sequential continuation
            self.seqcont()
        elif self.prob.cont_params["continuation"]["method"].lower() == "psa":
            # pseudo-arc length continuation
            self.psacont()

    def phase_condition(self):
        # parse and sort phase condition indices. range defined in json file is inclusive
        h_idx = []
        idx = self.prob.cont_params["continuation"]["phase_cond_index"]
        if idx:
            idx = idx.split(",")
            for i in range(len(idx)):
                if "-" in idx[i]:
                    idxrange = idx[i].split("-")
                    h_idx.extend(list(range(int(idxrange[0]), int(idxrange[1]) + 1)))
                else:
                    h_idx.append(int(idx[i]))
            h_idx = sorted(set(h_idx))

        # create phase condition matrix h
        self.nphase = len(h_idx)
        self.h = np.zeros((self.nphase, len(self.X0)))
        self.h[list(range(self.nphase)), h_idx] = 1.0

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

                [H, M, dHdt, pose_time, vel_time, self.energy0, cvg_zerof] = self.prob.zerofunction(self.T0, self.X0,
                                                                                        self.prob.cont_params)
                if not cvg_zerof:
                    raise Exception("Zero function failed.")

                # residual = np.linalg.norm(H) / np.linalg.norm(self.X0)
                residual = np.linalg.norm(H)
                print(f"{iter_firstpoint} \t {residual:.5e}")

                if residual < self.prob.cont_params["continuation"]["tol"]:
                    print("First point converged.")
                    print("\n^-_-^-_-^-_-^-_-^-_-^-_-^-_-^-_-^-_-^\n")
                    break

                iter_firstpoint += 1
                if not restart:
                    # update X0, T0
                    # Jacobian comprised of monodromy, phase condition, and orthogonality to linear solution
                    J = np.block([
                        [M, dHdt.reshape(-1, 1)],
                        [self.h, np.zeros((self.nphase, 1))],
                        [self.X0.reshape(-1, 1).T, np.zeros(1)]])
                    # +1 zero for orthogonality to linear solution
                    hx = np.matmul(self.h, self.X0)
                    H = np.vstack([H, hx.reshape(-1, 1), np.zeros(1)])
                    dxt = np.linalg.lstsq(J, -H, rcond=self.svd_rcond)[0]
                    self.T0 += dxt[-1, 0]
                    dx = dxt[:-1, 0]
                    self.X0 += dx

                elif restart and fixF:
                    # update X0 - ortho to restart solution?
                    ortho = True
                    if not ortho:
                        J = np.concatenate((M, self.h), axis=0)
                        hx = np.matmul(self.h, self.X0)
                        H = np.vstack([H, hx.reshape(-1, 1)])
                    else:
                        J = np.concatenate((M, self.h, self.X0.reshape(-1, 1).T), axis=0)
                        hx = np.matmul(self.h, self.X0)
                        H = np.vstack([H, hx.reshape(-1, 1), np.zeros(1)])
                    dx = np.linalg.lstsq(J, -H, rcond=self.svd_rcond)[0]
                    self.X0[:] += dx[:, 0]

            # find tangent: set one component to 1 and solve overdetermined system
            J = np.block([
                [M, dHdt.reshape(-1, 1)],
                [self.h, np.zeros((self.nphase, 1))],
                [np.zeros([1, len(self.X0)]), np.ones(1)]])
            Z = np.vstack([np.zeros((len(self.X0), 1)), np.zeros((self.nphase, 1)), np.ones(1)])
            self.tgt0 = np.linalg.lstsq(J, Z, rcond=self.svd_rcond)[0][:, 0]
            self.tgt0 /= np.linalg.norm(self.tgt0)

        elif restart and not fixF:
            [H, M, dHdt, pose_time, vel_time, energy0, cvg_zerof] = self.prob.zerofunction(self.T0, self.X0,
                                                                                                self.prob.cont_params)
            [H, M, dHdt, pose, energy, cvg] = self.prob.run_sim(self.T0, self.X0)
            # residual = np.linalg.norm(H) / np.linalg.norm(self.X0)
            residual = np.linalg.norm(H)
            print(f"{1} \t {residual:.5e}")
            print("First point is restarted solution.")
            print("\n^-_-^-_-^-_-^-_-^-_-^-_-^-_-^-_-^-_-^\n")

        # store solution in logger
        self.log.store(sol_X=self.X0.copy(), sol_T=self.T0.copy(), sol_tgt=self.tgt0.copy(), sol_pose_time=pose_time.copy(),
                       sol_vel_time=vel_time.copy(), sol_pose_base=self.pose_base0.copy(), sol_energy=self.energy0.copy())

        # update pose_base0 for next step
        self.pose_base0 = pose_time[:, 0]

    def seqcont(self):
        print("Sequential continuation started.")

        # start from first point solution
        X = self.X0.copy()
        T = self.T0.copy()

        # default values for first iteration
        itercont = 1  # point 0 is first point which is found already
        step = self.prob.cont_params["continuation"]["s0"]
        direction = self.prob.cont_params["continuation"]["dir"]

        # continuation loop
        while True:
            print("\n**************************************\n")
            print(f"Continuation point {itercont}, freq = {1 / T[0]:.3f}:")
            print(f"step: s = {step:.3e}.")

            if itercont > self.prob.cont_params["continuation"]["npts"]:
                print("Maximum number of continuation points reached.")
                break

            # increment period. seq doesn't turn so direction is fixed throughout
            # store in case we need to change back if not converged
            _T = T.copy()
            T -= step * direction

            if 1 / T > self.prob.cont_params["continuation"]["fmax"]:
                print("Maximum frequency reached.")
                break

            # shoot to find X
            itershoot = 0
            while True:
                # find residual
                [H, Mm0, dHdt, outputs, zerof_cvg] = self.prob.run_sim(T, X)

                if not zerof_cvg:
                    cvg = False
                    print("Zero function failed to converge.")
                    break

                residual = np.linalg.norm(H) / np.linalg.norm(X)
                print(f"{itershoot} \t {residual:.5e}")

                if (
                        residual < self.prob.cont_params["continuation"]["tol"]
                        and itershoot >= self.prob.cont_params["continuation"]["itermin"]
                ):
                    cvg = True
                    print("Solution converged.")
                    break

                if itershoot >= self.prob.cont_params["continuation"]["itermax"]:
                    cvg = False
                    print("Max number of iterations reached without convergence.")
                    break

                # correction
                itershoot += 1
                J = np.concatenate((Mm0, self.h), axis=0)
                hx = np.matmul(self.h, X)
                H = np.vstack([H, hx.reshape(-1, 1)])
                dx = np.linalg.lstsq(J, -H, rcond=self.svd_rcond)[0]
                X[:] += dx[:, 0]

            # adaptive step size for next point
            if (
                    itercont >= self.prob.cont_params["continuation"]["nadapt"]
                    or not zerof_cvg
            ):
                step = self.cont_step(step, itershoot, cvg)

            if cvg:
                itercont += 1
                # store solution in logger
                self.log.store(solX=X.copy(), solT=T.copy(), out=outputs.copy())
            else:
                # revert as convergence failed
                T = _T.copy()

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
        if direction == 1:
            # decreasing T (increasing F)
            stepsign = -np.sign(tgt[-1])
        elif direction == -1:
            # increasing T (decreasing F)
            stepsign = np.sign(tgt[-1])

        # continuation loop
        itercont = 1
        while True:
            print("\n**************************************\n")
            print(f"Continuation point {itercont}")
            print(f"Freq = {1/T:.2f} -- Energy = {energy:.2f}")
            print(f"Step = {stepsign * step:.3e}")
            print("Iter \t Residual")
            if itercont > self.prob.cont_params["continuation"]["npts"]:
                print("Maximum number of continuation points reached.")
                break

            # update config
            self.prob.updatefunction(pose_base)

            # prediction step along tangent (n.b. INC is 0 as pose_base is updated)
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
                # calculate residual
                [H, M, dHdt, pose_time, vel_time, energy_next, cvg_zerof] = self.prob.zerofunction(T_pred, X_pred,
                                                                                        self.prob.cont_params)
                if not cvg_zerof:
                    cvg_cont = False
                    print("Zero function failed to converge.")
                    break

                # residual = np.linalg.norm(H) / np.linalg.norm(X_pred)
                residual = np.linalg.norm(H)
                print(f"{itercorrect} \t {residual:.5e}")

                if (residual < self.prob.cont_params["continuation"]["tol"]
                        and itercorrect >= self.prob.cont_params["continuation"]["itermin"]):
                    cvg_cont = True
                    print("Solution converged.")
                    break
                if itercorrect >= self.prob.cont_params["continuation"]["itermax"]:
                    cvg_cont = False
                    print("Max number of iterations reached without convergence.")
                    break

                # apply corrections orthogonal to tangent
                itercorrect += 1
                J = np.block([
                    [M, dHdt.reshape(-1, 1)],
                    [self.h, np.zeros((self.nphase, 1))],
                    [tgt]])
                hx = np.matmul(self.h, X_pred)
                H = np.vstack([H, hx.reshape(-1, 1), np.zeros(1)])
                dxt = np.linalg.lstsq(J, -H, rcond=self.svd_rcond)[0]
                T_pred += dxt[-1, 0]
                dx = dxt[:-1, 0]
                X_pred += dx

            if cvg_cont:
                # find new tangent with final converged solution
                if frml == "peeters":
                    J = np.block([
                        [M, dHdt.reshape(-1, 1)],
                        [self.h, np.zeros((self.nphase, 1))],
                        [np.zeros([1, len(X_pred)]), np.ones(1)]])
                elif frml == "keller":  # not same as correction above as derivatives changed before "break"
                    J = np.block([
                        [M, dHdt.reshape(-1, 1)],
                        [self.h, np.zeros((self.nphase, 1))],
                        [tgt]])
                Z = np.vstack([np.zeros((len(X_pred), 1)), np.zeros((self.nphase, 1)), np.ones(1)])
                tgt_next = np.linalg.lstsq(J, Z, rcond=self.svd_rcond)[0][:, 0]
                tgt_next /= np.linalg.norm(tgt_next)

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
                step = self.cont_step(step, itercorrect, cvg_cont)

    def cont_step(self, step, itercorrect, cvg):
        if cvg:
            if itercorrect == self.prob.cont_params["continuation"]["iteropt"] or itercorrect == 0:
                step *= np.sqrt(2)
            else:
                step *= self.prob.cont_params["continuation"]["iteropt"] / itercorrect

            if step > self.prob.cont_params["continuation"]["smax"]:
                step = self.prob.cont_params["continuation"]["smax"]
                print("Maximum step size reached.")
            elif step < self.prob.cont_params["continuation"]["smin"]:
                step = self.prob.cont_params["continuation"]["smin"]
                print("Minimum step size reached.")
        else:
            print("Continuation step reducing.")
            step /= np.sqrt(2)
            if step < self.prob.cont_params["continuation"]["smin"]:
                raise Exception("Step size below smin, continuation cannot proceed.")
        return step
