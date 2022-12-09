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
        idx = self.prob.cont_params["phasecond"]["idx"]
        if idx:
            idx = idx.split(",")
            for i in range(len(idx)):
                if "-" in idx[i]:
                    idxrange = idx[i].split("-")
                    h_idx.extend(
                        list(range(int(idxrange[0]), int(idxrange[1]) + 1))
                    )
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
        restart = self.prob.cont_params["restart"]["file"]
        fixF = self.prob.cont_params["restart"]["fixF"]

        if not restart or (restart and fixF):
            iter = 0
            while True:
                iter += 1
                if iter > self.prob.cont_params["firstpoint"]["itermax"]:
                    raise Exception(
                        "Max number of iterations reached without convergence."
                    )

                [H, Mm0, dHdt, pose, outputs, zerof_cvg] = self.prob.zerofunction(self.T0, self.X0, self.prob.cont_params)
                if not zerof_cvg:
                    raise Exception("Zero function failed.")

                residual = np.linalg.norm(H) / np.linalg.norm(self.X0)
                print(f"{iter} \t {residual:.5e}")

                if residual < self.prob.cont_params["firstpoint"]["tol"]:
                    print("First point converged.")
                    print("\n^-_-^-_-^-_-^-_-^-_-^-_-^-_-^-_-^-_-^\n")
                    break

                if not restart:
                    # update X0, T0
                    # Jacobian comprised of monodromy, phase condition, and orthogonality to linear solution
                    J = np.block(
                        [
                            [Mm0, dHdt.reshape(-1, 1)],
                            [self.h, np.zeros((self.nphase, 1))],
                            [self.X0.reshape(-1, 1).T, np.zeros(1)],
                        ]
                    )
                    # +1 zero for orthogonality to linear solution
                    hx = np.matmul(self.h, self.X0)
                    H = np.vstack([H, hx.reshape(-1, 1), np.zeros(1)])
                    dxt = np.linalg.lstsq(J, -H, rcond=self.svd_cutoff)[0]
                    self.T0 += dxt[-1]
                    dx = dxt[:-1, :]
                    self.X0[:] += dx[:, 0]

                elif restart and fixF:
                    # update X0 - ortho to restart solution?
                    ortho = True
                    if not ortho:
                        J = np.concatenate((Mm0, self.h), axis=0)
                        hx = np.matmul(self.h, self.X0)
                        H = np.vstack([H, hx.reshape(-1, 1)])
                    else:
                        J = np.concatenate(
                            (Mm0, self.h, self.X0.reshape(-1, 1).T), axis=0
                        )
                        hx = np.matmul(self.h, self.X0)
                        H = np.vstack([H, hx.reshape(-1, 1), np.zeros(1)])
                    dx = np.linalg.lstsq(J, -H, rcond=self.svd_cutoff)[0]
                    self.X0[:] += dx[:, 0]

            # find tangent if converged
            # set one component to 1 and solve overdetermined system using pseudo inverse
            J = np.block(
                [
                    [Mm0, dHdt.reshape(-1, 1)],
                    [self.h, np.zeros((self.nphase, 1))],
                    [np.zeros([1, len(self.X0)]), np.ones(1)],
                ]
            )
            Z = np.vstack(
                [np.zeros((len(self.X0), 1)), np.zeros((self.nphase, 1)), np.ones(1)]
            )
            self.tgt0 = np.linalg.lstsq(J, Z, rcond=self.svd_cutoff)[0][:, 0]
            self.tgt0 /= np.linalg.norm(self.tgt0)

        elif restart and not fixF:
            [H, Mm0, dHdt, pose, outputs, zerof_cvg] = self.prob.run_sim(self.T0, self.X0)
            residual = np.linalg.norm(H) / np.linalg.norm(self.X0)
            print(f"{1} \t {residual:.5e}")
            print("First point is restarted solution.")
            print("\n^-_-^-_-^-_-^-_-^-_-^-_-^-_-^-_-^-_-^\n")

        # store solution in logger
        self.log.store(
            solX=self.X0.copy(),
            solT=self.T0.copy(),
            soltgt=self.tgt0.copy(),
            solpose=pose.copy(),
            out=outputs.copy(),
        )

    def seqcont(self):
        print("Sequential continuation started.")

        # start from first point solution
        X = self.X0.copy()
        T = self.T0.copy()

        # default values for first iteration
        itercont = 1  # point 0 is first point which is found already
        step = self.prob.cont_params["continuation"]["s0"]
        dir = self.prob.cont_params["continuation"]["dir"]

        # continuation loop
        while True:
            print("\n**************************************\n")
            print(f"Continuation point {itercont}, freq = {1/T[0]:.3f}:")
            print(f"step: s = {step:.3e}.")

            if itercont > self.prob.cont_params["continuation"]["npts"]:
                print("Maximum number of continuation points reached.")
                break

            # increment period. seq doesn't turn so direction is fixed throughout
            # store in case we need to change back if not converged
            _T = T.copy()
            T -= step * dir

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
                dx = np.linalg.lstsq(J, -H, rcond=self.svd_cutoff)[0]
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
        # formulation chosen
        frml = self.prob.cont_params["continuation"]["formula"].lower()
        print(f"{frml.title()} formulation.")
        if self.prob.cont_params["continuation"]["betacontrol"]:
            print("beta control is active.")

        # start from first point solution
        X = self.X0.copy()
        T = self.T0.copy()
        tgt = self.tgt0.copy()

        # default values for first iteration
        itercont = 1  # point 0 is first point which is found already
        step = self.prob.cont_params["continuation"]["s0"]

        # continuation direction
        dir = self.prob.cont_params["continuation"]["dir"]
        if dir == 1:
            # decreasing T (increasing F)
            stepsign = -np.sign(tgt[-1])
        elif dir == -1:
            # increasing T (decreasing F)
            stepsign = np.sign(tgt[-1])

        # continuation loop
        while True:
            print("\n**************************************\n")
            print(f"Continuation point {itercont}, freq = {1/T[0]:.3f}:")
            print(f"Step: s = {step:.3e}. sign = {stepsign}.")

            if itercont > self.prob.cont_params["continuation"]["npts"]:
                print("Maximum number of continuation points reached.")
                break

            # prediction
            # store in case we need to change back if not converged
            _X = X.copy()
            _T = T.copy()
            X += tgt[:-1] * step * stepsign
            T += tgt[-1] * step * stepsign

            if 1 / T > self.prob.cont_params["continuation"]["fmax"]:
                print("Maximum frequency reached.")
                break

            # correction
            itercorrect = 0
            while True:
                # find residual
                [H, Mm0, dHdt, pose, outputs, cvg] = self.prob.run_sim(T, X)

                if not cvg:
                    print("Zero function failed to converge.")
                    break

                residual = np.linalg.norm(H) / np.linalg.norm(X)
                residual_abs = np.linalg.norm(H)
                print(f"{itercorrect} \t {residual:.5e} \t {residual_abs:.5e}")
                residual = residual_abs
                if (
                    residual < self.prob.cont_params["continuation"]["tol"]
                    and itercorrect >= self.prob.cont_params["continuation"]["itermin"]
                ):
                    cvg = True
                    print("Solution converged.")
                    break

                if itercorrect >= self.prob.cont_params["continuation"]["itermax"]:
                    cvg = False
                    print("Max number of iterations reached without convergence.")
                    break

                # correction
                itercorrect += 1
                J = np.block(
                    [
                        [Mm0, dHdt.reshape(-1, 1)],
                        [self.h, np.zeros((self.nphase, 1))],
                        [tgt],
                    ]
                )
                hx = np.matmul(self.h, X)

                # we can't subtract the X in the Lie group setting. So for now just use the
                # peeters correction update
                # if frml == "keller":
                #     # all corrected solution lie on plane orthogonal to tangent
                #     N = (X - _X).T @ tgt[:-1] + (T - _T)[0] * tgt[-1] - step * stepsign
                #     H = np.vstack([H, hx.reshape(-1, 1), N])
                # elif frml == "peeters":
                #     # all delta corrections are orthogonal to tangent
                #     H = np.vstack([H, hx.reshape(-1, 1), np.zeros(1)])
                H = np.vstack([H, hx.reshape(-1, 1), np.zeros(1)])

                dxt = np.linalg.lstsq(J, -H, rcond=self.svd_cutoff)[0]
                T += dxt[-1]
                dx = dxt[:-1, :]
                X[:] += dx[:, 0]

            if cvg:
                # First, compare beta against betamax (if user requested)
                # proceed to next step if passed betacheck, otherwise cvg=False

                # update tangent using converged solution
                tgt_prev = tgt.copy()  # previous tangent
                if frml == "keller":
                    J = np.block(
                        [
                            [Mm0, dHdt.reshape(-1, 1)],
                            [self.h, np.zeros((self.nphase, 1))],
                            [tgt_prev],
                        ]
                    )
                elif frml == "peeters":
                    J = np.block(
                        [
                            [Mm0, dHdt.reshape(-1, 1)],
                            [self.h, np.zeros((self.nphase, 1))],
                            [np.zeros([1, len(self.X0)]), np.ones(1)],
                        ]
                    )
                Z = np.vstack(
                    [
                        np.zeros((len(self.X0), 1)),
                        np.zeros((self.nphase, 1)),
                        np.ones(1),
                    ]
                )
                tgt = np.linalg.lstsq(J, Z, rcond=self.svd_cutoff)[0][:, 0]
                tgt /= np.linalg.norm(tgt)

                # angle between tangents
                beta = np.array([np.degrees(np.arccos(tgt.T @ tgt_prev))])
                print(f"Beta = {beta[0]:.2f} deg")
                if (
                    self.prob.cont_params["continuation"]["betacontrol"]
                    and beta[0] > self.prob.cont_params["continuation"]["betamax"]
                ):
                    print(
                        "Beta exceeds maximum angle, roll back and reduce continuation step."
                    )
                    cvg = False
                    # revert the tangent
                    tgt = tgt_prev.copy()
                    pass

            if cvg:
                # only if cvg and also passed the beta check
                itercont += 1

                # obtain sign of next step
                if frml == "peeters":
                    stepsign = np.sign(stepsign * tgt.T @ tgt_prev)

                # store solution in logger
                self.log.store(
                    solX=X.copy(),
                    solT=T.copy(),
                    soltgt=tgt.copy(),
                    solpose=pose.copy(),
                    out=outputs.copy(),
                    beta=beta.copy(),
                )
            else:
                # revert as convergence failed
                X = _X.copy()
                T = _T.copy()

            # adaptive step size for next point
            itercontcheck = itercont > self.prob.cont_params["continuation"]["nadapt"]
            if itercontcheck or not cvg:
                step = self.cont_step(step, itercorrect, cvg)

    def cont_step(self, step, itercorrect, cvg):
        if cvg:
            step *= self.prob.cont_params["continuation"]["iteropt"] / itercorrect

            if step > self.prob.cont_params["continuation"]["smax"]:
                step = self.prob.cont_params["continuation"]["smax"]
                print("Maximum step size reached.")

            elif step < self.prob.cont_params["continuation"]["smin"]:
                step = self.prob.cont_params["continuation"]["smin"]
                print("Minimum step size reached.")

        else:
            print("Continuation step halving.")
            step /= 2
            if step < self.prob.cont_params["continuation"]["smin"]:
                raise Exception("Step size below smin, continuation cannot proceed.")
        return step
