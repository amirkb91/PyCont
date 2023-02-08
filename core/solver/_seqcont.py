import numpy as np
import scipy.linalg as spl

## ************* NEEDS REVISION ************* ##
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

            residual = spl.norm(H) / spl.norm(X)
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
            dx = spl.lstsq(J, -H, cond=None, check_finite=False, lapack_driver="gelsd")[0]
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