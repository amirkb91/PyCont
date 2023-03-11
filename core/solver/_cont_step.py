import numpy as np


def cont_step(self, step, itercorrect, cvg):
    if cvg:
        if itercorrect == self.prob.cont_params["continuation"]["iteropt"] or itercorrect == 0:
            step *= np.sqrt(2)
        else:
            step *= self.prob.cont_params["continuation"]["iteropt"] / itercorrect
        step = max(step, self.prob.cont_params["continuation"]["smin"])
        step = min(step, self.prob.cont_params["continuation"]["smax"])
    else:
        print("Continuation step reducing.")
        step /= 2
        if step < self.prob.cont_params["continuation"]["smin"]:
            raise Exception("Step size below smin, continuation cannot proceed.")
    return step
