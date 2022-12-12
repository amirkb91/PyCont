import h5py
import numpy as np


class StartingPoint:
    def __init__(self, prob):
        self.prob = prob
        self.X0 = None
        self.T0 = None
        self.pose0 = None
        self.vel0 = None
        self.tgt0 = None

    def restart(self):
        restartsol = h5py.File(self.prob.cont_params["restart"]["file"] + ".h5", "r+")
        X = restartsol["/X"]
        T = restartsol["/T"]
        tgt0 = restartsol["/Tangent"]
        self.X0 = X[:, self.prob.cont_params["restart"]["index"]]
        self.T0 = np.array([T[self.prob.cont_params["restart"]["index"]]])
        self.tgt0 = tgt0[:, self.prob.cont_params["restart"]["index"]]

        # If different frequency is specified
        if self.prob.cont_params["restart"]["fixF"]:
            self.T0 = np.array([1 / self.prob.cont_params["restart"]["F"]])

    def new_start(self, fxn):
        # User supplied function provides initial guess
        self.X0, self.T0 = fxn(self.prob.cont_params)
