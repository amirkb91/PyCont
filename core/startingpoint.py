import h5py
import numpy as np


class StartingPoint:
    def __init__(self, prob):
        self.prob = prob
        self.X0 = None
        self.T0 = None
        self.tgt0 = None

    def restart(self):
        restartsol = h5py.File(self.prob.parameters["restart"]["file"] + ".h5", "r+")
        X = restartsol["/X"]
        T = restartsol["/T"]
        tgt0 = restartsol["/Tangent"]
        self.X0 = X[:, self.prob.parameters["restart"]["index"]]
        self.T0 = np.array([T[self.prob.parameters["restart"]["index"]]])
        self.tgt0 = tgt0[:, self.prob.parameters["restart"]["index"]]

        # If different frequency is specified
        if self.prob.parameters["restart"]["fixF"]:
            self.T0 = np.array([1 / self.prob.parameters["restart"]["F"]])

    def initial_guess(self, guessfunction):
        # User supplied function provides initial guess
        self.X0, self.T0 = guessfunction()
