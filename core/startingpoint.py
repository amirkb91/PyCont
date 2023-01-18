import h5py
import numpy as np


class StartingPoint:
    def __init__(self, prob):
        self.prob = prob
        self.X0 = None
        self.T0 = None
        self.pose_base0 = None
        self.tgt0 = None
        self.energy0 = None

    def get_startingpoint(self):
        if self.prob.cont_params["first_point"]["restart"]["file_name"]:
            self.restart()
        else:
            self.new_start()

    def restart(self):
        # run icfunction to initialise config data
        [*_] = self.prob.icfunction(self.prob.cont_params)
        restartsol = h5py.File(self.prob.cont_params["first_point"]["restart"]["file_name"] + ".h5", "r+")
        index = self.prob.cont_params["first_point"]["restart"]["index"]
        self.X0 = restartsol["/X"][:, index]
        self.T0 = restartsol["/T"][index]
        self.pose_base0 = restartsol["/POSE_base"][:, index]
        self.tgt0 = restartsol["/Tangent"][:, index]
        self.energy0 = restartsol["/Energy"][index]

        # If different frequency is specified
        if self.prob.cont_params["first_point"]["restart"]["fixF"]:
            self.T0 = np.float64(1 / self.prob.cont_params["first_point"]["restart"]["F"])

    def new_start(self):
        # User supplied function provides initial guess
        self.X0, self.T0, self.pose_base0 = self.prob.icfunction(self.prob.cont_params)
