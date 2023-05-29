import h5py
import numpy as np


class StartingPoint:

    def __init__(self, prob):
        self.prob = prob
        self.X0 = None
        self.T0 = None
        self.pose0 = None
        self.tgt0 = None

    def get_startingpoint(self):
        if self.prob.cont_params["first_point"]["restart"]["file_name"]:
            self.restart()
        else:
            self.new_start()

    def restart(self):
        restartsol = h5py.File(self.prob.cont_params["first_point"]["restart"]["file_name"] + ".h5", "r+")
        index = self.prob.cont_params["first_point"]["restart"]["index"]
        self.T0 = restartsol["/T"][index]
        self.pose0 = restartsol["/Config/POSE"][:, index]
        vel = restartsol["/Config/VELOCITY"][:, index]
        self.tgt0 = restartsol["/Tangent"][:, index]

        # run icfunction to get dof data
        [*_] = self.prob.icfunction(self.prob.cont_params)
        dofdata = self.prob.doffunction()
        x = vel[dofdata["free_dof"]]
        self.X0 = np.concatenate([np.zeros(dofdata["ndof_free"]), x])

        # # If different frequency is specified
        # if self.prob.cont_params["first_point"]["restart"]["fixF"]:
        #     self.T0 = np.float64(1 / self.prob.cont_params["first_point"]["restart"]["F"])

    def new_start(self):
        # User supplied function provides initial guess
        self.X0, self.T0, self.pose0 = self.prob.icfunction(self.prob.cont_params)
