import h5py
import numpy as np


class StartingPoint:

    def __init__(self, prob):
        self.prob = prob
        self.X0 = None
        self.omega = None
        self.tau = None
        self.pose0 = None
        self.tgt0 = None

    def get_startingpoint(self):
        if self.prob.cont_params["first_point"]["restart"]["file_name"]:
            self.restart()
        else:
            self.new_start()

    def restart(self):
        restartsol = h5py.File(
            self.prob.cont_params["first_point"]["restart"]["file_name"] + ".h5", "r+")
        index = self.prob.cont_params["first_point"]["restart"]["index"]
        T = restartsol["/T"][index]
        self.pose0 = restartsol["/Config/POSE"][:, index]
        vel = restartsol["/Config/VELOCITY"][:, index]
        self.tgt0 = restartsol["/Tangent"][:, index]

        # run icfunction to get dof data and perform scaling if required
        _, T0, _ = self.prob.icfunction(self.prob.cont_params)
        dofdata = self.prob.doffunction()
        N = dofdata["ndof_free"]
        v = vel[dofdata["free_dof"]]
        self.X0 = np.concatenate([np.zeros(N), v])

        if self.prob.cont_params["shooting"]["scaling"] == True:
            self.omega = 1 / T0
            self.tau = self.omega * T
            self.X0[N:] *= 1 / self.omega  # scale velocities from X to Xtilde
        else:
            self.omega = 1.0
            self.tau = T

        # # If different frequency is specified
        # if self.prob.cont_params["first_point"]["restart"]["fixF"]:
        #     self.T0 = np.float64(1 / self.prob.cont_params["first_point"]["restart"]["F"])

    def new_start(self):
        # User supplied function provides initial guess
        self.X0, T0, self.pose0 = self.prob.icfunction(self.prob.cont_params)
        if self.prob.cont_params["shooting"]["scaling"] == True:
            self.omega = 1 / T0
            self.tau = 1.0
        else:
            self.omega = 1.0
            self.tau = T0
