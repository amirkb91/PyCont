import numpy as np
import h5py
import json
import matplotlib.pyplot as plt


class Logger:
    def __init__(self, prob):
        self.prob = prob
        self.store_index = 0
        self.sol_X = []
        self.sol_T = []
        self.sol_tgt = []
        self.sol_pose_time = []
        self.sol_vel_time = []
        self.sol_pose_base = []
        self.sol_energy = []
        self.sol_beta = []
        self.fig = None

    def store(self, **sol_data):
        self.store_index += 1
        for key, value in sol_data.items():
            if key == "sol_X":
                self.sol_X.append(value)
            elif key == "sol_T":
                self.sol_T.append(value)
            elif key == "sol_tgt":
                self.sol_tgt.append(value)
            elif key == "sol_pose_time":
                self.sol_pose_time.append(value)
            elif key == "sol_vel_time":
                self.sol_vel_time.append(value)
            elif key == "sol_pose_base":
                self.sol_pose_base.append(value)
            elif key == "sol_energy":
                self.sol_energy.append(value)
            elif key == "sol_beta":
                self.sol_beta.append(value)

        if self.store_index % self.prob.cont_params["Logger"]["save_frequency"] == 0:
            # save to disk and plot if required
            self.savetodisk()
            if self.prob.cont_params["continuation"]["plot"]:
                self.solplot()

    def savetodisk(self):
        savefile = h5py.File(self.prob.cont_params["Logger"]["file_name"] + ".h5", "w")
        savefile["/X"] = np.asarray(self.sol_X).T
        savefile["/T"] = np.asarray(self.sol_T).T
        savefile["/Tangent"] = np.asarray(self.sol_tgt).T
        savefile["/Config/POSE_time"] = np.transpose(np.asarray(self.sol_pose_time), (1, 2, 0))
        savefile["/Config/VELOCITY_time"] = np.transpose(np.asarray(self.sol_vel_time), (1, 2, 0))
        savefile["/POSE_base"] = np.asarray(self.sol_pose_base).T
        savefile["/Energy"] = np.asarray(self.sol_energy).T
        savefile["/beta"] = np.asarray(self.sol_beta).T
        savefile["/Parameters"] = json.dumps(self.prob.cont_params)
        savefile.close()

    def solplot(self):
        Energy = np.asarray(self.sol_energy)
        T = np.asarray(self.sol_T)
        beta = np.asarray(self.sol_beta)
        beta_xaxis = 10
        betaplot = False
        if self.prob.cont_params["continuation"]["method"] == "psa":
            betaplot = True

        if not self.fig:
            if betaplot:
                self.fig, (self.ax1, self.ax2) = plt.subplots(1, 2, figsize=(10, 7))
            else:
                self.fig, (self.ax1) = plt.subplots(1, 1, figsize=(6, 4.5))
            # Frequency energy plot
            self.ax1.grid()
            self.ax1.set_xscale("log")
            self.ax1.set_xlabel("Energy (J)")
            self.ax1.set_ylabel("Frequency (Hz)")
            self.ax1.ticklabel_format(useOffset=False, axis="y")
            # self.ax1.set_xlim(1e-6, 1e4)
            # self.ax1.set_ylim(self.prob.cont_params["continuation"]["fmin"],
            #                   self.prob.cont_params["continuation"]["fmax"])
            (self.line1,) = self.ax1.plot(Energy, 1 / T, marker=".", fillstyle="none")
            # beta plot
            if betaplot:
                self.ax2.grid()
                self.ax2.set_xlabel("Continuation Step")
                self.ax2.set_ylabel("beta (deg)")
                self.ax2.set_xlim(1, beta_xaxis)
                self.ax2.set_ylim(-5, 185)
                (self.line2,) = self.ax2.plot(range(1, len(beta) + 1), beta, marker=".", fillstyle="none", color="red")
            plt.pause(0.01)
        else:
            self.line1.set_data(Energy, 1 / T)
            self.ax1.relim()
            self.ax1.autoscale()
            if betaplot:
                self.line2.set_data(range(1, len(beta) + 1), beta)
                self.ax2.set_xlim(1, beta_xaxis * np.ceil(len(beta) / beta_xaxis))
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()
