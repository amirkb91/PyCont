import numpy as np
import h5py
import os
import json
import matplotlib.pyplot as plt


class Logger:
    def __init__(self, prob):
        self.prob = prob
        self.store_index = 0
        self.sol_X = []
        self.sol_T = []
        self.sol_tgt = []
        self.sol_pose = []
        self.sol_energy = []
        self.sol_beta = []
        self.fig = None
        self.foldername = os.path.basename(os.getcwd())

    def store(self, **sol_data):
        self.store_index += 1
        for key, value in sol_data.items():
            if key == "sol_X":
                self.sol_X.append(value[:, 0])
            elif key == "sol_T":
                self.sol_T.append(value)
            elif key == "sol_tgt":
                self.sol_tgt.append(value)
            elif key == "sol_pose":
                self.sol_pose.append(value)
            elif key == "sol_energy":
                self.sol_energy.append(value)
            elif key == "sol_beta":
                self.sol_beta.append(value)

        # save to disk and plot if required
        self.savetodisk()
        if self.solidx > 1 and self.prob.cont_params["plot"]:
            self.solplot()

    def savetodisk(self):
        savefile = h5py.File("continuation_sol.h5", "w")
        savefile["/X"] = np.asarray(self.sol_X).T
        savefile["/T"] = np.asarray(self.sol_T).T
        savefile["/Energy"] = np.asarray(self.sol_energy).T
        savefile["/Tangent"] = np.asarray(self.sol_tgt).T
        savefile["/POSE"] = np.transpose(np.asarray(self.sol_pose), (1, 2, 0))
        savefile["/Parameters"] = json.dumps(self.prob.cont_params)
        savefile.close()

    def solplot(self):
        beta_xaxis = 10
        if not self.fig:
            if self.beta.all():
                self.fig, (self.ax1, self.ax2) = plt.subplots(1, 2, figsize=(8, 6))
            else:
                self.fig, (self.ax1) = plt.subplots(1, 1, figsize=(4, 6))
            # Frequency energy plot
            self.ax1.grid()
            self.ax1.set_xscale("log")
            self.ax1.set_xlabel("Energy (J)")
            self.ax1.set_ylabel("Frequency (Hz)")
            self.ax1.ticklabel_format(useOffset=False, axis="y")
            # self.ax1.set_xlim(1e-10, 1e4)
            # self.ax1.set_ylim(
            #     self.prob.cont_params["continuation"]["fmin"],
            #     self.prob.cont_params["continuation"]["fmax"],
            # )
            (self.line1,) = self.ax1.plot(
                self.outputs["energy"], 1 / self.solT, marker=".", fillstyle="none"
            )
            # beta plot
            if self.beta.all():
                self.ax2.grid()
                self.ax2.set_xlabel("Continuation Step")
                self.ax2.set_ylabel("beta (deg)")
                self.ax2.set_xlim(1, beta_xaxis)
                self.ax2.set_ylim(0, 180)
                (self.line2,) = self.ax2.plot(
                    range(len(self.beta)),
                    self.beta,
                    marker=".",
                    fillstyle="none",
                    color="red",
                )

            plt.tight_layout()
            plt.pause(0.5)

        else:
            self.line1.set_data(self.outputs["energy"], 1 / self.solT)
            self.ax1.relim()
            self.ax1.autoscale()
            if self.beta.all():
                self.line2.set_data(range(len(self.beta)), self.beta)
                self.ax2.set_xlim(1, beta_xaxis * np.ceil(len(self.beta) / beta_xaxis))
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()
