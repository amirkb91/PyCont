import numpy as np
import h5py
import os
import json
import matplotlib.pyplot as plt

# **** WSL: import matplotlib fails unless XServer is running on Windows ****


class Logger:
    def __init__(self, prob):
        self.prob = prob
        self.solidx = 0
        self.solX = None
        self.solT = None
        self.soltgt = None
        self.solpose = None
        self.outputs = None
        self.beta = np.zeros(1)
        self.fig = None
        self.foldername = os.path.basename(os.getcwd())

    def store(self, **soldata):
        for key, value in soldata.items():

            if key == "solX":
                if self.solidx == 0:
                    self.solX = value.reshape(-1, 1)
                else:
                    self.solX = np.append(self.solX, value.reshape(-1, 1), axis=1)

            elif key == "solT":
                if self.solidx == 0:
                    self.solT = value
                else:
                    self.solT = np.append(self.solT, value)

            elif key == "soltgt":
                if self.solidx == 0:
                    self.soltgt = value.reshape(-1, 1)
                else:
                    self.soltgt = np.append(self.soltgt, value.reshape(-1, 1), axis=1)

            elif key == "solpose":
                if self.solidx == 0:
                    self.solpose = value
                else:
                    self.solpose = np.dstack((self.solpose, value))

            elif key == "out":
                if self.solidx == 0:
                    self.outputs = value
                else:
                    for _key, _value in value.items():
                        self.outputs[_key] = np.append(self.outputs[_key], _value)

            elif key == "beta":
                if self.solidx == 1:
                    self.beta = value
                else:
                    self.beta = np.append(self.beta, value)

        self.solidx += 1

        # save to disk and plot if required
        self.savetodisk()
        if self.solidx > 1 and self.prob.cont_params["plot"]:
            self.solplot()

    def savetodisk(self):
        savefile = h5py.File(self.foldername + "_contsol.h5", "w")
        savefile["/X"] = self.solX
        savefile["/T"] = self.solT
        savefile["/Energy"] = self.outputs["energy"]
        savefile["/Tangent"] = self.soltgt
        savefile["/POSE"] = self.solpose
        savefile["/Parameters"] = json.dumps(self.prob.cont_params)
        savefile.close()
        # np.savetxt("solXT.out", np.vstack([self.solX, self.solT]))
        # np.savetxt("energy.out", self.outputs["energy"])
        # np.savetxt("tgt.out", self.soltgt)

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
