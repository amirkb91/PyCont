import matplotlib.pyplot as plt
import sys
import h5py
import json
import numpy as np
import mplcursors
from beamcpp import BeamCpp


# show point data on figure
def show_annotation(sel):
    ind = int(sel.index)
    sel.annotation.set_text(f"index:{ind}")


SE = True
dat_output = False

node2plot = 21
if SE:
    pose_ind2plot = 4 * node2plot + 3  # Y disp
    normalise_freq = 41.82280070074808
else:
    pose_ind2plot = 3 * node2plot + 1  # Y disp
    normalise_freq = 53.1660
normalise_amp = 1.0  # beam thickness

files = sys.argv[1:]
for i, file in enumerate(files):
    if not file.endswith(".h5"):
        files[i] += ".h5"

plt.style.use("ggplot")
f, a = plt.subplots(figsize=(10, 7))
a.set(xlabel="F/\u03C9\u2099", ylabel="Normalised Amplitude")

# plot sols
line = []
for file in files:
    data = h5py.File(str(file), "r")
    pose_time = data["/Config_Time/POSE"][:]
    T = data["/T"][:]
    par = data["/Parameters"]
    par = json.loads(par[()])
    forced = par["continuation"]["forced"]

    n_solpoints = len(T)
    amp = np.zeros(n_solpoints)
    for i in range(n_solpoints):
        amp[i] = np.max(np.abs(pose_time[pose_ind2plot, :, i])) / normalise_amp

    if forced:
        stability = data["/Bifurcation/Stability"][:]
        stable_index = np.argwhere(np.diff(stability)).squeeze() + 1
        a.plot(
            1 / (T[stable_index] * normalise_freq),
            amp[stable_index],
            marker='o',
            linestyle="none",
            markerfacecolor="none",
            markeredgecolor="k"
        )

    line.append(
        a.plot(
            1 / (T * normalise_freq),
            amp,
            marker="none",
            fillstyle="none",
            label=file.split(".h5")[0]
        )
    )

    if dat_output:
        F_omega = 1 / (T * normalise_freq)
        np.savetxt(
            file.strip(".h5") + ".dat",
            np.concatenate([F_omega.reshape(-1, 1), amp.reshape(-1, 1)], axis=1)
        )

a.legend()

cursor = mplcursors.cursor(line[0], hover=False)
cursor.connect("add", show_annotation)
plt.draw()
plt.show()
