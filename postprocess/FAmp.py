import matplotlib.pyplot as plt
import sys
import h5py
import json
import numpy as np
import mplcursors


# show point data on figure
def show_annotation(sel):
    ind = int(sel.index)
    sel.annotation.set_text(f"index:{ind}")


node2plot = 30
normalise_freq = 2.62557
normalise_amp = 1.0
SE = True
dim = 3

if SE:
    iconfig = 6  # 2D: XY=2,3   ---   3D: XYZ=4,5,6
    pose_ind2plot = (3 * dim - 2) * node2plot + iconfig
else:
    iconfig = 6  # 2D: XY=0,1   ---   3D: XYZ=0,1,2
    pose_ind2plot = 3 * (dim - 1) * node2plot + iconfig

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

    if forced and "/Bifurcation/Stability" in data.keys():
        stability = data["/Bifurcation/Stability"][:]
        stable_index = np.argwhere(np.diff(stability)).squeeze() + 1
        a.plot(
            1 / (T[stable_index] * normalise_freq),
            amp[stable_index],
            marker="o",
            linestyle="none",
            markerfacecolor="none",
            markeredgecolor="k",
        )

    line.append(
        a.plot(
            1 / (T * normalise_freq),
            amp,
            marker="none",
            fillstyle="none",
            label=file.split(".h5")[0],
        )
    )

a.legend()

cursor = mplcursors.cursor(line[0], hover=False)
cursor.connect("add", show_annotation)
plt.draw()
plt.show()
