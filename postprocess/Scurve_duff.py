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


normalise_force = 1.0
normalise_amp = 1.0

files = sys.argv[1:]
for i, file in enumerate(files):
    if not file.endswith(".h5"):
        files[i] += ".h5"

plt.style.use("ggplot")
f, a = plt.subplots(figsize=(10, 7))
a.set(xlabel="Forcing Amplitude", ylabel="Max Position Amplitude")

# plot sols
line = []
for file in files:
    data = h5py.File(str(file), "r")
    pose_time = data["/Config_Time/POSE"][:]
    F = data["/Force_Amp"][:]
    par = data["/Parameters"]
    par = json.loads(par[()])
    forced = par["continuation"]["forced"]

    n_solpoints = len(F)
    amp = np.zeros(n_solpoints)
    for i in range(n_solpoints):
        amp[i] = np.max(np.abs(pose_time[0, :, i])) / normalise_amp

    if forced and "/Bifurcation/Stability" in data.keys():
        stability = data["/Bifurcation/Stability"][:]
        stable_index = np.argwhere(np.diff(stability)).squeeze() + 1
        a.plot(
            F[stable_index] / normalise_force,
            amp[stable_index],
            marker="o",
            linestyle="none",
            markerfacecolor="none",
            markeredgecolor="k",
        )

    line.append(
        a.plot(
            F / normalise_force,
            amp,
            marker="none",
            fillstyle="none",
            label=file.split(".h5")[0],
        )
    )

    data.close()

a.legend()

cursor = mplcursors.cursor(line[0], hover=False)
cursor.connect("add", show_annotation)
plt.draw()
plt.show()
