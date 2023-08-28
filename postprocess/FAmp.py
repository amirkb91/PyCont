import matplotlib.pyplot as plt
import sys
import h5py
import json
import numpy as np
import mplcursors
from examples.beam_2D.beamcpp import BeamCpp

# show point data on figure
def show_annotation(sel):
    ind = int(sel.index)
    sel.annotation.set_text(f"index:{ind}")

dat_output = False
pose_ind2plot = 63  # 63 for beam mid
normalise = 0.01  # beam thickness

files = sys.argv[1:]
for i, file in enumerate(files):
    if not file.endswith(".h5"):
        files[i] += ".h5"

plt.style.use("ggplot")
f, a = plt.subplots(figsize=(10, 7))
a.set(xlabel="Freq/Omega", ylabel="Normalised Amplitude")

# eig info from first file, all files should therefore be from same NNM
file1 = files[0]
data = h5py.File(str(file1), "r")
par = data["/Parameters"]
par = json.loads(par[()])
try:
    nnm = par["first_point"]["eig_start"]["NNM"]
except:
    nnm = par["first_point"]["_eig_start"]["NNM"]
# get pose0 and eig
[eig, frq, pose0] = BeamCpp.run_eig()
Omega = frq[nnm - 1, 0]
data.close()

# plot sols
line = []
for file in files:
    data = h5py.File(str(file), "r")
    pose_time = data["/Config_Time/POSE"][:]
    T = data["/T"][:]

    n_solpoints = len(T)
    amp = np.zeros(n_solpoints)
    for i in range(n_solpoints):
        amp[i] = np.max(np.abs(pose_time[pose_ind2plot, :, i] - pose0[pose_ind2plot])) / normalise

    line.append(
        a.plot(1 / (T * Omega), amp, marker=".", fillstyle="none", label=file.split(".h5")[0])
    )

    if dat_output:
        F_omega = 1 / (T * Omega)
        np.savetxt(
            file.strip(".h5") + ".dat",
            np.concatenate([F_omega.reshape(-1, 1), amp.reshape(-1, 1)], axis=1)
        )

a.legend()

cursor = mplcursors.cursor(line[0], hover=False)
cursor.connect("add", show_annotation)
plt.draw()
plt.show()
