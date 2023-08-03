import matplotlib.pyplot as plt
import sys
import h5py
import json
import numpy as np
from examples.beam_2D.beamcpp import BeamCpp

plt.style.use("ggplot")

# read solution file
pose_ind = int(input("POSE index to plot: "))  # 63 for beam mid

file = sys.argv[-1]
if not file.endswith(".h5"):
    file += ".h5"
data = h5py.File(str(file), "r")
pose_time = data["/Config_Time/POSE"][:]
T = data["/T"][:]
par = data["/Parameters"]
par = json.loads(par[()])
nnm = par["first_point"]["eig_start"]["NNM"]

# get pose0 and eig
[eig, frq, pose0] = BeamCpp.run_eig()
Omega = frq[nnm - 1, 0]

# calculate amplitude
# Need to find INC between pose0 and pose, for now just take difference (2D)
h = 0.01  # beam thickness
n_solpoints = len(T)
amp = np.zeros(n_solpoints)
for i in range(n_solpoints):
    amp[i] = np.max(np.abs(pose_time[pose_ind, :, i] - pose0[pose_ind])) / h

# figure properties
f, a = plt.subplots(figsize=(10, 7))
a.set(xlabel="Freq/Omega", ylabel="Amplitude/h")
a.plot(1 / (T * Omega), amp, marker=".", fillstyle="none", color="darkmagenta")

plt.draw()
plt.show()
