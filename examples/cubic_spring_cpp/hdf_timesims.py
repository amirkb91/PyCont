import h5py
import json
from alive_progress import alive_bar
import sys, shutil, os
import numpy as np
from springcpp import SpringCpp
""" Run time simulations for all NNM solutions and store """

# read solution file
file = sys.argv[1]
inplace = sys.argv[-1]
if inplace == "-i":
    inplace = True
else:
    inplace = False

if not file.endswith(".h5"):
    file += ".h5"
data = h5py.File(str(file), "r")
pose = data["/Config/POSE"][:]
vel = data["/Config/VELOCITY"][:]
T = data["/T"][:]
par = data["/Parameters"]
par = json.loads(par[()])
data.close()

# create new file to store time histories or append inplace
if inplace:
    new_file = file
else:
    new_file = file.strip(".h5") + "_withtime.h5"
    shutil.copy(file, new_file)

# run sims
SpringCpp.initialise(par)
SpringCpp.run_eig()  # To get nodal data in class
n_solpoints = len(T)
nsteps = par["shooting"]["single"]["nsteps_per_period"]
pose_time = np.zeros([np.shape(pose)[0], nsteps + 1, n_solpoints])
vel_time = np.zeros([np.shape(vel)[0], nsteps + 1, n_solpoints])
time = np.zeros([n_solpoints, nsteps + 1])

with alive_bar(n_solpoints) as bar:
    for i in range(n_solpoints):
        X = np.concatenate([pose[:, i], vel[:, i]])
        [pose_time[:, :, i], vel_time[:, :, i]] = SpringCpp.runsim_single(
            1.0, T[i], X, np.array([0, 0]), par, return_time=True
        )
        time[i, :] = np.linspace(0, T[i], nsteps + 1)
        bar()

# write to file
time_data = h5py.File(new_file, "a")
time_data["/Config_Time/POSE"] = pose_time
time_data["/Config_Time/VELOCITY"] = vel_time
time_data["/Config_Time/Time"] = time
time_data.close()
