import h5py
import json
from alive_progress import alive_bar
import sys, shutil, os
import numpy as np
from examples.beam_cpp.beamcpp import BeamCpp
""" Run time simulations for all NNM solutions and store """

new_nsteps = input("Specify new nsteps per period if desired: ")

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
BeamCpp.initialise(par)
BeamCpp.run_eig()  # To get nodal data in class
n_solpoints = len(T)
if new_nsteps:
    par["shooting"]["single"]["nsteps_per_period"] = int(new_nsteps)
nsteps = par["shooting"]["single"]["nsteps_per_period"]
pose_time = np.zeros([np.shape(pose)[0], nsteps + 1, n_solpoints])
vel_time = np.zeros([np.shape(vel)[0], nsteps + 1, n_solpoints])
time = np.zeros([n_solpoints, nsteps + 1])

with alive_bar(n_solpoints) as bar:
    for i in range(n_solpoints):
        x = vel[BeamCpp.free_dof, i]
        X = np.concatenate([np.zeros(BeamCpp.ndof_free), x])
        [pose_time[:, :, i], vel_time[:, :, i]] = BeamCpp.runsim_single(
            1.0, T[i], X, pose[:, i], par, return_time=True, sensoff=True
        )
        time[i, :] = np.linspace(0, T[i], nsteps + 1)
        bar()

# write to file
time_data = h5py.File(new_file, "a")
if "/Config_Time/POSE" in time_data.keys():
    del time_data["/Config_Time/POSE"]
if "/Config_Time/VELOCITY" in time_data.keys():
    del time_data["/Config_Time/VELOCITY"]
if "/Config_Time/Time" in time_data.keys():
    del time_data["/Config_Time/Time"]
time_data["/Config_Time/POSE"] = pose_time
time_data["/Config_Time/VELOCITY"] = vel_time
time_data["/Config_Time/Time"] = time
time_data.close()
