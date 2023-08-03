import h5py
import json
from alive_progress import alive_bar
import sys, shutil, os
import numpy as np
from examples.beam_2D.beamcpp import BeamCpp
''' Run time simulations for all NNM solutions and store '''

# read solution file
file = sys.argv[-1]
if not file.endswith(".h5"):
    file += ".h5"
data = h5py.File(str(file), "r")
pose = data["/Config/POSE"][:]
vel = data["/Config/VELOCITY"][:]
T = data["/T"][:]
par = data["/Parameters"]
par = json.loads(par[()])

# create new file to store time histories
new_file = file.strip(".h5") + "_withtime.h5"
if os.path.isfile(new_file):
    prompt = input("Time history file already exists, CTRL+C to stop, Enter to continue: \n")
shutil.copy(file, new_file)

# run sims
BeamCpp.run_eig(par)  # To get nodal data in class
n_solpoints = len(T)
nsteps = par["shooting"]["single"]["nsteps_per_period"]
pose_time = np.zeros([np.shape(pose)[0], nsteps + 1, n_solpoints])
vel_time = np.zeros([np.shape(vel)[0], nsteps + 1, n_solpoints])
time = np.zeros([n_solpoints, nsteps + 1])

with alive_bar(n_solpoints) as bar:
    for i in range(n_solpoints):
        x = vel[BeamCpp.free_dof, i]
        X = np.concatenate([np.zeros(BeamCpp.ndof_free), x])
        [pose_time[:, :, i], vel_time[:, :, i]] = BeamCpp.runsim_single(
            1.0, T[i], X, pose[:, i], par, return_time=True
        )
        time[i, :] = np.linspace(0, T[i], nsteps + 1)
        bar()

# write to file
time_data = h5py.File(new_file, "a")
time_data["/Config_Time/POSE"] = pose_time
time_data["/Config_Time/VELOCITY"] = vel_time
time_data["/Config_Time/Time"] = time
time_data.close()