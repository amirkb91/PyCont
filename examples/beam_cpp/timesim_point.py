import h5py
import json
import subprocess
import sys
import numpy as np
from examples.beam_cpp.beamcpp import BeamCpp

# inputs
solno = int(input("Solution Index: "))
nperiod = int(input("Number of periods: "))
nsteps = int(input("Steps per period: "))

# read solution file
file = sys.argv[-1]
if not file.endswith(".h5"):
    file += ".h5"
data = h5py.File(str(file), "r")
pose = data["/Config/POSE"][:, solno]
vel = data["/Config/VELOCITY"][:, solno]
T = data["/T"][solno]
par = data["/Parameters"]
par = json.loads(par[()])
method = par["shooting"]["method"]

# specifying T in initialise decouples force period from sim time (we don't care about sens in timesim)
BeamCpp.initialise(par, T)
BeamCpp.run_eig()  # To get nodal data in class

if method == "single":
    par["shooting"]["single"]["nperiod"] = nperiod
    par["shooting"]["single"]["nsteps_per_period"] = nsteps

    # run sim
    x = vel[BeamCpp.free_dof]
    X = np.concatenate([np.zeros(BeamCpp.ndof_free), x])
    BeamCpp.runsim_single(1.0, T, X, pose, par, sensitivity=False)

    # call plotbeam
    subprocess.run("cd " + BeamCpp.cpp_path + "&&" + "python3 plotbeam.py ", shell=True)

elif method == "multiple":
    npartition = par["shooting"]["multiple"]["npartition"]
    par["shooting"]["single"]["nsteps_per_period"] = nsteps

    delta_S = 1 / npartition
    pose = pose.reshape(BeamCpp.ndof_config, npartition, order="F")
    vel = vel.reshape(BeamCpp.ndof_all, npartition, order="F")
    pose_time = np.zeros([BeamCpp.ndof_config, nsteps + 1, npartition])
    vel_time = np.zeros([BeamCpp.ndof_all, nsteps + 1, npartition])
    time = np.zeros([nsteps + 1, npartition])

    for ipart in range(npartition):
        x = vel[BeamCpp.free_dof, ipart]
        X = np.concatenate([np.zeros(BeamCpp.ndof_free), x])
        [_, _, pose_time[:, :, ipart], vel_time[:, :, ipart], _, _] = BeamCpp.runsim_single(
            1.0, T * delta_S, X, pose[:, ipart], par, sensitivity=False, fulltime=True
        )

        partition_starttime = ipart * T * delta_S
        t_part = np.linspace(0, T * delta_S, nsteps + 1)
        time[:, ipart] = t_part + partition_starttime
    # write pose_time and vel_time to a new hdf file
    output_file = file.strip(".h5") + "_with_pointtime.h5"
    with h5py.File(output_file, "w") as f:
        f.create_dataset("/Config_Time/POSE", data=pose_time)
        f.create_dataset("/Config_Time/VELOCITY", data=vel_time)
        f.create_dataset("/Config_Time/Time", data=time)

    # one more time sime with initial conditions of first partition for whole orbit so we can plot
    x = vel[BeamCpp.free_dof, 0]
    X = np.concatenate([np.zeros(BeamCpp.ndof_free), x])
    BeamCpp.runsim_single(1.0, T, X, pose[:, 0], par, sensitivity=False)
    # call plotbeam
    subprocess.run("cd " + BeamCpp.cpp_path + "&&" + "python3 plotbeam.py ", shell=True)
