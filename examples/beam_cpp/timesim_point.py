import h5py
import json
import subprocess
import sys
import numpy as np
from beamcpp import BeamCpp

# inputs
solno = int(input("Solution Index: "))

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
    nsteps_par = par["shooting"]["single"]["nsteps_per_period"]
    nperiod = int(input("Number of periods: ") or 1)
    nsteps = int(input("Steps per period: ") or nsteps_par)
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
    nsteps = par["shooting"]["multiple"]["nsteps_per_partition"]
    delta_S = 1 / npartition

    pose = pose.reshape(BeamCpp.ndof_config, npartition, order="F")
    vel = vel.reshape(BeamCpp.ndof_all, npartition, order="F")
    time = np.array(
        [
            np.linspace(ipart * T * delta_S, (ipart + 1) * T * delta_S, nsteps + 1)
            for ipart in range(npartition)
        ]
    ).T

    x = vel[BeamCpp.free_dof, :]
    X = np.concatenate([np.zeros((BeamCpp.ndof_free, npartition)), x])
    X = X.flatten(order="F")
    [_, _, pose_time, vel_time, acc_time, energy, _] = BeamCpp.runsim_multiple(
        1.0, T, X, pose, par, fulltime=True
    )

    output_file = file.strip(".h5") + "_with_pointtime.h5"
    with h5py.File(output_file, "w") as f:
        f.create_dataset("/Config_Time/POSE", data=pose_time)
        f.create_dataset("/Config_Time/VELOCITY", data=vel_time)
        f.create_dataset("/Config_Time/ACCELERATION", data=acc_time)
        f.create_dataset("/Config_Time/Time", data=time)
        f.create_dataset("/Config_Time/Energy", data=energy)

    # copy solution over to C++ directory for plotting
    subprocess.run("cp " + output_file + " " + BeamCpp.cpp_path + "/beam_sim_mult.h5", shell=True)

    subprocess.run("cd " + BeamCpp.cpp_path + "&&" + "python3 plotbeam_mult.py ", shell=True)

    # # one more time sime with initial conditions of first partition for whole orbit so we can plot
    # x = vel[BeamCpp.free_dof, 0]
    # X = np.concatenate([np.zeros(BeamCpp.ndof_free), x])
    # BeamCpp.runsim_single(1.0, T, X, pose[:, 0], par, sensitivity=False)
    # # call plotbeam
    # subprocess.run("cd " + BeamCpp.cpp_path + "&&" + "python3 plotbeam.py ", shell=True)
