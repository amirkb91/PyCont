import h5py
import json
from alive_progress import alive_bar
import sys, shutil
import numpy as np
from beamcpp import BeamCpp
from postprocess.bifurcation import bifurcation_functions
from Frame import Frame

""" Run time simulations for all NNM solutions and store """

new_nsteps = input("Specify new nsteps per period if desired: ")
run_bif = input("Compute bifurcation functions? ")

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
[_, _, pose_ref] = BeamCpp.run_eig()  # To get nodal data in class
n_solpoints = len(T)
if new_nsteps:
    par["shooting"]["single"]["nsteps_per_period"] = int(new_nsteps)
nsteps = par["shooting"]["single"]["nsteps_per_period"]
pose_time = np.zeros([BeamCpp.ndof_config, nsteps + 1, n_solpoints])
vel_time = np.zeros([BeamCpp.ndof_all, nsteps + 1, n_solpoints])
inc_time = np.zeros([BeamCpp.ndof_all, nsteps + 1, n_solpoints])
time = np.zeros([n_solpoints, nsteps + 1])
if run_bif:
    print("\033[0;31m*** ENSURE apply_SE_correction = FALSE *** \033[0m\n")
    Floquet = np.zeros([2 * BeamCpp.ndof_free, n_solpoints], dtype=np.complex128)
    Stability = np.zeros(n_solpoints)
    Fold = np.zeros(n_solpoints)
    Flip = np.zeros(n_solpoints)
    Neimark_Sacker = np.zeros(n_solpoints)

with alive_bar(n_solpoints) as bar:
    for i in range(n_solpoints):
        x = vel[BeamCpp.free_dof, i]
        X = np.concatenate([np.zeros(BeamCpp.ndof_free), x])
        if run_bif:
            [_, J, pose_time[:, :, i], vel_time[:, :, i], _, _, _] = BeamCpp.runsim_single(
                1.0, T[i], X, pose[:, i], par, fulltime=True
            )
            M = J[:, :-1] + np.eye(2 * BeamCpp.ndof_free)
            bifurcation_out = bifurcation_functions(M)
            Floquet[:, i] = bifurcation_out[0]
            Stability[i] = bifurcation_out[1]
            Fold[i] = bifurcation_out[2]
            Flip[i] = bifurcation_out[3]
            Neimark_Sacker[i] = bifurcation_out[4]
        else:
            [_, _, pose_time[:, :, i], vel_time[:, :, i], _, _, _] = BeamCpp.runsim_single(
                1.0, T[i], X, pose[:, i], par, sensitivity=False, fulltime=True
            )
        time[i, :] = np.linspace(0, T[i], nsteps + 1)

        for j in range(nsteps + 1):
            for k in range(BeamCpp.nnodes_all):
                f = Frame.relative_frame(
                    BeamCpp.n_dim,
                    pose_ref[k * BeamCpp.config_per_node : (k + 1) * BeamCpp.config_per_node],
                    pose_time[
                        k * BeamCpp.config_per_node : (k + 1) * BeamCpp.config_per_node, j, i
                    ],
                )
                inc_time[k * BeamCpp.dof_per_node : (k + 1) * BeamCpp.dof_per_node, j, i] = (
                    Frame.get_parameters_from_frame(BeamCpp.n_dim, f)
                )
        bar()

# write to file
time_data = h5py.File(new_file, "a")
if "/Config_Time/POSE" in time_data.keys():
    del time_data["/Config_Time/POSE"]
if "/Config_Time/VELOCITY" in time_data.keys():
    del time_data["/Config_Time/VELOCITY"]
if "/Config_Time/INC" in time_data.keys():
    del time_data["/Config_Time/INC"]
if "/Config_Time/Time" in time_data.keys():
    del time_data["/Config_Time/Time"]
time_data["/Config_Time/POSE"] = pose_time
time_data["/Config_Time/VELOCITY"] = vel_time
time_data["/Config_Time/INC"] = inc_time
time_data["/Config_Time/Time"] = time
if run_bif:
    if "/Bifurcation/Floquet" in time_data.keys():
        del time_data["/Bifurcation/Floquet"]
    if "/Bifurcation/Stability" in time_data.keys():
        del time_data["/Bifurcation/Stability"]
    if "/Bifurcation/Fold" in time_data.keys():
        del time_data["/Bifurcation/Fold"]
    if "/Bifurcation/Flip" in time_data.keys():
        del time_data["/Bifurcation/Flip"]
    if "/Bifurcation/Neimark_Sacker" in time_data.keys():
        del time_data["/Bifurcation/Neimark_Sacker"]
    time_data["/Bifurcation/Floquet"] = Floquet
    time_data["/Bifurcation/Stability"] = Stability
    time_data["/Bifurcation/Fold"] = Fold
    time_data["/Bifurcation/Flip"] = Flip
    time_data["/Bifurcation/Neimark_Sacker"] = Neimark_Sacker
time_data.close()
