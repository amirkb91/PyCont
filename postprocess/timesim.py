import h5py
import json
import subprocess
import sys
from examples.beam_2D.beamcpp import BeamCpp

# inputs
solno = int(input("Solution Index: "))
nperiod = int(input("Number of periods: "))
nsteps = int(input("Steps per period: "))

# read solution file
file = sys.argv[-1]
data = h5py.File(str(file), "r")
pose_base = data["/POSE_base"][:, solno]
X = data["/X"][:, solno]
T = data["/T"][solno]

# read parameters from solution file and modify
par = data["/Parameters"]
par = json.loads(par[()])
par["shooting"]["nperiod"] = nperiod
par["shooting"]["nsteps_per_period"] = nsteps

# prescribe pose_base to ic file and run sim
BeamCpp.run_eig(par)  # To get nodal data
BeamCpp.config_update(pose_base)
BeamCpp.run_sim(T, X, par)

# call plotbeam
subprocess.run("cd " + BeamCpp.cpp_path + "&&" + "python3 plotbeam.py " + BeamCpp.simout_file + ".h5", shell=True)
