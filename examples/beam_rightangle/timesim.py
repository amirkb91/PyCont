import h5py
import json
import subprocess
import sys
from examples.beam_rightangle.beamcpp import BeamCpp

# inputs
solno = int(input("Solution Index: "))
nperiod = int(input("Number of periods: "))
nsteps = int(input("Steps per period: "))

# read solution file
file = sys.argv[-1]
if not file.endswith(".h5"):
    file += ".h5"
data = h5py.File(str(file), "r")
pose_base = data["/Config/POSE_base"][:, :, solno]
X = data["/X"][:, solno]
T = data["/T"][solno]

# read parameters from solution file and modify
par = data["/Parameters"]
par = json.loads(par[()])
par["shooting"]["single"]["nperiod"] = nperiod
par["shooting"]["single"]["nsteps_per_period"] = nsteps

# run sim
BeamCpp.run_eig(par)  # To get nodal data in class
BeamCpp.runsim_single(T, X, pose_base, par)

# call plotbeam
subprocess.run("cd " + BeamCpp.cpp_path + "&&" + "python3 plotbeam.py " + BeamCpp.simout_file + ".h5", shell=True)
