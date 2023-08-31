import h5py
import json
import subprocess
import sys
import numpy as np
from springcpp import SpringCpp

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
par["shooting"]["single"]["nperiod"] = nperiod
par["shooting"]["single"]["nsteps_per_period"] = nsteps
try:
    forced = par["continuation"]["forced"]
except:
    forced = False
data.close()

# run sim
SpringCpp.initialise(par)
SpringCpp.run_eig()  # To get nodal data in class
X = np.concatenate([pose, vel])
SpringCpp.runsim_single(1.0, T, X, np.array([0, 0]), par)

# call plotbeam
subprocess.run(
    "cd " + SpringCpp.cpp_path + "&&" + "python3 post_processing_simulation.py ", shell=True
)
