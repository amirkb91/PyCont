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
par["shooting"]["single"]["nperiod"] = nperiod
par["shooting"]["single"]["nsteps_per_period"] = nsteps
data.close()

# run sim
# specifying T in initialise decouples force period from sim time (we don't care about sens in timesim)
BeamCpp.initialise(par, T)
BeamCpp.run_eig()  # To get nodal data in class
x = vel[BeamCpp.free_dof]
X = np.concatenate([np.zeros(BeamCpp.ndof_free), x])
BeamCpp.runsim_single(1.0, T, X, pose, par, sensitivity=False)

# call plotbeam
subprocess.run("cd " + BeamCpp.cpp_path + "&&" + "python3 plotbeam.py ", shell=True)
