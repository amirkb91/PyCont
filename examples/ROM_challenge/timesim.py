import h5py
import json
import subprocess
import sys
import numpy as np
from examples.ROM_challenge.romchallenge import ROMChallenge

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
ROMChallenge.initialise(par, T)
ROMChallenge.run_eig()  # To get nodal data in class
x = vel[ROMChallenge.free_dof]
X = np.concatenate([np.zeros(ROMChallenge.ndof_free), x])
ROMChallenge.runsim_single(1.0, T, X, pose, par, sensitivity=False)
