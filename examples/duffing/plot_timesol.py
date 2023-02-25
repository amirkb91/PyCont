import h5py
import sys
import json
import numpy as np
import matplotlib.pyplot as plt

# inputs
solno = int(input("Solution Index: "))

# read solution file
file = sys.argv[-1]
if not file.endswith(".h5"):
    file += ".h5"
data = h5py.File(str(file), "r")
pos = data["/Config/POSE_time"][:, :, solno]
vel = data["/Config/VELOCITY_time"][:, :, solno]
T = data["/T"][solno]

# check parameters for time simulation data
par = data["/Parameters"]
par = json.loads(par[()])
nperiod = par["shooting"]["single"]["nperiod"]
nsteps = par["shooting"]["single"]["nsteps_per_period"]
t = np.linspace(0, T * nperiod, nsteps * nperiod)

# figure properties
f, (a1, a2, a3) = plt.subplots(1, 3, figsize=(10, 4))
f.suptitle(f"Frequency = {1 / T:.3f} Hz")
a1.set(xlabel="Time (s)", ylabel="Position (m)")
a1.plot(t, pos[0, :], '-', label="DoF 1")
a1.plot(t, pos[1, :], '--', label="DoF 2")
a1.legend()
a2.set(xlabel="Position DoF 1 (m)", ylabel="Position DoF 2 (m)")
a2.plot(pos[0, :], pos[1, :], '-')
a3.set(xlabel="Time (s)", ylabel="Velocity (m/s)")
a3.plot(t, vel[0, :], '-', label="DoF 1")
a3.plot(t, vel[1, :], '--', label="DoF 2")
a3.legend()

plt.tight_layout()
plt.show()
