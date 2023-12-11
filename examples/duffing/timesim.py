import h5py
import sys
import json
import numpy as np
from scipy.integrate import odeint
from duffing import Duffing
import matplotlib.pyplot as plt

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

# do time simulation
par = data["/Parameters"]
par = json.loads(par[()])
method = par["shooting"]["method"]
if method == "single":
    nperiod = par["shooting"]["single"]["nperiod"]
    nsteps = par["shooting"]["single"]["nsteps_per_period"]
    t = np.linspace(0, T * nperiod, nsteps * nperiod + 1)
    X = np.concatenate([pose, vel])
    timesol = np.array(odeint(Duffing.system_ode, X, t, rtol=1e-8, tfirst=True))
    pose_time = timesol[:, 0]
    vel_time = timesol[:, 1]
# elif method == "multiple":
#     npartition = par["shooting"]["multiple"]["npartition"]
#     nsteps = par["shooting"]["multiple"]["nsteps_per_partition"]
#     delta_S = 1 / npartition
#     pose = pose.reshape(2, npartition, order="F")
#     vel = vel.reshape(2, npartition, order="F")
#     print("Partition Timestamps (s):")
#     t = np.zeros([nsteps + 1, 1, npartition])
#     pose_time = np.zeros([nsteps + 1, 2, npartition])
#     vel_time = np.zeros([nsteps + 1, 2, npartition])
#     for ipart in range(npartition):
#         partition_starttime = ipart * T * delta_S
#         print(f"{partition_starttime:.3f}")
#         t_part = np.linspace(0, T * delta_S, nsteps + 1)
#         t[:, :, ipart] = t_part.reshape(-1, 1) + partition_starttime
#         X = np.concatenate([pose[:, ipart], vel[:, ipart]])
#         timesol = np.array(odeint(Duffing.system_ode, X, t_part, rtol=1e-8, tfirst=True))
#         pose_time[:, :, ipart] = timesol[:, :2]
#         vel_time[:, :, ipart] = timesol[:, 2:]

#     t = t.flatten(order="F")
#     pose_time = np.reshape(np.transpose(pose_time, (2, 0, 1)), (-1, 2))
#     vel_time = np.reshape(np.transpose(vel_time, (2, 0, 1)), (-1, 2))

# plot
f, (a1, a2, a3) = plt.subplots(1, 3, figsize=(10, 4))
f.suptitle(f"Frequency = {1 / T:.3f} Hz")
a1.set(xlabel="Time (s)", ylabel="Position (m)")
a1.plot(t, pose_time, "-")
a2.set(xlabel="Time (s)", ylabel="Velocity (m/s)")
a2.plot(t, vel_time, "-g")
a3.set(xlabel="Position (m)", ylabel="Velocity (m/s)")
a3.plot(pose_time, vel_time, "-k")

plt.tight_layout()
plt.show()
