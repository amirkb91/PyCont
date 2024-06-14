import h5py
import json
import numpy as np
from scipy.integrate import odeint
from cubic_spring import Cubic_Spring

# Period and residual values taken from debug as well
T_values = [4.493740101081286, 4.497073847490682, 4.49781822041693]
res_values = [1.527, 0.1843, 1.4310e-5]

# read pose and vel from txt files
pose_files = ["pose0_debug.txt", "pose1_debug.txt", "pose2_debug.txt"]
vel_files = ["vel0_debug.txt", "vel1_debug.txt", "vel2_debug.txt"]

pose_data = [np.loadtxt(file) for file in pose_files]
vel_data = [np.loadtxt(file) for file in vel_files]

# read solution file
file = "REF_NNM1_Mult_new.h5"
data = h5py.File(str(file), "r")

# do time simulation
par = data["/Parameters"]
par = json.loads(par[()])
npartition = par["shooting"]["multiple"]["npartition"]
nsteps = par["shooting"]["multiple"]["nsteps_per_partition"]
delta_S = 1 / npartition
t = np.zeros([nsteps + 1, 1, npartition])
pose_time = np.zeros([nsteps + 1, 2, npartition])
vel_time = np.zeros([nsteps + 1, 2, npartition])

for idx, (T, pose, vel) in enumerate(zip(T_values, pose_data, vel_data)):
    for ipart in range(npartition):
        partition_starttime = ipart * T * delta_S
        print(f"{partition_starttime:.3f}")
        t_part = np.linspace(0, T * delta_S, nsteps + 1)
        t[:, :, ipart] = t_part.reshape(-1, 1) + partition_starttime
        X = np.concatenate([pose[:, ipart], vel[:, ipart]])
        timesol = np.array(odeint(Cubic_Spring.model_ode, X, t_part, rtol=1e-8, tfirst=True))
        pose_time[:, :, ipart] = timesol[:, :2]
        vel_time[:, :, ipart] = timesol[:, 2:]
        np.savetxt(f"corr{idx}_partition_{ipart}_t.txt", t[:, 0, ipart])
        np.savetxt(f"corr{idx}_partition_{ipart}_pose_time.txt", pose_time[:, :, ipart])
        np.savetxt(f"corr{idx}_partition_{ipart}_vel_time.txt", vel_time[:, :, ipart])
