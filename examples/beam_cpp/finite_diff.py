from core.problem import Prob

from examples.beam_cpp.beamcpp import BeamCpp
from core.math.Frame import Frame

import numpy as np
import h5py

# Finite difference code to check sensitivity of multiple shooting periodicity condition

# Read ic0 files
with h5py.File("ic0.h5", "r") as f:
    pose0 = f["/dynamic_analysis/FEModel/POSE/MOTION"][:]
    vel0 = f["/dynamic_analysis/FEModel/VELOCITY/MOTION"][:]

# Read ic1 files
with h5py.File("ic1.h5", "r") as f:
    pose1 = f["/dynamic_analysis/FEModel/POSE/MOTION"][:]
    vel1 = f["/dynamic_analysis/FEModel/VELOCITY/MOTION"][:]

# Read eig
with h5py.File("beam_eig.h5", "r") as f:
    pose_ref = f["/eigen_analysis/POSE/MOTION"][:, 0]

# convert pose0 and vel0 to X0
p = np.array([])
for j in range(31):
    frame_0 = pose0[j * 4 : (j + 1) * 4]
    frame_ref = pose_ref[j * 4 : (j + 1) * 4]
    frame = Frame.relative_frame(2, frame_ref, frame_0)
    parameters = Frame.get_parameters_from_frame(2, frame)
    p = np.concatenate([p, parameters])
p = p[3:-3]
X0 = np.concatenate((p, vel0[3:-3]))

# convert pose1 and vel1 to X1
p = np.array([])
for j in range(31):
    frame_1 = pose1[j * 4 : (j + 1) * 4]
    frame_ref = pose_ref[j * 4 : (j + 1) * 4]
    frame = Frame.relative_frame(2, frame_ref, frame_1)
    parameters = Frame.get_parameters_from_frame(2, frame)
    p = np.concatenate([p, parameters])
p = p[3:-3]
X1 = np.concatenate((p, vel1[3:-3]))

# Period
T = 8e-3

# runsim_single from beamcpp.py
prob = Prob()
prob.read_contparams("contparameters.json")
BeamCpp.initialise(prob.cont_params, False, 1)
BeamCpp.run_eig()

BeamCpp.config_update(pose_ref)
cvg = BeamCpp.run_cpp(T, X0, 300, True)
with h5py.File(BeamCpp.cpp_path + BeamCpp.simout_file + ".h5", "r") as f:
    poseT = f["/dynamic_analysis/FEModel/POSE/MOTION"][:, -1]
    velT = f["/dynamic_analysis/FEModel/VELOCITY/MOTION"][:, -1]
    M = f["/Sensitivity/Monodromy"][:]
    M[: X0.size, : X0.size] += np.eye(X0.size)

H1 = BeamCpp.periodicity_INC_SE_local(pose1.reshape(-1, 1), poseT.reshape(-1, 1))
H2 = BeamCpp.periodicity_VEL_SE_local(H1, vel1.reshape(-1, 1), velT.reshape(-1, 1))
H = np.reshape(np.concatenate([H1, H2]), (-1, 1), order="F")

[dHdx0, dHdx1, dHdT] = BeamCpp.sens_periodicity_SEcorrection_local_multi(M, H1, velT, vel1)
np.savetxt("dHdX0.txt", dHdx0)
np.savetxt("dHdX1.txt", dHdx1)
np.savetxt("dHdT.txt", dHdT)

epsilon = 1e-5

# Central difference calculation for dHdx0
sens_X0 = np.zeros((H.size, H.size))
for i in range(H.size):
    X0_plus = X0.copy()
    X0_plus[i] += epsilon

    X0_minus = X0.copy()
    X0_minus[i] -= epsilon

    BeamCpp.config_update(pose_ref)
    cvg_plus = BeamCpp.run_cpp(T, X0_plus, 300, False)
    with h5py.File(BeamCpp.cpp_path + BeamCpp.simout_file + ".h5", "r") as f:
        poseT = f["/dynamic_analysis/FEModel/POSE/MOTION"][:, -1]
        velT = f["/dynamic_analysis/FEModel/VELOCITY/MOTION"][:, -1]

    H1 = BeamCpp.periodicity_INC_SE_local(pose1.reshape(-1, 1), poseT.reshape(-1, 1))
    H2 = BeamCpp.periodicity_VEL_SE_local(H1, vel1.reshape(-1, 1), velT.reshape(-1, 1))
    H_plus = np.reshape(np.concatenate([H1, H2]), (-1, 1), order="F")

    BeamCpp.config_update(pose_ref)
    cvg_minus = BeamCpp.run_cpp(T, X0_minus, 300, False)
    with h5py.File(BeamCpp.cpp_path + BeamCpp.simout_file + ".h5", "r") as f:
        poseT = f["/dynamic_analysis/FEModel/POSE/MOTION"][:, -1]
        velT = f["/dynamic_analysis/FEModel/VELOCITY/MOTION"][:, -1]
    H3 = BeamCpp.periodicity_INC_SE_local(pose1.reshape(-1, 1), poseT.reshape(-1, 1))
    H4 = BeamCpp.periodicity_VEL_SE_local(H3, vel1.reshape(-1, 1), velT.reshape(-1, 1))
    H_minus = np.reshape(np.concatenate([H3, H4]), (-1, 1), order="F")

    sens_X0[:, i] = (H_plus - H_minus).flatten() / (2 * epsilon)

np.savetxt("sens_X0.txt", sens_X0)

# Central difference calculation for dHdx1
sens_X1 = np.zeros((H.size, H.size))
for i in range(H.size):
    X1_plus = X1.copy()
    X1_plus[i] += epsilon

    pose1_plus = np.array([])
    p = np.concatenate([np.zeros(3), X1_plus[: X0.size // 2], np.zeros(3)])
    for j in range(31):
        parameters = p[j * 3 : (j + 1) * 3]
        frame = Frame.get_frame_from_parameters(2, parameters)
        pose1_plus = np.concatenate([pose1_plus, frame])

    X1_minus = X1.copy()
    X1_minus[i] -= epsilon

    pose1_minus = np.array([])
    p = np.concatenate([np.zeros(3), X1_minus[: X0.size // 2], np.zeros(3)])
    for j in range(31):
        parameters = p[j * 3 : (j + 1) * 3]
        frame = Frame.get_frame_from_parameters(2, parameters)
        pose1_minus = np.concatenate([pose1_minus, frame])

    BeamCpp.config_update(pose_ref)
    cvg_plus = BeamCpp.run_cpp(T, X0, 300, False)
    with h5py.File(BeamCpp.cpp_path + BeamCpp.simout_file + ".h5", "r") as f:
        poseT = f["/dynamic_analysis/FEModel/POSE/MOTION"][:, -1]
        velT = f["/dynamic_analysis/FEModel/VELOCITY/MOTION"][:, -1]

    H1 = BeamCpp.periodicity_INC_SE_local(pose1_plus.reshape(-1, 1), poseT.reshape(-1, 1))
    v1_ = np.concatenate([np.zeros(3), X1_plus[X0.size // 2 :], np.zeros(3)])
    H2 = BeamCpp.periodicity_VEL_SE_local(H1, v1_.reshape(-1, 1), velT.reshape(-1, 1))
    H_plus = np.reshape(np.concatenate([H1, H2]), (-1, 1), order="F")

    H3 = BeamCpp.periodicity_INC_SE_local(pose1_minus.reshape(-1, 1), poseT.reshape(-1, 1))
    v1_ = np.concatenate([np.zeros(3), X1_minus[X0.size // 2 :], np.zeros(3)])
    H4 = BeamCpp.periodicity_VEL_SE_local(H3, v1_.reshape(-1, 1), velT.reshape(-1, 1))
    H_minus = np.reshape(np.concatenate([H3, H4]), (-1, 1), order="F")

    sens_X1[:, i] = (H_plus - H_minus).flatten() / (2 * epsilon)

np.savetxt("sens_X1.txt", sens_X1)

# central difference calculation for period sensitivity
T_plus = T + epsilon
BeamCpp.config_update(pose_ref)
cvg_plus = BeamCpp.run_cpp(T_plus, X0, 300, False)
with h5py.File(BeamCpp.cpp_path + BeamCpp.simout_file + ".h5", "r") as f:
    pose_T = f["/dynamic_analysis/FEModel/POSE/MOTION"][:, -1]
    vel_T = f["/dynamic_analysis/FEModel/VELOCITY/MOTION"][:, -1]

H1 = BeamCpp.periodicity_INC_SE_local(pose1.reshape(-1, 1), pose_T.reshape(-1, 1))
H2 = BeamCpp.periodicity_VEL_SE_local(H1, vel1.reshape(-1, 1), vel_T.reshape(-1, 1))
H_plus = np.reshape(np.concatenate([H1, H2]), (-1, 1), order="F")

BeamCpp.config_update(pose_ref)
T_minus = T - epsilon
cvg_minus = BeamCpp.run_cpp(T_minus, X0, 300, False)
with h5py.File(BeamCpp.cpp_path + BeamCpp.simout_file + ".h5", "r") as f:
    pose_T = f["/dynamic_analysis/FEModel/POSE/MOTION"][:, -1]
    vel_T = f["/dynamic_analysis/FEModel/VELOCITY/MOTION"][:, -1]
H3 = BeamCpp.periodicity_INC_SE_local(pose1.reshape(-1, 1), pose_T.reshape(-1, 1))
H4 = BeamCpp.periodicity_VEL_SE_local(H3, vel1.reshape(-1, 1), vel_T.reshape(-1, 1))
H_minus = np.reshape(np.concatenate([H3, H4]), (-1, 1), order="F")

sens_period = (H_plus - H_minus).flatten() / (2 * epsilon)
np.savetxt("sens_T.txt", sens_period)
