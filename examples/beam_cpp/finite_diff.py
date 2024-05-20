from core.problem import Prob
from core.logger import Logger
from core.solver.continuation import ConX

from examples.beam_cpp.beamcpp import BeamCpp
from Frame import Frame

import numpy as np
import h5py

# read files X_IC and posebase_IC from the folder
X_IC = np.loadtxt("X_IC.dat")
posebase_IC = np.loadtxt("posebase_IC.dat")

# read 1 files from folder
pose_1 = np.loadtxt("pose_1.dat")
vel_1 = np.loadtxt("vel_1.dat")
# X_1 = np.loadtxt("X_1.dat")

# runsim_single from beamcpp.py
prob = Prob()
prob.read_contparams("contparameters.json")
BeamCpp.initialise(prob.cont_params, False, 1)
BeamCpp.run_eig()

T = 0.004223328338788308

BeamCpp.config_update(posebase_IC)
cvg = BeamCpp.run_cpp(T, X_IC, 300, True)
with h5py.File(BeamCpp.cpp_path + BeamCpp.simout_file + ".h5", "r") as f:
    pose_T = f["/dynamic_analysis/FEModel/POSE/MOTION"][:, -1]
    vel_T = f["/dynamic_analysis/FEModel/VELOCITY/MOTION"][:, -1]
    M = f["/Sensitivity/Monodromy"][:]
    M[:X_IC.size, :X_IC.size] += np.eye(X_IC.size)

H1 = BeamCpp.periodicity_INC_SE(pose_1.reshape(-1, 1), pose_T.reshape(-1, 1))
H2 = BeamCpp.periodicity_VEL_SE(H1, vel_1.reshape(-1, 1), vel_T.reshape(-1, 1))
H = np.reshape(np.concatenate([H1, H2]), (-1, 1), order="F")

[dHdx0, dHdx1, dHdT] = BeamCpp.sensitivity_periodicity_SE_correction_multi(M, H1, vel_T, vel_1)

epsilon = 1e-8

# Central difference calculation for dHdx0
sens_X0 = np.zeros((H.size, H.size))
for i in range(H.size):
    X_IC_plus = X_IC.copy()
    X_IC_plus[i] += epsilon

    X_IC_minus = X_IC.copy()
    X_IC_minus[i] -= epsilon

    BeamCpp.config_update(posebase_IC)
    cvg_plus = BeamCpp.run_cpp(T, X_IC_plus, 300, False)
    with h5py.File(BeamCpp.cpp_path + BeamCpp.simout_file + ".h5", "r") as f:
        pose_T = f["/dynamic_analysis/FEModel/POSE/MOTION"][:, -1]
        vel_T = f["/dynamic_analysis/FEModel/VELOCITY/MOTION"][:, -1]

    H1 = BeamCpp.periodicity_INC_SE(pose_1.reshape(-1, 1), pose_T.reshape(-1, 1))
    H2 = BeamCpp.periodicity_VEL_SE(H1, vel_1.reshape(-1, 1), vel_T.reshape(-1, 1))
    H_plus = np.reshape(np.concatenate([H1, H2]), (-1, 1), order="F")

    BeamCpp.config_update(posebase_IC)
    cvg_minus = BeamCpp.run_cpp(T, X_IC_minus, 300, False)
    with h5py.File(BeamCpp.cpp_path + BeamCpp.simout_file + ".h5", "r") as f:
        pose_T = f["/dynamic_analysis/FEModel/POSE/MOTION"][:, -1]
        vel_T = f["/dynamic_analysis/FEModel/VELOCITY/MOTION"][:, -1]
    H3 = BeamCpp.periodicity_INC_SE(pose_1.reshape(-1, 1), pose_T.reshape(-1, 1))
    H4 = BeamCpp.periodicity_VEL_SE(H3, vel_1.reshape(-1, 1), vel_T.reshape(-1, 1))
    H_minus = np.reshape(np.concatenate([H3, H4]), (-1, 1), order="F")

    sens_X0[:, i] = (H_plus - H_minus).flatten() / (2 * epsilon)

np.savetxt("sens_X0.txt", sens_X0)
np.savetxt("dHdX0.txt", dHdx0)

# central difference calculation for period sensitivity
# T_plus = T + epsilon
# BeamCpp.config_update(posebase_IC)
# cvg_plus = BeamCpp.run_cpp(T_plus, X_IC, 300, False)
# with h5py.File(BeamCpp.cpp_path + BeamCpp.simout_file + ".h5", "r") as f:
#     pose_T = f["/dynamic_analysis/FEModel/POSE/MOTION"][:,-1]
#     vel_T = f["/dynamic_analysis/FEModel/VELOCITY/MOTION"][:, -1]

# H1 = BeamCpp.periodicity_INC_SE(pose_1.reshape(-1,1), pose_T.reshape(-1,1))
# H2 = BeamCpp.periodicity_VEL_SE(H1, vel_1.reshape(-1,1), vel_T.reshape(-1,1))
# H_plus = np.reshape(np.concatenate([H1, H2]), (-1, 1), order="F")

# BeamCpp.config_update(posebase_IC)
# T_minus = T - epsilon
# cvg_minus = BeamCpp.run_cpp(T_minus, X_IC, 300, False)
# with h5py.File(BeamCpp.cpp_path + BeamCpp.simout_file + ".h5", "r") as f:
#     pose_T = f["/dynamic_analysis/FEModel/POSE/MOTION"][:,-1]
#     vel_T = f["/dynamic_analysis/FEModel/VELOCITY/MOTION"][:, -1]
# H3 = BeamCpp.periodicity_INC_SE(pose_1.reshape(-1,1), pose_T.reshape(-1,1))
# H4 = BeamCpp.periodicity_VEL_SE(H3, vel_1.reshape(-1,1), vel_T.reshape(-1,1))
# H_minus = np.reshape(np.concatenate([H3, H4]), (-1, 1), order="F")

# sens_period = (H_plus - H_minus).flatten() / (2 * epsilon)
# np.savetxt("sens_period.txt", sens_period)
# np.savetxt("dHdT.txt", dHdT)

# Central difference calculation for dHdx1
# BeamCpp.config_update(posebase_IC)
# cvg_plus = BeamCpp.run_cpp(T, X_IC, 300, False)
# with h5py.File(BeamCpp.cpp_path + BeamCpp.simout_file + ".h5", "r") as f:
#     pose_T = f["/dynamic_analysis/FEModel/POSE/MOTION"][:,-1]
#     vel_T = f["/dynamic_analysis/FEModel/VELOCITY/MOTION"][:, -1]

# p = np.array([])
# for j in range(31):
#     frame = pose_1[j*4:(j+1)*4]
#     parameters = Frame.get_parameters_from_frame(2, frame)
#     p = np.concatenate([p, parameters])
# p = p[3:-3]
# X_1 = np.concatenate((p, vel_1[3:-3]))

# epsilon = 1e-7
# sensitivity1 = np.zeros((H.size, H.size))
# for i in range(H.size):

#     X_1_plus = X_1.copy()
#     X_1_plus[i] += epsilon

#     p = np.concatenate([np.zeros(3), X_1_plus[:X_IC.size//2], np.zeros(3)])
#     pose1_plus = np.array([])
#     for j in range(31):
#         parameters = p[j*3:(j+1)*3]
#         frame = Frame.get_frame_from_parameters(2, parameters)
#         pose1_plus = np.concatenate([pose1_plus, frame])

#     H1 = BeamCpp.periodicity_INC_SE(pose1_plus.reshape(-1,1), pose_T.reshape(-1,1))
#     v1 = np.concatenate([np.zeros(3), X_1_plus[X_IC.size//2:], np.zeros(3)])
#     H2 = BeamCpp.periodicity_VEL_SE(H1, v1.reshape(-1,1), vel_T.reshape(-1,1))
#     H_plus = np.reshape(np.concatenate([H1, H2]), (-1, 1), order="F")

#     X_1_minus = X_1.copy()
#     X_1_minus[i] -= epsilon

#     p = np.concatenate([np.zeros(3), X_1_minus[:X_IC.size//2], np.zeros(3)])
#     pose1_minus = np.array([])
#     for j in range(31):
#         parameters = p[j*3:(j+1)*3]
#         frame = Frame.get_frame_from_parameters(2, parameters)
#         pose1_minus = np.concatenate([pose1_minus, frame])

#     H1_ = BeamCpp.periodicity_INC_SE(pose1_minus.reshape(-1,1), pose_T.reshape(-1,1))
#     v1_ = np.concatenate([np.zeros(3), X_1_minus[X_IC.size//2:], np.zeros(3)])
#     H2_ = BeamCpp.periodicity_VEL_SE(H1_, v1_.reshape(-1,1), vel_T.reshape(-1,1))
#     H_minus = np.reshape(np.concatenate([H1_, H2_]), (-1, 1), order="F")

#     sensitivity1[:, i] = (H_plus - H_minus).flatten() / (2 * epsilon)

# diff2 = np.abs(sensitivity1 - dHdx1) / np.abs(sensitivity1)
