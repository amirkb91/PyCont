import h5py
import numpy as np
import os

dim = 3
if dim == 2:
    number_of_q, offset_q = 2, 4
else:
    number_of_q, offset_q = 4, 7


def tilde(x: np.array):
    x = x.reshape((3, 1))[:, 0]
    x_tilde = np.array([[0., -x[2], x[1]], [x[2], 0., -x[0]], [-x[1], x[0], 0.]])
    return x_tilde


def rotate_v(q, x):
    R = np.eye(3) + 2. * q[0] * tilde(q[1:]) + 2. * tilde(q[1:]) @ tilde(q[1:])
    return R @ x


def composition_unit_quat(a_, b):
    a = np.copy(a_)
    a[0] *= b[0]
    a[0] -= np.dot(a[1:], b[1:])
    if dim == 2:
        a[1] = a_[0] * b[1] + b[0] * a_[1]
    else:
        a[1:] = a_[0] * b[1:] + b[0] * a_[1:] + np.cross(a_[1:], b[1:])
    return a


def unit_quat_to_angle(q, q_ref):
    q_rel = np.copy(q_ref)
    q_rel[1:] *= -1.
    q_rel = composition_unit_quat(q_rel, q)

    angle = 2. * np.arccos(np.clip(q_rel[0], -1.0, 1.0))
    if np.abs(angle) > 1.0e-8:
        sca = 1. / np.sin(0.5 * angle)
    elif np.abs(angle) > 1.0e-12:
        sca = 1. / (0.5 * angle)
    else:
        sca = 2.

    n = q_rel[1:] * sca
    if np.linalg.norm(n):
        n = n / np.linalg.norm(n)

    return rotate_v(q_ref, angle * n)


file_name = os.getcwd() + '/beam_rightangle_shooting.h5'
file_name_eig = os.getcwd() + '/beam_rightangle_eig.h5'

ind = -1
with h5py.File(file_name, mode='r', libver=('earliest', 'v112')) as h5_file:
    pose_motion = h5_file["continuation_shooting/FEModel/POSE/MOTION"][:, ind]
    freq = h5_file['continuation_shooting/period'][ind]

with h5py.File(file_name_eig, mode='r', libver=('earliest', 'v112')) as h5_file:
    pose_motion_ref = h5_file["eigen_analysis/POSE/MOTION"][:, 0]

n_config_SE = pose_motion.shape[0]
n_nodes = int(n_config_SE / (number_of_q + dim))
n_config_VK = n_nodes * (number_of_q - 1 + dim)

new_file_name = os.getcwd() + '/beam_rightangle_IC.h5'

with h5py.File(new_file_name, mode='w', libver=('earliest', 'v112')) as h5_file:
    group_analysis = h5_file.create_group("dynamic_analysis")
    group_FEModel = group_analysis.create_group("FEModel")
    group_POSE = group_FEModel.create_group("POSE")

    dataset_IC = group_POSE.create_dataset("MOTION", shape=(n_config_VK, 1), dtype=np.float64)
    for i in range(n_nodes):
        i_x_VK, i_x_SE = i * (number_of_q - 1 + dim), number_of_q + i * (number_of_q + dim)
        dataset_IC[i_x_VK:i_x_VK + dim, 0] = pose_motion[i_x_SE:i_x_SE + dim]

        i_r_VK, i_r_SE = i_x_VK + dim, i_x_SE - number_of_q
        dataset_IC[i_r_VK:i_r_VK + (number_of_q - 1), 0] \
            = unit_quat_to_angle(pose_motion[i_r_SE:i_r_SE + number_of_q], pose_motion_ref[i_r_SE:i_r_SE + number_of_q])
