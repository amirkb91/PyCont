import numpy as np
from .SO3 import UnitQuaternion, tilde


class Frame:
    @staticmethod
    def relative_frame(n_dim, frame_a, frame_b):
        # Calculate the relative frame from frame_a to frame_b
        # f = frame_a^-1 o frame_b
        if n_dim == 2:
            frame_a_q, frame_a_x = frame_a[:2], frame_a[2:]
            frame_b_q, frame_b_x = frame_b[:2], frame_b[2:]
        elif n_dim == 3:
            frame_a_q, frame_a_x = frame_a[:4], frame_a[4:]
            frame_b_q, frame_b_x = frame_b[:4], frame_b[4:]
        q = UnitQuaternion.relative_rotation(n_dim, frame_a_q, frame_b_q)
        x = UnitQuaternion.rotateT_vec(n_dim, frame_a_q, frame_b_x - frame_a_x)
        f = np.concatenate([q, x])
        return f

    @staticmethod
    def get_parameters_from_frame(n_dim, frame):
        # Calculate the vector parameters from the frame
        if n_dim == 2:
            p = np.zeros((3, ))
            p[2] = 2.0 * np.sign(frame[0]) * frame[1]
            p[:2] = frame[0] * frame[2:] + np.array([frame[3], -frame[2]]) * 0.5 * p[2]
        elif n_dim == 3:
            p = np.zeros((6, ))
            p[3:] = 2.0 * np.sign(frame[0]) * frame[1:4]
            Tm1T = frame[0] * np.eye(3) - tilde(0.5 * p[3:])
            p[:3] = np.matmul(Tm1T, frame[4:])
        return p

    @staticmethod
    def get_inverse_tangent_operator(n_dim, parameters):
        # Calculate the inverse tangent operator from the vector parameters
        parameters = parameters.flatten()
        if n_dim == 2:
            p0 = np.sqrt(1.0 - 0.25 * parameters[2]**2)
            Tm1 = p0 * np.eye(3)
            Tm1[0, 1] = -0.5 * parameters[2]
            Tm1[1, 0] = 0.5 * parameters[2]
            Tm1[0, 2] = 0.5 * parameters[1]
            Tm1[1, 2] = -0.5 * parameters[0]
        elif n_dim == 3:
            p0 = np.sqrt(1.0 - 0.25 * np.dot(parameters[3:], parameters[3:]))
            rho = np.dot(parameters[:3], parameters[3:])
            Tm1 = np.zeros((6, 6))
            Tm1[:3, :3] = Tm1[3:, 3:] = p0 * np.eye(3) + tilde(0.5 * parameters[3:])
            Tm1[:3, 3:] = -rho / (4.0 * p0) * np.eye(3) + tilde(0.5 * parameters[:3])
        return Tm1
