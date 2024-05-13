import numpy as np
from .SO3 import UnitQuaternion, tilde


class Frame:
    @classmethod
    def relative_frame(cls, n_dim, frame_a, frame_b):
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

    @classmethod
    def get_parameters_from_frame(cls, n_dim, frame):
        if n_dim == 2:
            p = np.zeros((3, ))
            p[2] = 2.0 * np.sign(frame[0]) * frame[1]
            p[:2] = frame[0] * frame[2:] + np.array([frame[3], -frame[2]]) * frame[1]
        elif n_dim == 3:
            p = np.zeros((6, ))
            p[3:] = 2.0 * np.sign(frame[0]) * frame[1:4]
            Tm1T = frame[0] * np.eye(3) - tilde(frame[1:4])
            p[:3] = np.matmul(Tm1T, frame[4:])
        return p
