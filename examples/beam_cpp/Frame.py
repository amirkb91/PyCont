import numpy as np
from .SO3 import UnitQuaternion, tilde


class Frame:
    @staticmethod
    def get_inverse(n_dim, frame):
        # Calculate the inverse of the frame
        if n_dim == 2:
            q, x = frame[:2], frame[2:]
            q_inv = np.array([q[0], -q[1]])
        elif n_dim == 3:
            q, x = frame[:4], frame[4:]
            q_inv = np.array([q[0], -q[1], -q[2], -q[3]])
        x_inv = -1 * UnitQuaternion.rotateT_vec(n_dim, q, x)
        frame_inv = np.concatenate([q_inv, x_inv])
        return frame_inv

    @staticmethod
    def composition(n_dim, frame_a, frame_b):
        # f = frame_a o frame_b
        if n_dim == 2:
            frame_a_q, frame_a_x = frame_a[:2], frame_a[2:]
            frame_b_q, frame_b_x = frame_b[:2], frame_b[2:]
        elif n_dim == 3:
            frame_a_q, frame_a_x = frame_a[:4], frame_a[4:]
            frame_b_q, frame_b_x = frame_b[:4], frame_b[4:]
        q = UnitQuaternion.composition(n_dim, frame_a_q, frame_b_q)
        x = frame_a_x + UnitQuaternion.rotate_vec(n_dim, frame_a_q, frame_b_x)
        f = np.concatenate([q, x])
        return f

    @staticmethod
    def get_adjoint(n_dim, frame):
        # Calculate the adjoint of the frame
        if n_dim == 2:
            q, x = frame[:2], frame[2:]
            Adj = np.eye(3)
            R = UnitQuaternion.get_rotation_matrix(n_dim, q)
            Adj[:2, :2] = R
            Adj[0, 2] = x[1]
            Adj[1, 2] = -x[0]
        elif n_dim == 3:
            q, x = frame[:4], frame[4:]
            R = UnitQuaternion.get_rotation_matrix(n_dim, q)
            Adj = np.block([[R, np.matmul(tilde(x), R)], [np.zeros((3, 3)), R]])
        return Adj

    @staticmethod
    def lie_bracket(n_dim, vec_a, vec_b):
        if n_dim == 2:
            bracket = np.array(
                [
                    vec_b[2] * vec_a[1] - vec_a[2] * vec_b[0],
                    vec_a[2] * vec_b[1] - vec_b[2] * vec_a[0],
                    0.0,
                ]
            )
        else:
            a_x = vec_a[:3]
            a_r = vec_a[3:]
            b_x = vec_b[:3]
            b_r = vec_b[3:]
            bracket = np.concatenate((np.cross(a_r, b_x) + np.cross(a_x, b_r), np.cross(a_r, b_r)))
        return bracket

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
            p = np.zeros((3,))
            p[2] = 2.0 * np.sign(frame[0]) * frame[1]
            p[:2] = frame[0] * frame[2:] + np.array([frame[3], -frame[2]]) * 0.5 * p[2]
        elif n_dim == 3:
            p = np.zeros((6,))
            p[3:] = 2.0 * np.sign(frame[0]) * frame[1:4]
            Tm1T = frame[0] * np.eye(3) - tilde(0.5 * p[3:])
            p[:3] = np.matmul(Tm1T, frame[4:])
        return p

    @staticmethod
    def get_frame_from_parameters(n_dim, parameters):
        if n_dim == 2:
            p0 = np.sqrt(1.0 - 0.25 * parameters[2] ** 2)
            q = np.array([p0, 0.5 * parameters[2]])
            p_tilde_over_2 = np.array([[0, -0.5 * parameters[2]], [0.5 * parameters[2], 0]])
            TT = 1.0 / p0 * (np.eye(2) + np.matmul(p_tilde_over_2, p_tilde_over_2)) + p_tilde_over_2
            x = np.matmul(TT, parameters[:2])
        elif n_dim == 3:
            p0 = np.sqrt(1.0 - 0.25 * np.dot(parameters[3:], parameters[3:]))
            q = np.concatenate([np.array([p0]), 0.5 * parameters[3:]])
            p_tilde_over_2 = tilde(0.5 * parameters[3:])
            TT = 1.0 / p0 * (np.eye(3) + np.matmul(p_tilde_over_2, p_tilde_over_2)) + p_tilde_over_2
            x = np.matmul(TT, parameters[:3])
        return np.concatenate([q, x])

    @classmethod
    def add_inc_to_frame_local(cls, n_dim, frame, inc):
        # out = frame o exp(inc)
        frame_from_inc = cls.get_frame_from_parameters(n_dim, inc)
        out = cls.composition(n_dim, frame, frame_from_inc)
        return out

    @classmethod
    def add_inc_to_frame_global(cls, n_dim, frame, inc):
        # out = exp(inc) o frame
        frame_from_inc = cls.get_frame_from_parameters(n_dim, inc)
        out = cls.composition(n_dim, frame_from_inc, frame)
        return out

    @staticmethod
    def get_inverse_tangent_operator(n_dim, parameters):
        # Calculate the inverse tangent operator from the vector parameters
        if n_dim == 2:
            p0 = np.sqrt(1.0 - 0.25 * parameters[2] ** 2)
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

    @staticmethod
    def get_derivative_inverse_tangent_operator(n_dim, parameters, direction):
        # Calculate the derivative inverse tangent operator from the vector parameters along a direction
        if n_dim == 2:
            DTinv = np.zeros((3, 3))
            p0 = np.sqrt(1.0 - 0.25 * parameters[2] ** 2)
            c = -0.25 / p0 * parameters[2]
            DTinv[0, 1] = 0.5 * direction[2]
            DTinv[1, 0] = -0.5 * direction[2]
            DTinv[2, 2] = -0.25 / p0 * direction[2] * parameters[2]
            DTinv[0, 2] = c * direction[0] - 0.5 * direction[1]
            DTinv[1, 2] = c * direction[1] + 0.5 * direction[0]

        elif n_dim == 3:
            p0 = np.sqrt(1.0 - 0.25 * np.dot(parameters[3:], parameters[3:]))
            rho = np.dot(parameters[:3], parameters[3:])
            pu, pw = parameters[:3].reshape((3, 1)), parameters[3:].reshape((3, 1))
            du, dw = direction[:3].reshape((3, 1)), direction[3:].reshape((3, 1))

            DTinv0u = tilde(-0.5 * du) - ((0.25 / p0) * (du @ pw.T))
            DTinv0u += ((-0.25 * rho * 0.25 / (p0**3)) * (dw @ pw.T)) - ((0.25 / p0) * (dw @ pu.T))
            DTinv0r = tilde(-0.5 * dw) - ((0.25 / p0) * (dw @ pw.T))

            DTinv = np.block([[DTinv0r, DTinv0u], [np.zeros((3, 3)), DTinv0r]])

        return DTinv
