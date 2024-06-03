import numpy as np


def tilde(x):
    x_tilde = np.zeros((3, 3))
    x_tilde[0, 1] = -x[2]
    x_tilde[0, 2] = x[1]
    x_tilde[1, 0] = x[2]
    x_tilde[1, 2] = -x[0]
    x_tilde[2, 0] = -x[1]
    x_tilde[2, 1] = x[0]
    return x_tilde


class UnitQuaternion:
    @staticmethod
    def relative_rotation(n_dim, q_a, q_b):
        if n_dim == 2:
            q = np.array([q_a[0] * q_b[0] + q_a[1] * q_b[1], q_a[0] * q_b[1] - q_a[1] * q_b[0]])
        elif n_dim == 3:
            q = np.array(
                [
                    q_a[0] * q_b[0] + q_a[1] * q_b[1] + q_a[2] * q_b[2] + q_a[3] * q_b[3],
                    q_a[0] * q_b[1] - q_a[1] * q_b[0] - q_a[2] * q_b[3] + q_a[3] * q_b[2],
                    q_a[0] * q_b[2] - q_a[2] * q_b[0] - q_a[3] * q_b[1] + q_a[1] * q_b[3],
                    q_a[0] * q_b[3] - q_a[3] * q_b[0] - q_a[1] * q_b[2] + q_a[2] * q_b[1],
                ]
            )
        return q

    @staticmethod
    def composition(n_dim, q_a, q_b):
        if n_dim == 2:
            q = np.array([q_a[0] * q_b[0] - q_a[1] * q_b[1], q_a[0] * q_b[1] + q_a[1] * q_b[0]])
        elif n_dim == 3:
            q = np.array(
                [
                    q_a[0] * q_b[0] - q_a[1] * q_b[1] - q_a[2] * q_b[2] - q_a[3] * q_b[3],
                    q_a[0] * q_b[1] + q_a[1] * q_b[0] + q_a[2] * q_b[3] - q_a[3] * q_b[2],
                    q_a[0] * q_b[2] + q_a[2] * q_b[0] + q_a[3] * q_b[1] - q_a[1] * q_b[3],
                    q_a[0] * q_b[3] + q_a[3] * q_b[0] + q_a[1] * q_b[2] - q_a[2] * q_b[1],
                ]
            )
        return q

    @staticmethod
    def get_rotation_matrix(n_dim, q):
        if n_dim == 2:
            s = -2 * q[1]
            c = 1 + q[1] * s
            s *= q[0]
            R = np.array([[c, -s], [s, c]])
        else:
            e0 = q[0]
            e = q[1:]
            e_tilde = tilde(e)
            R = np.eye(3) + 2.0 * (e0 * e_tilde + np.matmul(e_tilde, e_tilde))
        return R

    @staticmethod
    def rotate_vec(n_dim, q, vec):
        if n_dim == 2:
            s = 2 * q[1]
            c = 1 - s * q[1]
            s *= q[0]
            out = np.array([c * vec[0] - s * vec[1], s * vec[0] + c * vec[1]])
        else:
            e0e0 = q[0] * q[0]
            e1e2 = q[1] * q[2]
            e1e3 = q[1] * q[3]
            e0e2 = q[0] * q[2]
            e0e3 = q[0] * q[3]
            e2e3 = q[2] * q[3]
            e0e1 = q[0] * q[1]

            x = 2.0 * (
                (q[1] * q[1] + e0e0 - 0.5) * vec[0]
                + (e1e2 - e0e3) * vec[1]
                + (e0e2 + e1e3) * vec[2]
            )

            y = 2.0 * (
                (q[2] * q[2] + e0e0 - 0.5) * vec[1]
                + (e1e2 + e0e3) * vec[0]
                + (-e0e1 + e2e3) * vec[2]
            )

            z = 2.0 * (
                (q[3] * q[3] + e0e0 - 0.5) * vec[2]
                + (e1e3 - e0e2) * vec[0]
                + (e0e1 + e2e3) * vec[1]
            )

            out = np.array([x, y, z])

        return out

    @staticmethod
    def rotateT_vec(n_dim, q, vec):
        if n_dim == 2:
            s = -2 * q[1]
            c = 1 + q[1] * s
            s *= q[0]
            out = np.array([c * vec[0] - s * vec[1], s * vec[0] + c * vec[1]])
        else:
            e0e0 = q[0] * q[0]
            e1e2 = q[1] * q[2]
            e1e3 = q[1] * q[3]
            e0e2 = q[0] * q[2]
            e0e3 = q[0] * q[3]
            e2e3 = q[2] * q[3]
            e0e1 = q[0] * q[1]

            x = 2 * (
                (q[1] * q[1] + e0e0 - 0.5) * vec[0]
                + (e1e2 + e0e3) * vec[1]
                + (-e0e2 + e1e3) * vec[2]
            )
            y = 2 * (
                (q[2] * q[2] + e0e0 - 0.5) * vec[1]
                + (e1e2 - e0e3) * vec[0]
                + (e0e1 + e2e3) * vec[2]
            )
            z = 2 * (
                (q[3] * q[3] + e0e0 - 0.5) * vec[2]
                + (e1e3 + e0e2) * vec[0]
                + (-e0e1 + e2e3) * vec[1]
            )

            out = np.array([x, y, z])

        return out
