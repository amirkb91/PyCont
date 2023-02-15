import numpy as np
from scipy.integrate import odeint, solve_ivp
import scipy.linalg as spl
import matplotlib.pyplot as plt


def statespace(Z, t):
    # ODE of the mass spring system. Z\dot(t) = g(Z(t))
    M = np.eye(2)
    K = np.array([[2, -1], [-1, 2]])
    mu = 0.5

    X = Z[:2]
    Xdot = Z[2:]
    Minv = spl.inv(M)
    KX = K @ X
    fnl = np.array([mu * X[0] ** 3, 0])
    Zdot = np.concatenate((Xdot, -Minv @ (KX + fnl)))

    return Zdot


def statespace2(t, Z):
    # ODE of the mass spring system. Z\dot(t) = g(Z(t))
    M = np.eye(2)
    K = np.array([[2, -1], [-1, 2]])
    mu = 0.5

    X = Z[:2]
    Xdot = Z[2:]
    Minv = spl.inv(M)
    KX = K @ X
    fnl = np.array([mu * X[0] ** 3, 0])
    Zdot = np.concatenate((Xdot, -Minv @ (KX + fnl)))

    return Zdot


t = np.linspace(0, 2, 100)
Z0 = np.array([1, 1, 0, 0])
Z_all = odeint(statespace, Z0, t)

Z_all2 = solve_ivp(statespace2, [t[0], t[-1]], Z0, t_eval=t)
pp = Z_all2.y.T

plt.figure()
plt.plot(t, Z_all[:, :1], '-ro')
plt.plot(t, pp[:, :1], '-bo')
plt.show()
