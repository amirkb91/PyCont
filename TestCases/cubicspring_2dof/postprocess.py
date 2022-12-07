from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import mplcursors

# show point data on figure
def show_annotation(sel):
    ind = int(sel.index)
    x, y = sel.target
    sel.annotation.set_text(f"index:{ind} \nx:{x:.2f} y:{y:.2f}")


# load files
folder = ""
XT = np.loadtxt(folder + "solXT.out")
Energy = np.loadtxt(folder + "energy.out")
tan = np.loadtxt(folder + "tgt.out")
T = XT[-1, :]

# calculate beta
# beta = np.array([])
# for i in range(1, np.shape(tan)[1]):
#     beta = np.append(beta, np.rad2deg(np.arccos(tan[:, i - 1].T @ tan[:, i])))


# plot FEP and beta
f, (a1, a2) = plt.subplots(1, 2, figsize=(8, 6))

(l1,) = a1.plot(Energy, 1 / T, marker="o", fillstyle="none")
a1.set_xscale("log")
a1.grid()
a1.set_xlabel("Energy (J)")
a1.set_ylabel("Frequency (Hz)")
a1.set_xlim(1e-5, 1e3)
a1.set_ylim(0.1, 1)

# a2.plot(1 / T[1:], beta, marker="o", fillstyle="none", color="red")
# a2.grid()
# a2.set_xlabel("Frequency (Hz)")
# a2.set_ylabel("beta (deg)")
# a2.set_ylim(0, 90)

# --- secondary plot
# folder = "../../results/NNM1_abovetongue/"
# XT = np.loadtxt(folder + "solXT.out")
# Energy = np.loadtxt(folder + "energy.out")
# T = XT[-1, :]
# a1.plot(Energy, 1 / T, marker="*")

# Cursor and plt
# cursor = mplcursors.cursor(l1, hover=True)
# cursor.connect("add", show_annotation)
plt.draw()
plt.show()
