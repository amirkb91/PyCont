from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import mplcursors
import sys
import h5py

plt.style.use("ggplot")
files = sys.argv[1:]
if files[-1] == "y":
    saveflag = True
    files = files[:-1]
else:
    saveflag = False

# show point data on figure
def show_annotation(sel):
    ind = int(sel.index)
    sel.annotation.set_text(f"index:{ind}")


# figure properties
f, (a1, a2) = plt.subplots(1, 2, figsize=(8, 6))
a1.set(
    xlabel="Energy (J)",
    ylabel="Frequency (Hz)",
    # xlim=(1e-4, 1e2),
    # ylim=(50, 85),
    xscale="log",
)
a2.set(xlabel="Continuation step", ylabel="beta (deg)", ylim=(0, 180))

# plot sols
l = []
for file in files:
    # load data
    data = h5py.File(str(file), "r")
    T = data["/T"][:]
    Energy = data["/Energy"][:]

    # plot FEP
    l.append(
        a1.plot(Energy, 1 / T, marker=".", fillstyle="none", label=file.split(".h5")[0])
    )
    a1.plot(Energy[0], 1 / T[0], marker="x", fillstyle="full")

    # calculate beta
    try:
        tangent = data["/Tangent"][:]
        beta = np.array([])
        if np.shape(tangent)[1] > 1:  # psa solution file
            for i in range(1, np.shape(tangent)[1]):
                beta = np.append(
                    beta, np.rad2deg(np.arccos(tangent[:, i - 1].T @ tangent[:, i]))
                )
            a2.plot(range(len(beta)), beta, marker=".", fillstyle="none")
    except:
        pass

a1.legend()

# Cursor and plt
cursor = mplcursors.cursor(l[0], hover=False)
cursor.connect("add", show_annotation)
if saveflag:
    plt.savefig("FEP.pdf")
plt.draw()
plt.show()
