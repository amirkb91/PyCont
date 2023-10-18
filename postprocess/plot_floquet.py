import numpy as np
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider
import sys
import h5py

file = sys.argv[1]
if not file.endswith(".h5"):
    file += ".h5"

f, a = plt.subplots(figsize=(10, 7))
f.subplots_adjust(left=0.25)
a.set(xlabel="Real", ylabel="Imaginary")
a.axis("square")
a.grid("on")
a.plot(np.cos(np.linspace(0, 2 * np.pi, 1000)), np.sin(np.linspace(0, 2 * np.pi, 1000)), "-")
a.set_xlim([-1.5, 1.5])
a.set_ylim([-1.5, 1.5])

data = h5py.File(str(file), "r")
floquet = data["/Bifurcation/Floquet"][:]
n = np.shape(floquet)[1] - 1
(points,) = a.plot(floquet.real[:, 0], floquet.imag[:, 0], "o", markersize=5, color="orange")

ax = f.add_axes([0.1, 0.15, 0.0225, 0.63])
slider = Slider(
    ax=ax,
    label="Sol no.",
    valmin=0,
    valmax=n,
    valinit=0,
    valstep=1,
    orientation="vertical",
    color="orange",
)


def update(val):
    sol_no = slider.val
    points.set_xdata(floquet.real[:, sol_no])
    points.set_ydata(floquet.imag[:, sol_no])
    f.canvas.draw_idle()


slider.on_changed(update)
plt.show()
