import h5py
import os
import json

folder = "Results/Forced/Forced_offcentre_cclamped/"

files = os.listdir(folder)

for file in files:
    if "h5" in file:
        data = h5py.File(folder + file, "a")
        par = data["/Parameters"]
        par = json.loads(par[()])
        par["forcing"]["tau0"] = 0.0
        del data["/Parameters"]
        data["/Parameters"] = json.dumps(par)
        data.close()
