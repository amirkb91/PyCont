import h5py

data = h5py.File("nnm11.h5", "r")
solno = 12

pose = data["/Config/POSE"][:33, solno]
vel = data["/Config/VELOCITY"][:33, solno]

icdata = h5py.File("ic.h5", "w")
icdata["/dynamic_analysis/FEModel/POSE/MOTION"] = pose
icdata["/dynamic_analysis/FEModel/VELOCITY/MOTION"] = vel
icdata.close()
