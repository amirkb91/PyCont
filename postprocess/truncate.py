import sys
import h5py
import json

# inputs
file = sys.argv[-1]
solno = int(input("Solution Index: "))
truncate = input("Truncate (h)ead or (t)ail?: ")

file_t = file.strip(".h5") + "_trunc.h5"
data = h5py.File(file, "r")
data_t = h5py.File(file_t, "w")

if truncate == "t":
    data_t["/Energy"] = data["/Energy"][:solno + 1]
    data_t["/T"] = data["/T"][:solno + 1]
    data_t["/itercorrect"] = data["/itercorrect"][:solno + 1]
    data_t["/step"] = data["/step"][:solno + 1]
    data_t["/Tangent"] = data["/Tangent"][:, :solno + 1]
    data_t["/beta"] = data["/beta"][:solno + 1]
    data_t["/Parameters"] = json.dumps(json.loads(data["/Parameters"][()]))
    data_t["/Config/POSE"] = data["/Config/POSE"][:, :solno + 1]
    data_t["/Config/VELOCITY"] = data["/Config/VELOCITY"][:, :solno + 1]
elif truncate == "h":
    data_t["/Energy"] = data["/Energy"][solno:]
    data_t["/T"] = data["/T"][solno:]
    data_t["/itercorrect"] = data["/itercorrect"][solno:]
    data_t["/step"] = data["/step"][solno:]
    data_t["/Tangent"] = data["/Tangent"][:, solno:]
    data_t["/beta"] = data["/beta"][solno:]
    data_t["/Parameters"] = json.dumps(json.loads(data["/Parameters"][()]))
    data_t["/Config/POSE"] = data["/Config/POSE"][:, solno:]
    data_t["/Config/VELOCITY"] = data["/Config/VELOCITY"][:, solno:]
data.close()
data_t.close()
