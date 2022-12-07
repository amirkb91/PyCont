import h5py, json

file = "/home/akb110/Codes/PyCont/TestCases/cclamped/NNM2_VK_4.h5"
index = 29
truncate = "head"

file_t = file.strip(".h5") + "_trunc.h5"
data = h5py.File(file, "r")
data_t = h5py.File(file_t, "w")

if truncate == "tail":
    data_t["/Energy"] = data["/Energy"][: index + 1]
    data_t["/T"] = data["/T"][: index + 1]
    data_t["/X"] = data["/X"][:, : index + 1]
    data_t["/POSE"] = data["/POSE"][:, :, : index + 1]
    data_t["/Tangent"] = data["/Tangent"][:, : index + 1]
    data_t["/Parameters"] = json.dumps(json.loads(data["/Parameters"][()]))
elif truncate == "head":
    data_t["/Energy"] = data["/Energy"][index:]
    data_t["/T"] = data["/T"][index:]
    data_t["/X"] = data["/X"][:, index:]
    data_t["/POSE"] = data["/POSE"][:, :, index:]
    data_t["/Tangent"] = data["/Tangent"][:, index:]
    data_t["/Parameters"] = json.dumps(json.loads(data["/Parameters"][()]))

data.close()
data_t.close()
