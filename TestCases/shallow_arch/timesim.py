import numpy as np
import h5py
import json
import subprocess
import sys

"""
timeplot.py exists in beam_cpp directory.
run cpp simulation with IC from given continuation solution
call that script to get time plots for required continuation solution
"""

# inputs
solno = int(input("Solution Index: "))
nper = int(input("Number of periods: "))
ns = int(input("Steps per period: "))

# read solution file
file = sys.argv[-1]
data = h5py.File(str(file), "r")
X = data["/X"][:, solno]
T = data["/T"][solno]

# read parameters from solution file
par = data["/Parameters"]
par = json.loads(par[()])
rho = par["shooting"]["rho"]
rel_tol = par["shooting"]["rel_tol"]
beam_type = par["shooting"]["beam_type"]

# cpp run
# path and file names, open input file
case = "shallow_arch"
path2cpp = "/home/akb110/Codes/beam_cpp/examples/" + case + "/"
exec_cpp = "/home/akb110/Codes/beam_cpp/cmake-build-release/examples/" + case
paramfile = "parameters.json"
eigenfile = case + "_eig"
ouputfile = case + "_out"
inputdata = h5py.File(path2cpp + eigenfile + ".h5", "r+")
ndof_all = inputdata["number_of_dofs"][0][0]

# clear existing data in eigen file
if "Config/INC" in inputdata:
    del inputdata["Config/INC"]
del inputdata["Config/VELOCITY"]
inputdata.create_dataset("Config/INC", shape=(ndof_all,), dtype=np.dtype("float64"))
inputdata.create_dataset(
    "Config/VELOCITY", shape=(ndof_all,), dtype=np.dtype("float64")
)

# add X to input file (INC, VEL)
lenX_2 = len(X) // 2
inc = X[:lenX_2]
vel = X[lenX_2:]
# add zeros for bc nodes
inc = np.pad(inc, (6, 0), "constant")
vel = np.pad(vel, (6, 0), "constant")
# store in inputdata to use as initial condition
inputdata["/Config/INC"][:] = inc
inputdata["/Config/VELOCITY"][:] = vel
inputdata.close()

# edit C++ parameter file
cpp_parameter = json.load(open(path2cpp + paramfile))
cpp_parameter["output_file"] = ouputfile
cpp_parameter["input_state"] = eigenfile
cpp_parameter["beam_type"] = beam_type
cpp_parameter["Solver_parameters"]["number_of_steps"] = ns * nper
cpp_parameter["Solver_parameters"]["time"] = T * nper
cpp_parameter["Solver_parameters"]["rho"] = rho
cpp_parameter["Solver_parameters"]["rel_tol_res_forces"] = rel_tol
cpp_parameter["Sensitivity"]["pose"] = False
cpp_parameter["Sensitivity"]["velocity"] = False
cpp_parameter["Sensitivity"]["period"] = False
json.dump(cpp_parameter, open(path2cpp + "_" + paramfile, "w"), indent=2)

# run C++ sim
cpprun = subprocess.run(
    "cd " + path2cpp + "&&" + exec_cpp + " _" + paramfile,
    shell=True,
    stdout=open(path2cpp + "cpp.out", "w"),
    stderr=open(path2cpp + "cpp.err", "w"),
)

# call timeplot script
if beam_type == "VK":
    subprocess.run("cd " + path2cpp + "&&" + "python3 timeplot_VK.py", shell=True)
else:
    subprocess.run("cd " + path2cpp + "&&" + "python3 timeplot.py", shell=True)
