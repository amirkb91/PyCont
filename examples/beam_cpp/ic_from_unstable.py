import h5py
import json
import sys
import numpy as np
from beamcpp import BeamCpp
from postprocess.bifurcation import get_unstable_eigenvec
from Frame import Frame

""" Generate initial condition file to restart from by adding an unstable eigenvector """

scale = 2.0
print("\033[0;31m*** ENSURE apply_SE_correction = FALSE *** \033[0m\n")

solution_number = int(input("Enter solution number to restart from: "))

# read solution file
file = sys.argv[1]
if not file.endswith(".h5"):
    file += ".h5"
data = h5py.File(str(file), "r")
pose = data["/Config/POSE"][:, solution_number]
vel = data["/Config/VELOCITY"][:, solution_number]
T = data["/T"][solution_number]
par = data["/Parameters"]
par = json.loads(par[()])
data.close()

# run sim to get monodromy
BeamCpp.initialise(par)
BeamCpp.run_eig()  # To get nodal data in class

x = vel[BeamCpp.free_dof]
X = np.concatenate([np.zeros(BeamCpp.ndof_free), x])

[_, J, _, _, _, _, _] = BeamCpp.runsim_single(1.0, T, X, pose, par, fulltime=True)
M = J[:, :-1] + np.eye(2 * BeamCpp.ndof_free)

unstable_eigvec = get_unstable_eigenvec(M)
inc_eig = np.zeros(BeamCpp.ndof_all)
vel_eig = np.zeros(BeamCpp.ndof_all)
inc_eig[BeamCpp.free_dof] = unstable_eigvec[: BeamCpp.ndof_free].real
vel_eig[BeamCpp.free_dof] = unstable_eigvec[BeamCpp.ndof_free :].real


# add unstable eigenvector to solution
vel_ic = vel + scale * vel_eig
pose_ic = np.zeros_like(pose)

for i in range(BeamCpp.nnodes_all):
    f = Frame.get_frame_from_parameters(
        BeamCpp.n_dim,
        scale * inc_eig[i * BeamCpp.dof_per_node : (i + 1) * BeamCpp.dof_per_node],
    )
    pose_ic[i * BeamCpp.config_per_node : (i + 1) * BeamCpp.config_per_node] = Frame.composition(
        BeamCpp.n_dim, pose[i * BeamCpp.config_per_node : (i + 1) * BeamCpp.config_per_node], f
    )

# write to file to restart from
out_file = h5py.File("IC_unstable.h5", "w")
out_file["/Config/POSE"] = np.asarray(pose_ic.reshape(-1, 1))
out_file["/Config/VELOCITY"] = np.asarray(vel_ic.reshape(-1, 1))
out_file["/T"] = np.asarray(np.array([T]))
out_file["/Parameters"] = json.dumps(par)
out_file.close()


# out_file = h5py.File("IC_unstable_only.h5", "w")
# out_file["/dynamic_analysis/FEModel/INC/MOTION"] = np.asarray(scale*inc_eig.reshape(-1, 1))
# out_file["/dynamic_analysis/FEModel/VELOCITY/MOTION"] = np.asarray(0*scale*vel_eig.reshape(-1, 1))
# out_file.close()