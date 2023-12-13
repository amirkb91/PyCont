import h5py
import json
import subprocess
import numpy as np
from copy import deepcopy as dp
import os
import shutil
from concurrent.futures import ProcessPoolExecutor


class BeamCpp:
    # --------- Choose example case from mb_sef_cpp ---------#
    # cpp_example = "beam_2D"
    cpp_example = "beam_rightangle"

    if cpp_example == "beam_2D":
        # beam_2D (doubly clamped, arch, cantilever)
        cpp_path = "/home/akb110/Codes/mb_sef_cpp/examples/mybeam_2D/"
        cpp_exe = "/home/akb110/Codes/mb_sef_cpp/cmake-build-release/examples/mybeam_2D"
    elif cpp_example == "beam_rightangle":
        # beam_rightagnle
        cpp_path = "/home/akb110/Codes/mb_sef_cpp/examples/mybeam_rightangle/"
        cpp_exe = "/home/akb110/Codes/mb_sef_cpp/cmake-build-release/examples/mybeam_rightangle"
    # -------------------------------------------------------#

    n_core_parallel = 4
    cpp_modelfile = "model_def.json"
    cpp_paramfile_eig = "parameters_eig.json"
    cpp_paramfile_sim = "parameters_sim.json"
    model_def = json.load(open(cpp_path + cpp_modelfile))
    cpp_params_eig = json.load(open(cpp_path + cpp_paramfile_eig))
    cpp_params_sim = json.load(open(cpp_path + cpp_paramfile_sim))

    eig_file = cpp_params_eig["EigenSolverParameters"]["Logger"]["file_name"]
    simout_file = cpp_params_sim["TimeIntegrationSolverParameters"]["Logger"]["file_name"]
    ic_file = cpp_params_sim["TimeIntegrationSolverParameters"]["_initial_conditions"]["file_name"]
    model_name = model_def["ModelDef"]["model_name"]
    analysis_name = cpp_params_sim["TimeIntegrationSolverParameters"]["_initial_conditions"][
        "analysis_name"]

    free_dof = None
    fix_dof = None
    ndof_free = None
    ndof_fix = None
    ndof_all = None
    node_config = None
    ndof_config = None

    @classmethod
    def initialise(cls, cont_params, ForcePeriod=False):
        # prep para file and assign fixed values
        cls.cpp_params_sim["TimeIntegrationSolverParameters"]["rel_tol_res_forces"] = cont_params[
            "shooting"]["rel_tol"]
        cls.cpp_params_sim["TimeIntegrationSolverParameters"][
            "initial_conditions"] = cls.cpp_params_sim["TimeIntegrationSolverParameters"].pop(
                "_initial_conditions"
            )
        cls.cpp_params_sim["TimeIntegrationSolverParameters"].pop("_initial_conditions_eig")

        if not cont_params["continuation"]["forced"]:
            cls.model_def["ModelDef"]["amplitude"] = 0.0
            cls.model_def["ModelDef"]["tau0"] = 0.0
            cls.model_def["ModelDef"]["tau1"] = 0.0
            cls.cpp_params_sim["TimeIntegrationSolverParameters"]["rho"] = 1.0
        elif cont_params["continuation"]["forced"]:
            cls.model_def["ModelDef"]["amplitude"] = cont_params["forcing"]["amplitude"]
            cls.model_def["ModelDef"]["tau0"] = cont_params["forcing"]["tau0"]
            cls.model_def["ModelDef"]["tau1"] = cont_params["forcing"]["tau1"]
            cls.model_def["ModelDef"]["phase_ratio"] = cont_params["forcing"]["phase_ratio"]
            cls.model_def["ModelDef"]["def_period"] = 1.0
            cls.model_def["ModelDef"]["with_ratio"] = True
            cls.cpp_params_sim["TimeIntegrationSolverParameters"]["rho"] = cont_params["forcing"][
                "rho_GA"]

        if ForcePeriod:
            # if force period is specified, we want to de-couple force period from sim time
            # useful when running time simulation as post processing
            cls.model_def["ModelDef"]["def_period"] = ForcePeriod
            cls.model_def["ModelDef"]["with_ratio"] = False

        subprocess.run("cd " + cls.cpp_path + "&&" + "./clean_dir.sh", shell=True)
        json.dump(cls.model_def, open(cls.cpp_path + "_" + cls.cpp_modelfile, "w"), indent=2)

    @classmethod
    def run_eig(cls):
        subprocess.run(
            "cd " + cls.cpp_path + "&&" + cls.cpp_exe + " _" + cls.cpp_modelfile + " " +
            cls.cpp_paramfile_eig,
            shell=True,
            stdout=open(cls.cpp_path + "cpp.out", "w"),
            stderr=open(cls.cpp_path + "cpp.err", "w"),
        )
        eigdata = h5py.File(cls.cpp_path + cls.eig_file + ".h5", "r")
        eig = np.array(eigdata["/eigen_analysis/Eigenvectors/MOTION"])
        frq = np.array(eigdata["/eigen_analysis/Frequencies"])
        pose0 = np.array(eigdata["/eigen_analysis/POSE/MOTION"])[:, 0]
        cls.read_dofdata()

        return eig, frq, pose0

    @classmethod
    def runsim_single(
        cls, omega, tau, Xtilde, pose_base, cont_params, return_time=False, sensitivity=True
    ):
        nperiod = cont_params["shooting"]["single"]["nperiod"]
        nsteps = cont_params["shooting"]["single"]["nsteps_per_period"]
        N = cls.ndof_free

        T = tau / omega
        X = dp(Xtilde)
        X[N:] *= omega  # scale velocities from Xtilde to X

        cls.config_update(pose_base)
        cvg = cls.run_cpp(T * nperiod, X, nsteps * nperiod, sensitivity)
        if cvg:
            simdata = h5py.File(cls.cpp_path + cls.simout_file + ".h5", "r")
            energy = np.max(
                simdata["/dynamic_analysis/FEModel/energy"][:, :]
            )  # max energy important for forced response
            periodicity_inc = simdata["/dynamic_analysis/Periodicity/INC"][cls.free_dof]
            periodicity_vel = simdata["/dynamic_analysis/Periodicity/VELOCITY"][cls.free_dof]
            pose_time = simdata["/dynamic_analysis/FEModel/POSE/MOTION"][:]
            vel_time = simdata["/dynamic_analysis/FEModel/VELOCITY/MOTION"][:]
            # solution pose and vel taken from time 0 (initial values are those with inc and vel added)
            pose = pose_time[:, 0]
            vel = vel_time[:, 0]
            H = np.concatenate([periodicity_inc, periodicity_vel])
            if sensitivity:
                sens_H = simdata["/Sensitivity/Monodromy"][:]  # sens of periodicity function H wrt IC
                dHdtau = sens_H[:, -1] * nperiod * 1 / omega  # scale time derivative
                sens_H = np.delete(sens_H, -1, axis=1)
                sens_H[:, N:] *= omega  # scale velocity derivatives
                J = np.concatenate((sens_H, dHdtau.reshape(-1, 1)), axis=1)
                # if periodic solution is found, sens_H will equal Monodromy - eye.
                M = sens_H + np.eye(len(sens_H))
            else:
                J = M = None
            simdata.close()
        else:
            H = J = M = pose = vel = energy = None

        if not return_time:
            return H, J, M, pose, vel, energy, cvg
        elif return_time:
            return pose_time, vel_time

    @classmethod
    def runsim_multiple(cls, omega, tau, Xtilde, pose_base, cont_params):
        npartition = cont_params["shooting"]["multiple"]["npartition"]
        nsteps = cont_params["shooting"]["multiple"]["nsteps_per_partition"]
        N = cls.ndof_free
        twoN = 2 * N
        delta_S = 1 / npartition
        T = tau / omega

        # initialise
        J = np.zeros((npartition * twoN, npartition * twoN + 1))
        pose_time = np.zeros((cls.ndof_config, (nsteps + 1), npartition))
        vel_time = np.zeros((cls.ndof_all, (nsteps + 1), npartition))
        energy = np.zeros(npartition)
        cvg = [None] * npartition

        for ipart in range(npartition):
            # index values required for looping the partitions
            i = ipart * twoN
            i1 = (ipart + 1) * twoN
            j = (ipart + 1) % npartition * twoN
            j1 = ((ipart + 1) % npartition + 1) * twoN

            X = dp(Xtilde[i:i1])
            X[N:] *= omega  # scale velocities from Xtilde to X
            cls.config_update(pose_base[:, ipart])
            cvg[ipart] = cls.run_cpp(T * delta_S, X, nsteps)
            simdata = h5py.File(cls.cpp_path + cls.simout_file + ".h5", "r")
            M = simdata["/Sensitivity/Monodromy"][:]
            dHdtau = M[:, -1] * delta_S * 1 / omega  # scale time derivative
            M = np.delete(M, -1, axis=1)
            M[:, N:] *= omega  # scale velocity derivatives
            J[i:i1, i:i1] = M
            J[i:i1, j:j1] -= np.eye(twoN)
            J[i:i1, -1] = dHdtau
            pose_time[:, :, ipart] = simdata["/dynamic_analysis/FEModel/POSE/MOTION"][:]
            vel_time[:, :, ipart] = simdata["/dynamic_analysis/FEModel/VELOCITY/MOTION"][:]
            energy[ipart] = simdata["/dynamic_analysis/FEModel/energy"][:, -1][0]
            simdata.close()
        energy = np.mean(energy)
        cvg = all(cvg)

        # solution pose and vel taken from time 0
        pose = pose_time[:, 0, :]
        vel = vel_time[:, 0, :]

        # periodicity (to be put in seperate method)
        partition_order = (np.arange(npartition) + 1) % npartition
        H = np.array([])
        for ipart in range(npartition):
            h_pose = (
                pose_time[cls.free_dof, -1, ipart] -
                pose_time[cls.free_dof, 0, partition_order[ipart]]
            )
            h_vel = (
                vel_time[cls.free_dof, -1, ipart] -
                vel_time[cls.free_dof, 0, partition_order[ipart]]
            )
            H = np.append(H, np.concatenate([h_pose, h_vel]))
        H = H.reshape(-1, 1)

        return H, J, pose, vel, energy, cvg

    @classmethod
    def run_cpp(cls, T, X, nsteps, sensitivity):
        inc = np.zeros(cls.ndof_all)
        vel = np.zeros(cls.ndof_all)
        inc[cls.free_dof] = X[:cls.ndof_free]
        vel[cls.free_dof] = X[cls.ndof_free:]

        icdata = h5py.File(cls.cpp_path + cls.ic_file + ".h5", "a")
        icdata["/" + cls.analysis_name + "/FEModel/INC/MOTION"] = inc.reshape(-1, 1)
        icdata["/" + cls.analysis_name + "/FEModel/VELOCITY/MOTION"] = vel.reshape(-1, 1)
        icdata.close()

        cls.cpp_params_sim["TimeIntegrationSolverParameters"]["number_of_steps"] = nsteps
        cls.cpp_params_sim["TimeIntegrationSolverParameters"]["time"] = T
        if not sensitivity:
            cls.cpp_params_sim["TimeIntegrationSolverParameters"].pop("direct_sensitivity")
        if cls.n_core_parallel == 1:
            json.dump(cls.cpp_params_sim, open(cls.cpp_path + "_" + cls.cpp_paramfile_sim, "w"), indent=2)
            try:
                cpprun = subprocess.run(
                    "cd " + cls.cpp_path + "&&" + cls.cpp_exe + " _" + cls.cpp_modelfile + " _" +
                    cls.cpp_paramfile_sim,
                    shell=True,
                    # timeout=30,
                    stdout=open(cls.cpp_path + "cpp.out", "w"),
                    stderr=open(cls.cpp_path + "cpp.err", "w"),
                )
                cvg = not bool(cpprun.returncode)
            except subprocess.TimeoutExpired:
                print("C++ code timed out ------- ", end="")
                cvg = False
                os.remove(cls.cpp_path + cls.simout_file + ".h5")
        elif sensitivity and cls.n_core_parallel > 1:
            # Calculate the basic split size and the number of splits that need an extra column
            basic_split_size, extra_splits = divmod(cls.ndof_free, cls.n_core_parallel)
            start_indices = np.arange(cls.n_core_parallel) * basic_split_size + np.minimum(np.arange(cls.n_core_parallel), extra_splits)
            end_indices = np.roll(start_indices, -1)
            end_indices[-1] = cls.ndof_free  # Correct the last end index        
            split_indices = zip(start_indices, end_indices)        
        
            with ProcessPoolExecutor() as executor:
                convergence = list(executor.map(cls.run_cpp_parallel, zip(split_indices, range(1, cls.n_core_parallel + 1))))
            cvg = np.all(convergence)

        return cvg
    
    @classmethod
    def run_cpp_parallel(cls, combined_args):
        (start_index, end_index), run_id = combined_args

        requested_cols = np.array([start_index, end_index]).tolist() 
        suffix = f"_{run_id:03d}"
        cls.cpp_params_sim["TimeIntegrationSolverParameters"]["direct_sensitivity"]["requested_cols"] = requested_cols  
        cls.cpp_params_sim["TimeIntegrationSolverParameters"]["Logger"]["file_name"] = cls.simout_file + suffix 
        cls.cpp_params_sim["TimeIntegrationSolverParameters"]["initial_conditions"]["file_name"] = cls.ic_file + suffix 
        cls.model_def["ModelDef"]["model_name"] = cls.model_name + suffix
        json.dump(cls.cpp_params_sim, open(cls.cpp_path + "_" + cls.cpp_paramfile_sim.split('.')[0] + suffix + ".json", "w"), indent=2)
        json.dump(cls.model_def, open(cls.cpp_path + "_" + cls.cpp_modelfile.split('.')[0] + suffix + ".json", "w"), indent=2)
        shutil.copyfile(cls.cpp_path + cls.ic_file + ".h5", cls.cpp_path + cls.ic_file + suffix + ".h5")

        cpprun = subprocess.run("cd " + cls.cpp_path + "&&" + cls.cpp_exe + " _" + cls.cpp_modelfile.split('.')[0] + suffix + ".json" + " _" + 
                                cls.cpp_paramfile_sim.split('.')[0] + suffix + ".json",
                                shell=True,
                                stdout=open(cls.cpp_path + "cpp" + suffix + ".out", "w"),
                                stderr=open(cls.cpp_path + "cpp" + suffix + ".err", "w"),)
        return not bool(cpprun.returncode)

    @classmethod
    def partition_singleshooting_solution(cls, omega, tau, Xtilde, pose_base, cont_params):
        npartition = cont_params["shooting"]["multiple"]["npartition"]
        nsteps = cont_params["shooting"]["multiple"]["nsteps_per_partition"]
        rel_tol = cont_params["shooting"]["rel_tol"]
        N = cls.ndof_free
        slicing_index = nsteps * np.arange(npartition)

        T = tau / omega
        X = dp(Xtilde)
        X[N:] *= omega  # scale velocities from Xtilde to X

        cls.config_update(pose_base)
        # do time integration along whole orbit before slicing.
        # run nsteps per partition to ensure slicing is done at correct indices
        cvg = cls.run_cpp(T, X, nsteps * npartition, rel_tol)
        simdata = h5py.File(cls.cpp_path + cls.simout_file + ".h5", "r")
        pose_time = simdata["/dynamic_analysis/FEModel/POSE/MOTION"][:]
        vel_time = simdata["/dynamic_analysis/FEModel/VELOCITY/MOTION"][:]

        pose = pose_time[:, slicing_index]
        V = vel_time[cls.free_dof][:, slicing_index]
        # set inc to zero as solution stored in pose, keep velocity but scale first
        V *= 1 / omega
        X_out = np.concatenate((np.zeros((cls.ndof_free, npartition)), V))
        X_out = np.reshape(X_out, (-1), order="F")
        return X_out, pose

    @classmethod
    def config_update(cls, pose):
        # update beam configuration by writing initial conditions pose
        icdata = h5py.File(cls.cpp_path + cls.ic_file + ".h5", "w")
        icdata["/" + cls.analysis_name + "/FEModel/POSE/MOTION"] = pose.reshape(-1, 1)
        icdata.close()

    @classmethod
    def read_dofdata(cls):
        data = h5py.File(cls.cpp_path + cls.model_name + ".h5", "r")
        cls.free_dof = np.array(data["/FEModel/loc_dof_free/MOTION"])[:, 0]
        cls.fix_dof = np.array(data["/FEModel/loc_dof_fix/MOTION"])[:, 0]
        NodalFrame = list(data["/FEModel/Nodes_config/"].keys())[0]
        cls.node_config = np.array(data["/FEModel/Nodes_config/" + NodalFrame])[1:, :]
        cls.ndof_free = len(cls.free_dof)
        cls.ndof_fix = len(cls.fix_dof)
        cls.ndof_config = np.size(cls.node_config)
        cls.ndof_all = cls.ndof_free + cls.ndof_fix
        data.close()

    @classmethod
    def get_dofdata(cls):
        return {
            "free_dof": cls.free_dof,
            "ndof_all": cls.ndof_all,
            "ndof_fix": cls.ndof_fix,
            "ndof_free": cls.ndof_free,
            "node_config": cls.node_config,
            "ndof_config": cls.ndof_config,
        }

    # @classmethod
    # def periodicity(cls, pose, vel, target):
    #     if len(pose) == cls.ndof_all:
    #         # VK formulation
    #         posevel = np.concatenate([pose[cls.free_dof], vel[cls.free_dof]])
    #         H = posevel - target
    #     else:
    #         # SE formulation
    #         H = None
    #     return H
