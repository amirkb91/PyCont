import h5py
import json
import subprocess
import numpy as np
import sys


class BeamCpp:
    cppeig_exe = "/home/akb110/Codes/mb_sef_cpp/cmake-build-release/examples/mybeam_2D_eig"
    cppsim_exe = "/home/akb110/Codes/mb_sef_cpp/cmake-build-release/examples/mybeam_2D_sim"
    cpp_path = "/home/akb110/Codes/mb_sef_cpp/examples/mybeam_2D/"
    cpp_paramfile = "parameters.json"

    cpp_params = json.load(open(cpp_path + cpp_paramfile))
    model_name = cpp_params["ModelDef"]["model_name"]
    eig_file = cpp_params["EigenSolverParameters"]["Logger"]["file_name"]
    simout_file = cpp_params["TimeIntegrationSolverParameters"]["Logger"]["file_name"]
    ic_file = cpp_params["TimeIntegrationSolverParameters"]["_initial_conditions"]["file_name"]
    analysis_name = cpp_params["TimeIntegrationSolverParameters"]["_initial_conditions"][
        "analysis_name"]

    free_dof = None
    fix_dof = None
    ndof_free = None
    ndof_fix = None
    ndof_all = None
    node_config = None
    ndof_config = None

    @classmethod
    def run_eig(cls, cont_params):
        subprocess.run(
            "cd " + cls.cpp_path + "&&" + "./clean_dir.sh" + "&&" + cls.cppeig_exe + " " +
            cls.cpp_paramfile,
            shell=True,
            stdout=open(cls.cpp_path + "cpp.out", "w"),
            stderr=open(cls.cpp_path + "cpp.err", "w"))
        eigdata = h5py.File(cls.cpp_path + cls.eig_file + ".h5", "r")
        cls.read_dofdata()

        eig = np.array(eigdata["/eigen_analysis/Eigenvectors/MOTION"])
        frq = eigdata["/eigen_analysis/Frequencies"]
        # eig[np.abs(eig) < 1e-10] = 0.0
        nnm = cont_params["first_point"]["eig_start"]["NNM"]
        scale = cont_params["first_point"]["eig_start"]["scale"]
        x0 = scale * eig[:, nnm - 1]
        x0 = x0[cls.free_dof]
        v0 = np.zeros_like(x0)
        X0 = np.concatenate([x0, v0])
        T0 = 1 / frq[nnm - 1, 0]
        pose0 = np.array(eigdata["/eigen_analysis/POSE/MOTION"])[:, 0]

        return X0, T0, pose0

    @classmethod
    def runsim_single(cls, T, X, pose_base, cont_params):
        nperiod = cont_params["shooting"]["single"]["nperiod"]
        nsteps = cont_params["shooting"]["single"]["nsteps_per_period"]
        rel_tol = cont_params["shooting"]["rel_tol"]

        cls.config_update(pose_base)
        simdata, cvg = cls.run_cpp(T * nperiod, X, nsteps * nperiod, rel_tol)
        if cvg:
            energy = simdata["/dynamic_analysis/FEModel/energy"][:, -1][0]
            periodicity_inc = simdata["/dynamic_analysis/Periodicity/INC"][cls.free_dof]
            periodicity_vel = simdata["/dynamic_analysis/Periodicity/VELOCITY"][cls.free_dof]
            # solution pose and vel taken from time 0
            pose = simdata["/dynamic_analysis/FEModel/POSE/MOTION"][:, 0]
            vel = simdata["/dynamic_analysis/FEModel/VELOCITY/MOTION"][:, 0]
            H = np.concatenate([periodicity_inc, periodicity_vel])
            M = simdata["/Sensitivity/Monodromy"][:]
            dHdt = M[:, -1] * nperiod
            M = np.delete(M, -1, axis=1)
            M -= np.eye(len(M))
            J = np.concatenate((M, dHdt.reshape(-1, 1)), axis=1)
            simdata.close()
        else:
            H = J = pose = vel = energy = None

        return H, J, pose, vel, energy, cvg

    @classmethod
    def runsim_multiple(cls, T, X, pose_base, cont_params):
        npartition = cont_params["shooting"]["multiple"]["npartition"]
        nsteps = cont_params["shooting"]["multiple"]["nsteps_per_partition"]
        rel_tol = cont_params["shooting"]["rel_tol"]
        N = cls.ndof_free
        twoN = 2 * N
        delta_S = 1 / npartition

        # initialise
        J = np.zeros((npartition * twoN, npartition * twoN + 1))
        pose_time = np.zeros((cls.ndof_config, (nsteps + 1), npartition))
        vel_time = np.zeros((cls.ndof_all, (nsteps + 1), npartition))
        energy = np.zeros(npartition)
        cvg = [None] * npartition

        for ipart in range(npartition):
            i = ipart * twoN
            i1 = (ipart + 1) * twoN
            j = (ipart + 1) % npartition * twoN
            j1 = ((ipart + 1) % npartition + 1) * twoN

            cls.config_update(pose_base[:, ipart])
            simdata, cvg[ipart] = cls.run_cpp(T * delta_S, X[i:i1], nsteps, rel_tol)
            M = simdata["/Sensitivity/Monodromy"][:]
            dHdt = M[:, -1] * delta_S
            M = np.delete(M, -1, axis=1)
            J[i:i1, i:i1] = M
            J[i:i1, j:j1] -= np.eye(twoN)
            J[i:i1, -1] = dHdt
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
            h_pose = pose_time[cls.free_dof, -1, ipart] - pose_time[cls.free_dof, 0, partition_order[ipart]]
            h_vel = vel_time[cls.free_dof, -1, ipart] - vel_time[cls.free_dof, 0, partition_order[ipart]]
            H = np.append(H, np.concatenate([h_pose, h_vel]))
        H = H.reshape(-1, 1)

        return H, J, pose, vel, energy, cvg

    @classmethod
    def run_cpp(cls, T, X, nsteps, rel_tol):
        inc = np.zeros(cls.ndof_all)
        vel = np.zeros(cls.ndof_all)
        inc[cls.free_dof] = X[:cls.ndof_free]
        vel[cls.free_dof] = X[cls.ndof_free:]

        icdata = h5py.File(cls.cpp_path + cls.ic_file + ".h5", "a")
        icdata["/" + cls.analysis_name + "/FEModel/INC/MOTION"] = inc.reshape(-1, 1)
        icdata["/" + cls.analysis_name + "/FEModel/VELOCITY/MOTION"] = vel.reshape(-1, 1)
        icdata.close()

        cls.cpp_params["TimeIntegrationSolverParameters"]["number_of_steps"] = nsteps
        cls.cpp_params["TimeIntegrationSolverParameters"]["time"] = T
        cls.cpp_params["TimeIntegrationSolverParameters"]["rel_tol_res_forces"] = rel_tol
        cls.cpp_params["TimeIntegrationSolverParameters"]["initial_conditions"] = \
            cls.cpp_params["TimeIntegrationSolverParameters"]["_initial_conditions"]
        json.dump(cls.cpp_params, open(cls.cpp_path + "_" + cls.cpp_paramfile, "w"), indent=2)

        try:
            cpprun = subprocess.run(
                "cd " + cls.cpp_path + "&&" + cls.cppsim_exe + " _" + cls.cpp_paramfile,
                shell=True,
                timeout=300,
                stdout=open(cls.cpp_path + "cpp.out", "w"),
                stderr=open(cls.cpp_path + "cpp.err", "w"))
        except subprocess.TimeoutExpired:
            sys.exit("C++ code timed out.")
        simdata = h5py.File(cls.cpp_path + cls.simout_file + ".h5", "r")
        cvg = not bool(cpprun.returncode)
        return simdata, cvg

    @classmethod
    def partition_singleshooting_solution(cls, T, X, pose_base, cont_params):
        npartition = cont_params["shooting"]["multiple"]["npartition"]
        nsteps = cont_params["shooting"]["multiple"]["nsteps_per_partition"]
        rel_tol = cont_params["shooting"]["rel_tol"]
        slicing_index = nsteps * np.arange(npartition)

        cls.config_update(pose_base)
        # do time integration along whole orbit before slicing
        simdata, cvg = cls.run_cpp(T, X, nsteps * npartition, rel_tol)
        pose_time = simdata["/dynamic_analysis/FEModel/POSE/MOTION"][:]
        vel_time = simdata["/dynamic_analysis/FEModel/VELOCITY/MOTION"][:]

        pose = pose_time[:, slicing_index]
        V = vel_time[cls.free_dof][:, slicing_index]
        # set inc to zero as solution stored in pose, keep velocity
        X_out = np.concatenate((np.zeros((cls.ndof_free, npartition)), V))
        X_out = np.reshape(X_out, (-1), order='F')
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
            "ndof_config": cls.ndof_config
        }

    @classmethod
    def periodicity(cls, pose, vel, target):
        if len(pose) == cls.ndof_all:
            # VK formulation
            posevel = np.concatenate([pose[cls.free_dof], vel[cls.free_dof]])
            H = posevel - target
        else:
            # SE formulation
            H = None

        return H
