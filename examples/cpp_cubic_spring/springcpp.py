import h5py
import json
import subprocess
import numpy as np
from copy import deepcopy as dp
import os


class SpringCpp:
    cpp_path = "/home/akb110/Codes/mb_sef_cpp/examples/my_cubic_spring/"
    cpp_exe = "/home/akb110/Codes/mb_sef_cpp/cmake-build-release/examples/my_cubic_spring"
    cpp_modelfile = "model_def.json"
    cpp_paramfile_eig = "parameters_eig.json"
    cpp_paramfile_sim = "parameters_sim.json"

    model_def = json.load(open(cpp_path + cpp_modelfile))
    cpp_params_eig = json.load(open(cpp_path + cpp_paramfile_eig))
    cpp_params_sim = json.load(open(cpp_path + cpp_paramfile_sim))

    model_name = model_def["ModelDef"]["model_name"]
    eig_file = cpp_params_eig["EigenSolverParameters"]["Logger"]["file_name"]
    simout_file = cpp_params_sim["TimeIntegrationSolverParameters"]["Logger"]["file_name"]
    ic_file = cpp_params_sim["TimeIntegrationSolverParameters"]["_initial_conditions"]["file_name"]
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
    def run_eig(cls):
        subprocess.run(
            "cd " + cls.cpp_path + "&&" + "./clean_dir.sh" + "&&" + cls.cpp_exe + " " +
            cls.cpp_modelfile + " " + cls.cpp_paramfile_eig,
            shell=True,
            stdout=open(cls.cpp_path + "cpp.out", "w"),
            stderr=open(cls.cpp_path + "cpp.err", "w")
        )
        eigdata = h5py.File(cls.cpp_path + cls.eig_file + ".h5", "r")
        eig = np.array(eigdata["/eigen_analysis/Eigenvectors/MOTION"])
        frq = np.array(eigdata["/eigen_analysis/Frequencies"])
        pose0 = np.array(eigdata["/eigen_analysis/POSE/MOTION"])[:, 0]
        cls.read_dofdata()

        return eig, frq, pose0

    @classmethod
    def runsim_single(cls, omega, tau, Xtilde, pose_base, cont_params, return_time=False):
        nperiod = cont_params["shooting"]["single"]["nperiod"]
        nsteps = cont_params["shooting"]["single"]["nsteps_per_period"]
        rel_tol = cont_params["shooting"]["rel_tol"]
        N = cls.ndof_free

        T = tau / omega
        X = dp(Xtilde)
        X[N:] *= omega  # scale velocities from Xtilde to X

        cls.config_update(pose_base)
        cvg = cls.run_cpp(T * nperiod, X, nsteps * nperiod, rel_tol)
        if cvg:
            simdata = h5py.File(cls.cpp_path + cls.simout_file + ".h5", "r")
            energy = simdata["/dynamic_analysis/FEModel/energy"][:, -1][0]
            periodicity_inc = simdata["/dynamic_analysis/Periodicity/INC"][cls.free_dof]
            periodicity_vel = simdata["/dynamic_analysis/Periodicity/VELOCITY"][cls.free_dof]
            # solution pose and vel taken from time 0 (initial values are those with inc and vel added)
            pose = simdata["/dynamic_analysis/FEModel/POSE/MOTION"][:, 0]
            vel = simdata["/dynamic_analysis/FEModel/VELOCITY/MOTION"][:, 0]
            pose_time = simdata["/dynamic_analysis/FEModel/POSE/MOTION"][:]
            vel_time = simdata["/dynamic_analysis/FEModel/VELOCITY/MOTION"][:]
            H = np.concatenate([periodicity_inc, periodicity_vel])
            M = simdata["/Sensitivity/Monodromy"][:]
            dHdtau = M[:, -1] * nperiod * 1 / omega  # scale time derivative
            M = np.delete(M, -1, axis=1)
            M[:, N:] *= omega  # scale velocity derivatives
            M -= np.eye(len(M))
            J = np.concatenate((M, dHdtau.reshape(-1, 1)), axis=1)
            simdata.close()
        else:
            H = J = pose = vel = energy = None

        if not return_time:
            return H, J, pose, vel, energy, cvg
        elif return_time:
            return pose_time, vel_time

    @classmethod
    def runsim_forced(cls, omega, tau, Xtilde, pose_base, cont_params, return_time=False):
        nsteps = cont_params["shooting"]["single"]["nsteps_per_period"]
        rel_tol = cont_params["shooting"]["rel_tol"]
        ga_rho_forced = cont_params["shooting"]["ga_rho_forced"]
        amplitude = cont_params["forcing"]["amplitude"]
        phase_ratio = cont_params["forcing"]["phase_ratio"]
        damping = cont_params["forcing"]["damping"]
        N = cls.ndof_free

        T = tau / omega
        X = dp(Xtilde)

        cls.config_update(pose_base)
        cvg = cls.run_cpp_forced(
            T, X, amplitude, phase_ratio, damping, nsteps, rel_tol, ga_rho_forced
        )
        if cvg:
            simdata = h5py.File(cls.cpp_path + cls.simout_file + ".h5", "r")
            energy = np.max(simdata["/dynamic_analysis/FEModel/energy"][:])
            periodicity_inc = simdata["/dynamic_analysis/Periodicity/INC"][cls.free_dof]
            periodicity_vel = simdata["/dynamic_analysis/Periodicity/VELOCITY"][cls.free_dof]
            # solution pose and vel taken from time 0 (initial values are those with inc and vel added)
            pose = simdata["/dynamic_analysis/FEModel/POSE/MOTION"][:, 0]
            vel = simdata["/dynamic_analysis/FEModel/VELOCITY/MOTION"][:, 0]
            pose_time = simdata["/dynamic_analysis/FEModel/POSE/MOTION"][:]
            vel_time = simdata["/dynamic_analysis/FEModel/VELOCITY/MOTION"][:]
            H = np.concatenate([periodicity_inc, periodicity_vel])
            M = simdata["/Sensitivity/Monodromy"][:]
            dHdtau = M[:, -1]
            M = np.delete(M, -1, axis=1)
            M -= np.eye(len(M))
            J = np.concatenate((M, dHdtau.reshape(-1, 1)), axis=1)
            simdata.close()
        else:
            H = J = pose = vel = energy = None

        if not return_time:
            return H, J, pose, vel, energy, cvg
        elif return_time:
            return pose_time, vel_time

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

        cls.cpp_params_sim["TimeIntegrationSolverParameters"]["number_of_steps"] = nsteps
        cls.cpp_params_sim["TimeIntegrationSolverParameters"]["time"] = T
        cls.cpp_params_sim["TimeIntegrationSolverParameters"]["rel_tol_res_forces"] = rel_tol
        cls.cpp_params_sim["TimeIntegrationSolverParameters"]["initial_conditions"] = cls.cpp_params_sim[
            "TimeIntegrationSolverParameters"]["_initial_conditions"]
        json.dump(cls.cpp_params_sim, open(cls.cpp_path + "_" + cls.cpp_paramfile_sim, "w"), indent=2)

        try:
            cpprun = subprocess.run(
                "cd " + cls.cpp_path + "&&" + cls.cpp_exe + " " + cls.cpp_modelfile + " _" +
                cls.cpp_paramfile_sim,
                shell=True,
                timeout=10,
                stdout=open(cls.cpp_path + "cpp.out", "w"),
                stderr=open(cls.cpp_path + "cpp.err", "w")
            )
            cvg = not bool(cpprun.returncode)
        except subprocess.TimeoutExpired:
            print("C++ code timed out ------- ", end="")
            cvg = False
            os.remove(cls.cpp_path + cls.simout_file + ".h5")

        return cvg

    @classmethod
    def run_cpp_forced(cls, T, X, amplitude, phase_ratio, damping, nsteps, rel_tol, ga_rho_forced):
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
        cls.cpp_params_sim["TimeIntegrationSolverParameters"]["rel_tol_res_forces"] = rel_tol
        cls.cpp_params_sim["TimeIntegrationSolverParameters"]["rho"] = ga_rho_forced
        cls.cpp_params_sim["ForcingParameters"]["amplitude"] = amplitude
        cls.cpp_params_sim["ForcingParameters"]["phase_ratio"] = phase_ratio
        cls.cpp_params_sim["ModelDef"]["tau"] = damping
        cls.cpp_params_sim["TimeIntegrationSolverParameters"]["initial_conditions"] = \
            cls.cpp_params_sim["TimeIntegrationSolverParameters"]["_initial_conditions"]
        json.dump(cls.cpp_params_sim, open(cls.cpp_path + "_" + cls.cpp_paramfile_sim, "w"), indent=2)

        try:
            cpprun = subprocess.run(
                "cd " + cls.cpp_path + "&&" + cls.cppfrc_exe + " _" + cls.cpp_paramfile_sim,
                shell=True,
                timeout=10,
                stdout=open(cls.cpp_path + "cpp.out", "w"),
                stderr=open(cls.cpp_path + "cpp.err", "w")
            )
            cvg = not bool(cpprun.returncode)
        except subprocess.TimeoutExpired:
            print("C++ code timed out ------- ", end="")
            cvg = False
            os.remove(cls.cpp_path + cls.simout_file + ".h5")

        return cvg

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
        cls.fix_dof = np.array([])
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
