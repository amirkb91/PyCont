import h5py
import json
import subprocess
import numpy as np
import os


class SpringCpp:
    cpp_path = "/home/akb110/Codes/mb_sef_cpp/examples/mycubic_spring/"
    cpp_exe = "/home/akb110/Codes/mb_sef_cpp/cmake-build-release/examples/mycubic_spring"
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
        "analysis_name"
    ]

    free_dof = None
    fix_dof = None
    ndof_free = None
    ndof_fix = None
    ndof_all = None
    node_config = None
    ndof_config = None

    @classmethod
    def initialise(cls, cont_params):
        if not cont_params["continuation"]["forced"]:
            cls.model_def["ModelDef"]["amplitude"] = 0.0
            cls.model_def["ModelDef"]["damping_M"] = 0.0
            cls.cpp_params_sim["TimeIntegrationSolverParameters"]["rho"] = 1.0
        elif cont_params["continuation"]["forced"]:
            cls.model_def["ModelDef"]["amplitude"] = cont_params["forcing"]["amplitude"]
            cls.model_def["ModelDef"]["phase_ratio"] = cont_params["forcing"]["phase_ratio"]
            cls.model_def["ModelDef"]["damping_M"] = cont_params["forcing"]["damping"]
            cls.cpp_params_sim["TimeIntegrationSolverParameters"]["rho"] = cont_params["forcing"][
                "rho_GA"
            ]

        subprocess.run("cd " + cls.cpp_path + "&&" + "./clean_dir.sh", shell=True)
        json.dump(cls.model_def, open(cls.cpp_path + "_" + cls.cpp_modelfile, "w"), indent=2)

    @classmethod
    def run_eig(cls):
        subprocess.run(
            "cd "
            + cls.cpp_path
            + "&&"
            + cls.cpp_exe
            + " _"
            + cls.cpp_modelfile
            + " "
            + cls.cpp_paramfile_eig,
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
    def runsim_single(cls, omega, tau, Xtilde, pose_base, cont_params, return_time=False):
        nperiod = cont_params["shooting"]["single"]["nperiod"]
        nsteps = cont_params["shooting"]["single"]["nsteps_per_period"]
        rel_tol = cont_params["shooting"]["rel_tol"]
        N = cls.ndof_free

        T = tau / omega

        pose_sim = Xtilde[:N].copy() + pose_base.copy()
        vel_sim = Xtilde[N:] * omega  # scale velocities
        X = np.concatenate([pose_sim, vel_sim])

        cvg = cls.run_cpp(T * nperiod, X, nsteps * nperiod, rel_tol)
        if cvg:
            simdata = h5py.File(cls.cpp_path + cls.simout_file + ".h5", "r")
            energy = simdata["/dynamic_analysis/FEModel/energy"][:, -1][0]
            periodicity_inc = simdata["/dynamic_analysis/Periodicity/INC"][cls.free_dof]
            periodicity_vel = simdata["/dynamic_analysis/Periodicity/VELOCITY"][cls.free_dof]
            pose_time = simdata["/dynamic_analysis/FEModel/POSE/MOTION"][:]
            vel_time = simdata["/dynamic_analysis/FEModel/VELOCITY/MOTION"][:]
            # solution pose and vel taken from time 0
            pose = pose_time[:, 0]
            vel = vel_time[:, 0]
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
    def run_cpp(cls, T, X, nsteps, rel_tol):
        pose = np.zeros(cls.ndof_all)
        vel = np.zeros(cls.ndof_all)
        pose[cls.free_dof] = X[: cls.ndof_free]
        vel[cls.free_dof] = X[cls.ndof_free :]

        icdata = h5py.File(cls.cpp_path + cls.ic_file + ".h5", "w")
        icdata["/" + cls.analysis_name + "/FEModel/POSE/MOTION"] = pose.reshape(-1, 1)
        icdata["/" + cls.analysis_name + "/FEModel/VELOCITY/MOTION"] = vel.reshape(-1, 1)
        icdata.close()

        cls.cpp_params_sim["TimeIntegrationSolverParameters"]["number_of_steps"] = nsteps
        cls.cpp_params_sim["TimeIntegrationSolverParameters"]["time"] = T
        cls.cpp_params_sim["TimeIntegrationSolverParameters"]["rel_tol_res_forces"] = rel_tol
        cls.cpp_params_sim["TimeIntegrationSolverParameters"]["initial_conditions"] = (
            cls.cpp_params_sim["TimeIntegrationSolverParameters"]["_initial_conditions"]
        )
        json.dump(
            cls.cpp_params_sim, open(cls.cpp_path + "_" + cls.cpp_paramfile_sim, "w"), indent=2
        )

        try:
            cpprun = subprocess.run(
                "cd "
                + cls.cpp_path
                + "&&"
                + cls.cpp_exe
                + " _"
                + cls.cpp_modelfile
                + " _"
                + cls.cpp_paramfile_sim,
                shell=True,
                stdout=open(cls.cpp_path + "cpp.out", "w"),
                stderr=open(cls.cpp_path + "cpp.err", "w"),
            )
            cvg = not bool(cpprun.returncode)
        except subprocess.TimeoutExpired:
            print("C++ code timed out ------- ", end="")
            cvg = False
            os.remove(cls.cpp_path + cls.simout_file + ".h5")

        return cvg

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
            "ndof_config": cls.ndof_config,
        }
