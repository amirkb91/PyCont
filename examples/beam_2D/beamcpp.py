import h5py
import json
import subprocess
import numpy as np


class BeamCpp:
    cppeig_exe = "/home/akb110/Codes/mb_sef_cpp/cmake-build-release/examples/beam_eig"
    cppsim_exe = "/home/akb110/Codes/mb_sef_cpp/cmake-build-release/examples/beam_sim"
    cpp_path = "/home/akb110/Codes/mb_sef_cpp/examples/beam_2D/"
    cpp_paramfile = "parameters.json"

    cpp_params = json.load(open(cpp_path + cpp_paramfile))
    eig_file = cpp_params["EigenSolverParameters"]["Logger"]["file_name"]
    simout_file = cpp_params["TimeIntegrationSolverParameters"]["Logger"]["file_name"]
    ic_file = cpp_params["TimeIntegrationSolverParameters"]["_initial_conditions"]["file_name"]

    free_dof = None
    ndof_all = None
    ndof_fix = None
    ndof_free = None

    @classmethod
    def run_eig(cls, cont_params):
        subprocess.run("cd " + cls.cpp_path + "&&" + "./clean_dir.sh" + "&&" +
                       cls.cppeig_exe + " " + cls.cpp_paramfile + "&&" +
                       cls.cppsim_exe + " " + cls.cpp_paramfile,
                       shell=True,
                       stdout=open(cls.cpp_path + "cpp.out", "w"),
                       stderr=open(cls.cpp_path + "cpp.err", "w"))
        eigdata = h5py.File(cls.cpp_path + cls.eig_file + ".h5", "r")
        simdata = h5py.File(cls.cpp_path + cls.simout_file + ".h5", "r")
        cls.set_dofdata(simdata)
        simdata.close()

        eig = np.array(eigdata["/eigen_analysis/Eigenvectors"])
        frq = eigdata["/eigen_analysis/Frequencies"]
        eig[np.abs(eig) < 1e-10] = 0.0
        nnm = cont_params["first_point"]["eig_start"]["NNM"]
        scale = cont_params["first_point"]["eig_start"]["scale"]
        x0 = scale * eig[:, nnm - 1]
        x0 = x0[cls.free_dof]
        v0 = np.zeros_like(x0)
        X0 = np.concatenate([x0, v0])
        T0 = 1 / frq[nnm - 1, 0]
        pose_base0 = np.array(eigdata["/eigen_analysis/Config_ref"])[:, 0]

        return X0, T0, pose_base0

    @classmethod
    def runsim_multiple(cls, T, X, pose_base, cont_params):
        npartition = cont_params["shooting"]["multiple"]["npartition"]
        nsteps = cont_params["shooting"]["multiple"]["nsteps_per_partition"]
        rel_tol = cont_params["shooting"]["rel_tol"]
        delta_S = 1 / npartition
        timesol_partition_index_start = nsteps * np.arange(npartition) + np.arange(npartition)
        timesol_partition_index_end = timesol_partition_index_start - 1
        block_order = (np.arange(npartition) + 1) % npartition

        N = cls.ndof_free
        twoN = 2 * N
        J = np.zeros((npartition * twoN, npartition * twoN + 1))
        pose_time = np.zeros((np.shape(pose_base)[0], (nsteps + 1) * npartition))
        vel_time = np.zeros((cls.ndof_all, (nsteps + 1) * npartition))
        energy = np.zeros(npartition)
        cvg = [None] * npartition

        for ipart in range(npartition):
            i = ipart * twoN
            i1 = (ipart + 1) * twoN
            j = (ipart + 1) % npartition * twoN
            j1 = ((ipart + 1) % npartition + 1) * twoN
            p = ipart * (nsteps + 1)
            p1 = (ipart + 1) * (nsteps + 1)

            cls.config_update(pose_base[:, ipart])
            simdata, cvg[ipart] = cls.run_cpp(T * delta_S, X[i:i1], nsteps, rel_tol)
            energy[ipart] = simdata["/Model_0/energy"][:, -1][0]
            M = simdata["/Sensitivity/Monodromy"][:]
            dHdt = M[:, -1] * delta_S
            M = np.delete(M, -1, axis=1)
            pose_time[:, p:p1] = simdata["/Config/POSE"][:]
            vel_time[:, p:p1] = simdata["/Config/VELOCITY"][:]

            J[i:i1, i:i1] = M
            J[i:i1, j:j1] -= np.eye(twoN)
            J[i:i1, -1] = dHdt
            simdata.close()

        pose_base_new = pose_time[:, timesol_partition_index_start]
        energy = np.mean(energy)
        cvg = all(cvg)
        H1 = pose_time[cls.free_dof][:, timesol_partition_index_end[block_order]] - \
            pose_time[cls.free_dof][:, timesol_partition_index_start[block_order]]
        H2 = vel_time[cls.free_dof][:, timesol_partition_index_end[block_order]] - \
            vel_time[cls.free_dof][:, timesol_partition_index_start[block_order]]
        H = np.reshape(np.concatenate([H1, H2]), (-1, 1), order='F')

        return H, J, pose_time, vel_time, pose_base_new, energy, cvg

    @classmethod
    def runsim_single(cls, T, X, pose_base, cont_params):
        nperiod = cont_params["shooting"]["single"]["nperiod"]
        nsteps = cont_params["shooting"]["single"]["nsteps_per_period"]
        rel_tol = cont_params["shooting"]["rel_tol"]

        cls.config_update(pose_base)
        simdata, cvg = cls.run_cpp(T * nperiod, X, nsteps * nperiod, rel_tol)
        if cvg:
            energy = simdata["/Model_0/energy"][:, -1][0]
            periodicity_inc = simdata["/Periodicity/INC"][cls.ndof_fix:]
            periodicity_vel = simdata["/Periodicity/VELOCITY"][cls.ndof_fix:]
            pose_time = simdata["/Config/POSE"][:]
            pose_base_new = pose_time[:, 0]
            vel_time = simdata["/Config/VELOCITY"][:]
            H = np.concatenate([periodicity_inc, periodicity_vel])
            M = simdata["/Sensitivity/Monodromy"][:]
            dHdt = M[:, -1] * nperiod
            M = np.delete(M, -1, axis=1)
            M -= np.eye(len(M))
            J = np.concatenate((M, dHdt.reshape(-1, 1)), axis=1)
            simdata.close()
        else:
            H = J = pose_time = vel_time = pose_base_new = energy = None

        return H, J, pose_time, vel_time, pose_base_new, energy, cvg

    @classmethod
    def run_cpp(cls, T, X, nsteps, rel_tol):
        inc = X[:cls.ndof_free]
        vel = X[cls.ndof_free:]
        inc = np.pad(inc, (cls.ndof_fix, 0), "constant")
        vel = np.pad(vel, (cls.ndof_fix, 0), "constant")

        icdata = h5py.File(cls.cpp_path + cls.ic_file + ".h5", "a")
        if "/Config/INC" in icdata:
            del icdata["Config/INC"]
        if "/Config/VELOCITY" in icdata:
            del icdata["Config/VELOCITY"]
        icdata["/Config/INC"] = inc.reshape(-1, 1)
        icdata["/Config/VELOCITY"] = vel.reshape(-1, 1)
        icdata.close()

        cls.cpp_params["TimeIntegrationSolverParameters"]["number_of_steps"] = nsteps
        cls.cpp_params["TimeIntegrationSolverParameters"]["time"] = T
        cls.cpp_params["TimeIntegrationSolverParameters"]["rel_tol_res_forces"] = rel_tol
        cls.cpp_params["TimeIntegrationSolverParameters"]["initial_conditions"] = \
            cls.cpp_params["TimeIntegrationSolverParameters"]["_initial_conditions"]
        json.dump(cls.cpp_params, open(cls.cpp_path + "_" + cls.cpp_paramfile, "w"), indent=2)

        cpprun = subprocess.run("cd " + cls.cpp_path + "&&" + cls.cppsim_exe + " _" + cls.cpp_paramfile,
                                shell=True,
                                stdout=open(cls.cpp_path + "cpp.out", "w"),
                                stderr=open(cls.cpp_path + "cpp.err", "w"))
        simdata = h5py.File(cls.cpp_path + cls.simout_file + ".h5", "r")
        cvg = not bool(cpprun.returncode)
        return simdata, cvg

    @classmethod
    def partition_singleshooting_solution(cls, T, X, pose_base, cont_params):
        npartition = cont_params["shooting"]["multiple"]["npartition"]
        nsteps = cont_params["shooting"]["multiple"]["nsteps_per_partition"]
        rel_tol = cont_params["shooting"]["rel_tol"]
        cls.config_update(pose_base)
        # nsteps has to equal to total steps for multiple shooting so it can be sliced
        simdata, cvg = cls.run_cpp(T, X, nsteps * npartition, rel_tol)
        pose_time = simdata["/Config/POSE"][:]
        vel_time = simdata["/Config/VELOCITY"][:]
        slicing_index = nsteps * np.arange(npartition)
        V = vel_time[cls.free_dof][:, slicing_index]
        # update pose_base and set inc to zero
        pose_base = pose_time[:, slicing_index]
        X = np.concatenate((np.zeros((cls.ndof_free, npartition)), V))
        X = np.reshape(X, (-1), order='F')
        return X, pose_base

    @classmethod
    def config_update(cls, pose):
        # update beam configuration by writing initial conditions pose
        icdata = h5py.File(cls.cpp_path + cls.ic_file + ".h5", "w")
        icdata["/Config/POSE"] = pose.reshape(-1, 1)
        icdata.close()

    @classmethod
    def set_dofdata(cls, simdata):
        cls.free_dof = np.array(simdata["/Model_0/free_dof"])[:, 0]
        cls.ndof_all = simdata["/Model_0/number_of_dofs"][0][0]
        cls.ndof_fix = len(np.array(simdata["/Model_0/fix_dof"]))
        cls.ndof_free = len(cls.free_dof)

    @classmethod
    def get_dofdata(cls):
        return {"free_dof": cls.free_dof, "ndof_all": cls.ndof_all, "ndof_fix": cls.ndof_fix,
                "ndof_free": cls.ndof_free}

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
