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
        # Run eigen solver and store initial solution and run sim to get dof data
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

        # partition X0 and pose_base0 for multiple shooting
        if cont_params["shooting"]["method"] == "multiple":
            cls.config_update(pose_base0)
            simdata = cls.runsim_single_oneorbit(T0, X0, cont_params)
            X0, pose_base0 = cls.partition_singleshooting_solution(simdata, cont_params)

        return X0, T0, pose_base0

    @classmethod
    def runsim_single_oneorbit(cls, T, X, cont_params):
        nperiod = 1
        nsteps = cont_params["shooting"]["single"]["nsteps_per_period"]
        rel_tol = cont_params["shooting"]["rel_tol"]
        simdata = cls.run_cpp(T, X, nperiod, nsteps, rel_tol)
        return simdata

    @classmethod
    def run_cpp(cls, T, X, nperiod, nsteps, rel_tol):
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

        cls.cpp_params["TimeIntegrationSolverParameters"]["number_of_steps"] = nsteps * nperiod
        cls.cpp_params["TimeIntegrationSolverParameters"]["time"] = T * nperiod
        cls.cpp_params["TimeIntegrationSolverParameters"]["rel_tol_res_forces"] = rel_tol
        cls.cpp_params["TimeIntegrationSolverParameters"]["initial_conditions"] = \
            cls.cpp_params["TimeIntegrationSolverParameters"]["_initial_conditions"]
        json.dump(cls.cpp_params, open(cls.cpp_path + "_" + cls.cpp_paramfile, "w"), indent=2)

        subprocess.run("cd " + cls.cpp_path + "&&" + cls.cppsim_exe + " _" + cls.cpp_paramfile,
                       shell=True,
                       stdout=open(cls.cpp_path + "cpp.out", "w"),
                       stderr=open(cls.cpp_path + "cpp.err", "w"))
        simdata = h5py.File(cls.cpp_path + cls.simout_file + ".h5", "r")
        return simdata

    @classmethod
    def run_sim_single(cls):
        pass

    @classmethod
    def run_sim_multiple(cls, T, X, cont_params):
        # unpack run cont_params
        npartition = cont_params["shooting"]["npartition_multipleshooting"]
        nperiod = 1
        nsteps = cont_params["shooting"]["nsteps_per_period"]
        nsteps //= npartition
        rel_tol = cont_params["shooting"]["rel_tol"]

        # get INC and VEL from X
        inc = X[:cls.ndof_free]
        vel = X[cls.ndof_free:]
        inc = np.pad(inc, (cls.ndof_fix, 0), "constant")
        vel = np.pad(vel, (cls.ndof_fix, 0), "constant")

        # write initial conditions to ic_file
        icdata = h5py.File(cls.cpp_path + cls.ic_file + ".h5", "a")
        if "/Config/INC" in icdata:
            del icdata["Config/INC"]
        if "/Config/VELOCITY" in icdata:
            del icdata["Config/VELOCITY"]
        icdata["/Config/INC"] = inc.reshape(-1, 1)
        icdata["/Config/VELOCITY"] = vel.reshape(-1, 1)
        icdata.close()

        # edit C++ parameter file
        cls.cpp_params["TimeIntegrationSolverParameters"]["number_of_steps"] = nsteps * nperiod
        cls.cpp_params["TimeIntegrationSolverParameters"]["time"] = T * nperiod
        cls.cpp_params["TimeIntegrationSolverParameters"]["rel_tol_res_forces"] = rel_tol
        cls.cpp_params["TimeIntegrationSolverParameters"]["initial_conditions"] = \
            cls.cpp_params["TimeIntegrationSolverParameters"]["_initial_conditions"]

        run_twice = False
        while True:
            json.dump(cls.cpp_params, open(cls.cpp_path + "_" + cls.cpp_paramfile, "w"), indent=2)
            cpprun = subprocess.run(
                "cd " + cls.cpp_path + "&&" + cls.cppsim_exe + " _" + cls.cpp_paramfile,
                shell=True,
                stdout=open(cls.cpp_path + "cpp.out", "w"),
                stderr=open(cls.cpp_path + "cpp.err", "w"),
            )
            if cpprun.returncode == 0:
                cvg = True
                break
            else:
                cvg = False
                if run_twice:
                    print(f"Time Sim failed - Running with {fine_factor}x points")
                    cls.cpp_params["TimeIntegrationSolverParameters"]["number_of_steps"] = nsteps_fine * nperiod
                    run_twice = False
                else:
                    break

        if cvg:
            simdata = h5py.File(cls.cpp_path + cls.simout_file + ".h5", "r")
            energy = simdata["/Model_0/energy"][:, -1][0]
            periodicity_inc = simdata["/Periodicity/INC"][cls.ndof_fix:]
            periodicity_vel = simdata["/Periodicity/VELOCITY"][cls.ndof_fix:]
            M = simdata["/Sensitivity/Monodromy"][:]
            dHdt = nperiod * M[:, -1]
            M = np.delete(M, -1, axis=1)
            pose = simdata["/Config/POSE"][:]
            vel = simdata["/Config/VELOCITY"][:]
            if target is None:
                H = np.concatenate([periodicity_inc, periodicity_vel])
            else:
                # multiple shooting: periodicity calculated with respect to target
                # take final pose and vel to compare with target
                H = cls.periodicity(pose[:, -1], vel[:, -1], target)

            simdata.close()
        else:
            H = M = dHdt = pose = vel = energy = None

        return H, M, dHdt, pose, vel, energy, cvg

    @classmethod
    def run_sim(cls, T, X, cont_params, mult=False, target=None):
        # unpack run cont_params
        npartition = cont_params["shooting"]["npartition_multipleshooting"]
        nperiod = cont_params["shooting"]["nperiod_singleshooting"]
        nsteps = cont_params["shooting"]["nsteps_per_period"]
        rel_tol = cont_params["shooting"]["rel_tol"]
        fine_factor = 2
        # modify nsteps if multiple shooting
        if mult:
            nsteps //= npartition
            nperiod = 1
        nsteps_fine = nsteps * fine_factor

        # get INC and VEL from X
        inc = X[:cls.ndof_free]
        vel = X[cls.ndof_free:]
        inc = np.pad(inc, (cls.ndof_fix, 0), "constant")
        vel = np.pad(vel, (cls.ndof_fix, 0), "constant")

        # write initial conditions to ic_file
        icdata = h5py.File(cls.cpp_path + cls.ic_file + ".h5", "a")
        if "/Config/INC" in icdata:
            del icdata["Config/INC"]
        if "/Config/VELOCITY" in icdata:
            del icdata["Config/VELOCITY"]
        icdata["/Config/INC"] = inc.reshape(-1, 1)
        icdata["/Config/VELOCITY"] = vel.reshape(-1, 1)
        icdata.close()

        # edit C++ parameter file
        cls.cpp_params["TimeIntegrationSolverParameters"]["number_of_steps"] = nsteps * nperiod
        cls.cpp_params["TimeIntegrationSolverParameters"]["time"] = T * nperiod
        cls.cpp_params["TimeIntegrationSolverParameters"]["rel_tol_res_forces"] = rel_tol
        cls.cpp_params["TimeIntegrationSolverParameters"]["initial_conditions"] = \
            cls.cpp_params["TimeIntegrationSolverParameters"]["_initial_conditions"]

        run_twice = False
        while True:
            json.dump(cls.cpp_params, open(cls.cpp_path + "_" + cls.cpp_paramfile, "w"), indent=2)
            cpprun = subprocess.run(
                "cd " + cls.cpp_path + "&&" + cls.cppsim_exe + " _" + cls.cpp_paramfile,
                shell=True,
                stdout=open(cls.cpp_path + "cpp.out", "w"),
                stderr=open(cls.cpp_path + "cpp.err", "w"),
            )
            if cpprun.returncode == 0:
                cvg = True
                break
            else:
                cvg = False
                if run_twice:
                    print(f"Time Sim failed - Running with {fine_factor}x points")
                    cls.cpp_params["TimeIntegrationSolverParameters"]["number_of_steps"] = nsteps_fine * nperiod
                    run_twice = False
                else:
                    break

        if cvg:
            simdata = h5py.File(cls.cpp_path + cls.simout_file + ".h5", "r")
            energy = simdata["/Model_0/energy"][:, -1][0]
            periodicity_inc = simdata["/Periodicity/INC"][cls.ndof_fix:]
            periodicity_vel = simdata["/Periodicity/VELOCITY"][cls.ndof_fix:]
            M = simdata["/Sensitivity/Monodromy"][:]
            dHdt = nperiod * M[:, -1]
            M = np.delete(M, -1, axis=1)
            pose = simdata["/Config/POSE"][:]
            vel = simdata["/Config/VELOCITY"][:]
            if target is None:
                H = np.concatenate([periodicity_inc, periodicity_vel])
            else:
                # multiple shooting: periodicity calculated with respect to target
                # take final pose and vel to compare with target
                H = cls.periodicity(pose[:, -1], vel[:, -1], target)

            simdata.close()
        else:
            H = M = dHdt = pose = vel = energy = None

        return H, M, dHdt, pose, vel, energy, cvg

    @classmethod
    def config_update(cls, pose):
        # update beam configuration by writing initial conditions pose
        icdata = h5py.File(cls.cpp_path + cls.ic_file + ".h5", "w")
        icdata["/Config/POSE"] = pose.reshape(-1, 1)
        icdata.close()

    @classmethod
    def partition_singleshooting_solution(cls, simdata, cont_params):
        pose_time = simdata["/Config/POSE"][:]
        vel_time = simdata["/Config/VELOCITY"][:]
        npartition = cont_params["shooting"]["multiple"]["npartition"]
        nsteps = cont_params["shooting"]["multiple"]["nsteps_per_partition"]
        timesol_partition_index = nsteps * np.arange(npartition)
        V = vel_time[cls.free_dof][:, timesol_partition_index]
        # INC is zero as pose_base is also updated
        X = np.concatenate((np.zeros((cls.ndof_free, npartition)), V))
        X = np.reshape(X, (-1), order='F')
        pose_base = pose_time[:, timesol_partition_index]
        return X, pose_base

    @classmethod
    def get_dofdata(cls):
        return {"free_dof": cls.free_dof, "ndof_all": cls.ndof_all, "ndof_fix": cls.ndof_fix,
                "ndof_free": cls.ndof_free}

    @classmethod
    def set_dofdata(cls, simdata):
        cls.free_dof = np.array(simdata["/Model_0/free_dof"])[:, 0]
        cls.ndof_all = simdata["/Model_0/number_of_dofs"][0][0]
        cls.ndof_fix = len(np.array(simdata["/Model_0/fix_dof"]))
        cls.ndof_free = len(cls.free_dof)

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
