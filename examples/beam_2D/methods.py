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

    @classmethod
    def run_eig(cls, cont_params):
        # Run eigen solver and store initial solution
        # Run sim solver to get config data
        cpprun = subprocess.run(
            "cd " + cls.cpp_path + "&&" + cls.cppeig_exe + " " + cls.cpp_paramfile +
            "&&" + cls.cppsim_exe + " " + cls.cpp_paramfile,
            shell=True,
            stdout=open(cls.cpp_path + "cpp.out", "w"),
            stderr=open(cls.cpp_path + "cpp.err", "w"),
        )

        eigdata = h5py.File(cls.cpp_path + cls.eig_file + ".h5", "r")
        simdata = h5py.File(cls.cpp_path + cls.simout_file + ".h5", "r")

        cls.free_dof = np.array(simdata["/Model_0/free_dof"])
        cls.ndof_all = simdata["/Model_0/number_of_dofs"][0][0]
        cls.ndof_fix = len(np.array(simdata["/Model_0/fix_dof"]))

        eig = np.array(eigdata["/eigen_analysis/Eigenvectors"])
        frq = eigdata["/eigen_analysis/Frequencies"]
        eig[np.abs(eig) < 1e-10] = 0.0
        nnm = cont_params["continuation"]["NNM"]
        scale = cont_params["continuation"]["eig_scale"]
        x0 = scale * eig[:, nnm - 1]
        x0 = x0[cls.free_dof][:, 0]
        v0 = np.zeros_like(x0)
        X0 = np.concatenate([x0, v0])
        T0 = 1 / frq[nnm - 1, 0]

        return X0, T0

    @classmethod
    def run_sim(cls, T, X, par):
        # unpack run cont_params
        ns = par["shooting"]["npts"]
        ns_fail = par["shooting"]["npts_fail"]
        nperiod = par["shooting"]["nperiod"]
        rel_tol = par["shooting"]["rel_tol"]

        # get INC and VEL from X
        lenX_2 = len(X) // 2
        inc = X[:lenX_2]
        vel = X[lenX_2:]
        inc = np.pad(inc, (cls.ndof_fix, 0), "constant")
        vel = np.pad(vel, (cls.ndof_fix, 0), "constant")

        # write initial conditions to ic_file
        icdata = h5py.File(cls.cpp_path + cls.ic_file + ".h5", "w")
        icdata["/Config/INC"] = inc.reshape(-1, 1)
        icdata["/Config/VELOCITY"] = vel.reshape(-1, 1)
        icdata.close()

        # edit C++ parameter file
        cls.cpp_params["TimeIntegrationSolverParameters"]["number_of_steps"] = ns * nperiod
        cls.cpp_params["TimeIntegrationSolverParameters"]["time"] = T * nperiod
        cls.cpp_params["TimeIntegrationSolverParameters"]["rel_tol_res_forces"] = rel_tol
        cls.cpp_params["TimeIntegrationSolverParameters"]["initial_conditions"] = \
            cls.cpp_params["TimeIntegrationSolverParameters"]["_initial_conditions"]

        runtwice = False
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
                if runtwice:
                    break
                cvg = False
                cls.cpp_params["TimeIntegrationSolverParameters"]["number_of_steps"] = ns_fail * nperiod
                runtwice = True

        if cvg:
            simdata = h5py.File(cls.cpp_path + cls.simout_file + ".h5", "r")
            energy = simdata["/Model_0/energy"][:, -1][0]
            periodicity_inc = simdata["/Periodicity/INC"][cls.ndof_fix:]
            periodicity_vel = simdata["/Periodicity/VELOCITY"][cls.ndof_fix:]
            H = np.concatenate([periodicity_inc, periodicity_vel])
            M = simdata["/Sensitivity/Monodromy"][:]
            dHdt = nperiod * M[:, -1]
            M = np.delete(M, -1, axis=1)
            M = M - np.eye(len(M))

            pose = simdata["/Config/POSE"][:]
            vel = simdata["/Config/VELOCITY"][:]

            simdata.close()

        else:
            H = M = dHdt = pose = vel = energy = None

        return H, M, dHdt, pose, vel, energy, cvg
