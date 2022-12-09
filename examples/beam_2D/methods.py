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
    sim_file = cpp_params["TimeIntegrationSolverParameters"]["Logger"]["file_name"]

    @classmethod
    def run_eig(cls, cont_params):
        # Run eigen solver and store initial solution
        # Run sim solver to get model data
        cpprun = subprocess.run(
            "cd " + cls.cpp_path + "&&" + cls.cppeig_exe + " " + cls.cpp_paramfile +
            "&&" + cls.cppsim_exe + " " + cls.cpp_paramfile,
            shell=True,
            stdout=open(cls.cpp_path + "cpp.out", "w"),
            stderr=open(cls.cpp_path + "cpp.err", "w"),
        )

        eigdata = h5py.File(cls.cpp_path + cls.eig_file + ".h5", "r")
        simdata = h5py.File(cls.cpp_path + cls.sim_file + ".h5", "r")
        eig = np.array(eigdata["/eigen_analysis/Eigenvectors"])
        frq = eigdata["/eigen_analysis/Frequencies"]
        eig[np.abs(eig) < 1e-10] = 0.0
        free_dof = np.array(simdata["/Model_0/free_dof"])

        nnm = cont_params["continuation"]["NNM"]
        scale = cont_params["continuation"]["eig_scale"]
        x0 = scale * eig[:, nnm - 1]
        x0 = x0[free_dof]
        v0 = np.zeros_like(x0)
        X0 = np.concatenate([x0, v0])
        T0 = 1 / frq[nnm - 1]

        return X0, T0

    def run_sim(self, T, X, par):
        """
        User-defined zero function
        Here the function is the C++ simulation
        """
        # unpack run cont_params
        ns = par["shooting"]["npts"]
        ns_fail = par["shooting"]["npts_fail"]
        nperiod = par["shooting"]["nperiod"]
        rho = par["shooting"]["rho"]
        rel_tol = par["shooting"]["rel_tol"]
        beam_type = par["shooting"]["beam_type"]

        # path and file names, open eigen file to write INC and VEL
        path2cpp = "/home/akb110/Codes/beam_cpp/examples/cclamped/"
        exec_cpp = "/home/akb110/Codes/beam_cpp/cmake-build-release/examples/cclamped"
        paramfile = "cont_params.json"
        eigenfile = "cclamped_eig"
        ouputfile = "cclamped_out"
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

        # add X to eigen file (INC, VEL)
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
        cpp_parameter["Solver_parameters"]["number_of_steps"] = ns * nperiod
        cpp_parameter["Solver_parameters"]["time"] = T[0] * nperiod
        cpp_parameter["Solver_parameters"]["rho"] = rho
        cpp_parameter["Solver_parameters"]["rel_tol_res_forces"] = rel_tol
        cpp_parameter["Sensitivity"]["pose"] = True
        cpp_parameter["Sensitivity"]["velocity"] = True
        cpp_parameter["Sensitivity"]["period"] = True
        json.dump(cpp_parameter, open(path2cpp + "_" + paramfile, "w"), indent=2)

        # run C++ sim
        cpprun = subprocess.run(
            "cd " + path2cpp + "&&" + exec_cpp + " _" + paramfile,
            shell=True,
            stdout=open(path2cpp + "cpp.out", "w"),
            stderr=open(path2cpp + "cpp.err", "w"),
        )

        # if C++ hasn't converged, try one more time with more samples, else reduce continuation step
        cvg = True
        if cpprun.returncode != 0:
            print("C++ high samlpe.")
            cpp_parameter["Solver_parameters"]["number_of_steps"] = ns_fail
            json.dump(cpp_parameter, open(path2cpp + "_" + paramfile, "w"), indent=2)
            cpprun2 = subprocess.run(
                "cd " + path2cpp + "&&" + exec_cpp + " _" + paramfile,
                shell=True,
                stdout=open(path2cpp + "cpp.out", "w"),
                stderr=open(path2cpp + "cpp.err", "w"),
            )
            if cpprun2.returncode != 0:
                cvg = False

        if cvg:
            # read periodicity and gradients
            outputdata = h5py.File(path2cpp + ouputfile + ".h5", "r")
            inc_out = outputdata["/Periodicity/INC"]
            vel_out = outputdata["/Periodicity/VELOCITY"]
            H = np.concatenate([inc_out[6:], vel_out[6:]])
            Mm0 = outputdata["/Sensitivity/Monodromy_m0"][:]
            dHdt = nperiod * np.concatenate([Mm0[6:ndof_all, -1], Mm0[ndof_all + 6:, -1]])
            Mm0 = np.delete(Mm0, -1, axis=1)  # stored inside dHdt so delete
            Mm0 = np.block(
                [
                    [Mm0[6:ndof_all, 6:ndof_all], Mm0[6:ndof_all, ndof_all + 6:]],
                    [Mm0[ndof_all + 6:, 6:ndof_all], Mm0[ndof_all + 6:, ndof_all + 6:]],
                ]
            )

            # additional user requested outputs
            pose = outputdata["/Config/POSE"][:]
            energy = outputdata["/energy"][:, -1][0]
            outputs = {"energy": np.array([energy])}

            outputdata.close()

        else:
            H = Mm0 = dHdt = pose = outputs = None

        return H, Mm0, dHdt, pose, outputs, cvg
