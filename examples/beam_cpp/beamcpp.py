import h5py
import json
import subprocess
import numpy as np
from copy import deepcopy as dp
import shutil
from concurrent.futures import ProcessPoolExecutor
from Frame import Frame


class BeamCpp:
    # --------- Choose example case from mb_sef_cpp ---------#
    cpp_example = "beam_2D"  # doubly clamped, arch, cantilever
    # cpp_example = "beam_rightangle"  # right angle beam
    # cpp_example = "beam_boxwing"  # Benchmark Aircraft
    # cpp_example = "beam_vertcant"  # 2D vertical cantilever beam
    # cpp_example = "beam_3Dcant"  # 3D cantilever beam
    # cpp_example = "beam_bristolwing" # 15% deflection Bristol wing

    if cpp_example == "beam_2D":
        folder, exe = "mybeam_2D/", "mybeam_2D"
    elif cpp_example == "beam_rightangle":
        folder, exe = "mybeam_rightangle/", "mybeam_rightangle"
    elif cpp_example == "beam_boxwing":
        folder, exe = "mybeam_boxwing/", "mybeam_boxwing"
    elif cpp_example == "beam_vertcant":
        folder, exe = "mybeam_vertcant/", "mybeam_vertcant_nosupport"
        cpp_def_period = 2.0
    elif cpp_example == "beam_3Dcant":
        folder, exe = "mybeam_3Dcant/", "mybeam_3Dcant"
    elif cpp_example == "beam_bristolwing":
        folder, exe = "mybeam_bristolwing/", "mybeam_bristolwing"
    # -------------------------------------------------------#

    cpp_path = "/home/akb110/Codes/mb_sef_cpp/examples/" + folder
    cpp_exe = "/home/akb110/Codes/mb_sef_cpp/cmake-build-release/examples/" + exe
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
    config_per_node = None
    nnodes_all = None
    dof_per_node = None
    nnodes_free = None
    n_dim = None
    if "cpp_def_period" not in vars():
        cpp_def_period = 1.0

    @classmethod
    def initialise(cls, cont_params, ForcePeriod=False, nprocs=1):
        # prep para file and assign fixed values
        cls.cpp_params_sim["TimeIntegrationSolverParameters"]["rel_tol_res_forces"] = cont_params[
            "shooting"
        ]["rel_tol"]
        cls.cpp_params_sim["TimeIntegrationSolverParameters"]["initial_conditions"] = (
            cls.cpp_params_sim["TimeIntegrationSolverParameters"].pop("_initial_conditions")
        )
        cls.cpp_params_sim["TimeIntegrationSolverParameters"].pop("_initial_conditions_eig")

        if not cont_params["continuation"]["forced"]:
            cls.model_def["ModelDef"]["amplitude"] = 0.0
            cls.model_def["ModelDef"]["tau0"] = 0.0
            cls.model_def["ModelDef"]["tau1"] = 0.0
            cls.cpp_params_sim["TimeIntegrationSolverParameters"]["rho"] = 0.9
        elif cont_params["continuation"]["forced"]:
            cls.model_def["ModelDef"]["amplitude"] = cont_params["forcing"]["amplitude"]
            cls.model_def["ModelDef"]["tau0"] = cont_params["forcing"]["tau0"]
            cls.model_def["ModelDef"]["tau1"] = cont_params["forcing"]["tau1"]
            cls.model_def["ModelDef"]["phase_ratio"] = cont_params["forcing"]["phase_ratio"]
            cls.model_def["ModelDef"]["def_period"] = cls.cpp_def_period
            cls.model_def["ModelDef"]["with_ratio"] = True
            cls.cpp_params_sim["TimeIntegrationSolverParameters"]["rho"] = cont_params["forcing"][
                "rho_GA"
            ]

        if ForcePeriod:
            # if force period is specified, we want to de-couple force period from sim time
            # useful when running time simulation as post processing
            cls.model_def["ModelDef"]["def_period"] = ForcePeriod / cls.cpp_def_period
            cls.model_def["ModelDef"]["with_ratio"] = False

        subprocess.run("cd " + cls.cpp_path + "&&" + "./clean_dir.sh", shell=True)
        json.dump(cls.model_def, open(cls.cpp_path + "_" + cls.cpp_modelfile, "w"), indent=2)

        cls.nprocs = nprocs  # number of CPU cores
        cls.column_split = False  # split Jacobian columns in C++ code

    @classmethod
    def run_eig(cls):
        cmd = "cd " + cls.cpp_path + "&&" + cls.cpp_exe + " _" + cls.cpp_modelfile + " " + cls.cpp_paramfile_eig  # fmt: skip
        subprocess.run(
            cmd,
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
        cls, omega, tau, Xtilde, pose_base, cont_params, sensitivity=True, fulltime=False
    ):
        nperiod = cont_params["shooting"]["single"]["nperiod"]
        nsteps = cont_params["shooting"]["single"]["nsteps_per_period"]
        N = cls.ndof_free

        T = tau / omega
        X = Xtilde.copy()
        X[N:] *= omega  # scale velocities from Xtilde to X
        H = J = pose = vel = energy = None

        cls.config_update(pose_base)
        cvg = cls.run_cpp(T * nperiod, X, nsteps * nperiod, sensitivity)

        if cvg:
            simdata = h5py.File(cls.cpp_path + cls.simout_file + ".h5", "r")
            energy = np.max(simdata["/dynamic_analysis/FEModel/energy"][:, :])
            periodicity_inc = simdata["/dynamic_analysis/Periodicity/INC"][cls.free_dof]
            periodicity_vel = simdata["/dynamic_analysis/Periodicity/VELOCITY"][cls.free_dof]
            pose_time = simdata["/dynamic_analysis/FEModel/POSE/MOTION"][:]
            vel_time = simdata["/dynamic_analysis/FEModel/VELOCITY/MOTION"][:]
            acc_time = simdata["/dynamic_analysis/FEModel/ACCELERATION/MOTION"][:]
            pose = pose_time[:, 0]
            vel = vel_time[:, 0]
            H = np.concatenate([periodicity_inc, periodicity_vel])
            if sensitivity:
                # scale velocity and time derivatives with omega and nperiod
                J = simdata["/Sensitivity/Monodromy"][:, :]
                J[:, N:] *= omega
                J[:, -1] *= nperiod / omega
            simdata.close()

        # elif cvg and cls.column_split:
        #     for i in range(1, cls.nprocs + 1):
        #         suffix = f"_{i:03d}"
        #         simdata = h5py.File(cls.cpp_path + cls.simout_file + suffix + ".h5", "r")
        #         sens_H_i = simdata["/Sensitivity/Monodromy"][:, :-1]
        #         _N = np.shape(sens_H_i)[1] // 2
        #         sens_H_i[:, _N:] *= omega
        #         if i == 1:
        #             energy = np.max(simdata["/dynamic_analysis/FEModel/energy"][:, :])
        #             periodicity_inc = simdata["/dynamic_analysis/Periodicity/INC"][cls.free_dof]
        #             periodicity_vel = simdata["/dynamic_analysis/Periodicity/VELOCITY"][
        #                 cls.free_dof
        #             ]
        #             pose_time = simdata["/dynamic_analysis/FEModel/POSE/MOTION"][:]
        #             vel_time = simdata["/dynamic_analysis/FEModel/VELOCITY/MOTION"][:]
        #             pose = pose_time[:, 0]
        #             vel = vel_time[:, 0]
        #             H = np.concatenate([periodicity_inc, periodicity_vel])
        #             dHdtau = simdata["/Sensitivity/Monodromy"][:, -1] * nperiod * 1 / omega
        #             sens_H_inc = sens_H_vel = np.empty((2 * cls.ndof_free, 0))
        #         sens_H_inc = np.concatenate((sens_H_inc, sens_H_i[:, :_N]), axis=1)
        #         sens_H_vel = np.concatenate((sens_H_vel, sens_H_i[:, _N:]), axis=1)
        #         simdata.close()
        #     sens_H = np.concatenate((sens_H_inc, sens_H_vel), axis=1)
        #     J = np.concatenate((sens_H, dHdtau.reshape(-1, 1)), axis=1)

        if not fulltime:
            return H, J, pose, vel, energy, cvg
        else:
            return H, J, pose_time, vel_time, acc_time, energy, cvg

    @classmethod
    def runsim_multiple(
        cls, omega, tau, Xtilde, pose_base, cont_params, sensitivity=True, fulltime=False
    ):
        # multiple shooting sensitivity SE correction has to be done in Python
        cls.cpp_params_sim["TimeIntegrationSolverParameters"]["direct_sensitivity"][
            "apply_SE_correction"
        ] = False
        npartition = cont_params["shooting"]["multiple"]["npartition"]
        nsteps = cont_params["shooting"]["multiple"]["nsteps_per_partition"]
        N = cls.ndof_free
        twoN = 2 * N
        delta_S = 1 / npartition
        T = tau / omega  # scaling is not implemented in this method, omega will not be used

        # Precomputations
        partition_extremeties = np.arange(npartition + 1) * (nsteps + 1)
        indices_start = partition_extremeties[:npartition]
        indices_end = partition_extremeties[1:] - 1
        block_order = (np.arange(npartition) + 1) % npartition

        # Initialisations
        J = np.zeros((npartition * twoN, npartition * twoN + 1))
        M_all = np.empty((npartition, twoN, twoN + 1))
        pose_time = np.zeros((cls.ndof_config, (nsteps + 1) * npartition))
        vel_time = np.zeros((cls.ndof_all, (nsteps + 1) * npartition))
        acc_time = np.zeros((cls.ndof_all, (nsteps + 1) * npartition))
        energy_time = np.zeros((nsteps + 1, npartition))
        energy = 0
        cvg = [None] * npartition
        H = pose = vel = None

        # Run all partitions, in parallel if user specifies
        if cls.nprocs > 1:
            with ProcessPoolExecutor(max_workers=cls.nprocs) as executor:
                results = list(
                    executor.map(
                        cls.run_partition,
                        range(npartition),
                        [Xtilde] * npartition,
                        pose_base.T,
                        [T] * npartition,
                        [delta_S] * npartition,
                        [nsteps] * npartition,
                        [sensitivity] * npartition,
                    )
                )
        else:
            results = []
            for ipart in range(npartition):
                res = cls.run_partition(
                    ipart, Xtilde, pose_base[:, ipart], T, delta_S, nsteps, sensitivity
                )
                results.append(res)

        # Collect results
        for res in results:
            ipart, cvg_ipart, E, pose_seg, vel_seg, acc_seg, energy_time_seg, M = res
            p0, p1 = partition_extremeties[ipart], partition_extremeties[ipart + 1]
            cvg[ipart] = cvg_ipart
            energy = max(energy, E)
            pose_time[:, p0:p1] = pose_seg
            vel_time[:, p0:p1] = vel_seg
            acc_time[:, p0:p1] = acc_seg
            energy_time[:, ipart] = energy_time_seg
            M_all[ipart] = M
        cvg = all(cvg)

        # Compute periodicity and Jacobian
        if cvg:
            if cls.SEbeam:
                H1 = cls.periodicity_INC_SE_local(
                    pose_time[:, indices_start[block_order]],
                    pose_time[:, indices_end],
                )
                H2 = cls.periodicity_VEL_SE_local(
                    H1,
                    vel_time[:, indices_start[block_order]],
                    vel_time[:, indices_end],
                )
                # Jacobian computation for SE beam
                for ipart in range(npartition):
                    i0, i1 = ipart * twoN, (ipart + 1) * twoN
                    j0, j1 = (ipart + 1) % npartition * twoN, ((ipart + 1) % npartition + 1) * twoN
                    VT_thispart = vel_time[:, indices_end[block_order][ipart]]
                    V0_nextpart = vel_time[:, indices_start[block_order][ipart]]
                    [dHdx0, dHdx1, dHdT] = cls.sens_periodicity_SEcorrection_local_multi(
                        M_all[ipart], H1[:, ipart], VT_thispart, V0_nextpart
                    )
                    J[i0:i1, i0:i1] = dHdx0
                    J[i0:i1, j0:j1] += dHdx1
                    J[i0:i1, -1] = dHdT * delta_S  # scale time derivative
            else:
                H1 = cls.periodicity_INC_linear(
                    pose_time[:, indices_start[block_order]],
                    pose_time[:, indices_end],
                )
                H2 = cls.periodicity_VEL_linear(
                    vel_time[:, indices_start[block_order]],
                    vel_time[:, indices_end],
                )
                # Jacobian computation for VK beam
                for ipart in range(npartition):
                    i0, i1 = ipart * twoN, (ipart + 1) * twoN
                    j0, j1 = (ipart + 1) % npartition * twoN, (
                        ((ipart + 1) % npartition) + 1
                    ) * twoN
                    J[i0:i1, i0:i1] = M_all[ipart, :, :-1]
                    J[i0:i1, j0:j1] += -np.eye(twoN)
                    J[i0:i1, -1] = M_all[ipart, :, -1] * delta_S  # scale time derivative

            H = np.reshape(np.concatenate([H1, H2]), (-1, 1), order="F")
            # solution pose and vel at time 0 for each partition
            pose = pose_time[:, indices_start]
            vel = vel_time[:, indices_start]

        if not fulltime:
            return H, J, pose, vel, energy, cvg
        else:
            p3d = np.transpose(pose_time.reshape(cls.ndof_config, npartition, nsteps + 1), (0, 2, 1))  # fmt: skip
            v3d = np.transpose(vel_time.reshape(cls.ndof_all, npartition, nsteps + 1), (0, 2, 1))
            a3d = np.transpose(acc_time.reshape(cls.ndof_all, npartition, nsteps + 1), (0, 2, 1))
            return H, J, p3d, v3d, a3d, energy_time, cvg

    @classmethod
    def run_partition(cls, ipart, Xtilde, pose_base_column, T, delta_S, nsteps, sensitivity):
        N = cls.ndof_free
        twoN = 2 * N
        i0, i1 = ipart * twoN, (ipart + 1) * twoN
        X = Xtilde[i0:i1].copy()

        cls.config_update(pose_base_column, run_id=ipart + 1)
        cvg_ipart = cls.run_cpp(T * delta_S, X, nsteps, sensitivity, run_id=ipart + 1)

        # Read the output
        suffix = f"_{ipart+1:03d}"
        filepath = cls.cpp_path + cls.simout_file + suffix + ".h5"
        with h5py.File(filepath, "r") as simdata:
            energy_arr = simdata["/dynamic_analysis/FEModel/energy"][:]
            E = np.max(energy_arr)
            pose_seg = simdata["/dynamic_analysis/FEModel/POSE/MOTION"][:]
            vel_seg = simdata["/dynamic_analysis/FEModel/VELOCITY/MOTION"][:]
            acc_seg = simdata["/dynamic_analysis/FEModel/ACCELERATION/MOTION"][:]
            energy_time_seg = simdata["/dynamic_analysis/FEModel/energy"][:]
            M = simdata["/Sensitivity/Monodromy"][:]
            # cpp gives M - I when apply_SE_correction=false
            M[:twoN, :twoN] += np.eye(twoN)

        return ipart, cvg_ipart, E, pose_seg, vel_seg, acc_seg, energy_time_seg, M

    @classmethod
    def run_cpp(cls, T, X, nsteps, sensitivity, run_id=None):
        inc = np.zeros(cls.ndof_all)
        vel = np.zeros(cls.ndof_all)
        inc[cls.free_dof] = X[: cls.ndof_free]
        vel[cls.free_dof] = X[cls.ndof_free :]

        suffix = f"_{run_id:03d}" if run_id is not None else ""
        icdata = h5py.File(cls.cpp_path + cls.ic_file + suffix + ".h5", "a")
        icdata["/" + cls.analysis_name + "/FEModel/INC/MOTION"] = inc.reshape(-1, 1)
        icdata["/" + cls.analysis_name + "/FEModel/VELOCITY/MOTION"] = vel.reshape(-1, 1)
        icdata.close()

        cls.cpp_params_sim["TimeIntegrationSolverParameters"]["number_of_steps"] = nsteps
        cls.cpp_params_sim["TimeIntegrationSolverParameters"]["time"] = T

        # Note that multiple shooting first point does single shooting, so check run_id
        if run_id is None:
            cpp_params_sim = dp(cls.cpp_params_sim)  # don't want pop to permanently pop dict
            if not sensitivity:
                cpp_params_sim["TimeIntegrationSolverParameters"].pop("direct_sensitivity")
            json.dump(
                cpp_params_sim, open(cls.cpp_path + "_" + cls.cpp_paramfile_sim, "w"), indent=2
            )
            cmd = "cd " + cls.cpp_path + "&&" + cls.cpp_exe + " _" + cls.cpp_modelfile + " _" + cls.cpp_paramfile_sim  # fmt: skip
            cpprun = subprocess.run(
                cmd,
                shell=True,
                stdout=open(cls.cpp_path + "cpp.out", "w"),
                stderr=open(cls.cpp_path + "cpp.err", "w"),
            )
            cvg = not bool(cpprun.returncode)

        elif run_id is not None:
            cpp_params_sim = dp(cls.cpp_params_sim)
            if not sensitivity:
                cpp_params_sim["TimeIntegrationSolverParameters"].pop("direct_sensitivity")

            cpp_params_sim["TimeIntegrationSolverParameters"]["Logger"]["file_name"] = (
                cls.simout_file + suffix
            )
            cpp_params_sim["TimeIntegrationSolverParameters"]["initial_conditions"]["file_name"] = (
                cls.ic_file + suffix
            )
            json.dump(
                cpp_params_sim,
                open(
                    cls.cpp_path + "_" + cls.cpp_paramfile_sim.split(".")[0] + suffix + ".json", "w"
                ),
                indent=2,
            )

            cmd = "cd " + cls.cpp_path + "&&" + cls.cpp_exe + " _" + cls.cpp_modelfile + " _" + cls.cpp_paramfile_sim.split(".")[0] + suffix + ".json"  # fmt: skip
            cpprun = subprocess.run(
                cmd,
                shell=True,
                stdout=open(cls.cpp_path + "cpp" + suffix + ".out", "w"),
                stderr=open(cls.cpp_path + "cpp" + suffix + ".err", "w"),
            )
            cvg = not bool(cpprun.returncode)

        # elif cls.column_split:
        #     # Calculate the basic split size and the number of splits that need an extra column
        #     basic_split_size, extra_splits = divmod(cls.ndof_free, cls.nprocs)
        #     start_indices = np.arange(cls.nprocs) * basic_split_size + np.minimum(
        #         np.arange(cls.nprocs), extra_splits
        #     )
        #     end_indices = np.roll(start_indices, -1)
        #     end_indices[-1] = cls.ndof_free  # Correct the last end index
        #     split_indices = zip(start_indices, end_indices)

        #     with ProcessPoolExecutor(max_workers=cls.nprocs) as executor:
        #         convergence = list(
        #             executor.map(
        #                 cls.run_cpp_column_split, zip(split_indices, range(1, cls.nprocs + 1))
        #             )
        #         )
        #     cvg = np.all(convergence)

        return cvg

    # @classmethod
    # def run_cpp_column_split(cls, combined_args):
    #     (start_index, end_index), run_id = combined_args

    #     requested_cols = np.array([start_index, end_index]).tolist()
    #     suffix = f"_{run_id:03d}"
    #     cls.cpp_params_sim["TimeIntegrationSolverParameters"]["direct_sensitivity"][
    #         "requested_cols"
    #     ] = requested_cols
    #     cls.cpp_params_sim["TimeIntegrationSolverParameters"]["Logger"]["file_name"] = (
    #         cls.simout_file + suffix
    #     )
    #     cls.cpp_params_sim["TimeIntegrationSolverParameters"]["initial_conditions"]["file_name"] = (
    #         cls.ic_file + suffix
    #     )
    #     json.dump(
    #         cls.cpp_params_sim,
    #         open(cls.cpp_path + "_" + cls.cpp_paramfile_sim.split(".")[0] + suffix + ".json", "w"),
    #         indent=2,
    #     )
    #     shutil.copyfile(
    #         cls.cpp_path + cls.ic_file + ".h5", cls.cpp_path + cls.ic_file + suffix + ".h5"
    #     )

    #     cmd = "cd " + cls.cpp_path + "&&" + cls.cpp_exe + " _" + cls.cpp_modelfile + " _" + cls.cpp_paramfile_sim.split(".")[0] + suffix + ".json"  # fmt: skip
    #     cpprun = subprocess.run(
    #         cmd,
    #         shell=True,
    #         stdout=open(cls.cpp_path + "cpp" + suffix + ".out", "w"),
    #         stderr=open(cls.cpp_path + "cpp" + suffix + ".err", "w"),
    #     )
    #     return not bool(cpprun.returncode)

    @classmethod
    def partition_singleshooting_solution(cls, tau, Xtilde, pose_base, cont_params):
        npartition = cont_params["shooting"]["multiple"]["npartition"]
        nsteps = cont_params["shooting"]["multiple"]["nsteps_per_partition"]
        N = cls.ndof_free
        slicing_index = nsteps * np.arange(npartition)

        T = tau
        X = Xtilde.copy()
        cls.config_update(pose_base)

        # do time integration along whole orbit before slicing
        cls.run_cpp(T, X, nsteps * npartition, True)
        simdata = h5py.File(cls.cpp_path + cls.simout_file + ".h5", "r")
        pose_time = simdata["/dynamic_analysis/FEModel/POSE/MOTION"][:]
        vel_time = simdata["/dynamic_analysis/FEModel/VELOCITY/MOTION"][:]
        pose = pose_time[:, slicing_index]
        vel = vel_time[cls.free_dof][:, slicing_index]
        # set inc to zero as solution stored in pose
        Xsol = np.concatenate((np.zeros((N, npartition)), vel))
        Xsol = np.reshape(Xsol, (-1), order="F")
        return Xsol, pose

    @classmethod
    def config_update(cls, pose, run_id=None):
        # update beam configuration by writing initial conditions pose
        suffix = f"_{run_id:03d}" if run_id is not None else ""
        icdata = h5py.File(cls.cpp_path + cls.ic_file + suffix + ".h5", "w")
        icdata["/" + cls.analysis_name + "/FEModel/POSE/MOTION"] = pose.reshape(-1, 1)
        icdata.close()

    @classmethod
    def read_dofdata(cls):
        with h5py.File(cls.cpp_path + cls.model_name + ".h5", "r") as data:
            cls.free_dof = np.array(data["/FEModel/loc_dof_free/MOTION"])[:, 0]
            cls.fix_dof = np.array(data["/FEModel/loc_dof_fix/MOTION"])[:, 0]
            NodalFrame = list(data["/FEModel/Nodes_config/"].keys())[0]
            cls.node_config = np.array(data["/FEModel/Nodes_config/" + NodalFrame])[1:, :]
        cls.ndof_free = cls.free_dof.size
        cls.ndof_fix = cls.fix_dof.size
        cls.ndof_config = cls.node_config.size
        cls.ndof_all = cls.ndof_free + cls.ndof_fix
        cls.config_per_node = cls.node_config.shape[0]
        cls.nnodes_all = cls.node_config.shape[1]
        cls.dof_per_node = cls.ndof_all // cls.nnodes_all
        cls.nnodes_free = cls.ndof_free // cls.dof_per_node
        cls.n_dim = 2 if cls.dof_per_node == 3 else 3
        cls.SEbeam = cls.config_per_node != cls.dof_per_node

    @classmethod
    def get_dofdata(cls):
        return {
            "free_dof": cls.free_dof,
            "ndof_all": cls.ndof_all,
            "ndof_fix": cls.ndof_fix,
            "ndof_free": cls.ndof_free,
            "node_config": cls.node_config,
            "ndof_config": cls.ndof_config,
            "config_per_node": cls.config_per_node,
            "nnodes_all": cls.nnodes_all,
            "dof_per_node": cls.dof_per_node,
            "n_dim": cls.n_dim,
            "SEbeam": cls.SEbeam,
        }

    @classmethod
    def periodicity_INC_linear(cls, pose_a, pose_b):
        return (pose_b - pose_a)[cls.free_dof, :]

    @classmethod
    def periodicity_VEL_linear(cls, vel_a, vel_b):
        return (vel_b - vel_a)[cls.free_dof, :]

    @classmethod
    def periodicity_INC_SE_local(cls, pose_a, pose_b):
        # inc = pose_a^-1 o pose_b
        # loop over all partitions, inputs have to be 2D arrays
        periodicity_inc = np.zeros((cls.ndof_all, pose_a.shape[1]))
        for i in range(pose_a.shape[1]):
            for j in range(cls.nnodes_all):
                f = Frame.relative_frame(
                    cls.n_dim,
                    pose_a[j * cls.config_per_node : (j + 1) * cls.config_per_node, i],
                    pose_b[j * cls.config_per_node : (j + 1) * cls.config_per_node, i],
                )
                p = Frame.get_parameters_from_frame(cls.n_dim, f)
                periodicity_inc[j * cls.dof_per_node : (j + 1) * cls.dof_per_node, i] = p
        return periodicity_inc[cls.free_dof, :]

    @classmethod
    def periodicity_VEL_SE_local(cls, inc, vel_a, vel_b):
        # loop over all partitions, inputs have to be 2D arrays
        # inc only has free dofs already, so slice vel_a and vel_b
        # input inc should be periodicity_inc_local
        vel_a = vel_a[cls.free_dof, :]
        vel_b = vel_b[cls.free_dof, :]
        periodicity_vel = np.zeros((cls.ndof_free, vel_a.shape[1]))
        for i in range(vel_a.shape[1]):
            for j in range(cls.nnodes_free):
                v = (
                    -Frame.get_inverse_tangent_operator(
                        cls.n_dim, -inc[j * cls.dof_per_node : (j + 1) * cls.dof_per_node, i]
                    )
                    @ vel_a[j * cls.dof_per_node : (j + 1) * cls.dof_per_node, i]
                    + Frame.get_inverse_tangent_operator(
                        cls.n_dim, inc[j * cls.dof_per_node : (j + 1) * cls.dof_per_node, i]
                    )
                    @ vel_b[j * cls.dof_per_node : (j + 1) * cls.dof_per_node, i]
                )
                periodicity_vel[j * cls.dof_per_node : (j + 1) * cls.dof_per_node, i] = v
        return periodicity_vel

    @classmethod
    def periodicity_INC_SE_global(cls, pose_a, pose_b):
        # inc = pose_b o pose_a^-1
        # loop over all partitions, inputs have to be 2D arrays
        periodicity_inc = np.zeros((cls.ndof_all, pose_a.shape[1]))
        for i in range(pose_a.shape[1]):
            for j in range(cls.nnodes_all):
                finv = Frame.get_inverse(
                    cls.n_dim, pose_a[j * cls.config_per_node : (j + 1) * cls.config_per_node, i]
                )
                f = Frame.composition(
                    cls.n_dim,
                    pose_b[j * cls.config_per_node : (j + 1) * cls.config_per_node, i],
                    finv,
                )
                p = Frame.get_parameters_from_frame(cls.n_dim, f)
                periodicity_inc[j * cls.dof_per_node : (j + 1) * cls.dof_per_node, i] = p
        return periodicity_inc[cls.free_dof, :]

    @classmethod
    def periodicity_VEL_SE_global(cls, inc, pose_a, vel_a, vel_b):
        # loop over all partitions, inputs have to be 2D arrays
        # inc only has free dofs already, so slice vel_a and vel_b
        # input inc should be periodicity_inc_global
        vel_a = vel_a[cls.free_dof, :]
        vel_b = vel_b[cls.free_dof, :]
        periodicity_vel = np.zeros((cls.ndof_free, vel_a.shape[1]))
        for i in range(vel_a.shape[1]):
            for j in range(cls.nnodes_free):
                v = (
                    Frame.get_inverse_tangent_operator(
                        cls.n_dim, inc[j * cls.dof_per_node : (j + 1) * cls.dof_per_node, i]
                    )
                    @ Frame.get_adjoint(
                        cls.n_dim,
                        pose_a[j * cls.config_per_node : (j + 1) * cls.config_per_node, i],
                    )
                    @ (
                        vel_b[j * cls.dof_per_node : (j + 1) * cls.dof_per_node, i]
                        - vel_a[j * cls.dof_per_node : (j + 1) * cls.dof_per_node, i]
                    )
                )
                periodicity_vel[j * cls.dof_per_node : (j + 1) * cls.dof_per_node, i] = v
        return periodicity_vel

    @classmethod
    def sens_periodicity_SEcorrection_local(cls, monodromy, periodicity_inc, vel_time):
        """
        Compute the sensitivity of the periodicity condition with SE correction.
        **** This applies to the LOCAL periodcity formulations. ****

        When SE correction is set to false in C++, the monodromy matrix is returned as M - I
        If we add I we get the actual monodromy.
        This method then corrects the monodromy to account for the SE contributions.
        It gives the same result as the C++ code when apply_SE_correction == true.
        """
        N = cls.ndof_free
        nodes = cls.nnodes_free
        dpn = cls.dof_per_node
        n_dim = cls.n_dim
        vel_time = vel_time[cls.free_dof]
        sens_SE = np.zeros((2 * N, 2 * N + 1))

        # Precompute indices
        indices_N = np.array([(i * dpn, (i + 1) * dpn) for i in range(nodes)])
        indices_2N = np.array([(i * dpn, (i + 1) * dpn) for i in range(nodes * 2)])
        indices_comb = np.array(
            [(i * dpn, (i + 1) * dpn, N + i * dpn, N + (i + 1) * dpn) for i in range(nodes)]
        )

        # Precompute inverse tangent operators
        T_pos_neg = np.array(
            [
                (
                    Frame.get_inverse_tangent_operator(
                        n_dim, periodicity_inc[i * dpn : (i + 1) * dpn, 0]
                    ),
                    Frame.get_inverse_tangent_operator(
                        n_dim, -periodicity_inc[i * dpn : (i + 1) * dpn, 0]
                    ),
                )
                for i in range(nodes)
            ]
        )

        # POSE sensitivities
        for i, (i0, i1) in enumerate(indices_N):
            T_pos, T_neg = T_pos_neg[i]
            for j0, j1 in indices_2N:
                sens_SE[i0:i1, j0:j1] = T_pos @ monodromy[i0:i1, j0:j1]
            # add T_neg on the diagonal elements only
            sens_SE[i0:i1, i0:i1] += -T_neg

            # POSE sensitivities wrt period
            sens_SE[i0:i1, -1] = T_pos @ monodromy[i0:i1, -1]

        # VEL sensitivities:
        for i, (i0, i1, I0, I1) in enumerate(indices_comb):
            psi = periodicity_inc[i0:i1, 0]
            v0 = vel_time[i0:i1, 0]
            vT = vel_time[i0:i1, -1]
            T_pos, T_neg = T_pos_neg[i]
            DT_pos = Frame.get_derivative_inverse_tangent_operator(n_dim, psi, vT)
            DT_neg = Frame.get_derivative_inverse_tangent_operator(n_dim, -psi, v0)
            for j0, j1 in indices_2N:
                pose_sens = sens_SE[i0:i1, j0:j1]
                sens_SE[I0:I1, j0:j1] = (
                    (T_pos @ monodromy[I0:I1, j0:j1]) + DT_neg @ pose_sens + DT_pos @ pose_sens
                )
            # add T_neg on the diagonal elements only
            sens_SE[I0:I1, I0:I1] += -T_neg

            # VEL sensitivities wrt period
            pose_sens = sens_SE[i0:i1, -1]
            sens_SE[I0:I1, -1] = (
                T_pos @ monodromy[I0:I1, -1] + DT_neg @ pose_sens + DT_pos @ pose_sens
            )

        return sens_SE

    @classmethod
    def sens_periodicity_SEcorrection_local_multi(
        cls, monodromy, periodicity_inc, VT_thispart, V0_nextpart
    ):
        """
        Compute the sensitivity of the periodicity condition with SE correction for multiple shooting.
        **** This applies to the LOCAL periodcity formulations. ****

        Running multiple shooting with a single partition gives the same result as single shooting.
        Since single shooting sensitivities are computed in C++ with apply_SE_correction, this means
        this method is implemented correctly and is therefore verified. Both flavours of velocity tested.
        """
        N = cls.ndof_free
        nodes = cls.nnodes_free
        dpn = cls.dof_per_node
        n_dim = cls.n_dim
        VT_thispart = VT_thispart[cls.free_dof]
        V0_nextpart = V0_nextpart[cls.free_dof]
        dHdx0 = np.zeros((2 * N, 2 * N))
        dHdx1 = np.zeros((2 * N, 2 * N))
        dHdT = np.zeros((2 * N))

        # Precompute indices
        indices_N = np.array([(i * dpn, (i + 1) * dpn) for i in range(nodes)])
        indices_2N = np.array([(i * dpn, (i + 1) * dpn) for i in range(nodes * 2)])
        indices_comb = np.array(
            [(i * dpn, (i + 1) * dpn, N + i * dpn, N + (i + 1) * dpn) for i in range(nodes)]
        )

        # Precompute inverse tangent operators
        T_pos_neg = np.array(
            [
                (
                    Frame.get_inverse_tangent_operator(
                        n_dim, periodicity_inc[i * dpn : (i + 1) * dpn]
                    ),
                    Frame.get_inverse_tangent_operator(
                        n_dim, -periodicity_inc[i * dpn : (i + 1) * dpn]
                    ),
                )
                for i in range(nodes)
            ]
        )

        # POSE sensitivities
        for i, (i0, i1) in enumerate(indices_N):
            T_pos, T_neg = T_pos_neg[i]
            for j0, j1 in indices_2N:
                dHdx0[i0:i1, j0:j1] = T_pos @ monodromy[i0:i1, j0:j1]

            dHdx1[i0:i1, i0:i1] = -T_neg

            # POSE sensitivities wrt period
            dHdT[i0:i1] = T_pos @ monodromy[i0:i1, -1]

        # VEL sensitivities:
        for i, (i0, i1, I0, I1) in enumerate(indices_comb):
            psi = periodicity_inc[i0:i1]
            v0 = V0_nextpart[i0:i1]
            vT = VT_thispart[i0:i1]
            T_pos, T_neg = T_pos_neg[i]
            DT_pos = Frame.get_derivative_inverse_tangent_operator(n_dim, psi, vT)
            DT_neg = Frame.get_derivative_inverse_tangent_operator(n_dim, -psi, v0)
            for j0, j1 in indices_2N:
                pose_sens_x0 = dHdx0[i0:i1, j0:j1]
                pose_sens_x1 = dHdx1[i0:i1, j0:j1]
                dHdx0[I0:I1, j0:j1] = (
                    (T_pos @ monodromy[I0:I1, j0:j1])
                    + DT_pos @ pose_sens_x0
                    + DT_neg @ pose_sens_x0
                )
                dHdx1[I0:I1, j0:j1] = DT_pos @ pose_sens_x1 + DT_neg @ pose_sens_x1
            dHdx1[I0:I1, I0:I1] = -T_neg

            # VEL sensitivities wrt period
            pose_sens = dHdT[i0:i1]
            dHdT[I0:I1] = T_pos @ monodromy[I0:I1, -1] + DT_neg @ pose_sens + DT_pos @ pose_sens

        # **** If you choose to do VEL periodicity linearly, use the following code instead ****
        # VEL sensitivities:
        # dHdx0[N:,:] = monodromy[N:, :-1]
        # dHdT[N:] = monodromy[N:, -1]
        # dHdx1[N:, N:] = -np.eye(N)

        return dHdx0, dHdx1, dHdT
