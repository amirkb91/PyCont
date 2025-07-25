import os
import json


class Prob:
    def __init__(self):
        self.cont_params = None
        self.doffunction = None
        self.icfunction = None
        self.zerofunction = None
        self.zerofunction_firstpoint = None
        self.partitionfunction = None

    def read_contparams(self, cont_paramfile):
        with open(cont_paramfile) as f:
            data = json.load(f)

        # Default parameters
        defaults = {
            "continuation": {
                "forced": False,
                "method": "psa",
                "tangent": "peeters",
                "dir": 1,
                "npts": 500,
                "tol": 0.001,
                "itermin": 1,
                "iteropt": 3,
                "itermax": 5,
                "iterjac": 1,
                "nadapt": 2,
                "s0": 1e-3,
                "smin": 1e-6,
                "smax": 1e-1,
                "betacontrol": False,
                "betamax": 20,
                "fmin": 20,
                "fmax": 100,
                "Emax": 1e5,
                "phase_index_unforced": "allvel",
            },
            "shooting": {
                "method": "single",
                "scaling": False,
                "rel_tol": 1e-08,
                "single": {"nperiod": 1, "nsteps_per_period": 200},
                "multiple": {"npartition": 3, "nsteps_per_partition": 100},
            },
            "forcing": {
                "amplitude": 1,
                "phase_ratio": 0.5,
                "tau0": 1e-4,
                "tau1": 1e-4,
                "rho_GA": 0.95,
                "starting_freq_scale": 0.7,
            },
            "first_point": {
                "from_eig": True,
                "itermax": 30,
                "eig_start": {"NNM": 1, "scale": 0.01},
                "restart": {
                    "file_name": "",
                    "index": 50,
                    "recompute_tangent": False,
                    "fixF": False,
                    "F": 60,
                },
            },
            "Logger": {"plot": True, "file_name": "temp"},
        }

        self.cont_params = self.fill_defaults(data, defaults)

    def fill_defaults(self, data, defaults):
        for key, value in defaults.items():
            if key not in data:
                data[key] = value
            elif isinstance(value, dict):
                data[key] = self.fill_defaults(data[key], value)
        return data

    def add_doffunction(self, fxn):
        self.doffunction = fxn

    def add_icfunction(self, fxn):
        self.icfunction = fxn

    def add_zerofunction(self, fxn, fxn2=None):
        self.zerofunction = fxn
        if not fxn2:
            self.zerofunction_firstpoint = fxn
        else:
            self.zerofunction_firstpoint = fxn2

    def add_zerofunction_firstpoint(self, fxn):
        self.zerofunction_firstpoint = fxn

    def add_partitionfunction(self, fxn):
        self.partitionfunction = fxn
