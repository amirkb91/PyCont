import os
import json


class Prob:
    def __init__(self):
        self.cont_params = None
        self.zerofunction = None
        self.icfunction = None
        self.updatefunction = None
        self.doffunction = None

    def read_contparams(self, cont_paramfile):
        if os.path.exists(cont_paramfile):
            self.cont_params = json.load(open(cont_paramfile))
        else:
            raise Exception("Continuation parameter file does not exist!")

    def add_zerofunction(self, fxn):
        self.zerofunction = fxn

    def add_icfunction(self, fxn):
        self.icfunction = fxn

    def add_updatefunction(self, fxn):
        self.updatefunction = fxn

    def add_doffunction(self, fxn):
        self.doffunction = fxn
