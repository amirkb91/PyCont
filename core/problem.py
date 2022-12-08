import os
import json


class Prob:
    def __init__(self):
        self.parameters = None
        self.zerofunction = None

    def read_parameters(self, contparameter_file):
        if os.path.exists(contparameter_file):
            self.parameters = json.load(open(contparameter_file))
        else:
            raise Exception("Continuation parameter file does not exist!")

    def add_zerofunction(self, fxn):
        self.zerofunction = fxn
