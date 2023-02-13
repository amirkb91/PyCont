from ._phase_condition import phase_condition
from ._first_point import first_point
from ._seqcont import seqcont
from ._psacont import psacont


class ConX:
    def __init__(self, prob, start, log):
        self.h = None
        self.nphase = None
        self.prob = prob
        self.X0 = start.X0
        self.T0 = start.T0
        self.pose_base0 = start.pose_base0
        self.tgt0 = None
        self.energy0 = None
        self.pose_time0 = None
        self.vel_time0 = None
        self.log = log

    def solve(self):
        # calculate phase condition matrix h
        phase_condition(self)
        # correct starting solution
        first_point(self)
        if self.prob.cont_params["continuation"]["method"] == "seq":
            # sequential continuation
            seqcont(self)
        elif self.prob.cont_params["continuation"]["method"] == "psa":
            psacont(self)
