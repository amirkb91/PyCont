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
        self.energy0 = start.energy0
        self.tgt0 = start.tgt0
        self.log = log

    def solve(self):
        # calculate phase condition matrix h
        phase_condition(self)

        # compute first point of the branch
        first_point(self)

        if self.prob.cont_params["continuation"]["method"].lower() == "seq":
            # sequential continuation
            seqcont(self)
        elif self.prob.cont_params["continuation"]["method"].lower() == "psa":
            # pseudo-arc length continuation
            psacont(self)
