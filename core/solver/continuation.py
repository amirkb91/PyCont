from ._phase_condition import phase_condition
from ._first_point import first_point
from ._seqcont import seqcont
from ._psacont import psacont
from ._psacont_mult import psacont_mult


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
        elif self.prob.cont_params["continuation"]["method"] == "psa" and \
                self.prob.cont_params["shooting"]["method"] == "single":
            # pseudo-arc length continuation single shooting
            psacont(self)
        elif self.prob.cont_params["continuation"]["method"] == "psa" and \
                self.prob.cont_params["shooting"]["method"] == "multiple":
            # pseudo-arc length continuation multiple shooting
            psacont_mult(self)
