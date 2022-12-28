from core.problem import Prob
from core.logger import Logger
from core.solver.continuation import ConX
from core.startingpoint import StartingPoint

from beamcpp import BeamCpp

# Problem
prob = Prob()
prob.read_contparams("contparameters.json")
prob.add_zerofunction(BeamCpp.run_sim)
prob.add_icfunction(BeamCpp.run_eig)
prob.add_updatefunction(BeamCpp.config_update)

# Continuation starting point
start = StartingPoint(prob)
start.get_startingpoint()

# Logger
log = Logger(prob)

# Solve continuation on problem
con = ConX(prob, start, log)
con.solve()
