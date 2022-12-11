from core.problem import Prob
from core.logger import Logger
from core.continuation import ConX
from core.startingpoint import StartingPoint

from methods import BeamCpp

# Problem object
prob = Prob()
prob.read_contparams("contparameters.json")
prob.add_zerofunction(BeamCpp.run_sim)
prob.add_updatefunction(BeamCpp.config_update)

# Continuation starting point object
start = StartingPoint(prob)
if prob.cont_params["restart"]["file"]:
    start.restart()
else:
    start.new_start(BeamCpp.run_eig)

# Logger object
log = Logger(prob)

# Solve continuation on problem
con = ConX(prob, start, log)
con.solve()
