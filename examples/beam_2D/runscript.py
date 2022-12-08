from core.problem import Prob
from core.logger import Logger
from core.continuation import ConX
from core.startingpoint import StartingPoint

from methods import beam_2d_sim, beam_2d_eig

# Problem object
prob = Prob()
prob.read_parameters("contparameters.json")
prob.add_zerofunction(beam_2d_sim)

# Continuation starting point object
start = StartingPoint(prob)
if prob.parameters["restart"]["file"]:
    start.restart()
else:
    start.clean_start(beam_2d_eig)

# Logger object
log = Logger(prob)

# Solve continuation on problem
con = ConX(prob, start, log)
con.solve()
