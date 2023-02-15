from core.problem import Prob
from core.logger import Logger
from core.solver.continuation import ConX
from core.startingpoint import StartingPoint

from duffing import Duffing

# Problem
prob = Prob()
prob.read_contparams("contparameters.json")
prob.add_icfunction(Duffing.eigen_solve)
prob.add_zerofunction(Duffing.time_solve)
prob.add_zerofunction_firstpoint(Duffing.time_solve)

# Continuation starting point
start = StartingPoint(prob)
start.get_startingpoint()

# Logger
log = Logger(prob)

# Solve continuation on problem
con = ConX(prob, start, log)
con.solve()
