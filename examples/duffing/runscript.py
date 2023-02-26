from core.problem import Prob
from core.logger import Logger
from core.solver.continuation import ConX
from core.startingpoint import StartingPoint

from duffing import Duffing

# Problem
prob = Prob()
prob.read_contparams("contparameters.json")
prob.add_icfunction(Duffing.eigen_solve)
prob.add_doffunction(Duffing.get_fe_data)
if prob.cont_params["shooting"]["method"] == "single":
    prob.add_zerofunction(Duffing.time_solve)
    prob.zerofunction_firstpoint = prob.zerofunction
elif prob.cont_params["shooting"]["method"] == "multiple":
    prob.add_zerofunction(Duffing.time_solve_multiple)
    prob.add_zerofunction_firstpoint(Duffing.time_solve)
    prob.add_partitionfunction(Duffing.partition_singleshooting_solution)

# Continuation starting point
start = StartingPoint(prob)
start.get_startingpoint()

# Logger
log = Logger(prob)

# Solve continuation on problem
con = ConX(prob, start, log)
con.solve()
