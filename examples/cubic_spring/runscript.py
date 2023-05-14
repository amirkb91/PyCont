from core.problem import Prob
from core.logger import Logger
from core.solver.continuation import ConX
from core.startingpoint import StartingPoint

from cubic_spring import Cubic_Spring

# Problem
prob = Prob()
prob.read_contparams("contparameters.json")
prob.add_icfunction(Cubic_Spring.eigen_solve)
prob.add_doffunction(Cubic_Spring.get_fe_data)
if prob.cont_params["shooting"]["method"] == "single":
    prob.add_zerofunction(Cubic_Spring.time_solve)
    prob.zerofunction_firstpoint = prob.zerofunction
elif prob.cont_params["shooting"]["method"] == "multiple":
    prob.add_zerofunction(Cubic_Spring.time_solve_multiple)
    prob.add_zerofunction_firstpoint(Cubic_Spring.time_solve)
    prob.add_partitionfunction(Cubic_Spring.partition_singleshooting_solution)

# Continuation starting point
start = StartingPoint(prob)
start.get_startingpoint()

# Logger
log = Logger(prob)

# Solve continuation on problem
con = ConX(prob, start, log)
con.solve()
