from core.problem import Prob
from core.logger import Logger
from core.solver.continuation import ConX
from core.startingpoint import StartingPoint

from springcpp import SpringCpp

# Problem
prob = Prob()
prob.read_contparams("contparameters.json")
prob.add_doffunction(SpringCpp.get_dofdata)
prob.add_icfunction(SpringCpp.run_eig)
if prob.cont_params["continuation"]["forced"]:
    prob.add_zerofunction(SpringCpp.runsim_forced)
else:
    prob.add_zerofunction(SpringCpp.runsim_single)
prob.zerofunction_firstpoint = prob.zerofunction

# Continuation starting point
start = StartingPoint(prob)
start.get_startingpoint()

# Logger
log = Logger(prob)

# Solve continuation on problem
con = ConX(prob, start, log)
con.solve()
