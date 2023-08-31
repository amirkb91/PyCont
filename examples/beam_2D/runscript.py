from core.problem import Prob
from core.logger import Logger
from core.solver.continuation import ConX
from core.startingpoint import StartingPoint

from beamcpp import BeamCpp

# Problem
prob = Prob()
prob.read_contparams("contparameters.json")
prob.add_doffunction(BeamCpp.get_dofdata)
prob.add_icfunction(BeamCpp.run_eig)
if prob.cont_params["shooting"]["method"] == "single":
    prob.add_zerofunction(BeamCpp.runsim_single)
elif prob.cont_params["shooting"]["method"] == "multiple":
    prob.add_zerofunction(BeamCpp.runsim_multiple, BeamCpp.runsim_single)    
    prob.add_partitionfunction(BeamCpp.partition_singleshooting_solution)

# Initialise class based on continuation parameters
BeamCpp.initialise(prob.cont_params)

# Continuation starting point
start = StartingPoint(prob)
start.get_startingpoint()

# Logger
log = Logger(prob)

# Solve continuation on problem
con = ConX(prob, start, log)
con.solve()
