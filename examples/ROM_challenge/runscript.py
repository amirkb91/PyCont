from core.problem import Prob
from core.logger import Logger
from core.solver.continuation import ConX
from core.startingpoint import StartingPoint

from examples.ROM_challenge.romchallenge import ROMChallenge

# Problem
prob = Prob()
prob.read_contparams("contparameters.json")
prob.add_doffunction(ROMChallenge.get_dofdata)
prob.add_icfunction(ROMChallenge.run_eig)
if prob.cont_params["shooting"]["method"] == "single":
    prob.add_zerofunction(ROMChallenge.runsim_single)
elif prob.cont_params["shooting"]["method"] == "multiple":
    prob.add_zerofunction(ROMChallenge.runsim_multiple, ROMChallenge.runsim_single)
    prob.add_partitionfunction(ROMChallenge.partition_singleshooting_solution)

# Initialise class based on continuation parameters
ROMChallenge.initialise(prob.cont_params, False, 4)

# Continuation starting point
start = StartingPoint(prob)
start.get_startingpoint()

# Logger
log = Logger(prob)

# Solve continuation on problem
con = ConX(prob, start, log)
con.solve()
