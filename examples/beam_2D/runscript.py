from core.problem import Prob
from core.logger import Logger
from core.continuation import ConX
from core.startingpoint import StartingPoint
from methods import zerofunction, initialguess_eig

# Problem object
prob = Prob()
prob.read_parameters("contparameters.json")
prob.add_zerofunction(zerofunction)

# Continuation starting point object
start = StartingPoint(prob)
if prob.parameters["restart"]["file"]:
    start.restart()
else:
    start.initial_guess(initialguess_eig)

# Logger object
log = Logger(prob)

# Solve continuation on problem
con = ConX(prob, start, log)
con.solve()
