# this script tests class problem
import lsm_classes
import matplotlib.pyplot as plt
import numpy as np
import numpy.random 
from lsm_classes import PyLSMSolver


lsm_solver = lsm_classes.PyLSMSolver(161,81,0.5)
(bpts_xy, af, segmentlength) = lsm_solver.discretize()

movelimit = 0.5
nbpts = bpts_xy.shape[0]
bpts_sens = np.random.rand(nbpts)
print('-1')
lsm_solver.preprocess(movelimit, bpts_sens)
print('0')
(ub,lb) = lsm_solver.get_bounds()
print('a')
lambdas = np.zeros(2)
lambdas[0] = -1.2
lambdas[1] = 2
print('b')
displacement = lsm_solver.computeDisplacements(lambdas)
