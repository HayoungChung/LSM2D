from mdao_module_v3 import LSM_M2DO, FEA_M2DO, Sensitivity_M2DO, Optimiser_pre, Optimiser_post
from openmdao.api import Component, Group, IndepVarComp, Problem, ScipyOptimizer
import numpy as np
from slsm_Module import *
import matplotlib.pyplot as plt
from Opt_Module import *



lxy = [160, 80]
exy = [160, 80]

AllowedAreaFraction = 0.01
moveLimit = 0.9
maxTime = 200
minArea = 0.5 
VoidMaterial = 1e-6

# initialize the classes 
fea = FEA_M2DO(lxy,exy)
lsm = LSM_M2DO(exy)

lambdas = vector__double__()
lambdas.append(0)
lambdas.append(0)

Multipliers = vector__double__()
Multipliers.append(0)

sens = Sensitivity_M2DO(fea.CLinElasticity, np.zeros([len(lsm.boundary.points),2]))
meshArea = exy[0]*exy[1]
minArea = 0.5

constraintDistances = vector__double__()
constraintDistances.append(0)

for ii in range(0,1):
    constraintDistances[0] = meshArea*minArea - lsm.boundary.area

    boundary_deepcopy = slsm.vector_BoundaryPoint_()
    boundary_deepcopy.extend(lsm.boundary.points)
    # boundary_deepcopy = copy.deepcopy(lsm.boundary)
    
    b = optim.Optimise(boundary_deepcopy, constraintDistances, lambdas, Multipliers, 0.9, False) # bpt: call-by-refs
    
    optim_pre  = Optimiser_pre(b, boundary_deepcopy)
    optim_post = Optimiser_post(b, boundary_deepcopy)

    # ======== GROUPING ===================================
    group_1 = Group() # group at the outermost layer
    group_1.add('fea',fea)
    group_1.add('lsm',lsm)
    group_1.add('sens',sens)
    group_1.add('optim_pre',optim_pre)
    group_1.add('optim_post',optim_post)

    # FEA ~ LSM
    group_1.connect('lsm.AreaFraction', 'fea.AreaFraction')
    
    # FEA ~ Sens 
    group_1.connect('fea.Field', 'sens.Field')

    # Sens ~ LSM
    # group_1.connect('lsm.AreaFraction', 'sens.AreaFraction')
    group_1.connect('lsm.BoundaryPointCoords', 'sens.BoundaryPointCoords')

    # Sens ~ iotim
    group_1.connect('sens.BoundaryPointSensitivities', 'optim_pre.BoundaryPointSensitivities')

    # LSM ~ optim_out
    group_1.connect('optim_post.timestep','lsm.timeStep')
    group_1.connect('optim_post.BoundaryPoints','lsm.BoundaryPoints')
    
    #group_1.connect('optim_pre.lambdas','optim_post.lambdas')
    

    # ======== PROBLEM DEF ===================================
    top = Problem()
    top.root = group_1
    top.root.add('p1', IndepVarComp('lambdas',np.zeros(2)))
    top.root.connect('p1.lambdas', 'optim_pre.lambdas')
    top.root.connect('p1.lambdas', 'optim_post.lambdas')

    top.root.deriv_options['type'] = 'fd'

    top.driver = ScipyOptimizer()
    top.driver.options['optimizer'] = 'SLSQP'
    # top.driver.options['tol'] = 1e-8
    # top.driver.options['maxiter'] = 1000

    top.driver.add_desvar('p1.lambdas',\
            lower = np.array(optim_pre.get_lowerLimits()), \
            upper = np.array(optim_pre.get_upperLimits())) # must be changed

    top.driver.add_objective('optim_pre.F_obj')
    top.driver.add_constraint('optim_pre.F_cons',upper=0.0)

    #top.check_setup()
    top.setup()
    top.run()
    # top.run_once()
    # =========================================================
    # ======== 1st Loop for time integration ==================
    # =========================================================
    lambdas_x = top['p1.lambdas']
    res = top['optim_pre.F_obj']
    cons = top['optim_pre.F_cons']
    
    # BoundaryPoints = lsm.get_BoundaryCoords()
    # plt.plot(BoundaryPoints[:,0],BoundaryPoints[:,1],'o')
    # plt.show()
