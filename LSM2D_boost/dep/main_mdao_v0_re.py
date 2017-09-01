from mdao_module_v1 import LSM_M2DO, FEA_M2DO, Sensitivity_M2DO, Optim_pre_#, Optim_post_
from openmdao.api import Component, Group, IndepVarComp, Problem, ScipyOptimizer
import numpy as np
from slsm_Module import *

lxy = [160,80]
exy = [160,80]

# initialize the classes 
lsm = LSM_M2DO(exy)
fea = FEA_M2DO(lxy,exy)
ElemArea = np.zeros(exy[0]*exy[1])
for ii in range(0,exy[0]*exy[1]):
    ElemArea[ii] = lsm.mesh.elements[ii].area
    
fea.set_AreaFraction(ElemArea) # initial set
sens = Sensitivity_M2DO(fea.CLinElasticity,len(lsm.boundary.points))
optim_pre = Optim_pre_(0.9,lsm.boundary.points)
#optim_post = Optim_post_(0.9)

# ======== GROUPING ===================================
group_1 = Group() # group at the outermost layer
group_1.add('fea',fea)
group_1.add('lsm',lsm)
group_1.add('sens',sens)
group_1.add('optim_pre',optim_pre)
#group_1.add('optim_post',optim_post)

# FEA ~ LSM
group_1.connect('lsm.AreaFraction', 'fea.AreaFraction')

# FEA ~ Sens 
group_1.connect('fea.Field', 'sens.Field')

# Sens ~ LSM
group_1.connect('lsm.BoundaryPointCoords', 'sens.BoundaryPointCoords')
group_1.connect('lsm.AreaFraction', 'sens.AreaFraction')

# Sens ~ Optim_pre_
group_1.connect('sens.BoundaryPointSensitivities', 'optim_pre.BoundaryPointSensitivities')

# LSM ~ Optim_pre / _post
group_1.connect('lsm.BoundaryPoints', 'optim_pre.BoundaryPoints')
#group_1.connect('lsm.BoundaryPoints', 'optim_post.BoundaryPoints')

group_1.connect('lsm.ConstraintDistance', 'optim_pre.ConstraintDistance')
#group_1.connect('lsm.ConstraintDistance', 'optim_post.ConstraintDistance')

#group_1.connect('optim_post.BoundaryPointVelocities', 'lsm.BoundaryPointVelocities')
#group_1.connect('optim_post.timeStep', 'lsm.timeStep')

# btw optims
#group_1.connect('optim_pre.lambdas','optim_post.lambdas')
#group_1.connect('optim_post.lambdas_rescaled','optim_post.lambdas')
#group_1.connect('optim_pre.F_obj','optim_post.res_fun')

## promote problem functions
#group_1.add('optim_pre', optim_pre, promotes = ['F_obj'])
#group_1.add('optim_pre', optim_pre, promotes = ['F_cons'])
#group_1.add('optim_pre', optim_pre, promotes = ['lambdas'])

# QUESTION: inserting optim_post into the main python script?
# since here every modules are updating: 
# only velocity and lambdas are updated... 
# so let that be initinalized.. 

# ======== PROBLEM DEF ===================================
top = Problem()
top.root = group_1

top.root.add('p1',IndepVarComp('lambdas',np.zeros(2)))
top.root.connect('p1.lambdas','optim_pre.lambdas')

top.root.deriv_options['type'] = 'fd'

top.driver = ScipyOptimizer()
top.driver.options['optimizer'] = 'SLSQP'
top.driver.options['tol'] = 1e-8
top.driver.options['maxiter'] = 1000

top.driver.add_desvar('p1.lambdas',\
            lower = np.array([-238057, -1165]), \
            upper = np.array([17320147, 1159])) # must be changed

top.driver.add_objective('optim_pre.F_obj')
top.driver.add_constraint('optim_pre.F_cons',upper=0.0)

#top.check_setup()
top.setup()
top.run()
# =========================================================
# ======== 1st Loop for time integration ==================
# =========================================================
