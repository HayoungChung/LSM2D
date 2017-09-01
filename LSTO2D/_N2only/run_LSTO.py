import numpy as np

from openmdao.api import Problem, view_model, ScipyOptimizer

from fem2d import PyFEMSolver
from forces import get_forces

from lsm2d_SLP_Group import LSM2D_slpGroup 
# from pySensitivity_vec import computeSensitivity
# from lsm2d import discretize, get_boundaryPoints, get_sensitivity, updates

# FEM Mesh
num_nodes_x = 81
num_nodes_y = 41
 
# LSM Mesh
num_param_x = num_nodes_x
num_param_y = num_nodes_y

length_x = num_nodes_x-1
length_y = num_nodes_y-1

num_dofs = num_nodes_x * num_nodes_y * 2 + num_nodes_y 

# FEA properties
E = 1.
nu = 0.3
f = -1.

num_gpts = 2

# LSM properties 
radius = 5
movelimit = 0.5

# FEA initialize 
fem_solver = PyFEMSolver(num_nodes_x, num_nodes_y, length_x, length_y, E, nu)
forces = get_forces(num_nodes_x, num_nodes_y, f=f)

# LSM iniitliaze
# lsm_solver = PyLSMSolver(num_nodes_x, num_nodes_y) # iniital swbpts_yisscheese config 
# ID_fix = lsm_solver.get_fixedID()

# HJ loop
max_loop = 100
# for iHJ in range(0,maX_loop):
#     # 0. discretize
#     bpts_xy, areafraction = lsm_solver.discretize()






# irs = np.zeros(num_dofs)
# jcs = np.zeros(num_dofs)
# data = np.zeros(num_dofs)
# fem_solver.get_stiffness(irs, jcs, data) # should be constant
# Kmat = scipy..coo_sparse(data,(irs,jcs))
# ilu = ...
# y_dof = ilu.solve(forces)
# u_dof = y_dof(1::num_dofs)

# LSM






# fixedGpts_sens = fem_solver.get_sensitivity(u_dof, num_gpts)

# LSM
# lsm_module = pyLSMsolver(num_nodes_x, num_nodes_y)
# lsm_mesh = lsm_module.get_mesh()
bpts_x = bpts_y = np.zeros(10)
fixedGpts_x = fixedGpts_y = np.zeros(10)
fixedGpts_sens = np.zeros(10)


model = LSM2D_slpGroup(
    # lsm_module = lsm_module,
    # lsm_mesh = lsm_mesh,
    bpts_x = bpts_x, bpts_y = bpts_y, # boundary points
    fixedGpts_x = fixedGpts_x, fixedGpts_y = fixedGpts_y, # Gausspoints at fixed grid
    fixedGpts_sens = fixedGpts_sens, # sensitivity at Gausspoints at fixed grid 
    radius = radius, # radius for least square methods
    movelimit = movelimit,
)

prob = Problem(model)

prob.driver = ScipyOptimizer()
prob.driver.options['optimizer'] = 'SLSQP'
prob.driver.options['tol'] = 1e-9
prob.driver.options['disp'] = True

prob.setup()

view_model(prob)

# prob.run_driver()
# lambdas = model.lambdas