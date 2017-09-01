import numpy as np
import scipy.sparse 
import scipy.sparse.linalg

from fem2d import PyFEMSolver
from forces import get_forces

#from plot import get_mesh, plot_solution
from lsm_classes import PyLSMSolver

# QUICKFIX before fem2d is implemented
from optim_refact_v3 import Cantilever

from ls_sensitivity import _LeastSquare

# FEM Mesh
num_nodes_x = 161
num_nodes_y = 81
num_nodes = num_nodes_x * num_nodes_y
num_elems = (num_nodes_x-1) * (num_nodes_y-1)

# FEA properties
E = 1.
nu = 0.3
f = -1.

# LSM Mesh
num_param_x = num_nodes_x
num_param_y = num_nodes_y

length_x = num_nodes_x-1
length_y = num_nodes_y-1

num_dofs = num_nodes_x * num_nodes_y * 2 + num_nodes_y * 2

# LSM initialize (swisscheese config)
lsm_solver = PyLSMSolver(num_nodes_x, num_nodes_y, 0.5) 
(bpts_xy, areafraction, segmentLength) = lsm_solver.discretize()
#areafraction = np.ones(num_elems)
## fem2d ==================================
# FEA initialize 
fem_solver = PyFEMSolver(num_nodes_x, num_nodes_y, length_x, length_y, E, nu)
forces = get_forces(num_nodes_x, num_nodes_y, f=f)
rhs = np.zeros(num_dofs)
rhs[:num_nodes_x*num_nodes_y*2] = forces

num_sparse = num_elems * 64 * 4 + 2 * 2 * num_nodes_y
irs = np.zeros(num_sparse, dtype=np.int32)
jcs = np.zeros(num_sparse, dtype=np.int32)
data = np.zeros(num_sparse)

fem_solver.get_stiffness_matrix_LSTO(areafraction, data, irs, jcs)

# 1. get stiffness matrix & solve (LU decomposition)
mtx = scipy.sparse.csc_matrix((data, (irs, jcs)), shape=(num_dofs, num_dofs))
lumtx = scipy.sparse.linalg.splu(mtx)
_dofs = lumtx.solve(rhs)
u_fem2d = _dofs[:num_nodes*2]

# plot_solution(orig_nodes, deflected_nodes=u_fem2d)
# QUICKFIX ===============================
cantilever = Cantilever(True)
cantilever.set_fea()

Order_gpts = 2
num_gpts = num_elems*Order_gpts**2

u_cant = cantilever.get_u(areafraction)
# ==========================================
(IntegrationPointSensitivities, BoundarySensitivities) = cantilever.get_sens(bpts_xy)

xpos = np.ones(num_elems*4)
ypos = np.ones(num_elems*4)
sens_compl = np.ones(num_elems*4)

fem_solver.get_sensitivity_LSTO(u_fem2d, xpos, ypos, sens_compl)
xpos = xpos.reshape((num_elems,Order_gpts**2))
ypos = ypos.reshape((num_elems,Order_gpts**2))
fixedGpts_xy = np.zeros((num_elems,Order_gpts**2,2))
fixedGpts_xy[:,:,0] = xpos
fixedGpts_xy[:,:,1] = ypos
fixedGpts_sens = sens_compl.reshape((num_elems,Order_gpts**2))

leastsquare = _LeastSquare(bpts_xy, fixedGpts_xy,
                fixedGpts_sens, areafraction, 2)
bpts_sens = leastsquare.get_sens_compliance()



 

'''
# LSM properties 
radius = 5
movelimit = 0.5

# sensitivity calculation 
fixedGpts_sens = IntegrationPointSensitivities.flatten(order = 'C')
'''
