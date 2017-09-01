import numpy as np 
import scipy.sparse
import scipy.sparse.linalg

from forces import get_forces
from fem2d import PyFEMSolver
from ma57_pyWrap import PyMA57Solver

num_nodes_x = 81
num_nodes_y = 41

length_x = 2
length_y = 1

E = 1.
nu = 0.3
f = -1.

state_size = 2 * num_nodes_x * num_nodes_y + 2 * num_nodes_y # dof + lag.multiplier
num_nodes = num_nodes_x * num_nodes_y

fem_solver = PyFEMSolver(num_nodes_x, num_nodes_y, length_x, length_y, E, nu)
# nodes = get_mesh(num_nodes_x, num_nodes_y, length_x, length_y)
forces = get_forces(num_nodes_x, num_nodes_y, f=f)
rhs = np.zeros(state_size)
rhs[:disp_size] = forces

size = (num_nodes_x - 1) * (num_nodes_y - 1) * 64 * 4 + 2 * 2 * num_nodes_y
data = np.zeros(size)
rows = np.zeros(size, np.int32)
cols = np.zeros(size, np.int32)

fem_solver.get_stiffness_matrix(np.ones(num_nodes), data, rows, cols)

# sparse solver
mtx = scipy.sparse.csc_matrix((data, (rows, cols)), shape=(state_size, state_size))
ilu = scipy.sparse.linalg.spilu(mtx, drop_tol=1e-14)

def PC(object):
    def __init__(self, ilu):
        self.ilu = ilu
    def __call__(self, rhs):
        return self.ilu.solve(rhs, 'N')
   

pc_op = scipy.sparse.linalg.LinearOperator((size, size), matvec = PC(ilu))
sol[:] = scipy.sparse.linalg.gmres(mtx, rhs, x0=sol, M=pc_op, tol=1e-10, restart=200,)[0]

# ma57 direct solver
ma57 = PyMA57Solver(False, False)
ma57.pyCompute(size, mtx.nonzero(), rows, cols, data)
ma57.pySolve(rhs)

