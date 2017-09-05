#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 13:09:56 2017

@author: hac210
"""

import numpy as np
import scipy.sparse 
import scipy.sparse.linalg

from openmdao.api import Problem, view_model, ScipyOptimizer

from fem2d import PyFEMSolver
from forces import get_forces

from lsm_classes import PyLSMSolver

from ls_sensitivity import _LeastSquare
import scipy.optimize as sp_optim


import matplotlib.pyplot as plt 
# FEM Mesh
num_nodes_x = 161
num_nodes_y = 81
num_nodes = num_nodes_x * num_nodes_y
num_elems = (num_nodes_x-1) * (num_nodes_y-1)

# LSM Mesh
num_param_x = num_nodes_x
num_param_y = num_nodes_y

length_x = num_nodes_x-1
length_y = num_nodes_y-1

num_dofs = num_nodes_x * num_nodes_y * 2 + num_nodes_y * 2

# FEA properties
E = 1.
nu = 0.3
f = -1.
fem_solver = PyFEMSolver(num_nodes_x, num_nodes_y, length_x, length_y, E, nu) #, True)
forces = get_forces(num_nodes_x, num_nodes_y, f=f)
rhs = np.zeros(num_dofs)
rhs[:num_nodes_x*num_nodes_y*2] = forces

Order_gpts = 2
num_gpts = num_elems*Order_gpts**2
 

#num_gpts = 2

# LSM properties 
radius = 2
movelimit = 0.1

# LSM initialize (swisscheese config)
lsm_solver = PyLSMSolver(num_nodes_x, num_nodes_y, 0.5) 

# HJ loop
max_loop = 1
for i_HJ in range(0,max_loop):
    # 0. discretize
    (bpts_xy, areafraction, segmentLength) = lsm_solver.discretize()
    
    # plt.plot(bpts_xy[:,0],bpts_xy[:,1],'o')
    # plt.show()

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
    
    xpos = np.ones(num_elems*4)
    ypos = np.ones(num_elems*4)
    sens_compl = np.ones(num_elems*4)

    fem_solver.get_sensitivity_LSTO(u_fem2d, xpos, ypos, sens_compl) # verified
    xpos = xpos.reshape((num_elems,Order_gpts**2))
    ypos = ypos.reshape((num_elems,Order_gpts**2))
    fixedGpts_xy = np.zeros((num_elems,Order_gpts**2,2))
    fixedGpts_xy[:,:,0] = xpos
    fixedGpts_xy[:,:,1] = ypos
    fixedGpts_sens = sens_compl.reshape((num_elems,Order_gpts**2))

    leastsquare = _LeastSquare(bpts_xy, fixedGpts_xy,
                fixedGpts_sens, areafraction, radius)
    bpts_sens = leastsquare.get_sens_compliance()

    lsm_solver.preprocess(movelimit, bpts_sens)
    # print (lsm_solver.get_constraintDistances())

    (ub, lb) = lsm_solver.get_bounds()
    
    lambdas = np.ones(2)
    
    displacement = lsm_solver.computeDisplacements(lambdas)
    displacement_np = np.asarray(displacement)
    # (obj, obj2) = lsm_solver.computeFunction(displacement_np, 0)
    # (con, con2) = lsm_solver.computeFunction(displacement_np, 1)

    # scipy minimizer
    def objF(x):
        # displacement = lsm_solver.computeDisplacements(x)
        # displacement_np = np.asarray(displacement)
        # return  lsm_solver.computeFunction(displacement_np, 0)[0]
        return lsm_solver.callback(x,0)
        
    def conF(x):
        # displacement = lsm_solver.computeDisplacements(x)
        # displacement_np = np.asarray(displacement)
        # return  lsm_solver.computeFunction(displacement_np, 0)[1]
        return lsm_solver.callback(x,1)
    
    def objF_dp(x):
        # displacement = lsm_solver.computeDisplacements(x)
        # displacement_np = np.asarray(displacement)
        # return  lsm_solver.computeFunction(displacement_np, 0)[0]
        return lsm_solver.callback_dp(x,0)[0]
        
    def conF_dp(x):
        # displacement = lsm_solver.computeDisplacements(x)
        # displacement_np = np.asarray(displacement)
        # return  lsm_solver.computeFunction(displacement_np, 0)[1]
        return lsm_solver.callback_dp(x,1)[0]

    def objF_nocallback(x):
        displacement = lsm_solver.computeDisplacements(x)
        displacement_np = np.asarray(displacement)
        return  lsm_solver.computeFunction(displacement_np, 0)[0]
        
    def conF_nocallback(x):
        displacement = lsm_solver.computeDisplacements(x)
        displacement_np = np.asarray(displacement)
        return  lsm_solver.computeFunction(displacement_np, 1)[1]

    ## comparison
    lambda0 = np.ones(2)*0.0
    print( objF(lambda0), conF(lambda0))
    # print( lsm_solver.callback(lambda0, 0), lsm_solver.callback(lambda0,1))
    print( objF_nocallback(lambda0), conF_nocallback(lambda0))

    lambda0 = np.ones(2)*0.1
    print( objF(lambda0), conF(lambda0))
    print( objF_nocallback(lambda0), conF_nocallback(lambda0))


    displacement = lsm_solver.computeDisplacements(lambda0)
    disp1 = lsm_solver.callback_dp(lambda0,0)[1]
    disp2 = lsm_solver.callback_dp(lambda0,1)[1]

    cons = ({'type' : 'eq', 'fun': lambda x: conF_nocallback(x)})
    res = sp_optim.minimize(objF_nocallback, np.zeros(2), method='SLSQP', options={'disp': True}, \
                   bounds=((lb[0],ub[0]),(lb[1],ub[1])),\
                   constraints = cons)
    ## ## ##


    print(lb, ub)
    lambdas[0] = res.x[0]
    lambdas[1] = res.x[1]
    print(lambdas)
