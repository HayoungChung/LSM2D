from Cantilever import Cantilever
import numpy as np 
import matplotlib.pyplot as plt
from openmdao.api import Component, Group, IndepVarComp, Problem, ScipyOptimizer

from slsm_Module import *
from Opt_Module import *

# tofix
from main_mdao_v1 import Optim_pre

CC = Cantilever(True)
nELEM = 160 * 80
nNODE = 161 * 81
AllowedAreaFraction = 1e-3
VoidMaterial = 1e-6 # if it is zero: nan is obtained
CC.set_fea()
meshArea = 160*80
minArea = 0.5

# # testing LinElasticity ================================
# u0 = CC.get_u(np.ones(nELEM))
# u0_3 = np.reshape(u0,[2,nNODE],order='F').transpose();
# scale = 1./100.
# plt.scatter(CC.CMesh.Nodes[:,0]+u0_3[:,0]*scale,CC.CMesh.Nodes[:,1]+u0_3[:,1]*scale)
# plt.show()

# LSM
CC.set_slsm()
# # 1.1. Volume fraction
# AreaFraction = np.zeros(nELEM)
# for i in range(0,nELEM):
#     if CC.mesh.elements[i].area < AllowedAreaFraction:
#         AreaFraction[i] = VoidMaterial
#     else:
#         AreaFraction[i] = CC.mesh.elements[i].area

# u1 = CC.get_u(AreaFraction)
# u1_3 = np.reshape(u1,[2,nNODE],order='F').transpose();
# scale = 1./100.
# plt.scatter(CC.CMesh.Nodes[:,0]+u1_3[:,0]*scale,CC.CMesh.Nodes[:,1]+u1_3[:,1]*scale)
# plt.show()

# nested part
lambdas = vector__double__()
lambdas.append(0)
lambdas.append(0)

Multipliers = vector__double__()
Multipliers.append(0)

class HJ_solver(Component):
    def __init__(self, *args):
        super(HJ_solver, self).__init__()
        self.add_param('AreaFraction', val = np.zeros(nELEM))
        self.add_output('Compliance', val = float(0))
        self.add_output('Area', val = float(0))
        self.add_output('AreaFraction_1', val = np.zeros(nELEM))
        # self.add_output('BoundaryPoints', val = vector__double__()())

    def update(self, timeStep):
        CC.levelSet.computeVelocities(CC.boundary.points)
        CC.levelSet.ComputeGradients()        
        isReinitialised = CC.levelSet.update(timeStep)
        CC.levelSet.reinitialise()

    def discretise(self):
        CC.boundary.discretise(False)
        CC.boundary.computeAreaFractions()
        CC.boundary.ComputeNormalVectors()

        AreaFraction = np.zeros(nELEM)
        for ii in range(0,nELEM):
            if CC.mesh.elements[ii].area < 0.01: #AllowedAreaFraction
                AreaFraction[ii] = 1e-6 #VoidMaterial
            else:
                AreaFraction[ii] = CC.mesh.elements[ii].area

        return AreaFraction
        
    def solve_nonlinear(self, params, unknowns, resides):
        # 1. import levelset
        AreaFraction = params['AreaFraction']
        
        # # 1.1. discretise
        # CC.boundary.discretise(False)
        # CC.boundary.computeAreaFractions();
        # CC.boundary.ComputeNormalVectors();

        # # 2. solve elasticity problem
        # AreaFraction = np.zeros(nELEM)
        # for ii in range(0,nELEM):
        #     if CC.mesh.elements[ii].area < AllowedAreaFraction:
        #         AreaFraction[ii] = VoidMaterial
        #     else:
        #         AreaFraction[ii] = CC.mesh.elements[i].area
        
        u_tmp = CC.get_u(AreaFraction) 
        
        # ## check field ------
        # u_tmp3 = np.reshape(u_tmp,[2,nNODE],order='F').transpose();
        # scale = 1./100.
        # plt.scatter(CC.CMesh.Nodes[:,0]+u_tmp3[:,0]*scale,CC.CMesh.Nodes[:,1]+u_tmp3[:,1]*scale)
        # plt.show() 
        # ----------------------
                
        # 2.1. Sensitivity setup
        nBpts = len(CC.boundary.points)
        BoundaryPoints = np.zeros([nBpts,2])   
    
        for ii in range(0,nBpts):
            BoundaryPoints[ii,0] = CC.boundary.points[ii].coord.x
            BoundaryPoints[ii,1] = CC.boundary.points[ii].coord.y #VERIFIED

        CC.set_sens()
        BoundarySensitivities = CC.get_sens(BoundaryPoints)
        for ii in range(0,nBpts):
            CC.boundary.points[ii].sensitivities[0] = BoundarySensitivities[ii] #VERIFIED
            CC.boundary.points[ii].sensitivities[1] = -1.0

        # 3. nested update  ===============================================================
        constraintDistances = vector__double__()
        constraintDistances.append(meshArea*minArea - CC.boundary.area)
        
        b = optim.Optimise(CC.boundary.points, constraintDistances, lambdas, Multipliers, CC.levelSet.moveLimit, False)
        b.preprocess()

        optim_pre = Optim_pre(b)

        group_1 = Group()
        group_1.add('optim_pre', optim_pre)
        
        top = Problem()
        top.root = group_1

        top.root.add('p1',IndepVarComp('lambdas',np.array([lambdas[0],lambdas[1]])))
        top.root.connect('p1.lambdas','optim_pre.lambdas')

        top.root.deriv_options['type'] = 'fd'

        top.driver = ScipyOptimizer()
        top.driver.options['optimizer'] = 'SLSQP'
        # top.driver.options['tol'] = 1e-8
        # top.driver.options['maxiter'] = 1000

        top.driver.add_desvar('p1.lambdas',\
                lower = np.array([b.negativeLambdaLimits[0], b.negativeLambdaLimits[1]]), \
                upper = np.array([b.positiveLambdaLimits[0], b.positiveLambdaLimits[1]])) # must be changed

        top.driver.add_objective('optim_pre.F_obj')
        top.driver.add_constraint('optim_pre.F_cons',upper=0.0)    # scipy minimizer
        
        top.setup()
        top.run()   

        # 3.1. post-process
        lambdas_x = top['p1.lambdas']
        res = top['optim_pre.F_obj']
        cons = top['optim_pre.F_cons']
        
        lambdas[0] = lambdas_x[0]
        lambdas[1] = lambdas_x[1]
        
        res2 = b.postprocess(res)       
        # ====================================================================================

        # 6. Update
        timeStep = -lambdas[0]
        # time += timeStep
        self.update(timeStep)
        
        CC.boundary.discretise(False)
        CC.boundary.computeAreaFractions()
        CC.boundary.ComputeNormalVectors()

        # Unknowns
        unknowns['Compliance'] = CC.CLinElasticity.get_compliance(u_tmp)
        unknowns['Area'] = CC.boundary.area/meshArea

        AreaFraction1 = np.zeros(nELEM)
        for ii in range(0,nELEM):
            if CC.mesh.elements[ii].area < 1e-6: #AllowedAreaFraction
                AreaFraction1[ii] = 1e-6 #VoidMaterial
            else:
                AreaFraction1[ii] = CC.mesh.elements[ii].area
        unknowns['AreaFraction_1'] = AreaFraction1