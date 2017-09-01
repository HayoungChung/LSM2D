# mymain.py

import numpy as np
import scipy.optimize as sp_optim
import LinE as LinE # module of linE of my own
import Sensitivity_vec as Sens
from slsm_Module import *
import time as clock_t
#from ctypes import *
import matplotlib.pyplot as plt
from Opt_Module import *
from openmdao.api import Component, Group, IndepVarComp, Problem, ScipyOptimizer

class Optim_pre(Component):
    # called for every iterations
    def __init__(self, optimClass): # boundarypoints, constraintDistances, Multipliers, moveLimit):
        super(Optim_pre,self).__init__()
        self.optimise = optimClass        
        self.add_param('lambdas', val = np.zeros(2))
        self.add_output('F_obj', val = float(0))
        self.add_output('F_cons', val = float(0))
        
        # self.boundarypoints = boundarypoints
        # self.constraintDistances = constraintDistances
        # self.Multipliers = Multipliers
        # self.moveLimit = moveLimit

    def solve_nonlinear(self, params, unknowns, resides):
        lambdas = params['lambdas']
        q_desvar = vector__double__()
        q_desvar.append(lambdas[0])
        q_desvar.append(lambdas[1])
        #b = optim.Optimise(self.boundarypoints, self.constraintDistances, \
                        # q_desvar, self.Multipliers, self.moveLimit, False)
        # self.optimise.preprocess()
        # self.negLimit = [self.optimise.negativeLambdaLimits[0],self.optimise.negativeLambdaLimits[1]]
        # self.posLimit = [self.optimise.positiveLambdaLimits[0],self.optimise.positiveLambdaLimits[1]]

        unknowns['F_obj'] = self.optimise.callback(q_desvar,vector__double__(),0)
        unknowns['F_cons'] = self.optimise.callback(q_desvar,vector__double__(),1)

# 0. setup
AllowedAreaFraction = 0.01
moveLimit = 0.9
maxTime = 200
minArea = 0.5 
VoidMaterial = 1e-6
 
# 1. Define domains
lxy = [160,80]
exy = [160,80]
CMesh = LinE.FEAMeshQ4(lxy,exy)
mesh  = slsm.Mesh(exy[0],exy[1],False)
meshArea = exy[0]*exy[1]

Holes = slsm.vector__Holes__()
Holes.append(slsm.Hole(16, 14, 5))
Holes.append(slsm.Hole(32, 27, 5))
Holes.append(slsm.Hole(48, 14, 5))
Holes.append(slsm.Hole(64, 27, 5))
Holes.append(slsm.Hole(80, 14, 5))
Holes.append(slsm.Hole(96, 27, 5))
Holes.append(slsm.Hole(112, 14, 5))
Holes.append(slsm.Hole(128, 27, 5))
Holes.append(slsm.Hole(144, 14, 5))
Holes.append(slsm.Hole(16, 40, 5))
Holes.append(slsm.Hole(32, 53, 5))
Holes.append(slsm.Hole(48, 40, 5))
Holes.append(slsm.Hole(64, 53, 5))
Holes.append(slsm.Hole(80, 40, 5))
Holes.append(slsm.Hole(96, 53, 5))
Holes.append(slsm.Hole(112, 40, 5))
Holes.append(slsm.Hole(128, 53, 5))
Holes.append(slsm.Hole(144, 40, 5))
Holes.append(slsm.Hole(16, 66, 5))
Holes.append(slsm.Hole(48, 66, 5))
Holes.append(slsm.Hole(80, 66, 5))
Holes.append(slsm.Hole(112, 66, 5))
Holes.append(slsm.Hole(144, 66, 5))

levelSet = slsm.LevelSet(mesh,Holes,moveLimit,6,False) # moveLimit is defined elsewhere
levelSet.reinitialise()

# 1.1. Volume fraction
for i in range(0,CMesh.nELEM):
    if mesh.elements[i].area < VoidMaterial:
        CMesh.AreaFraction[i] = VoidMaterial
    else:
        CMesh.AreaFraction[i] = mesh.elements[i].area

boundary = slsm.Boundary(levelSet)
boundary.discretise(False)
boundary.computeAreaFractions()
boundary.ComputeNormalVectors()

io = slsm.InputOutput()

# 2. Lienar Elasticity
E = 1.0 
v = 0.3
thickness = 1.0

Cijkl = LinE.LinearElasticMaterial.get_Cijkl_E_v(E,v)
xtip1 = CMesh.get_NodeID([lxy[0],int(lxy[1]/2)],1e-3,1e-3)
xtip2 = CMesh.get_NodeID([lxy[0],int(lxy[1]/2)-1],1e-3,1e-3)
xtip3 = CMesh.get_NodeID([lxy[0],int(lxy[1]/2)+1],1e-3,1e-3)

BC_force1 = CMesh.get_dof('y',xtip1)
BC_force2 = CMesh.get_dof('y',xtip2)
BC_force3 = CMesh.get_dof('y',xtip3)

CLinElasticity = LinE.LinearElasticity(CMesh,Cijkl)


CSensitivities = Sens.ElasticitySensitivities(CLinElasticity)

for ii in range(0,CMesh.nELEM):
    if (mesh.elements [ii].area < VoidMaterial):
        CMesh.AreaFraction[ii] = VoidMaterial
    else:
        CMesh.AreaFraction[ii] = mesh.elements [ii].area #VERIFIED
        
# 4. 1st nest for boundary evolution

Radius = 2
Weights = 1
AllowedAreaFraction = 0.01
WeightFlag = 5
nReinit = 0

meshArea = exy[0]*exy[1]
Max_Iter = 200

# results 
times = [0]
time = 0
Objectives = [None]
areas = [None]

lambdas = vector__double__()
lambdas.append(0)
lambdas.append(0)

#lambdas2 = vector__double__()
#lambdas2.append(0)
#lambdas2.append(0)

Multipliers = vector__double__()
Multipliers.append(0)
#%%
t = 0.0
for iIter in range(0,Max_Iter):
    #%%
  # =========== LINEAR ELASTICITY MODULE ===============================================
#    t = clock_t.time()
    CLinElasticity.Assembly() 
#    print("Assembly takes %d sec\n",clock_t.time()-t)
    # BCs 
    Xlo_id = CMesh.get_NodeID([0,0],1e-3,np.inf)
    BC_fixed = CMesh.get_dof('xy',Xlo_id)

    CLinElasticity.Apply_BC(BC_fixed)

    CLinElasticity.set_F()
    CLinElasticity.set_F(BC_force1,-5)
    CLinElasticity.set_F(BC_force2,-2.5)
    CLinElasticity.set_F(BC_force3,-2.5)

#    t = clock_t.time()
    u = CLinElasticity.solve()
#    compliance = CLinElasticity.get_compliance(u) #VERIFIED
#    print("Linear solver takes %d sec\n",clock_t.time()-t)
    #%%
    nBpts = len(boundary.points)
    BoundaryPoints = np.zeros([nBpts,2])
    
    
    for ii in range(0,nBpts):
        BoundaryPoints[ii,0] = boundary.points[ii].coord.x
        BoundaryPoints[ii,1] = boundary.points[ii].coord.y #VERIFIED
        
    plt.plot(BoundaryPoints[:,0],BoundaryPoints[:,1],'o')
    plt.show()
#    t = clock_t.time()
    BoundarySensitivities = CSensitivities.Compliance(BoundaryPoints, Weights, Radius, WeightFlag, AllowedAreaFraction)
#    print("bpt sensitivity calculation takes %d sec\n",clock_t.time()-t)

    for ii in range(0,nBpts):
        boundary.points[ii].sensitivities[0] = BoundarySensitivities[ii] #VERIFIED
        boundary.points[ii].sensitivities[1] = -1.0
    
    #%%
    constraintDistances = vector__double__()
    constraintDistances.append(meshArea*minArea - boundary.area)
    
    ## should be changed to MDAO =============================
#    print("Attempt to optimise\n")    
#    timeStep = float()
#    optimise = slsm.Optimise(boundary.points, constraintDistances, lambdas, Multipliers, timeStep)
#    optimise.solve()
#    
#    t = clock_t.time()
#    a = optim.Optimise(boundary.points, constraintDistances, lambdas2, Multipliers, levelSet.moveLimit, False)
#    a.preprocess()
#    obj = a.nlopt_activate()
#    obj2 = a.postprocess(obj)
#    print("optimization takes %d sec\n",clock_t.time()-t)
#    
#    t = clock_t.time()
    # b = optim.Optimise(boundary.points, constraintDistances, lambdas, Multipliers, levelSet.moveLimit, False)
    # b.preprocess()    

    b = optim.Optimise(boundary.points, constraintDistances, lambdas, Multipliers, levelSet.moveLimit, False)
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
        #    print("optimization takes %d sec\n",clock_t.time()-t)
    
#    t = clock_t.time()
#    b = slsm_Module()
#    timeStep = b.wrapper_optimise(boundary.points, constraintDistances, lambdas2, Multipliers, levelSet.moveLimit)
#    print("optimization takes %d sec\n",clock_t.time()-t)

    #%%
    lambdas_x = top['p1.lambdas']
    res = top['optim_pre.F_obj']
    cons = top['optim_pre.F_cons']
    
    lambdas[0] = lambdas_x[0]
    lambdas[1] = lambdas_x[1]
    
    res2 = b.postprocess(res)
    
    timeStep = 1 # -lambdas[0]

#    t = clock_t.time()
    levelSet.computeVelocities(boundary.points)
    levelSet.ComputeGradients()
    
    isReinitialised = levelSet.update(timeStep)
    
    if isReinitialised == False:
        if nReinit == 20:
            levelSet.reinitialise()
        nReinit += 1
    else:
        nReinit = 0
    
    nReinit += 1

    boundary.discretise(False)
    boundary.computeAreaFractions()


    # export Area fraction to FEA
    # same as 1.1. Volume fraction
    for ii in range(0,CMesh.nELEM):
        if (mesh.elements [ii].area < VoidMaterial):
            CMesh.AreaFraction[ii] = VoidMaterial
        else:
            CMesh.AreaFraction[ii] = mesh.elements [ii].area
            
    boundary.ComputeNormalVectors()
    time += timeStep
    
    times.append(time)
#    print("update takes %d sec\n",clock_t.time()-t)
    compliance = CLinElasticity.get_compliance(u)
    print(iIter, compliance, boundary.area/meshArea)
    
#    io.saveLevelSetVTK(times.size(), levelSet); #TOFIX
#    io.saveBoundarySegmentsTXT(times.size(), boundary); #TOFIX
    
    Objectives.append(compliance)
    areas.append(boundary.area/meshArea)
    

U2 = CLinElasticity.Field.reshape(CMesh._dpn,int(CMesh.nNODE),order='F').transpose()
node_f = CMesh.Nodes + U2

idE = CMesh.Elements[:,[0,1,2,3,0]].flatten(order='C').astype(int)
xorder = node_f[idE,0].reshape(int(len(idE)/(CMesh._npe+1)),CMesh._npe+1)
yorder = node_f[idE,1].reshape(int(len(idE)/(CMesh._npe+1)),CMesh._npe+1)
plt.plot(xorder.transpose(),yorder.transpose(),'b-')
