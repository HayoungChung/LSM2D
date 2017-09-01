# mymain.py
# many parts are from mymain_v2.py
from openmdao.api import Component, Group, IndepVarComp, Problem
import copy
import numpy as np
import scipy.optimize as sp_optim
import LinE as LinE # module of linE of my own
import Sensitivity_vec as Sens
from slsm_Module import slsm, vector__double__
from Opt_Module import optim
# import time as clock_t
#from ctypes import *
import matplotlib.pyplot as plt

class FEA_M2DO(Component):
    def __init__(self,lxy,exy,*args):
        super(FEA_M2DO,self).__init__()
        self.VoidMaterial = 1e-6

        # Mesh Generation ==============================================
        self.CMesh = LinE.FEAMeshQ4(lxy,exy)        
        self._nELEM = exy[1]*exy[0]
        # BC (moved to init to save somputation time)===================
        # 0. clamped end
        Xlo_id = self.CMesh.get_NodeID([0,0],1e-3,np.inf)
        self.BC_fixed = self.CMesh.get_dof('xy',Xlo_id)
        # 1. Force
        xtip1 = self.CMesh.get_NodeID([lxy[0],int(lxy[1]/2)],1e-3,1e-3)
        xtip2 = self.CMesh.get_NodeID([lxy[0],int(lxy[1]/2)-1],1e-3,1e-3)
        xtip3 = self.CMesh.get_NodeID([lxy[0],int(lxy[1]/2)+1],1e-3,1e-3)

        self.BC_force1 = self.CMesh.get_dof('y',xtip1)
        self.BC_force2 = self.CMesh.get_dof('y',xtip2)
        self.BC_force3 = self.CMesh.get_dof('y',xtip3)
        
        # Material set ===================================================
        E = 1.0 ; v = 0.3# ; thickness = 1.0
        Cijkl = LinE.LinearElasticMaterial.get_Cijkl_E_v(E,v)

        # LinElasticity set ==============================================
        self.CLinElasticity = LinE.LinearElasticity(self.CMesh,Cijkl)

        # parameters and outputs =========================================
        # from slsm
        self.add_param('AreaFraction',val = np.zeros(self.CMesh.nELEM))
        # to sensitivity
        self.add_output('Field',val = np.zeros((self.CMesh.nNODE*self.CMesh._dpn,1)))
        
    def set_AreaFraction(self,ElementArea):
        if len(ElementArea) != self._nELEM:
            print('numbers of ElementArea do not match\n')
            quit()
        for ii in range(0,self._nELEM):
            if ElementArea[ii].area < self.VoidMaterial:
                self.CMesh.AreaFraction[ii] = self.VoidMaterial
            else:
                self.CMesh.AreaFraction[ii] = ElementArea[ii].area

    def set_LinElasticity(self):
        # LinearElasticity Module

        self.CLinElasticity.Assembly()                 
        self.CLinElasticity.Apply_BC(self.BC_fixed)

        self.CLinElasticity.set_F()
        self.CLinElasticity.set_F(self.BC_force1,-5)
        self.CLinElasticity.set_F(self.BC_force2,-2.5)
        self.CLinElasticity.set_F(self.BC_force3,-2.5)
    
    def solve_nonlinear(self,params,unknowns,resids):
        ElementArea = params['AreaFraction']
        self.set_AreaFraction(ElementArea)

        self.set_LinElasticity()
        U = self.CLinElasticity.solve()
        unknowns['Field'] = U

class Sensitivity_M2DO(Component):
    def __init__(self, CLinElasticity,*args):
        super(Sensitivity_M2DO, self).__init__(*args)
        # shallowcopy is used since Sensitivity module does not change Elasticity object
        self.CSensitivity = Sens.ElasticitySensitivities(CLinElasticity)

        # from linearElasticity
        self.add_param('Field',val = np.zeros((CLinElasticity.CMesh.nNODE*CLinElasticity.CMesh._dpn,1)))
        # from slsm
        self.add_param('BoundaryPointCoords', shape = 2)
        self.add_param('AreaFraction',val = np.zeros(CLinElasticity.CMesh.nELEM))

        # to optimization
        self.add_output('BoundaryPointSensitivities',shape = 2)
    
    def solve_nonlinear(self,params,unknowns,resids):
        Field = params['Field']
        if Field != self.CSensitivity.CLinElasticity.Field:
            print('given field does not coincide with CLinElasticity.Field\n')
            quit()
        BoundaryPointCoords = params['BoundaryPointCoords']
        BoundaryPointSensitivities = self.CSensitivities.Compliances\
                (BoundaryPointCoords, 1, 2, 5, 0.01)
        unknowns['BoundaryPointSensitivities'] = BoundaryPointSensitivities
    
class LSM_M2DO(Component):
    def __init__(self, exy, *args):
        super(LSM_M2DO, self).__init__(*args)
        
        # lsm mesh ========================================================
        self.mesh = slsm.Mesh(exy[0],exy[1],False)
        self.meshArea = exy[0]*exy[1]

        # levelset mesh ===================================================
        Holes = self.get_Holes()
        self.levelSet = slsm.LevelSet(self.mesh,Holes,0.9,6,False)
        self.levelSet.reinitialise()
        
        # boundary ========================================================
        self.boundary = slsm.Boundary(self.levelSet)
        self.boundary.discretise(False)
        self.boundary.ComputeNormalVectors()
        self.boundary.computeAreaFractions()
        
        # counter for reinitialization
        self.nReinit = 0
        
        nBpts = self.boundary.points.__len__()

        # to FEA
        self.add_output('AreaFraction',val = np.zeros(exy[0]*exy[1]))
        # to Sens
        self.add_output('BoundaryPointCoords',shape = 2)
        # from Optimization_post (required for update)
        self.add_param('BoundaryPointVelocities',val = np.zeros((nBpts,2)))
        self.add_param('timeStep',val = float(0)) # -lambdas[0]

        # to Optimization_pre                
        self.add_output('ConstraintDistance',val = float(0.0))
        self.add_output('BoundaryPoints', shape=2)
        #self.add_output('BoundaryPointCoords',shape = 2)
        # self.add_output('BoundaryLengths',val = np.zeros(nBpts))
        # self.add_output('BoundaryLimits',val = np.zeros((nBpts,2)))

    def get_Holes(self):
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

        return Holes

    def pre_optimization(self):
        nPoints = boundary.points.__len__()
        displacements = np.zeros()

    def get_AreaFrac_deepcopy(self):
        AreaFrac_deep = np.zeros(self.mesh.elements.__len__())
        for ii in range(0,self.mesh.elements.__len__()):
            AreaFrac_deep = self.mesh.elements[ii].area
        return AreaFrac_deep

    def get_constraintDistance(self):
        constraintDistance = vector__double__()
        # minArea = 0.5
        constraintDistance.append(self.meshArea*0.5 - self.boundary.area)
        return constraintDistance        

    def get_BoundaryCoords(self):
        nBpts = len(self.boundary.points)
        BoundaryPoints = np.zeros([nBpts,2])
    
        for ii in range(0,nBpts):
            BoundaryPoints[ii,0] = self.boundary.points[ii].coord.x
            BoundaryPoints[ii,1] = self.boundary.points[ii].coord.y #VERIFIED
        
        return BoundaryPoints
    
        
    def solve_nonlinear(self, params, unknowns, resids):
        
        BoundaryPointVelocities = params['BoundaryPointVelocities']
        # wrapping to struct data
        for ii in range(0,self.boundary.points.__len__()):
            self.boundary.points[ii].velocity = BoundaryPointVelocities[ii]

        timeStep = params['timeStep']
        self.levelSet.computeVelocities(boundary.points)
        self.levelSet.ComputeGradients()
        isReinitialise = self.levelSet.update(timeStep)

        if isReinitialise == False:
            if nReinit == 20:
                self.levelSet.reinitialise()
            #nReinit += 1
        else:
            nReinit = 0
        
        self.boundary.discretise(False)
        self.boundary.computeAreaFractions()
        
        unknowns['AreaFraction'] = self.get_AreaFrac_deepcopy()        
        unknowns['ConstraintDistance'] = self.get_constraintDistance()
        unknowns['BoundaryPoints'] = copy.deepcopy(self.boundary.points)
        unknowns['BoundaryPointCoords'] = self.get_BoundaryCoords() # to sens
        
class Optim_pre_(Component):
    # def __init__(self, BoundaryPoints, constraintDistance, lambdas, Multipliers, moveLimit, *args):
    def __init__(self,moveLimit, *args):
        super(Optim_pre_, self).__init__(*args)
        # self.BoundaryPoints = BoundaryPoints
        # self.ConstraintDistance = ConstraintDistance
        # self.lambdas = lambdas
        self.Multipliers = vector__double__()
        self.Multipliers.append(0)
        self.moveLimit = moveLimit

        self.add_param('BoundaryPoints', shape=2)
        self.add_param('BoundaryPointSensitivities', shape =2)
        self.add_param('ConstraintDistance', val=float(0))
        self.add_param('lambdas', val = np.zeros(2))
        self.add_output('F_obj',val=float(0))
        self.add_output('F_cons',val=float(0))

    def solve_nonlinear(self,params, unknowns,resids):
        boundarypoints = params['BoundaryPoints']
        boundarypointsSensitivites = params['BoundaryPointSensitivities']
        for ii in range(0,boundarypoints.__len__()):
            boundarypoints[ii].velocity = boundarypointsSensitivites[ii]
        
        lambdas = params['lambdas']
        q_desvar = vector__double__()
        q_desvar.extend(lambdas[0])
        q_desvar.append(lambdas[1])
        constraintDistance = params['ConstraintDistance']
        b = optim.Optimise(boundarypoints, constraintDistance, q_desvar, self.Multipliers, self.moveLimit, False)
        b.preprocess()    
        m = vector__double__()
        m.extend(lambdas[0])
        m.append(lambdas[1])

        # unknowns['negLimit'] = b.negativeLambdaLimits
        # unknowns['posLimit'] = b.positiveLambdaLimits
        self.negLimit = b.negativeLambdaLimits
        self.posLimit = b.positiveLambdaLimits
        unknowns['F_obj'] = b.callback(m,vector__double__(),0)
        unknowns['F_cons'] = b.callback(m,vector__double__(),1)

class Optim_post_(Component):
    def __init__(self, moveLimit, *args):
        super(Optim_post_, self).__init__(*args)
        # self.BoundaryPoints = BoundaryPoints
        # self.constraintDistance = constraintDistance
        # self.lambdas = lambdas
        self.Multipliers = vector__double__()
        self.Multipliers.append(0)       
        self.moveLimit = moveLimit

        # from slsm
        self.add_param('BoundaryPoints', shape=2)
        self.add_param('ConstraintDistance', val=float(0))
        # from optimization
        self.add_param('lambdas', val = np.zeros(2))
        self.add_param('res_fun',val = float(0))
        # to slsm
        self.add_output('BoundaryPointVelocities',shape=2)
        self.add_output('lambdas_rescaled',val=np.zeros(2))
        self.add_output('timeStep', val = float(0) )
    def get_boundaryVels(self):
        BoundaryPointVelocities = np.zeros(boundarypoints.__len__())
        for ii in range(0,boundarypoints.__len__()):
            BoundaryPointVelocities[ii] = boundarypoints[ii].velocity
        return BoundaryPointVelocities

    def solve_nonlinear(self,params, unknowns,resids):
        boundarypoints = params['BoundaryPoints']
        ConstraintDistance = params['ConstraintDistance']
        lambdas = params['lambdas']
        res_fun = params['res_fun']
        b = optim.Optimise(boundarypoints, ConstraintDistance, lambdas, self.Multipliers, self.moveLimit, False)
        b.postprocess(res_fun)    
        unknowns['BoundaryPointVelocities'] = self.get_boundaryVels()
        unknowns['timeStep'] = -np.abs(lambdas[0])
        unknowns['lambdas_rescaled'] = [lambdas[0],lambdas[1]]
