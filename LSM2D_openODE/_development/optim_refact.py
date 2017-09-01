import numpy as np
import LinE as LinE
import Sensitivity_vec as Sens
from slsm_Module import *
from Opt_Module import *
import matplotlib.pyplot as plt


class OptimRefact(object):
    def __init__(self, *args):
        super(OptimRefact, self).__init__()

        self.counter = 0

        self.cant = Cantilever(True)
        self.cant.set_fea()
        self.cant.set_slsm()

        self.phi = np.zeros(self.cant.mesh.nNodes)
        for ii in range(0,self.cant.mesh.nNodes):
            self.phi[ii] = self.cant.levelSet.signedDistance[ii]

        Multipliers = vector__double__()
        Multipliers.append(0)

        self.Multipliers = Multipliers

        lambdas = vector__double__()
        lambdas.append(0)
        lambdas.append(0)
        self.lambdas = lambdas
       
        self.constraintDistances = vector__double__()
        self.constraintDistances.append(0)
    def get_phis(self, phi):
        self.phi = phi
        u = self.cant.get_u(self.phi)

        nBpts = len(self.cant.boundary.points)
        BoundaryPoints = np.zeros((nBpts, 2))
        for ii in range(0, nBpts):
            BoundaryPoints[ii, 0] = self.cant.boundary.points[ii].coord.x
            BoundaryPoints[ii, 1] = self.cant.boundary.points[ii].coord.y
        bndSensitivities = self.cant.get_sens(BoundaryPoints)

        for ii in range(0, nBpts):
            self.cant.boundary.points[ii].sensitivities[0] = bndSensitivities[ii]
            self.cant.boundary.points[ii].sensitivities[1] = -1.0

        # output: phi, dphi_dt = from optim output
        self.constraintDistances[0] = (self.cant.lxy[0]*self.cant.lxy[1] * 0.5 - self.cant.boundary.area)
        print("constraints = %d"  % self.constraintDistances[0])

        myOptim = optim.Optimise(self.cant.levelSet, self.cant.boundary.points, \
               self.constraintDistances, self.lambdas, self.Multipliers, 0.9, False)
        # myOptim.preprocess()

        self.cant.levelSet.reinitialise()

        myOptim.solve()

        myOptim.compute_dphi_dt()	
        dphi_dt = np.zeros(self.cant.mesh.nNodes)

        for ii in range(0,self.cant.mesh.nNodes):
            self.phi[ii] = myOptim.phi[ii]
            dphi_dt[ii] = myOptim.dphi_dt[ii]

	    # num_boundary_pts = self.cant.mesh.nNodes
        plt.clf()
        plt.plot(BoundaryPoints[:,0], BoundaryPoints[:,1], 'o')
        plt.savefig('plots/Bpts_%i.png' % self.counter)

        plt.clf()
        plt.scatter(BoundaryPoints[:,0],BoundaryPoints[:,1], 20, bndSensitivities, marker = 'o')
        plt.savefig('plots/Bsens_%i.png' % self.counter)

        plt.clf()
        plt.scatter(self.cant.CMesh.Nodes[:,0], self.cant.CMesh.Nodes[:,1], 20, dphi_dt, marker = 'o')
        plt.savefig('plots/dphi_dt_%i.png' % self.counter)

        plt.clf()
        plt.scatter(self.cant.CMesh.Nodes[:,0], self.cant.CMesh.Nodes[:,1], 20, self.phi, marker = 'o')
        plt.savefig('plots/phi_%i.png' % self.counter)
        self.counter += 1
        print("counter %i" % self.counter)
        # plt.show()
        return dphi_dt


class Cantilever(object):
    def __init__(self, isHoles, *args):
        super(Cantilever, self).__init__()
        self.lxy = lxy = [160, 80]
        self.exy = exy = [160, 80]
        # nELEM = 160*80
        # levelset option
        # AllowedAreaFraction = 0.01
        # moveLimit = 0.9
        # VoidMaterial = 1e-6

        self.CMesh = LinE.FEAMeshQ4(lxy, exy)
        self.mesh  = slsm.Mesh(exy[0], exy[1], False)
        # print("lxy = " + str(lxy))

        if isHoles:
            self.__Holes = self.Holes()
        else:
            self.__Holes = vector_Holes()

    def Holes(self):
        # print ("creating holes")
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

    def set_fea(self, *args):
        lxy = self.lxy
        # print ("setup fea")

        E = 1.0
        v = 0.3
        thickness = 1.0
        # print("(E,v,h) = " + str(E) + "," + str(v) + "," + str(thickness))

        Cijkl = LinE.LinearElasticMaterial.get_Cijkl_E_v(E, v)
        xtip1 = self.CMesh.get_NodeID([lxy[0], int(lxy[1]/2)], 1e-3, 1e-3)
        xtip2 = self.CMesh.get_NodeID([lxy[0], int(lxy[1]/2)-1], 1e-3, 1e-3)
        xtip3 = self.CMesh.get_NodeID([lxy[0], int(lxy[1]/2)+1], 1e-3, 1e-3)

        self.__BC_force1 = self.CMesh.get_dof('y', xtip1)
        self.__BC_force2 = self.CMesh.get_dof('y', xtip2)
        self.__BC_force3 = self.CMesh.get_dof('y', xtip3)

        self.CLinElasticity = LinE.LinearElasticity(self.CMesh, Cijkl)
        self.CSensitivities = Sens.ElasticitySensitivities(self.CLinElasticity)

    def set_slsm(self, *args):
        # print ("setting slsm")
        # initialize levelset and boundary discretization
        self.levelSet = slsm.LevelSet(self.mesh, self.__Holes, 0.9, \
        6, False) # moveLimit is defined elsewhere
        self.levelSet.reinitialise()

        self.boundary = slsm.Boundary(self.levelSet)
        self.boundary.discretise(False)
        self.boundary.computeAreaFractions()
        self.boundary.ComputeNormalVectors()

        # self.initAreaFraction = np.zeros(nELEM)
        # for ii in range(0,nELEM):
        #     if self.mesh.elements[ii].area < 1e-6: #AllowedAreaFraction
        #         self.initAreaFraction[ii] = 1e-6 #VoidMaterial
        #     else:
        #         self.initAreaFraction[ii] = self.mesh.elements[ii].area

    def get_u(self, phi):
        self.phi = phi # grid of lsmMesh
        phi_to_discretise = vector__double__()
        # update phi
        for ii in range(0, len(phi)):
            self.levelSet.signedDistance[ii] = phi[ii]
        self.boundary.discretise(False)
        self.boundary.computeAreaFractions()
        self.boundary.ComputeNormalVectors()

        nELEM = self.CMesh.nELEM
        AreaFraction = np.zeros(nELEM)
        for ii in range(0,nELEM):
            if self.mesh.elements[ii].area < 1e-6: #AllowedAreaFraction
                AreaFraction[ii] = 1e-6 #VoidMaterial
            else:
                AreaFraction[ii] = self.mesh.elements[ii].area

         # print ("Assembling ...")
        self.CLinElasticity.Assembly() 
        # print ("... done")
        Xlo_id = self.CMesh.get_NodeID([0,0],1e-3,np.inf)
        BC_fixed = self.CMesh.get_dof('xy',Xlo_id)

        self.CLinElasticity.Apply_BC(BC_fixed)

        self.CLinElasticity.set_F()
        self.CLinElasticity.set_F(self.__BC_force1,-5)
        self.CLinElasticity.set_F(self.__BC_force2,-2.5)
        self.CLinElasticity.set_F(self.__BC_force3,-2.5)

        # print ("computing displacement")
        u = self.CLinElasticity.solve()

        NodesF = self.CMesh.Nodes + u.reshape(2,self.CMesh.nNODE,order='F').transpose()
        idE = self.CMesh.Elements[:,[0,1,2,3,0]].flatten(order='C').astype(int)
        xorder = self.CMesh.Nodes[idE,0].reshape(int(len(idE)/(self.CMesh._npe+1)),self.CMesh._npe+1)
        yorder = self.CMesh.Nodes[idE,1].reshape(int(len(idE)/(self.CMesh._npe+1)),self.CMesh._npe+1)
        plt.plot(xorder.transpose(),yorder.transpose())
        plt.savefig("DEFORM.png")

        
        return u    
        
    # def set_sens(self, *args):
    #     self.CSensitivities = Sens.ElasticitySensitivities(self.CLinElasticity)
    #         
    def get_sens(self, BoundaryPoints, *args):
        # sensitivity option
        Radius = 2
        Weights = 1
        AllowedAreaFraction = 0.01
        WeightFlag = 5
        BoundarySensitivities = self.CSensitivities.Compliance(BoundaryPoints, Weights, Radius, WeightFlag, AllowedAreaFraction)
        return BoundarySensitivities



        
        
