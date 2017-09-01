import numpy as np 
import LinE as LinE
from Opt_Module import *
import Sensitivity_vec as Sens
from slsm_Module import *
import matplotlib.pyplot as plt

lxy = [160,80];
exy = [160,80];
nELEM = 160*80
# levelset option
AllowedAreaFraction = 0.01
moveLimit = 0.9
VoidMaterial = 1e-6

# sensitivity option
Radius = 2
Weights = 1
AllowedAreaFraction = 0.01
WeightFlag = 5

class Cantilever(object):
    def __init__(self, isHoles, *args):
        super(Cantilever, self).__init__()
        self.CMesh = LinE.FEAMeshQ4(lxy,exy)
        self.mesh  = slsm.Mesh(exy[0],exy[1],False);
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
        
        # print ("setup fea")

        E = 1.0 
        v = 0.3
        thickness = 1.0
        # print("(E,v,h) = " + str(E) + "," + str(v) + "," + str(thickness))

        Cijkl = LinE.LinearElasticMaterial.get_Cijkl_E_v(E,v)
        xtip1 = self.CMesh.get_NodeID([lxy[0],int(lxy[1]/2)],1e-3,1e-3)
        xtip2 = self.CMesh.get_NodeID([lxy[0],int(lxy[1]/2)-1],1e-3,1e-3)
        xtip3 = self.CMesh.get_NodeID([lxy[0],int(lxy[1]/2)+1],1e-3,1e-3)

        self.__BC_force1 = self.CMesh.get_dof('y',xtip1)
        self.__BC_force2 = self.CMesh.get_dof('y',xtip2)
        self.__BC_force3 = self.CMesh.get_dof('y',xtip3)

        self.CLinElasticity = LinE.LinearElasticity(self.CMesh,Cijkl)
        self.CSensitivities = Sens.ElasticitySensitivities(self.CLinElasticity)
    
    def set_slsm(self, *args):       
        # print ("setting slsm") 
        # initialize levelset and boundary discretization
        self.levelSet = slsm.LevelSet(self.mesh,self.__Holes,moveLimit,6,False) # moveLimit is defined elsewhere
        self.levelSet.reinitialise()

        self.boundary = slsm.Boundary(self.levelSet)
        self.boundary.discretise(False)
        self.boundary.computeAreaFractions()
        self.boundary.ComputeNormalVectors()

        self.initAreaFraction = np.zeros(nELEM)
        for ii in range(0,nELEM):
            if self.mesh.elements[ii].area < 1e-6: #AllowedAreaFraction
                self.initAreaFraction[ii] = 1e-6 #VoidMaterial
            else:
                self.initAreaFraction[ii] = self.mesh.elements[ii].area
    
    def get_u(self, elemArea):        
        for ii in range(0,self.CMesh.nELEM):
            if (elemArea[ii] < VoidMaterial):
                self.CMesh.AreaFraction[ii] = VoidMaterial
    
            else:
                self.CMesh.AreaFraction[ii] = elemArea[ii] 

        # plt.plot(self.CMesh.AreaFraction)
        # plt.show()

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
        
        return u    
        
    def set_sens(self, *args):
        self.CSensitivities = Sens.ElasticitySensitivities(self.CLinElasticity)
    def get_sens(self, BoundaryPoints, *args):
        BoundarySensitivities = self.CSensitivities.Compliance\
                    (BoundaryPoints, Weights, Radius, WeightFlag, AllowedAreaFraction)
        return BoundarySensitivities

        
        
