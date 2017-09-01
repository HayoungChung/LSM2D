from libcpp.vector cimport vector
from cpython cimport array
import numpy as np
cimport numpy as np

cdef extern from "LS_toolbox.h":
    cdef cppclass Levelset2D:
        Levelset2D(int, int, double) except +
        void discretize()
        # void get_area_fractions(double* areafraction)
        # double get_num_boundary_coords()
        void get_phi(double* phi)
        void set_signedDistance(double* phi)
        void get_boundary_coords(double* x, double* y)
        void set_sensitivities(double* FEAsensitivities)
        void reinitialize()
        void update()
        void compute_velocity()
        compute_bndVel(double* lambdas)

        # WIP: TODO
        void get_bound(double movelimit, double* bound_upper, double bound_lower) 
        void set_lambdas(double* lambdas)
        void recover_scale(double* scales)
        void get_segLength(double* segLength)
        void get_lambdaLimits(double* negLambdaLim, double* posLambdaLim)
        
        # MAKE IT DEPRECATED
        void get_delphi(double* delphi)
        
cdef class PyLSMSolver:
    cdef Levelset2D *thisptr
    def __cinit__(self, int num_nodes_x, int num_nodes_y, double maxArea):
        self.thisptr = new Levelset2D(num_nodes_x,num_nodes_y,maxArea)
        self.nNode = num_nodes_x * num_nodes_y
        self.nElem = (num_nodes_x - 1) * (num_nodes_y - 1)
        self.nLambdas = 2
    def __dealloc__(self):
        del self.thisptr
    def discretize(self):
        self.discretize()
        areafraction = np.ndarray[self.nElem]
        self.thisptr.get_area_fractions(&areafraction[0])
        nBpts = self.get_num_boundary_coords()
        posx = np.ndarray[nBpts]
        posy = np.ndarray[nBpts]
        self.get_boundary_coords(&posx[0], &posy[0])
        segLength = np.ndarray[nBpts]
        self.get_segLength(&segLength[0])
                
        bpts_xy = np.zeros((nBpts,2))
        bpts_xy[:,0] = posx
        bpts_xy[:,1] = posy
        return (bpts_xy, areafraction, segLength, )

    def set_sensitivities(self,np.ndarray[double] FEAsensitivities):
        self.thisptr.set_sensitivities(&FEAsensitivities[0])
    
    def get_lambdaLimits(self):
        negLambdaLim = np.ndarray[2]
        posLambdaLim = np.ndarray[2]
        return (negLambdaLim, posLambdaLim, )      

    def get_area_fractions(self, np.ndarray[double] areafraction):
        self.thisptr.get_area_fractions(&areafraction[0])
    def get_num_boundary_coords(self):
        return self.thisptr.get_num_boundary_coords()
    def get_phi(self, np.ndarray[double] phi):
        self.thisptr.get_phi(&phi[0])
    def get_boundary_coords(self, np.ndarray[double] x, np.ndarray[double] y):
        self.thisptr.get_boundary_coords(&x[0], &y[0])
    def get_delphi(self,np.ndarray[double] delphi):
        self.thisptr.get_delphi(&delphi[0])
    def set_signedDistance(self,np.ndarray[double] phi):
        self.thisptr.set_signedDistance(&phi[0])
    def reinitialize(self):
        self.thisptr.reinitialize()
    def update(self):
        self.thisptr.update()
    
