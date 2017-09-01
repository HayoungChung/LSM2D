from libcpp.vector cimport vector
from cpython cimport array
import numpy as np
cimport numpy as np

cdef extern from "LS_toolbox.h":
    cdef cppclass Levelset2D:
        Levelset2D(int, int, double) except +
        void get_area_fractions(double* areafraction)
        double get_num_boundary_coords()
        void get_phi(double* phi)
        void set_signedDistance(double* phi)
        void get_boundary_coords(double* x, double* y)
        void get_delphi(double* delphi)
        void set_sensitivities(double* FEAsensitivities)
        void reinitialize()
        void update()

cdef class pyLevelset:
    cdef Levelset2D *thisptr
    def __cinit__(self, int num_nodes_x, int num_nodes_y, double maxArea):
        self.thisptr = new Levelset2D(num_nodes_x,num_nodes_y,maxArea)
    def __dealloc__(self):
        del self.thisptr
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
    def set_sensitivities(self,np.ndarray[double] FEAsensitivities):
        self.thisptr.set_sensitivities(&FEAsensitivities[0])
    def set_signedDistance(self,np.ndarray[double] phi):
        self.thisptr.set_signedDistance(&phi[0])
    def reinitialize(self):
        self.thisptr.reinitialize()
    def update(self):
        self.thisptr.update()
