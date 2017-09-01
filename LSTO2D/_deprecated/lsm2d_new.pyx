from libcpp.vector cimport vector
from cpython cimport array
import numpy as np
cimport numpy as np

cdef extern from "LS_toolbox.h":
    cdef cppclass Levelset2D:
        Levelset2D(int, int) except +
        void discretize()
        # void get_area_fractions(double* areafraction)
        # double get_num_boundary_coords()

        
        void get_phi(double* phi)
        void set_signedDistance(double* phi)
        void get_boundary_coords(double* x, double* y)
        void get_delphi(double* delphi)
        void set_sensitivities(double* FEAsensitivities)
        void reinitialize()
        void update()
        void set_lambdas(double*)
        void recover_scales(double*)
        void compute_velocity()
        void update()

        void get_bound(double, double*, double*)        

cdef class PyLSMSolver:
    cdef Levelset2D *thisptr
    cdef __cinit__(self, int num_nodes_x, int num_nodes_y)
        self.thisptr = new Levelset2D(num_nods_x, num_nodes_y)
    def __dealloc__(self):
        del self.thisptr
    def set_lambdas(self, np.ndarray[double] lambdas):
        self.thisptr.set_lambdas(&lambdas[0])
    def recover_scales(self, np.ndarray[double] scales):
        self.thispts.recover_scales(&scales[0])
    def compute_velocity(self):
        self.thispts.compute_velocity()
    def update(self):
        self.thispts.update()
    def get_bound(self, double movelimit, np.ndarray[double] bound_upper, np.ndarray[double] bound_lower):
        self.thispts.get_bound(movelimit, &bound_upper[0], &bound_lower[0])   


'''