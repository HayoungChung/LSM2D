from libcpp.vector cimport vector
from libcpp cimport bool
from cpython cimport array
import numpy as np
cimport numpy as np


cdef extern from "fem_solver.h":
    cdef cppclass FEMSolver:
        FEMSolver(int, int, double, double, double, double) except +
        void get_stiffness_matrix(double * multipliers, double * data, int * rows, int * cols, bool isNodal)
        # void get_stiffness_matrix(double* multipliers, double* data, int* rows, int* cols)
        void get_stiffness_matrix(double * data, int * rows, int * cols)
        void get_stiffness_matrix_derivs(double * states, double * data, int * rows, int * cols, bool isNodal)
        void get_stiffness_matrix_derivs(double * states, double * data, int * rows, int * cols)
        void set_area_fractions(double * areafraction)
        void get_sensitivity_LSTO(double * u, double * xpos, double * ypos, double * sens)

cdef class PyFEMSolver:

    cdef FEMSolver * thisptr

    def __cinit__(self, int num_nodes_x, int num_nodes_y, double length_x, double length_y,
                  double E, double nu,
                  ):
        self.thisptr = new FEMSolver(num_nodes_x, num_nodes_y, length_x, length_y, E, nu)

    def __dealloc__(self):
        del self.thisptr

    # def get_stiffness_matrix(
    #         self, np.ndarray[double] multipliers,
    #         np.ndarray[double] data, np.ndarray[int] rows, np.ndarray[int] cols):
    #     self.thisptr.get_stiffness_matrix( & multipliers[0], & data[0], & rows[0], & cols[0], True)

    def get_stiffness_matrix(
            self, np.ndarray[double] multipliers,
            np.ndarray[double] data, np.ndarray[int] rows, np.ndarray[int] cols, bool isNodal=True):
        self.thisptr.get_stiffness_matrix( & multipliers[0], & data[0], & rows[0], & cols[0], isNodal)

    def get_stiffness_matrix(
            self, np.ndarray[double] data, np.ndarray[int] rows, np.ndarray[int] cols):
        self.thisptr.get_stiffness_matrix(& data[0], & rows[0], & cols[0])

    # def get_stiffness_matrix_derivs(
    #         self, np.ndarray[double] states,
    #         np.ndarray[double] data, np.ndarray[int] rows, np.ndarray[int] cols):
    #     self.thisptr.get_stiffness_matrix_derivs( & states[0], & data[0], & rows[0], & cols[0], True)

    def get_stiffness_matrix_derivs(
            self, np.ndarray[double] states,
            np.ndarray[double] data, np.ndarray[int] rows, np.ndarray[int] cols, bool isNodal=True):
        self.thisptr.get_stiffness_matrix_derivs( & states[0], & data[0], & rows[0], & cols[0], isNodal)

    def get_stiffness_matrix_LSTO(
            self, np.ndarray[double] areafraction, np.ndarray[double] data,
            np.ndarray[int] rows, np.ndarray[int] cols):
        self.thisptr.set_area_fractions(& areafraction[0])
        self.thisptr.get_stiffness_matrix(& data[0], & rows[0], & cols[0])

    def get_sensitivity_LSTO(self, np.ndarray[double] u, np.ndarray[double] xpos,
                             np.ndarray[double] ypos, np.ndarray[double] sens):
        self.thisptr.get_sensitivity_LSTO( & u[0], & xpos[0], & ypos[0], & sens[0])
