from libcpp.vector cimport vector
from libcpp cimport bool
from cpython cimport array
import numpy as np
cimport numpy as np


cdef extern from "ma57_solver.h":
  cdef cppclass MA57Solver:
    MA57Solver(bool, bool) except +
    void pyCompute (int n, int nz, int* irn, int* jrn, double* a)
    void pySolve(double* b)

cdef class PyMA57Solver:

    cdef MA57Solver *thisptr
    def __cinit__(self, bool use_metis, bool print_output):
        self.thisptr = new MA57Solver(use_metis, print_output)
    def __dealloc__(self):
        del self.thisptr
    def pyCompute(self, int n, int nz, np.ndarray[int] irn, np.ndarray[int] jrn, np.ndarray[double] a):
        self.thisptr.pyCompute(n, nz, &irn[0], &jrn[0], &a[0])
    def pySolve(self, np.ndarray[double] b):
        self.thisptr.pySolve(&b[0])