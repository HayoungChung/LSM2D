from libcpp.vector cimport vector
from cpython cimport array
# import cython.operator.dereference
import numpy as np
cimport numpy as np

cdef extern from "common.h":
    cdef struct Coord:
        double x 
        double y
    cdef struct BoundaryPoint:
        Coord coord
        double velocity
        vector[double] sensitivities
        double length

cdef extern from "mesh.h":
    cdef struct Element:
        Coord coord
        double area 

    cdef cppclass Mesh:
        Mesh(int, int, bool) except +
        const unsigned int nNodes
        vector[Element] elements
        
cdef extern from "level_set.h":
    cdef cppclass LevelSet:
        LevelSet(Mesh*) except +
        void reinitialise()
        void update(double)
        void computeVelocities(vector[BoundaryPoint]&)
        void computeGradients()
        
        vector[double] signedDistance
        vector[double] velocity


cdef extern from "boundary.h":
    cdef cppclass Boundary:
        Boundary(LevelSet*) except +
        void discretise()
        void computeAreaFractions()
        void computeNormalVectors()
        vector[BoundaryPoint] points
        unsigned int nPoints
        double area

cdef extern from "optimise_noNLOPT.h":
    cdef cppclass Optimise:
        Optimise(vector[BoundaryPoint]&, vector[double]&, vector[double]&, double&, double)

        void computeScaleFactors()
        void computeLambdaLimits()
        void computeConstraintDistances()
        void computeGradients(vector[double]&, vector[double]&, unsigned int)
        void computeDisplacements(vector[double]&)
        # double computeFunction(int)
        double computeFunction(vector[double]&, int)
        double rescaleDisplacements()

        vector[double] scaleFactors
        vector[double] negativeLambdaLimits
        vector[double] positiveLambdaLimits
        vector[unsigned int] indexMap
        vector[double] displacements


cdef class PyLSMSolver:
    cdef Mesh *meshptr # as we don't have nullary constructor, a pointer must be used
    cdef LevelSet *levelsetptr
    cdef Boundary *boundaryptr
    cdef Optimise *optimiseptr
    cdef int nNode 
    cdef int nElem 
    cdef double targetArea
    cdef int nBpts
    # cdef vector[BoundaryPoint]

    def __cinit__(self, int num_nodes_x, int num_nodes_y, double minArea):
        self.meshptr = new Mesh(num_nodes_x-1,num_nodes_y-1,False)
        self.levelsetptr = new LevelSet(self.meshptr)
        self.boundaryptr = new Boundary(self.levelsetptr)
        self.boundaryptr.discretise()
        
        meshArea = (num_nodes_x-1)*(num_nodes_y-1)
        self.nNode = num_nodes_x * num_nodes_y
        self.nElem = (num_nodes_x-1) * (num_nodes_y-1)
        self.targetArea = meshArea*minArea

    def __dealloc__(self):
        del self.meshptr  
        del self.levelsetptr  
        del self.boundaryptr
        del self.optimiseptr  

    def discretize(self):
        self.boundaryptr.discretise()
        self.boundaryptr.computeAreaFractions()
        areafraction = np.ndarray(self.nElem)
        self.nBpts = self.boundaryptr.nPoints
        for ii in range(0,self.nElem):
            areafraction[ii] = self.meshptr.elements[ii].area
            if areafraction[ii] < 1e-3:
                areafraction[ii] = 1e-3

        bpts_xy = np.ndarray((self.nBpts,2))
        segLength = np.zeros(self.nBpts)

        for ii in range(0, self.boundaryptr.nPoints):
            bpts_xy[ii,0] = self.boundaryptr.points[ii].coord.x
            bpts_xy[ii,1] = self.boundaryptr.points[ii].coord.y
            segLength[ii] = self.boundaryptr.points[ii].length

        return (bpts_xy, areafraction, segLength) 
        
    # this is a temporary fix =======================================
    def preprocess(self, double movelimit, np.ndarray[double] BptsSensitivity):
        cdef vector[double] lambdas
        lambdas.push_back(0.0)
        lambdas.push_back(0.0)
        
        for ii in range(0, self.boundaryptr.nPoints):
            self.boundaryptr.points[ii].sensitivities[0] = BptsSensitivity[ii]
            self.boundaryptr.points[ii].sensitivities[1] = -1.0 
        
        cdef vector[double] constraintDistances 
        constraintDistances.push_back(self.targetArea-self.boundaryptr.area)

        self.optimiseptr = new Optimise(
            self.boundaryptr.points, constraintDistances, lambdas, np.abs(lambdas[0]), movelimit)

        self.optimiseptr.computeConstraintDistances();

        # // Initializing index map making all constraints to be active.
        # // Compute the scale factors for the objective and constraints.
        self.optimiseptr.computeScaleFactors()

        # // Compute the lambda limits.
        self.optimiseptr.computeLambdaLimits()
        
        # // Compute scaled constraint change distances and remove inactive
        # // inequality constraints.
        # self.optimiseptr.computeConstraintDistances()

    def get_bounds(self):        
        return (self.optimiseptr.positiveLambdaLimits, self.optimiseptr.negativeLambdaLimits)

    def computeDisplacements(self, np.ndarray[double] lambdas):
        cdef vector[double] lambda_v
        lambda_v.push_back(lambdas[0])
        lambda_v.push_back(lambdas[1])
        self.optimiseptr.computeDisplacements(lambda_v)
        return self.optimiseptr.displacements

    def computeFunction(self,np.ndarray[double] displacement, int index):
        return self.optimiseptr.computeFunction(displacement, index)
        
    def computeGradients(self, np.ndarray[double] lambdas,unsigned int index):
        ndvs = lambdas.shape[0]
        gradient = np.ndarray(ndvs,dtype=float) 
        self.optimiseptr.computeGradients(lambdas, gradient, index)
        return gradient

    def postprocess(self,np.ndarray[double] lambdas):
        # // Compute the optimum displacement vector.
        cdef vector[double] lambda_v
        lambda_v.push_back(lambdas[0])
        lambda_v.push_back(lambdas[1])
        self.optimiseptr.computeDisplacements(lambda_v)
        
        # return gradient        

        # # // Rescale the displacements and lambda values (if necessary).
        self.optimiseptr.rescaleDisplacements()
        # # // Calculate the unscaled lambda values.
        for ii in range(0,2):
            lambda_v[ii] *= self.optimiseptr.scaleFactors[ii];

        # # // Effective time step.
        timeStep = -(lambda_v[0])
        # # // Calculate boundary point velocities.
        for ii in range(0, self.nBpts):
            self.boundaryptr.points[ii].velocity = self.optimiseptr.displacements[ii] / timeStep
        # print("?")
        # return self.boundaryptr.points # ERROR: pointer issue
        # ======================================================================================
        
    def computeVelocities(self):
        self.levelsetptr.computeVelocities(self.boundaryptr.points)
        self.levelsetptr.computeGradients()
        return self.levelsetptr.velocity

    def reinitialise(self):
        self.levelsetptr.reinitialise()
        return self.levelsetptr.signedDistance

    def update(self,double t):
        self.levelsetptr.update(t)
        return self.levelsetptr.signedDistance

