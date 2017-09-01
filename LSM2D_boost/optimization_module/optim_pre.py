class optim_pre(object):
    def __init__(self, *args):
        super(optim_pre, self).__init__(*args)
        pass
    
    def computeScaleFactors(self): 
        '''
        for effective optimization
        '''

        pass

    def computeLambdaLimits(self): 
        '''
        computes lamba limits by minimum displacements that violates CFL condition 
        for each function
        '''
        pass

    def computeConstraintDistances(self):
        '''
        if constraints are violated: scale constraint distances
        if constraints are satisfied: remove constraints
        '''

        pass


        