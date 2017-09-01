class optim_post(object):
    def __init__(self, *args):
        super(optim_post, self).__init__(*args)
    
    def computeDisplacements(self,lambdas):
        '''
        displacement = scalefactor * lambda * boundarySensitivity * length
        '''
        pass

    def rescaleDisplacement(self):
        '''
        check for CFL violation
        '''

        pass

    def update_velocity(self):
        '''
        boudarypoints.Velocity = displacement / timestep 
        '''
        pass

    
