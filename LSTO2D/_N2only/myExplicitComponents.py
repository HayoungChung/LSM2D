import numpy as np
from openmdao.api import ExplicitComponent

# from lsm_module import PyLSM2D

class ComplSensComp(ExplicitComponent):
    def initialize(self):
        self.metadata.declare('bpts_x', type_=np.ndarray, required=True)
        self.metadata.declare('bpts_y', type_=np.ndarray, required=True)
        self.metadata.declare('fixedGpts_x', type_=np.ndarray, required=True)
        self.metadata.declare('fixedGpts_y', type_=np.ndarray, required=True)
        self.metadata.declare('radius', type_=(int, float), required=True)
            
    def setup(self):
        bpts_x = self.metadata['bpts_x']
        bpts_y = self.metadata['bpts_y']
        fixedGpts_x = self.metadata['fixedGpts_x']
        fixedGpts_y = self.metadata['fixedGpts_y']
        radius = self.metadata['radius']

        num = fixedGpts_x.shape[0]
        self.add_input('fixedGpts_sens', shape=num, val=0.0)
        self.add_output('bpts_sens', shape=num, val=0.0)

        ''' 
        WIP:this part computes a mtx for least square 
        # function of radius, bpts ...
        # self.declare_partials()
        '''
        
    def compute(self, inputs, outputs):
        outputs['bpts_sens'] = self.mtx.dot(inputs['fixedGpts_sens'])

class NormComp(ExplicitComponent):
    def initialize(self):
        self.metadata.declare('num_bpts', type_=int, required=True)
        self.metadata.declare('num_dvs', type_=int, required=True)
        pass
    def setup(self):             
        num_bpts = self.metadata['num_bpts']
        num_dvs = self.metadata['num_dvs']
        self.add_input('compliance', shape = num_bpts)
        self.add_output('sens', shape = (num_bpts, num_dvs))
        ''' 
        WIP: normalize 
        # self.declare_partials()
        '''
    def compute(self, inputs, outputs):
        outputs['sens'] = self.mtx.dot(inputs['compliance'])

class DispComp(ExplicitComponent):
    def initialize(self):
        self.metadata.declare('num_bpts', type_=int, required=True)
        self.metadata.declare('num_dvs', type_=int, required=True)
        # self.metadata.declare('lsm_module', type_=PyLSM2D, required=True)
        # self.metadata.declare('lsm_mesh', type_=np.ndarray, required=True)
        pass
    def setup(self):
        num_bpts = self.metadata['num_bpts']
        num_dvs = self.metadata['num_dvs']
        # lsm_module = self.metadata['lsm_module']
        # lsm_mesh = self.metadata['lsm_mesh']
    
        self.add_input('sens', shape = (num_bpts,num_dvs))
        self.add_input('lambdas', shape = num_dvs)
        self.add_output('distance', shape = num_bpts)

        '''
        WIP: get length of the segment
        '''
        self.seglength = 0
    def compute(self, inputs, outputs):
        outputs['distance'] = inputs['sens'][0]*inputs['lambdas'][0]*self.seglength + \
                            inputs['sens'][1]*inputs['lambdas'][1]*self.seglength 
        

        