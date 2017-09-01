from openmdao.api import ExplicitComponent
from openode.api import ODEFunction
import numpy as np

class ODEComp(ExplicitComponent):

    def initialize(self):
        self.metadata.declare('optim', required=True)
        self.metadata.declare('num', default=1, type_=int)

    def setup(self):
        optim = self.metadata['optim']
        num = self.metadata['num']
        num_boundary_pts = len(optim.phi)
        
        self.add_input('phi', shape=(num, num_boundary_pts))
        self.add_output('dphi_dt', shape=(num, num_boundary_pts))

        self.declare_partials('*', '*', dependent=False)

    def compute(self, inputs, outputs):
        optim = self.metadata['optim']
        num = self.metadata['num']
        # dphi_dt = np.zeros(len(optim.phi))
        dphi_dt = self.metadata['optim'].get_delphi(inputs['phi'])
        outputs['dphi_dt'] = dphi_dt.reshape(1,len(dphi_dt))

class MyODEFunction(ODEFunction):
    
    def __init__(self, system_init_kwargs):
        super(MyODEFunction, self).__init__()

        self.set_system(ODEComp, system_init_kwargs=system_init_kwargs)
        self.declare_state('phi', shape=(len(system_init_kwargs['optim'].phi),), rate_target='dphi_dt', state_targets=['phi'])
