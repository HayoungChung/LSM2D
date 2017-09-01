import numpy as np

from openmdao.api import ExplicitComponent


class DispComp(ExplicitComponent):

    def initialize(self):
        self.metadata.declare('num_nodes_x', type_=int, required=True)
        self.metadata.declare('num_nodes_y', type_=int, required=True)

    def setup(self):
        num_nodes_x = self.metadata['num_nodes_x']
        num_nodes_y = self.metadata['num_nodes_y']

        state_size = 2 * num_nodes_x * num_nodes_y + 2 * num_nodes_y
        disp_size = 2 * num_nodes_x * num_nodes_y

        self.add_input('states', shape=state_size)
        self.add_output('disp', shape=disp_size)

        self.declare_partials('disp', 'states',
            val=np.ones(disp_size), rows=np.arange(disp_size), cols=np.arange(disp_size))

    def compute(self, inputs, outputs):
        num_nodes_x = self.metadata['num_nodes_x']
        num_nodes_y = self.metadata['num_nodes_y']

        disp_size = 2 * num_nodes_x * num_nodes_y

        outputs['disp'] = inputs['states'][:disp_size]
