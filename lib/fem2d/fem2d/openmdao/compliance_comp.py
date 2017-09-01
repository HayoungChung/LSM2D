import numpy as np

from openmdao.api import ExplicitComponent


class ComplianceComp(ExplicitComponent):

    def initialize(self):
        self.metadata.declare('num_nodes_x', type_=int, required=True)
        self.metadata.declare('num_nodes_y', type_=int, required=True)

    def setup(self):
        num_nodes_x = self.metadata['num_nodes_x']
        num_nodes_y = self.metadata['num_nodes_y']

        disp_size = 2 * num_nodes_x * num_nodes_y

        self.add_input('disp', shape=disp_size)
        self.add_input('forces', shape=disp_size)
        self.add_output('compliance')

    def compute(self, inputs, outputs):
        outputs['compliance'] = np.dot(inputs['disp'], inputs['forces'])

    def compute_partials(self, inputs, partials):
        partials['compliance', 'disp'][0, :] = inputs['forces']
        partials['compliance', 'forces'][0, :] = inputs['disp']
