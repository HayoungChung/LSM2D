import numpy as np

from openmdao.api import ExplicitComponent


class WeightComp(ExplicitComponent):

    def initialize(self):
        self.metadata.declare('num', type_=int, required=True)
        self.metadata.declare('num_nodes_y', type_=int, required=True)

    def setup(self):
        num = self.metadata['num']

        self.add_input('x', shape=num)
        self.add_output('weight')

        derivs = np.ones((1, num)) / num
        self.declare_partials('weight', 'x', val=derivs)

    def compute(self, inputs, outputs):

        outputs['weight'] = sum(inputs['x']) / self.metadata['num']
