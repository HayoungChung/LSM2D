import numpy as np

from openmdao.api import ExplicitComponent


class PenalizationComp(ExplicitComponent):

    def initialize(self):
        self.metadata.declare('p', type_=(int, float), required=True)
        self.metadata.declare('num', type_=int, required=True)

    def setup(self):
        num = self.metadata['num']

        self.add_input('x', shape=num)
        self.add_output('y', shape=num)

        arange = np.arange(num)
        self.declare_partials('y', 'x', rows=arange, cols=arange)

    def compute(self, inputs, outputs):
        p = self.metadata['p']

        outputs['y'] = inputs['x'] ** p

    def compute_partials(self, inputs, partials):
        p = self.metadata['p']

        partials['y', 'x'] = p * inputs['x'] ** (p - 1)
