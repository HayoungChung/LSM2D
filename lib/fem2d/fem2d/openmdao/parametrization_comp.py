import numpy as np

from openmdao.api import ExplicitComponent


class ParametrizationComp(ExplicitComponent):

    def initialize(self):
        self.metadata.declare('num_rows', type_=int, required=True)
        self.metadata.declare('num_cols', type_=int, required=True)
        self.metadata.declare('mtx', required=True)

    def setup(self):
        num_rows = self.metadata['num_rows']
        num_cols = self.metadata['num_cols']
        mtx = self.metadata['mtx']

        self.add_input('x', shape=num_cols)
        self.add_output('y', shape=num_rows)
        self.declare_partials('y', 'x', val=mtx)

    def compute(self, inputs, outputs):
        outputs['y'] = self.metadata['mtx'].dot(inputs['x'])
