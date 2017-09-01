import numpy as np

from openmdao.api import ExplicitComponent


class ObjectiveComp(ExplicitComponent):

    def initialize(self):
        self.metadata.declare('w', type_=(int, float), required=True)

    def setup(self):
        w = self.metadata['w']

        self.add_input('weight')
        self.add_input('compliance')
        self.add_output('objective')
        self.declare_partials('objective', 'weight', val=w)
        self.declare_partials('objective', 'compliance', val=1-w)

    def compute(self, inputs, outputs):
        w = self.metadata['w']

        outputs['objective'] = w * inputs['weight'] + (1 - w) * inputs['compliance']
