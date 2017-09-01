import numpy as np 

from openmdao.api import Group, IndepVarComp
from myExplicitComponents import ComplSensComp, NormComp, DispComp
# from lsm_mododule import PyLSM2D

class LSM2D_slpGroup(Group):
    def initialize(self):
        # self.metadata.declare('lsm_module', type_=PyLSM2D, required=True)
        # self.metadata.declare('lsm_mesh', type_=np.ndarray, required=True)
        self.metadata.declare('bpts_x', type_=np.ndarray, required=True)
        self.metadata.declare('bpts_y', type_=np.ndarray, required=True)
        self.metadata.declare('fixedGpts_x', type_=np.ndarray, required=True)
        self.metadata.declare('fixedGpts_y', type_=np.ndarray, required=True)
        self.metadata.declare('fixedGpts_sens', type_=np.ndarray, required=True)
        self.metadata.declare('radius', type_=(int, float), required=True)
        self.metadata.declare('movelimit', type_=float, required=True)

    def setup(self):
        # lsm_module = self.metadata['lsm_module']
        bpts_x = self.metadata['bpts_x']
        bpts_y = self.metadata['bpts_y']
        fixedGpts_sens = self.metadata['fixedGpts_sens']
        radius = self.metadata['radius']
        movelimit = self.metadata['movelimit']
        fixedGpts_x = self.metadata['fixedGpts_x']
        fixedGpts_y = self.metadata['fixedGpts_y']

        num_dvs = 2 # number of lambdas
        num_bpts = bpts_x.shape[0]

        bound_upper = np.ones(num_dvs)
        bound_lower = np.zeros(num_dvs)
        # moveBound = lsm_module.get_bound(lsm_mesh, lambdas,bound_upper,bound_lower,movelimit)        
        
        # inputs (IndepVarComp: component)
        comp = IndepVarComp()
        comp.add_output('lambdas', val = 0.0, shape = num_dvs)
        comp.add_output('fixedGpts_sens', val = fixedGpts_sens)
        comp.add_design_var('lambdas', lower=bound_lower, upper=bound_upper)
        self.add_subsystem('inputs_comp', comp)
        self.connect('inputs_comp.fixedGpts_sens', 'compl_sens_comp.fixedGpts_sens')
        self.connect('inputs_comp.lambdas', 'disp_comp.lambdas')

        # compliance sensitivity at each boundary points
        comp = ComplSensComp(bpts_x = bpts_x, bpts_y = bpts_y , 
            fixedGpts_x = fixedGpts_x, fixedGpts_y = fixedGpts_y, radius = radius)
        self.add_subsystem('compl_sens_comp', comp)
        self.connect('compl_sens_comp.bpts_sens', 'norm_comp.compliance')
        
        # normalize the sensitivity
        comp = NormComp(num_bpts = num_bpts, num_dvs = num_dvs)
        self.add_subsystem('norm_comp', comp)
        self.connect('norm_comp.sens', 'disp_comp.sens')

        # displacement (z) calculation
        comp = DispComp(num_bpts = num_bpts, num_dvs = num_dvs) #lsm_module = lsm_module, lsm_mesh = lsm_mesh)
        self.add_subsystem('disp_comp', comp)
        self.add_objective('disp_comp.distance')
