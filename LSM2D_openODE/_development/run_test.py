from openmdao.api import Problem
from openode.api import ExplicitTMIntegrator, ForwardEuler
from ODEfunction import MyODEFunction
import numpy as np
from optim_refact_v3 import OptimRefact


num = 2

optim = OptimRefact()
system_init_kwargs = {'optim': optim}
initial_phi = optim.phi
optim.dphi_dt = optim.get_delphi(initial_phi)
ode_function = MyODEFunction(system_init_kwargs)

phi2 = optim.dphi_dt + initial_phi
# opti.get_delphi(phi2)

# integrator = ExplicitTMIntegrator(
#     ode_function=ode_function, 
#     scheme=ForwardEuler(),
#     initial_conditions={'phi': initial_phi}, time_spacing=np.arange(num),
#     start_time=0.,
#     end_time=1.,
# )

# prob = Problem(integrator)
# prob.setup()
# prob.run_model()
