from HJ_solver import HJ_solver, CC
from Cantilever import Cantilever
from openmdao.api import Component, Group, IndepVarComp, Problem, ScipyOptimizer
import numpy as np 
# import numpy as np
import matplotlib.pyplot as plt

nELEM = 160*80
AreaFraction = np.ones(nELEM)

Group1 = Group()
HJ = HJ_solver()

Group1.add('HJ',HJ)
top = Problem()

top.root = Group1
AreaFraction = CC.initAreaFraction
Group1.add('Indep_var', IndepVarComp('AreaFraction',AreaFraction))
Group1.connect('Indep_var.AreaFraction','HJ.AreaFraction')

for ii in range(0,10):
    top.setup()
    top.run()

    AA = top['HJ.AreaFraction']
    Comp = top['HJ.Compliance']
    Area = top['HJ.Area']

    nBpts = len(CC.boundary.points)
    BoundaryPoints = np.zeros([nBpts,2])
    for ii in range(0,nBpts):
        BoundaryPoints[ii,0] = CC.boundary.points[ii].coord.x
        BoundaryPoints[ii,1] = CC.boundary.points[ii].coord.y #VERIFIED
        
    plt.plot(BoundaryPoints[:,0],BoundaryPoints[:,1],'o')
    plt.show()

    AA1 = top['HJ.AreaFraction_1']

    # plt.plot(AA1-AA)
    # plt.show()