from openmdao.api import Component, Group, IndepVarComp, Problem
import numpy as np
import Sensitivity_vec as Sens 
import matplotlib.pyplot as plt

from Cantilever import Cantilever
# todo 
import HJ_solver as HJ # import nested_solver as Vel

cantilever = Cantilever(True)
nELEM = 160 * 80
nNODE = 161 * 81

HJ = HJ(cantilever); # components

Group_HJ = Group()
Group_HJ.add('HJ',HJ);
Group_HJ.add('indep_var',IndepVarComp('AreaFraction',np.zeros(nELEM)))
Group_HJ.connect('indep_var.AreaFraction','HJ.CLinE.AreaFraction')


top = Problem()
top.root = Group_HJ


MAX_ITER = 100
for iITER in range(0,MAX_ITER):
    top.setup()
    top.run_once()

    compliance = Group_HJ.HJ.compliance
    area = Group_HJ.HJ.boundary.area
    boundarypoints = Group_HJ.HJ.boundarypoints

    nBpts = len(boundarypoints)
    BoundaryPoints = np.zeros([nBpts,2])

    for ii in range(0,nBpts):
        BoundaryPoints[ii,0] = boundarypoints[ii].coord.x
        BoundaryPoints[ii,1] = boundarypoints[ii].coord.y #VERIFIED
        
    plt.plot(BoundaryPoints[:,0],BoundaryPoints[:,1],'o')
    plt.show()



    