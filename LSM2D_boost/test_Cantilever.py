from Cantilever import Cantilever
import numpy as np
import matplotlib.pyplot as plt

CC = Cantilever(True)
nELEM = 160 * 80
nNODE = 161 * 81
AllowedAreaFraction = 1e-3
VoidMaterial = 1e-6 # if it is zero: nan is obtained

# testing LINE
CC.set_fea()
u0 = CC.get_u(np.ones(nELEM))
u0_3 = np.reshape(u0,[2,nNODE],order='F').transpose();
scale = 1./100.
plt.scatter(CC.CMesh.Nodes[:,0]+u0_3[:,0]*scale,CC.CMesh.Nodes[:,1]+u0_3[:,1]*scale)
plt.show()

# LSM
CC.set_slsm()
# 1.1. Volume fraction
AreaFraction = np.zeros(nELEM)
for i in range(0,nELEM):
    if CC.mesh.elements[i].area < AllowedAreaFraction:
        AreaFraction[i] = VoidMaterial
    else:
        AreaFraction[i] = CC.mesh.elements[i].area

u1 = CC.get_u(AreaFraction)
u1_3 = np.reshape(u1,[2,nNODE],order='F').transpose();
plt.scatter(CC.CMesh.Nodes[:,0]+u1_3[:,0]*scale,CC.CMesh.Nodes[:,1]+u1_3[:,1]*scale)
plt.show()
