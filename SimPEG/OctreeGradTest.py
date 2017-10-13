# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 09:31:58 2017

@author: DominiqueFournier
"""

import numpy as np
import scipy.sparse as sp
import unittest
from SimPEG import Mesh, Maps, Models, Utils, PF, Regularization, Directives
from SimPEG import InvProblem, Optimization, Inversion, DataMisfit
import inspect
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

mesh = Mesh.TreeMesh([8,8])
mesh.refine(2)

def refine(cell):
    level=3
    xyz = cell.center + mesh.x0
    if np.all([xyz[0]>.6, xyz[0]<0.7,
               xyz[1]>.6, xyz[1]<0.7]):
        return level
    return 0

mesh.refine(refine)
mesh.plotGrid()
plt.show()

mesh.number()

# TODO: Preallocate!
I, J, V = [], [], []
PM = [-1, 1]*mesh.dim  # plus / minus

# TODO total number of faces?
offset = [0]*2 + [mesh.ntFx]*2 + [mesh.ntFx+mesh.ntFy]*2

for ii, ind in enumerate(mesh._sortedCells):

    p = mesh._pointer(ind)
    w = mesh._levelWidth(p[-1])

    if mesh.dim == 2:
        faces = [
                 mesh._fx2i[mesh._index([p[0], p[1], p[2]])],
                 mesh._fx2i[mesh._index([p[0] + w, p[1], p[2]])],
                 mesh._fy2i[mesh._index([p[0], p[1], p[2]])],
                 mesh._fy2i[mesh._index([p[0], p[1] + w, p[2]])]
                ]
    elif mesh.dim == 3:
        faces = [
                 mesh._fx2i[mesh._index([p[0], p[1], p[2], p[3]])],
                 mesh._fx2i[mesh._index([p[0] + w, p[1], p[2], p[3]])],
                 mesh._fy2i[mesh._index([p[0], p[1], p[2], p[3]])],
                 mesh._fy2i[mesh._index([p[0], p[1] + w, p[2], p[3]])],
                 mesh._fz2i[mesh._index([p[0], p[1], p[2], p[3]])],
                 mesh._fz2i[mesh._index([p[0], p[1], p[2] + w, p[3]])]
                ]
    for off, pm, face in zip(offset, PM, faces):
        I += [ii]
        J += [face + off]
        V += [pm]

D = sp.csr_matrix((V, (I, J)), shape=(mesh.nC, mesh.ntF))
R = mesh._deflationMatrix('F', large=True)
S = sp.hstack([sp.identity(mesh.ntFx),
                           sp.csr_matrix((mesh.ntFx, mesh.ntFy))])
Gx = ((D*R).T)

bc = Utils.mkvc(Gx.sum(1) == 0)
Gx=mesh._deflationMatrix('Fx', large=True).T * (S * mesh._deflationMatrix('F', large=True) *Gx)

#Gx_hang = []
#for ii in range(Gx.shape[0]):
#    if Gx[ii,:].count_nonzero() > 2:
#        row = Gx[ii,:].copy()
#        
#        ind = np.where(row.todense()!=0)[1]
#        
#        # The largest absolute value is the cell
#        ind1 = np.argmax(np.abs(row))
#        ind = ind[ind!=ind1]
#        
#        
#        for jj in range(len(ind)):
#            
#            if jj == 0:
#                Gx[ii,ind1] = np.sign(row[0,ind1])*1
#                Gx[ii,ind] = 0
#                Gx[ii,ind[jj]] = np.sign(row[0,ind[jj]])*1
#
#            else:
#                temp = row.copy()
#                temp[0,ind1] = np.sign(row[0,ind1])*1
#                temp[0,ind] = 0
#                temp[0,ind[jj]] = np.sign(row[0,ind[jj]])*1
#        
#                Gx_hang += [temp]
#        
#Gx = sp.vstack([Gx,sp.vstack(Gx_hang)])