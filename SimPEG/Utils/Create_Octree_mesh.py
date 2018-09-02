# -*- coding: utf-8 -*-
"""
Created on Wed May  9 13:20:56 2018

@author: craigm


xyzLocs are the points of topo and obs
h is the core discretization. right now it has to be [x, x, x], (cubes)
padDist is how far out you want to go for the padding
padCore is the number of cells in the three finest octree levels
beyond the xyzLocs
gridLoc changes where the top of topography is, either in the center ('CC') or the top ('N') of the octree mesh
the total mesh width and depth dimensions should be 2^n like, 2,4,8,16,32,64,128,256

"""
from SimPEG import Mesh, Utils
import SimPEG.PF as PF
import numpy as np
import os

work_dir = "C:/Users/DominiqueFournier/ownCloud/Research/Synthetic/Block_Gaussian_topo/"
out_dir = "C:/Users/DominiqueFournier/ownCloud/Research/Synthetic/Block_Gaussian_topo/SimPEG_Susc_Inv/"
input_file = "SimPEG_MAG.inp"

os.chdir(work_dir)
# %%
# Read in the input file which included all parameters at once
# (mesh, topo, model, survey, inv param, etc.)
print('Loading data')
driver = PF.MagneticsDriver.MagneticsDriver_Inv(work_dir + input_file)


# Access the mesh and survey information
meshInput = driver.mesh
survey = driver.survey
actv = driver.activeCells

xyzLocs = survey.srcField.rxList[0].locs.copy()

topo = None
if driver.topofile is not None:
    topo = np.genfromtxt(driver.basePath + driver.topofile,
                         skip_header=1)
#    xyzLocs = np.r_[xyzLocs, topo]

if isinstance(meshInput, Mesh.TensorMesh):
    # Define an octree mesh based on the provided tensor
    h = np.r_[meshInput.hx.min(), meshInput.hy.min(), meshInput.hz.min()]
    coreX = meshInput.hx == h[0]
    coreY = meshInput.hy == h[1]
    coreZ = meshInput.hz == h[2]

    padx = meshInput.hx[~coreX].sum()
    pady = meshInput.hy[~coreY].sum()
    padz = meshInput.hz[~coreZ].sum()

    padDist = np.r_[np.c_[padx, padx], np.c_[pady, pady], np.c_[padz, padz]]

    print("Creating TreeMesh. Please standby...")
    mesh = Utils.modelutils.meshBuilder(topo, h, padDist,
                                        meshGlobal=meshInput,
                                        meshType='TREE',
                                        verticalAlignment='top')

    mesh = Utils.modelutils.refineTree(mesh, topo, dtype='surface',
                                       nCpad=[0, 2, 2], finalize=False)

    mesh = Utils.modelutils.refineTree(mesh, xyzLocs, dtype='surface',
                                       nCpad=[10, 5, 0], finalize=True)

    # mesh = Utils.modelutils.refineTree(mesh, xyzLocs, dtype='surface',
    #                                   nCpad=[0, 10, 0], finalize=True)


else:
    mesh = Mesh.TreeMesh.readUBC(driver.basePath + driver.mshfile)

if driver.topofile is not None:
    # Find the active cells
    actv = Utils.surface2ind_topo(mesh, topo)
else:
    actv = mesh.gridCC[:, 2] <= meshInput.vectorNz[-1]

nC = int(np.sum(actv))

# write the mesh and a file of active cells
Mesh.TreeMesh.writeUBC(
    mesh,
    out_dir + 'OctreeMesh_Domtest_' + str(int(h[0])) + '_' + str(int(h[2])) + 'm.msh',
    models={out_dir + 'ActiveOctree_Domtest_' + str(int(h[0])) + '_' + str(int(h[2])) + 'm.dat': actv}
)

# write to VTK for paraview or VisIT.
# Mesh.TreeMesh.writeVTK(mesh, out_dir + 'OctreeMesh_' + str(int(h[0])) + 'm_3_3_3.vtk')
#                       models={out_dir + 'ActiveOctree_' + str(int(h[0])) + 'm_3_3_3_vtk.dat': actv})
