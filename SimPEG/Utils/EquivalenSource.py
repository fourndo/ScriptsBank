# -*- coding: utf-8 -*-
"""
Created on Wed May  9 13:20:56 2018

@author: fourndo@gmail.com


Run an equivalent source inversion

"""
from SimPEG import (
    Mesh, Utils, Maps, Regularization, Regularization,
    DataMisfit, Inversion, InvProblem, Directives, Optimization,
    )
import SimPEG.PF as PF
import numpy as np
import os
import json
from scipy.spatial import Delaunay
from scipy.interpolate import NearestNDInterpolator

workDir = "C:\\Users\\DominiqueFournier\\Dropbox\\Projects\\Synthetic\\Nut_Cracker\\"
outDir = "EquivalentSource\\"
inputFile = "SimPEG_ES.json"
padLen = 100
maxRAM = 1.0
n_cpu = 8

octreeObs = [2, 0]  # Octree levels below observation points
octreeTopo = [0, 1]

meshType = 'TREE'

os.system('mkdir ' + workDir + outDir)

with open(workDir + inputFile, 'r') as f:
    driver = json.load(f)


# Deal with the data
if driver["observed"][0] == 'GRAV':

    survey = Utils.io_utils.readUBCgravityObservations(workDir + driver["observed"][1])

elif driver["observed"][0] == 'MAG':

    survey = Utils.io_utils.readUBCmagneticsObservations(workDir + driver["observed"][1])

else:
    assert False, "Equivalent source only implemented for 'observed' 'GRAV' | 'MAG' "


if "topography" in list(driver.keys()):
    topo = np.genfromtxt(workDir + driver["topography"],
                         skip_header=1)
else:
    topo = False

if "targetChi" in list(driver.keys()):
    targetChi = driver["targetChi"]
else:
    targetChi = 1

# print(list(driver.keys()))
# if "forward" in list(driver.keys()):
#     if driver["forward"][0] == "DRAPE":
#         # Define an octree mesh based on the data

# else:
#     forward = False
locs = survey.srcField.rxList[0].locs
tri = Delaunay(locs[:, :2])
p1, p2 = tri.points[tri.vertices[:, 0], :], tri.points[tri.vertices[:, 1], :]
dl = np.min(np.sum((p1-p2)**2., axis=1)**0.5)

# Compute distance of obs above topo
F = NearestNDInterpolator(topo[:, :2], topo[:, 2])
dz = locs[:, 2] - F(locs[:, 0], locs[:, 1])


h = np.r_[dl.min(), dl.min(), dz.min()]
print(h)
# LOOP THROUGH TILES
surveyMask = np.ones(survey.nD, dtype='bool')
# Going through all problems:
# 1- Pair the survey and problem
# 2- Add up sensitivity weights
# 3- Add to the ComboMisfit

# Create first mesh outside the parallel process
padDist = np.r_[np.c_[padLen, padLen], np.c_[padLen, padLen], np.c_[padLen, 0]]

print("Creating Global Octree")
mesh = Utils.modelutils.meshBuilder(
        topo, h, padDist, meshType='TREE',
        verticalAlignment='center'
    )

if topo is not None:
    mesh = Utils.modelutils.refineTree(
        mesh, topo, dtype='surface',
        nCpad=octreeTopo, finalize=False
    )

mesh = Utils.modelutils.refineTree(
    mesh, locs, dtype='surface',
    nCpad=octreeObs, finalize=True
)

mesh.vectorCCx

surf = Utils.modelutils.activeTopoLayer(mesh, topo)


Mesh.TreeMesh.writeUBC(
      mesh, workDir + outDir + 'OctreeMeshGlobal.msh',
      models={workDir + outDir + 'ActiveSurface.act': surf}
    )


# Get the layer of cells directyl below topo
#surf = Utils.actIndFull2layer(mesh, active)
nC = int(surf.sum())  # Number of active cells
print(nC)
# Create active map to go from reduce set to full
surfMap = Maps.InjectActiveCells(mesh, surf, -100)

# Create identity map
idenMap = Maps.IdentityMap(nP=nC)

# Create static map
if driver["observed"][0] == 'GRAV':
    prob = PF.Gravity.GravityIntegral(
        mesh, rhoMap=idenMap, actInd=surf,
        Jpath=workDir+outDir+"sensitivity.zarr", equiSourceLayer=True
        )
elif driver["observed"][0] == 'MAG':
    prob = PF.Magnetics.MagneticIntegral(mesh, chiMap=idenMap, actInd=surf, equiSourceLayer=True)

# Pair the survey and problem
survey.pair(prob)

wr = prob.getJtJdiag(np.zeros(nC))
wr = (wr/np.max(wr))
wr = wr**0.5

# Create a regularization function, in this case l2l2
reg = Regularization.Sparse(mesh, indActive=surf, mapping=idenMap)
reg.mref = np.zeros(nC)
reg.alpha_z = 0
reg.cell_weights = wr

# Specify how the optimization will proceed, set susceptibility bounds to inf
opt = Optimization.ProjectedGNCG(maxIter=25, lower=-np.inf,
                                 upper=np.inf, maxIterLS=20,
                                 maxIterCG=20, tolCG=1e-3)

# Define misfit function (obs-calc)
dmis = DataMisfit.l2_DataMisfit(survey)
dmis.W = 1./survey.std

# Create the default L2 inverse problem from the above objects
invProb = InvProblem.BaseInvProblem(dmis, reg, opt)

# Specify how the initial beta is found
betaest = Directives.BetaEstimate_ByEig()

# Beta schedule for inversion
betaSchedule = Directives.BetaSchedule(coolingFactor=2., coolingRate=1)

# Target misfit to stop the inversion,
# try to fit as much as possible of the signal, we don't want to lose anything
targetMisfit = Directives.TargetMisfit(chifact=targetChi)

# Put all the parts together
inv = Inversion.BaseInversion(invProb,
                              directiveList=[betaest, betaSchedule, targetMisfit])

# Run the equivalent source inversion
mstart = np.zeros(nC)
mrec = inv.run(mstart)

# Ouput result
Mesh.TreeMesh.writeUBC(
      mesh, workDir + outDir + 'OctreeMeshGlobal.msh',
      models={workDir + outDir + 'EquivalentSource.mod': surfMap * mrec}
    )

if driver["observed"][0] == 'GRAV':

    survey = Utils.io_utils.writeUBCgravityObservations(workDir + outDir + 'Predicted.dat', survey, invProb.dpred)

elif driver["observed"][0] == 'MAG':

    survey = Utils.io_utils.writeUBCmagneticsObservations(workDir + outDir + 'Predicted.dat', survey, invProb.dpred)


# h = np.r_[meshInput.hx.min(), meshInput.hy.min(), meshInput.hz.min()]
# coreX = meshInput.hx == h[0]
# coreY = meshInput.hy == h[1]
# coreZ = meshInput.hz == h[2]

# padx = meshInput.hx[~coreX].sum()
# pady = meshInput.hy[~coreY].sum()
# padz = meshInput.hz[~coreZ].sum()

# padDist = np.r_[np.c_[padx, padx], np.c_[pady, pady], np.c_[padz, padz]]

# print("Creating TreeMesh. Please standby...")
# mesh = Utils.modelutils.meshBuilder(topo, h, padDist,
#                                     meshGlobal=meshInput,
#                                     meshType='TREE',
#                                     verticalAlignment='top')

# mesh = Utils.modelutils.refineTree(mesh, topo, dtype='surface',
#                                    nCpad=[0, 2, 2], finalize=False)

# mesh = Utils.modelutils.refineTree(mesh, xyzLocs, dtype='surface',
#                                    nCpad=[10, 5, 0], finalize=True)

# # mesh = Utils.modelutils.refineTree(mesh, xyzLocs, dtype='surface',
# #                                   nCpad=[0, 10, 0], finalize=True)


# if driver.topofile is not None:
#     # Find the active cells
#     actv = Utils.surface2ind_topo(mesh, topo)
# else:
#     actv = mesh.gridCC[:, 2] <= meshInput.vectorNz[-1]

# nC = int(np.sum(actv))

# # write the mesh and a file of active cells
# Mesh.TreeMesh.writeUBC(
#     mesh,
#     out_dir + 'OctreeMesh_Domtest_' + str(int(h[0])) + '_' + str(int(h[2])) + 'm.msh',
#     models={out_dir + 'ActiveOctree_Domtest_' + str(int(h[0])) + '_' + str(int(h[2])) + 'm.dat': actv}
# )

# # write to VTK for paraview or VisIT.
# # Mesh.TreeMesh.writeVTK(mesh, out_dir + 'OctreeMesh_' + str(int(h[0])) + 'm_3_3_3.vtk')
# #                       models={out_dir + 'ActiveOctree_' + str(int(h[0])) + 'm_3_3_3_vtk.dat': actv})
