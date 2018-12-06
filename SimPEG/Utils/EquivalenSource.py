# -*- coding: utf-8 -*-
"""
Created on Wed May  9 13:20:56 2018

@author: fourndo@gmail.com


Run an equivalent source inversion

"""
from SimPEG import Mesh, Utils
import SimPEG.PF as PF
import numpy as np
import os
import json


workDir = "C:\\Users\\DominiqueFournier\\Dropbox\\Projects\\Synthetic\\Nut_Cracker\\"
inputFile = "SimPEG_ES.json"

os.chdir(workDir)

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
    topo = np.genfromtxt(workDir + driver["observed"],
                         skip_header=1)
else:
    topo = False

print(list(driver.keys()))
if "forward" in list(driver.keys()):
    if driver["forward"][0] == "DRAPE":
        # Define an octree mesh based on the data
        locs = survey.srcField.rxList[0].locs
        dx = locs[:-1, :] - locs[1:, :]
        dl = np.sum(dx**2., axis=1)**0.5

else:
    fwr = False


print(dl.min(), np.mean(dl), np.median(dl))


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
