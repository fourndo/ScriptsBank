# -*- coding: utf-8 -*-
"""
Created on Sun Nov 13 08:31:12 2016

@author: dominiquef
"""

import scipy.io
from SimPEG import Maps, Mesh, mkvc, Utils
import numpy.matlib as npm
import numpy as np
import scipy.sparse as sp
import SimPEG.PF as PF
import SimPEG.EM.Static.DC as DC
import SimPEG.EM.FDEM as FDEM
import SimPEG.EM.TDEM as TDEM
from SimPEG.EM.Static.Utils import drapeTopotoLoc, gen_DCIPsurvey, writeUBC_DCobs, readUBC_DC3Dobs
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import time
import gc
from test_TDEM_forward_Analytic import halfSpaceProblemAnaDiff
from shutil import copyfile
import os
import glob
import fileinput
import re


work_dir = 'C:\\LC\\Private\\dominiquef\\Projects\\4414_Minsim\\Modeling\\Forward'
topofile = 'ROT_DEM_30m.topo'


inp_dir = {'DC': '\\FWR_DC',
           'TDEM':'\\FWR_TDEM\\OUT_DIR',}

meshfile = 'Mesh_20m.msh'
mesh_global = Mesh.TensorMesh.readUBC(work_dir + '\\' + meshfile)

modelfile = {'MAG':'Mesh_20m_Susc\\Randn_Std_model.sus','DC':'Mesh_20m_Cond\\Randn_Std_model.con','Gz':'Mesh_20m_Dens\\Density_Randn_Std.den'}
nullcell = 'nullcell.dat'

topo = np.genfromtxt(work_dir + '\\' + topofile,
                     skip_header=1)

actvmod = mesh_global.readModelUBC(work_dir + '\\' + nullcell)
actv = np.ones(mesh_global.nC,dtype=bool)
actv[actvmod==0] = False

obsfile = 'Mesh_20m_Susc\\Mag_grid_50m.dat'
topofile = 'Mesh_20m_Susc\\ROT_DEM_30m.topo'

dType = 'DC'

ndv = {'MAG': 0, 'Gz': 0, 'DC': 1e-8, 'FDEM': 1e-8, 'TDEM': 1e-8}

npadxy = {'DC': 10, 'TDEM': 7, }
npadz = {'DC': 10, 'TDEM': 10, }
expf = 1.3

offTile = 500
lenTile = 500

# Original grid from "DataSimulator"
gCx = np.asarray(range(9))*lenTile + offTile + mesh_global.x0[0]
gCy = np.asarray(range(23))*lenTile + offTile + mesh_global.x0[1]

X, Y = np.meshgrid(gCx, gCy)
X, Y = mkvc(X), mkvc(Y)

dx = [20., 40.]

# Specify which stations to use
tID = [36, 39]

core = {'DC': np.r_[np.ones(80)*dx[0]],
        'TDEM': np.r_[np.ones(5)*dx[1], [33, 26],
                      np.ones(10)*dx[0], [26, 33],
                      np.ones(5)*dx[1]]}

if dType == 'DC':

    # # Fix bug in UBC export
    # os.chdir(work_dir + inp_dir[dType])
    # source_dir = os.getcwd()
    # target_dir = "clean"
    # source_files = [fname for fname in glob.glob(os.path.join(source_dir, "*.obs"))]

    # # check if target directory exists... if not, create it.
    # if not os.path.exists(target_dir):
    #      os.makedirs(target_dir)

    # for source_file in source_files:
    #     target_file = os.path.join(target_dir,os.path.basename(source_file))
    #     with open(source_file, 'r') as sfile:
    #         with open(target_file, 'w') as tfile:
    #             lines = sfile.readlines()
    #             # do the replacement in the second line.
    #             # (remember that arrays are zero indexed)
    #             lines[2] = re.sub("E\+0212100", 'E+02 12100', lines[2])
    #             tfile.writelines(lines)

    # print("DONE")

    os.chdir(work_dir + inp_dir[dType] + '\\' + target_dir)
    DcData = []
    files = glob.glob("*Tile*")
    txLoc = []

    print('Loading data')
    for file in files:

        obj = readUBC_DC3Dobs(file)

        DcData.append(obj['DCsurvey'])

        txLoc.append(obj['DCsurvey'].srcList[0].loc[0])

    print('Load complete!')
    txLoc = np.asarray(txLoc)

    # Grab the pole location for supplied grid ID
    endl = np.c_[np.r_[X[tID]], np.r_[Y[tID]]]

    # Create grid of survey electrodes
    grid_x1, grid_y1 = np.meshgrid(np.linspace(X[tID[0]]-1000., X[tID[0]]+1000., 40.),
                                   np.linspace(Y[tID[0]]+475., Y[tID[0]]+975., 10.))

    grid_x2, grid_y2 = np.meshgrid(np.linspace(X[tID[0]]-1000., X[tID[0]]+1000., 40.),
                                   np.linspace(Y[tID[0]]+525., Y[tID[0]]+1025., 10.))

    # Find the closest grids from source locations
    tid = np.argmin(np.abs(endl[0, 0] - txLoc[:, 0])+np.abs(endl[0, 1] - txLoc[:, 1]))

    # Interpolate the potentials ar suvey points from pre-computed grids
    P1_phi1 = scipy.interpolate.griddata(DcData[tid].srcList[0].rxList[0].locs[0][:, 0:2],
                                         DcData[tid].dobs,
                                         (grid_x1, grid_y1), method='linear')

    P2_phi1 = scipy.interpolate.griddata(DcData[tid].srcList[0].rxList[0].locs[0][:, 0:2],
                                         DcData[tid].dobs,
                                         (grid_x2, grid_y2), method='linear')

    tid = np.argmin(np.abs(endl[1, 0] - txLoc[:, 0])+np.abs(endl[1, 1] - txLoc[:, 1]))
    P1_phi2 = scipy.interpolate.griddata(DcData[tid].srcList[0].rxList[0].locs[0][:, 0:2],
                                         DcData[tid].dobs,
                                         (grid_x1, grid_y1), method='linear')

    P2_phi2 = scipy.interpolate.griddata(DcData[tid].srcList[0].rxList[0].locs[0][:, 0:2],
                                         DcData[tid].dobs,
                                         (grid_x2, grid_y2), method='linear')

    # Reshape the survey points into colums
    locRx1 = np.c_[mkvc(grid_x1), mkvc(grid_y1)]
    locRx2 = np.c_[mkvc(grid_x2), mkvc(grid_y2)]

    Pmid = (locRx1 + locRx2)/2.

    # Number of points
    nD = locRx1.shape[0]

    # Get distance between each electrodes
    rC1P1 = np.sqrt(np.sum((npm.repmat(endl[0, :], nD, 1) - locRx1)**2, axis=1))
    rC2P1 = np.sqrt(np.sum((npm.repmat(endl[1, :], nD, 1) - locRx1)**2, axis=1))
    rC1P2 = np.sqrt(np.sum((npm.repmat(endl[0, :], nD, 1) - locRx2)**2, axis=1))
    rC2P2 = np.sqrt(np.sum((npm.repmat(endl[1, :], nD, 1) - locRx2)**2, axis=1))

    # Compute potential difference
    volt = ((P1_phi1 - P1_phi2) - (P2_phi1 - P2_phi2))

    # Compute apparent resistivity
    rho = np.abs(mkvc(volt)) * np.pi * 2. / (1/rC1P1 - 1/rC2P1 - 1/rC1P2 + 1/rC2P2)

    # PLOT
    plt.figure()
    ax_prim = plt.subplot(1, 1, 1)
    plt.gca().set_aspect('equal', adjustable='box')

    plt.contourf(grid_x1[0, :],
                 grid_y1[:, 0],
                 rho.reshape(volt.shape, order='F'))
    plt.scatter(Pmid[:, 0], Pmid[:, 1], c=rho)
    plt.colorbar(ax=ax_prim)
    plt.title('Tiled')
    plt.show()

    # Plot grid
    fig, axs = plt.figure(figsize=(6, 10)), plt.subplot(1, 1, 1)
    axs.scatter(mkvc(X), mkvc(Y), c='r', s=20, edgecolor=None)
    axs.set_aspect('equal')
    for ii in range(len(X)):

        axs.text(X[ii], Y[ii], str(ii),
                 verticalalignment='bottom', horizontalalignment='center',
                 color='k', fontsize=10)

    axs.set_xlim(X.min(), X.max())
    axs.set_ylim(Y.min(), Y.max())

    d_out = np.c_[locRx1, locRx2, rho]
