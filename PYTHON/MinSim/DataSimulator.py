# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 13:31:29 2016

@author: dominiquef
"""

#from __future__ import unicode_literals

import scipy.io
from SimPEG import Maps, Mesh, mkvc, Utils
import numpy.matlib as npm
import numpy as np
import scipy.sparse as sp
import SimPEG.PF as PF
import SimPEG.EM.Static.DC as DC
import SimPEG.EM.FDEM as FDEM
import SimPEG.EM.TDEM as TDEM
from SimPEG.EM.Static.Utils import drapeTopotoLoc
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import time
import gc
from test_TDEM_forward_Analytic import halfSpaceProblemAnaDiff
from shutil import copyfile
import os
# from pymatsolver import BicgJacobiSolver as Solver

# %% Input parameters
workDir = 'C:\\Egnyte\\Private\\dominiquef\\Projects\\4414_Minsim\\Modeling\\Forward'
meshfile = 'Mesh_20m.msh'
locfile = 'Mag_grid_50m.dat'
modelfile = {'MAG': 'Mesh_20m_Susc\\Susc_Constrained.sus',
             'DC': 'Mesh_20m_Cond\\Cond_Median_EDT.con',
             'Gz': 'Mesh_20m_Dens\\Dens_Constrained.den'}

outDirDict = {'MAG': "\\FWR_MAG\\MAG_Tile",
              'DC': "\\FWR_DC\\DC_Tile",
              'Gz': "\\FWR_Gz\\Gz_Tile",
              'ATEM': "\\FWR_TDEM\\INP_DIR\\",
              'GTEM': "\\FWR_GTEM\\INP_DIR\\"}

nullcell = 'nullcell.dat'

obsfile = 'Mesh_20m_Susc\\Mag_grid_50m.dat'
topofile = 'ROT_DEM_30m.topo'

dtype = ['GTEM']

outDir = outDirDict[dtype[0]]

gc.collect()

# %% Load files
mesh_global = Mesh.TensorMesh.readUBC(workDir + '\\' + meshfile)

model = {'MAG': mesh_global.readModelUBC(workDir + '\\' + modelfile['MAG']),
         'Gz': mesh_global.readModelUBC(workDir + '\\' + modelfile['Gz']),
         'DC': mesh_global.readModelUBC(workDir + '\\' + modelfile['DC'])}

topo = np.genfromtxt(workDir + '\\' + topofile,
                     skip_header=1)

actvmod = mesh_global.readModelUBC(workDir + '\\' + nullcell)
actv = np.ones(mesh_global.nC, dtype=bool)
actv[actvmod == 0] = False

srvy = PF.MagneticsDriver.MagneticsDriver_Inv()
srvy.basePath = workDir
# data = srvy.readMagneticsObservations('\\'+obsfile)

# locs = data.srcField.rxList[0].locs

offTile = 500
lenTile = 500

dx = [20., 40.]
#
nCx = int(np.floor((500/dx[0]))-1)

core = {'DC': np.r_[np.ones(110)*dx[0]],
        'ATEM': np.r_[np.ones(5)*dx[1], [33, 26], np.ones(10)*dx[0],
                      [26, 33], np.ones(5)*dx[1]],
        'GTEM': np.r_[dx[1], [33, 26], np.ones(nCx)*dx[0], np.ones(nCx)*dx[0],
                      [26, 33], dx[1]],
        'Gz': np.r_[np.ones(15)*dx[1], [20, 15], np.ones(nCx)*dx[0],
                    [15, 20], np.ones(15)*dx[1]],
        'MAG': np.r_[np.ones(15)*dx[1], [20, 15], np.ones(nCx)*dx[0],
                     [15, 20], np.ones(15)*dx[1]]}

ndv = {'MAG': 0, 'Gz': 0, 'DC': 1e-8, 'FDEM': 1e-8, 'ATEM': 1e-8, 'GTEM': 1e-8}

expf = 1.3
npadxy = {'DC': 10, 'ATEM': 7, 'GTEM': 7, 'Gz': 7, 'MAG': 7}
npadz = {'DC': 10, 'ATEM': 10, 'GTEM': 10, 'Gz': 10, 'MAG': 10}

gCx = np.asarray(range(9))*lenTile + offTile + mesh_global.x0[0]
gCy = np.asarray(range(23))*lenTile + offTile + mesh_global.x0[1]

X, Y = np.meshgrid(gCx, gCy)

indT = np.zeros(X.shape, dtype=bool)
indT[:, :] = True
indT[:, :] = True

X, Y, indT = mkvc(X), mkvc(Y), mkvc(indT)

zvec = np.r_[[33], np.ones(25)*30, [23, 18, 15, 12],
             np.ones(15)*10, [12, 15, 18, 23],
             30*expf**(np.asarray(range(npadz[dtype[0]])))]


# Generate grid survey for given parameters
def grid_survey(spacing, width, x0=(0, 0, 0), topo=None):
    """
        grid_survey(spacing,width)
        Generate concentric grid surveys centered at the origin
        :grid: List of grid spacing
        :width: Grid width

        return: rxLoc an n-by-3 receiver locations

    """
    rxLoc = []

    if len(spacing) != len(width):
        raise 'Number of elements in spacing different than grid width'

    ii = -1
    for dx in spacing:

        ii += 1

        # Define survey grids centered at [0,0]

        # Number of grid points that fits in the width
        nC = int(width[ii]/dx)

        rxVec = -width[ii]/2. + dx/2. + np.asarray(range(nC))*dx

        rxGridx, rxGridy = np.meshgrid(rxVec, rxVec)

        rxGridx += x0[0]
        rxGridy += x0[1]

        if topo is not None:

            rxGridz = scipy.interpolate.griddata(topo[:, :2], topo[:, 2],
                                                 (rxGridx, rxGridy),
                                                 method='linear') + x0[2]

        else:
            rxGridz = np.zeros_like(rxGridx) + x0[2]

        # Remove points if already inside inner grid
        if ii > 0:
            indx = np.logical_and(np.logical_and(rxGridy < rxLoc[:, 1].max(),
                                                 rxGridy > rxLoc[:, 1].min()),
                                  np.logical_and(rxGridx > rxLoc[:, 0].min(),
                                  rxGridx < rxLoc[:, 0].max()))

            indx = indx.reshape(rxGridx.shape)
            rxGridx = rxGridx[indx == 0]
            rxGridy = rxGridy[indx == 0]
            rxGridz = rxGridz[indx == 0]

            rxLoc = np.vstack([rxLoc, np.c_[mkvc(rxGridx), mkvc(rxGridy),
                                            mkvc(rxGridz)]])

        else:
            rxLoc = np.c_[mkvc(rxGridx), mkvc(rxGridy), mkvc(rxGridz)]

    return rxLoc

# Loop through the tiles and create meshs
meshes = []
# fig, axs = plt.figure(figsize=(6, 10)), plt.subplot(111)
# axs.scatter(locs[:, 0], locs[:, 1], c=data.dobs, s=10, edgecolor=None)

# axs.scatter(mkvc(X), mkvc(Y), c='r', s=20, edgecolor=None)
# axs.set_aspect('equal')
# for ii in range(len(X)):
#     if indT[ii]:
#         axs.text(X[ii], Y[ii], str(ii),
#                  verticalalignment='center', horizontalalignment='center',
#                  color='k', fontsize=10)
#    axs.add_patch(
#    patches.Rectangle(
#        (mtemp.x0[0],mtemp.x0[1]),   # (x,y)
#        np.sum(mtemp.hx),          # width
#        np.sum(mtemp.hy), fill = False   ))
# axs.set_xlim(X.min(), X.max())
# axs.set_ylim(Y.min(), Y.max())


# %% Loop through the meshes and forward model
DcData = []
MagData = []
GzData = []
FdemData = []

os.chdir(workDir)

tiles = range(len(X[:]))

# Initialize ckdTree and store after the first time it is called
ckdTree = None

TxLoop = []
for tID in tiles:

    # Create active cells for sub mesh
    if not indT[tID]:
        continue

    # INTERPOLATE THE MODEL TO LOCAL MESH
    for dID in dtype:

        start_time = time.time()

        padxy = np.r_[dx[1]*expf**(np.asarray(range(npadxy[dID]))+1)]
        padz = np.r_[dx[1]*expf**(np.asarray(range(npadz[dID]))+1)]

        hx = np.r_[padxy[::-1], core[dID], padxy]
        hy = np.r_[padxy[::-1], core[dID], padxy]
        hz = np.r_[padz, zvec]

        mesh = Mesh.TensorMesh([hx, hy, hz], 'CC0')

        mesh._x0 = (mesh.x0[0] + X[tID], mesh.x0[1]+Y[tID],
                    mesh.x0[2]-(np.sum(hz[:(npadz[dID]+20)])))

        # Extract model from global to local mesh
        actv2 = Utils.surface2ind_topo(mesh, topo, 'N')

        if ckdTree is None:
            print("Creating ckTree... this might take a while")
            ckdTree = Maps.Mesh2MeshTopo([mesh_global, mesh],
                                         [actv, actv2], nIterpPts=12)

        P = Maps.Mesh2MeshTopo([mesh_global, mesh],
                               [actv, actv2], nIterpPts=12, tree=ckdTree.tree)

        if dID == 'FDEM':
            model_Tile = P*(model['DC'][actv])

        elif dID == 'ATEM':
            model_Tile = P*(model['DC'][actv])

        elif dID == 'GTEM':
            model_Tile = P*(model['DC'][actv])
        else:
            model_Tile = P*(model[dID][actv])

        m = np.ones(mesh.nC)*ndv[dID]
        m[actv2] = model_Tile
        # mtemp = m.reshape(mesh.vnC, order='F')
        # mtemp[:npadxy[dID],:,:] = ndv[dID]
        # mtemp[-npadxy[dID]:,:,:] = ndv[dID]
        # mtemp[:,:npadxy[dID],:] = ndv[dID]
        # mtemp[:,-npadxy[dID]:,:] = ndv[dID]
        # mtemp[:,:,:npadz[dID]] = ndv[dID]

#        m = mkvc(mtemp)
        print("Model Interpolated  %s seconds " % (time.time()-start_time))

        if dID == 'MAG':

            # OPTION TO WRITE OUT
            # Mesh.TensorMesh.writeUBC(mesh, work_dir + '\\MeshTile1.msh')
            # Mesh.TensorMesh.writeModelUBC(mesh,
            #                               workDir + '\\MeshTile1.mod', m)

            # Create grid of points
            spacing = (20, 40)
            width = (500, 1500)
            rxLoc = grid_survey(spacing, width,
                                x0=(X[tID], Y[tID], 60.), topo=topo)

            # PDE Formulation
            m_map = PF.BaseMag.BaseMagMap(mesh)
            prob = PF.Magnetics.Problem3D_DiffSecondary(mesh, muMap=m_map)

            survey = PF.BaseMag.BaseMagSurvey()
            survey.setBackgroundField(79, 31, 60000)
            survey.rxLoc = rxLoc
            start_time = time.time()
            prob.pair(survey)
            u = prob.fields(m)
            data = survey.projectFields(u)

            print("Solve PDE --- %s seconds ---" % (time.time() - start_time))

#            data = np.c_[rxLoc,data]
#            MagData.append(data)
            # Plot out

            fid = open(workDir + outDir + str(tID) + '.obs', 'w')
            np.savetxt(fid, np.c_[rxLoc, data],
                       fmt='%e', delimiter=' ', newline='\n')
            fid.close()

            del prob, u, survey, m_map, data

        elif dID == 'Gz':

            # Create grid of points
            spacing = (20, 40)
            width = (500, 1500)
            locXyz = grid_survey(spacing, width,
                                 x0=(X[tID], Y[tID], 60.), topo=topo)

            # Mesh.TensorMesh.writeUBC(mesh, workDir + outDir + 'MeshTile.msh')
            # Mesh.TensorMesh.writeModelUBC(mesh,
            #                               workDir + outDir + 'MeshTile.den',
            #                               m)

            m_map = PF.BaseGrav.BaseGravMap(mesh)
            prob = PF.Gravity.Problem3D_Diff(mesh, rhoMap=m_map)

            rxLoc = PF.BaseGrav.RxObs(locXyz)
            srcField = PF.BaseGrav.SrcField([rxLoc])
            survey = PF.BaseGrav.LinearSurvey(srcField)

            prob.pair(survey)
            u = prob.fields(m)

            gg = survey.projectFields(u)
            data = gg['gz']

            print("Solve PDE --- %s seconds ---" % (time.time() - start_time))

            fid = open(workDir + outDir + str(tID) + '.obs', 'w')
            np.savetxt(fid, np.c_[locXyz, data],
                       fmt='%e', delimiter=' ', newline='\n')
            fid.close()

            del prob, u, survey, m_map, data

        elif dID == 'DC':

            Mesh.TensorMesh.writeUBC(mesh, workDir + outDir + 'MeshTile.msh')
            Mesh.TensorMesh.writeModelUBC(mesh, workDir + outDir + 'MeshTile.mod', m)
            stype = 'gradient'

            # Get the cell in center of mesh
            tx = mkvc(drapeTopotoLoc(mesh, topo,
                                     np.c_[X[tID], Y[tID]], airind=actv2))

            # Create grid of points
            spacing = (20)
            width = (2200)
            rxLoc = grid_survey([spacing], [width],
                                x0=(X[tID]-5., Y[tID]-5., 0.))

            rxLocDrapped = drapeTopotoLoc(mesh, topo,
                                          rxLoc[:, :2], airind=actv2)

            # Write locations out
            fid = open(workDir + outDir + 'DCObslocs.dat', 'w')
            np.savetxt(fid, np.r_[tx, tx,
                                  rxLocDrapped.shape[0]].reshape((1, 7)),
                       fmt='%e', delimiter=' ', newline='\n')

            np.savetxt(fid, np.c_[rxLocDrapped, rxLocDrapped],
                       fmt='%e', delimiter=' ', newline='\n')
            fid.close()

            # Run forward model
            os.system('dcipf3d dcipf3d.inp')

            copyfile(workDir + outDir + 'dc3d.dat', workDir + outDir + 'Tile_' + str(tID) + '.obs')

            print("Solve Tile %i --- %s seconds ---" % (tID, time.time() - start_time))

        elif dID == 'FDEM':

            # OPTION TO WRITE OUT
            Mesh.TensorMesh.writeUBC(mesh, workDir + '\\MeshTile1.msh')
            Mesh.TensorMesh.writeModelUBC(mesh, workDir + '\\MeshTile1.mod', m)

            start_time = time.time()

            rxLoc = grid_survey([1], [1], x0=(X[tID], Y[tID], 60.), topo=topo)

            bzi = FDEM.Rx.Point_bSecondary(rxLoc, 'z', 'real')
            bzr = FDEM.Rx.Point_bSecondary(rxLoc, 'z', 'imag')

            freqs = np.asarray([400.])
            srcLoc = rxLoc
            srcLoc[0][2] += 20.
            srcList = [
                FDEM.Src.MagDipole([bzr, bzi], freq, srcLoc, orientation='Z')
                for freq in freqs
            ]

            mapping = Maps.IdentityMap(mesh)
            surveyFD = FDEM.Survey(srcList)
            prbFD = FDEM.Problem3D_b(mesh, sigmaMap=mapping, Solver=Solver)
            prbFD.pair(surveyFD)
            std = 0.0
            surveyFD.makeSyntheticData(m, std, force=True)

            print("Solve Tile %i --- %s seconds ---" % (tID,
                                                        time.time()-start_time))

            np.savetxt(fid, np.r_[srcLoc[0], surveyFD.dobs].reshape((1, 13)),
                       fmt='%e', delimiter=' ', newline='\n')

        elif dID == 'ATEM':

            # OPTION TO WRITE OUT
            Mesh.TensorMesh.writeUBC(mesh, workDir + outDir + 'Tile_' + str(tID) + '.msh')
            Mesh.TensorMesh.writeModelUBC(mesh, workDir + outDir + 'Tile_' + str(tID) + '.con', m)

            # halfSpaceProblemAnaDiff
            start_time = time.time()

            rxLoc = grid_survey([1], [1], x0=(X[tID], Y[tID], 30.), topo=topo)

            srcLoc = rxLoc

            times = np.loadtxt(workDir + '\\VTEM_times.dat')

            print("Solve Tile %i --- %s seconds ---" % (tID, time.time()-start_time))

            fid = open(workDir + outDir + 'Trx_Tile_' + str(tID) + '.loc', 'w')
            fid.write('N_TRX 1\n\n')
            fid.write('TRX_LOOP\n')
            np.savetxt(fid, np.r_[srcLoc[0], 13., 0, 0].reshape((1,6)), fmt='%e',delimiter=' ',newline='\n\n')
            fid.write('N_RECV 1\n')

            np.savetxt(fid, srcLoc[0].reshape((1,3)), fmt='%e',delimiter=' ',newline='\n\n')

            fid.close()

        elif dID == 'GTEM':

            # OPTION TO WRITE OUT
            Mesh.TensorMesh.writeUBC(mesh,
                                     workDir + outDir + 'Tile_' + str(tID) + '.msh')
            Mesh.TensorMesh.writeModelUBC(mesh,
                                          workDir + outDir + 'Tile_' + str(tID) + '.con', m)

#            halfSpaceProblemAnaDiff
            start_time = time.time()

            rxLoc = grid_survey([10], [1000],
                                x0=(X[tID], Y[tID], 1.), topo=topo)


            print("Solve Tile %i in %s seconds" % (tID, time.time()-start_time))

#            dtemp =
            fid = open(workDir + outDir + 'Trx_Tile_' + str(tID) + '.loc', 'wb')
#            fid.write('IGNORE -9.9999 \n\n')
            line = 'N_TRX 1\n\n' + 'TRX_ORIG\n' + '5\n'
            
            fid.write(line.encode('utf-8'))


            # Define points for square loop centered on grid
            xOrig = X[tID] + lenTile * np.r_[-0.5, 0.5]
            yOrig = Y[tID] + lenTile * np.r_[-0.5, 0.5]

            xOrig, yOrig = np.meshgrid(xOrig, yOrig)

            # Flip the last two nodes for clockwise ordering
            yOrig = np.c_[yOrig[:, 0], yOrig[::-1, 1]]

            zOrig = scipy.interpolate.griddata(topo[:, :2], topo[:, 2],
                                               (xOrig, yOrig),
                                               method='linear') + 1.

            xOrig, yOrig, zOrig = mkvc(xOrig), mkvc(yOrig), mkvc(zOrig)

            TxLoop.append(np.c_[mkvc(xOrig), mkvc(yOrig), mkvc(zOrig)])

            for ii in range(4):
                np.savetxt(fid, np.c_[xOrig[ii], yOrig[ii], zOrig[ii]],
                           fmt='%e', delimiter=' ', newline='\n')
            # repeat first node to close the loop
            np.savetxt(fid, np.c_[xOrig[0], yOrig[0], zOrig[0]],
                       fmt='%e', delimiter=' ', newline='\n')

            line = 'N_RECV ' + str(rxLoc.shape[0]) + '\n'
            fid.write(line.encode('utf-8'))
            
#            fid.write('N_TIME ' + str(len(times)) + '\n\n')

            np.savetxt(fid, rxLoc, fmt='%e', delimiter=' ', newline='\n')

            fid.close()

    gc.collect()



if dtype == 'GTEM':

    # %% Write out GOcad polyline
    fid = open(workDir + '\\FWR_GTEM\\GroundLoops.pl', 'w')
    fid.write('GOCAD PLine 1\n')
    fid.write('HEADER {name:GTEM_Loops}\n')
    fid.write('PROPERTIES ID\n')
    fid.write('PROPERTY_CLASSES id\n')
    fid.write('PROPERTY_CLASS_HEADER id {\n')
    fid.write('name:ID\n')
    fid.write('}\n')

    count_vrx = 0
    count_tx = -1
    for tx in TxLoop:

        fid.write('ILINE\n')

        count_tx += 1
        for ii in range(4):
            count_vrx += 1
            fid.write('PVRTX ')
            np.savetxt(fid, np.r_[count_vrx, tx[ii,:], count_tx].reshape((1,5)), fmt='%i',delimiter=' ')

        count_vrx -= 3

        for ii in range(3):
            fid.write('SEG ')
            np.savetxt(fid, np.r_[count_vrx+ii,count_vrx+ii+1].reshape((1,2)), fmt='%i',delimiter=' ')

        fid.write('SEG ')
        np.savetxt(fid, np.r_[count_vrx+ii+1,count_vrx].reshape((1,2)), fmt='%i',delimiter=' ')

        count_vrx += 4

    fid.write('END')
    fid.close()
