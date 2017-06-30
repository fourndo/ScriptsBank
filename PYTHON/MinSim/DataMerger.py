# -*- coding: utf-8 -*-
"""
Created on Mon Nov 07 08:03:02 2016

@author: dominiquef
"""

import numpy as np
import scipy.sparse as sp
from scipy.interpolate import griddata
from SimPEG import mkvc, Mesh
import SimPEG.PF as PF
import SimPEG.EM.Static.DC as DC
import SimPEG.EM.FDEM as FDEM
from SimPEG.EM.Static.Utils import readUBC_DC3Dobs
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import time
import gc
import os
import glob
import re

work_dir = 'C:\\Egnyte\\Private\\dominiquef\\Projects\\4414_Minsim\\Modeling\\Forward\\'

out_dir = 'C:\\Egnyte\\Private\\dominiquef\\Projects\\4414_Minsim\\Modeling\\Forward\\Compilation\\DC\\'

local_dir = {'MAG':'FWR_MAG\\Airborne\\',
             'Gz':'FWR_Gz\\Ground\\Complex\\',
             'TDEM':'FWR_TDEM\\OUT_DIR\\Grid65p5\\',
             'DC': '\\FWR_DC\\clean',
             'GDEM':'FWR_GTEM\\OUT_DIR\\',
             'Gxx':'FWR_Gz\\Airborne\\FWR_Gxx\\',
             'Gxy':'FWR_Gz\\Airborne\\FWR_Gxy\\',
             'Gxz':'FWR_Gz\\Airborne\\FWR_Gxz\\',
             'Gyy':'FWR_Gz\\Airborne\\FWR_Gyy\\',
             'Gyz':'FWR_Gz\\Airborne\\FWR_Gyz\\',
             'Gzz':'FWR_Gz\\Airborne\\FWR_Gzz\\'}
#work_dir = 'C:\\LC\\Private\\dominiquef\\Projects\\4414_Minsim\\Modeling\\VPmg\\FWR_Gz\\V2'

#work_dir = 'C:\\LC\\Private\\dominiquef\\Projects\\4414_Minsim\\Modeling\\VPmg\\
meshfile = 'C:\\Egnyte\\Private\\dominiquef\\Projects\\4414_Minsim\\Modeling\\Forward\\Mesh_20m.msh'
mesh = Mesh.TensorMesh.readUBC(meshfile)

ore = np.loadtxt(work_dir + 'Ore_Lenses.xyz')
dtype = 'MAG'

outfile = 'Airborne_MAG.dat'
#outfile = 'Ground_'+ dtype + '_Complex.dat'
label = 'mGal'#'TMI (nT)'#'$10^{-4}$ mGal/m'#'$\partial B_z/\partial z$'
vmin, vmax = -300, 1200#-20, 20#
xmin, xmax = 313520, 317450
ymin, ymax = 6065430, 6076085
#%% Progress function
def progress(iter, prog, final):
    """
    progress(iter,prog,final)

    Function measuring the progress of a process and print to screen the %.
    Useful to estimate the remaining runtime of a large problem.

    Created on Dec, 20th 2015

    @author: founrdo@gmail
    """
    arg = np.floor(float(iter)/float(final)*100.)

    if arg > prog:

        print("Done " + str(arg) + " %")
        prog = arg

    return prog

def read_h3d_pred(predFile):
    """
    read_h3d_pred(predFile)

    Function to read in h3d pred file

    Created on Dec, 12th 2016

    @author fourndo@gmail

    """
    sfile = open(predFile,'r')
    lines = sfile.readlines()

    obs = []
    ii = -1
    for line in lines[1:]:

        temp = np.array(line.split(), dtype=float)

        if np.any(temp):
            ii += 1
            obs.append(temp)

    fid.close()

    obs = np.asarray(obs)

    return np.asarray(obs)

#%% Body starts here
os.chdir(work_dir+local_dir[dtype])

if np.all([dtype != 'TDEM',dtype != 'GDEM', dtype != 'DC']):
    AllData = []
    for file in glob.glob("*Tile*"):
        AllData.append(np.loadtxt(file))



    # Create data CC grid
    offTile = 500
    dx = 30

    ndx = int((mesh.vectorNx[-1] - mesh.x0[0] - 2*offTile )/dx)
    ndy = int((mesh.vectorNy[-1] - mesh.x0[1] - 2*offTile )/dx)

    gCx = np.cumsum(np.ones(ndx)*dx) + offTile + mesh.x0[0] - dx
    gCy = np.cumsum(np.ones(ndy)*dx) + offTile + mesh.x0[1] - dx

    X, Y = np.meshgrid(gCx, gCy)
    X, Y = mkvc(X), mkvc(Y)

    # Assign tile index to grid and interpolate values and weights
    dataGrid = np.zeros_like(X)
    wghtGrid = np.zeros_like(X)

    count = -1
    ii = 0


    for tile in AllData:

        logicX = np.logical_and(np.min(tile[:, 0]) < X, np.max(tile[:, 0]) > X)
        logicY = np.logical_and(np.min(tile[:, 1]) < Y, np.max(tile[:, 1]) > Y)

        # Find points inside current tile
        inDex = np.logical_and(logicX, logicY)

        # Assume rectangular tiles and reshape the data locs
        rowx, rowy = np.unique(X[inDex]), np.unique(Y[inDex])

        # Reshape the grid
        gridx, gridy = np.meshgrid(rowx, rowy)

        # Radial weights
        r = ((gridx - np.median(tile[:, 0]))**2. +
             (gridy - np.median(tile[:, 1]))**2.) + 1e-1

        wght = 1./r**3.


        # Resample the grid
        data = griddata(tile[:, 0:2], tile[:, 3],
                        (gridx, gridy))

        # Append result to the full grid
        dataGrid[inDex] += mkvc((data*wght))

        # Append weights to full grid
        wghtGrid[inDex] += mkvc(wght)

        # Show progress
        count = progress(ii, count, len(AllData))

        ii += 1

    dataGrid /= wghtGrid

#    scale = (mkvc(Y)-np.min(Y))
#    dataGrid += (scale/np.max(scale))* 150.
#    dataGrid /= 0.91

    fid = open(out_dir + outfile, 'wb')
    np.savetxt(fid, np.c_[X,Y,dataGrid], fmt='%e',delimiter=' ',newline='\n')
    fid.close()

    #%% Plot result
    dataGrid = dataGrid.reshape(len(gCx),len(gCy))
    fig, axs = plt.figure(figsize=(5,9)), plt.subplot()
    plt.contourf(gCx,gCy,dataGrid.T,100, cmap = 'jet', vmin=vmin, vmax=vmax)
    plt.clim(vmin,vmax)
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    plt.colorbar(orientation='vertical', shrink=0.5, label=label)

    plt.contour(gCx,gCy,dataGrid.T,10,colors='k')
    plt.scatter(ore[:,0],ore[:,1],c='k', s = 1, linewidth=0)

    axs.set_aspect('equal')
    axs.set_title(outfile)
    plt.savefig(out_dir + outfile[:-4] + '.png')

    #%% Pplot true data

    dataObs = np.loadtxt(work_dir+'\\Observed_Gridded50m.dat')
    fig, axs = plt.figure(figsize=(5,9)), plt.subplot()

    gCx = np.unique(dataObs[:,0])
    gCy = np.unique(dataObs[:,1])

    dataGrid = dataObs[:,-1].reshape((len(gCx), len(gCy)), order='F')
    plt.contourf(gCx,gCy,dataGrid.T,100, cmap = 'jet', vmin=vmin, vmax=vmax)
    plt.clim(vmin,vmax)
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    plt.colorbar(orientation='vertical', shrink=0.5, label=label)

    plt.contour(gCx,gCy,dataGrid.T,10,colors='k')
    plt.scatter(ore[:,0],ore[:,1],c='k', s = 1, linewidth=0)

    axs.set_aspect('equal')
    axs.set_title('Observed')
    plt.savefig(out_dir + 'Observed.png')

#%% Deal with TDEM data
elif dtype == 'TDEM':
    AllData = []
    tID = []
    for file in glob.glob("*Tile*"):

        with open(file,'r') as sfile:
            lines = sfile.readlines()
            tID.append(int(re.search('\d+',file).group()))
            stn = []
            for line in lines:

                if line[0]!='%':

                    temp = map(float,line.split())

                    if len(temp)>0:
                        stn.append(temp)


        AllData.append(np.asarray(stn))

    AllData = np.asarray(AllData)

    fig = plt.figure(figsize=(10,10))

    offTile = 500
    dx = 250

    ndx = int((mesh.vectorNx[-1] - mesh.x0[0] - 2*offTile )/dx)
    ndy = int((mesh.vectorNy[-1] - mesh.x0[1] - 2*offTile )/dx)

    gCx = np.cumsum(np.ones(ndx)*dx) + offTile + mesh.x0[0] - dx
    gCy = np.cumsum(np.ones(ndy)*dx) + offTile + mesh.x0[1] - dx

    X, Y = np.meshgrid(gCx, gCy)

    for pp in range(AllData.shape[1]):

        axs = plt.subplot(4,5,pp+1)

        data = griddata(AllData[:,pp,0:2], AllData[:,pp,-1],
                        (X, Y))

        p = plt.contourf(gCx,gCy,np.log10(data),100, cmap='jet')
        plt.scatter(ore[:,0],ore[:,1],c=ore[:,2], s = 1, cmap='RdBu', linewidth=0)
        axs.set_xticklabels([])
        axs.set_yticklabels([])
#        p.clim([np.min(AllData[:,pp,-1]),np.max(AllData[:,pp,-1])*0.5])
        axs.set_aspect('equal')
        axs.set_title(str(AllData[0,pp,3]*1e+3) + ' ms')



    plt.tight_layout()
    axs.set_title(outfile)
    plt.savefig(out_dir + outfile[:-4] + '.png')

    fid = open(out_dir  + outfile, 'w')
    fid.write('X Y ')
    for tt in range(AllData.shape[1]): fid.write('t' + str(tt) +'_' + str(int(AllData[0,tt,3]*10**6.)) + 'us ')
    fid.write('\n')
    np.savetxt(fid, np.c_[AllData[:,0,0:2],AllData[:,:,-1]], fmt='%e',delimiter=' ',newline='\n')
    fid.close()

elif dtype == 'GDEM':

    # Re-organize data in to X,Y, t1, t2,... for each tile
    os.chdir(work_dir + local_dir[dtype])

    files = glob.glob("*Tile*")

    tID = []
    txLoc = []
    for file in files :


        obj = read_h3d_pred(file)

        times = np.unique(obj[:,3])
        ntimes = len(times)
        txLoc.append(np.c_[np.mean(obj[::ntimes,0]),np.mean(obj[::ntimes,1])])
        tID.append(int(re.search('\d+',file).group()))

        # Stack all the times
        dOut = np.c_[obj[::ntimes,0],obj[::ntimes,1]]

        for tt in range(ntimes):
            dOut = np.c_[dOut,-obj[tt::ntimes,-1]]


        fid = open(out_dir + '\\Tile_' + str(tID[-1]) + '_XYd.obs', 'w')
        fid.write('X Y ')
        for tt in range(ntimes): fid.write('t' + str(tt) +'_' + str(int(times[tt]*10**6.)) + 'us ')
        fid.write('\n')
        np.savetxt(fid, dOut, fmt='%e',delimiter=' ',newline='\n')
        fid.close()


    tID = np.asarray(tID)
    txLoc = np.squeeze(np.asarray(txLoc))
    inds = np.argsort(tID)

    fid = open(out_dir + '\\GridLocs.xyz', 'w')
    np.savetxt(fid,np.c_[tID[inds],txLoc[inds]],fmt='%i',delimiter=' ',newline='\n')
    fid.close()



elif dtype == 'DC':

    # Open all the potential grids in UBC format and convert to X,Y,phi
    os.chdir(work_dir + local_dir[dtype])
    DcData = []
    files = glob.glob("*Tile*")
    txLoc = []
    tID = []
    for file in files :


        obj = readUBC_DC3Dobs(file)

        DcData.append(obj['DCsurvey'])

        txLoc.append(obj['DCsurvey'].srcList[0].loc[0])
        tID.append(int(re.search('\d+',file).group()))

    tID = np.asarray(tID)
    txLoc = np.asarray(txLoc)
    inds = np.argsort(tID)

    fid = open(out_dir + '\\TxLocs.xyz', 'w')
    np.savetxt(fid,np.c_[tID[inds],txLoc[inds]],fmt='%i',delimiter=' ',newline='\n')
    fid.close()

    ii = -1
    for survey in DcData:

        ii += 1
        locXY = survey.srcList[0].rxList[0].locs[0][:,:2]
        phi = survey.dobs

        fid = open(out_dir + '\\Tile_' + str(tID[ii]) + '_XYphi.obs', 'w')
        np.savetxt(fid, np.c_[locXY,phi], fmt='%e',delimiter=' ',newline='\n')
        fid.close()

