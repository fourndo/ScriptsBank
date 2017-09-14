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
from scipy import interpolate
import numpy.matlib as npm
import time
import gc
import os
import glob
import re

work_dir = 'C:\\Egnyte\\Private\\dominiquef\\Projects\\4414_Minsim\\Modeling\\Forward\\'

out_dir = 'C:\\Egnyte\\Private\\dominiquef\\Projects\\4414_Minsim\\Modeling\\Forward\\Compilation\\'

local_dir = {'MAG':'FWR_MAG\\Ground\\',
             'Gz':'FWR_Gz\\Ground\\FWR_Gz\\',
             'TDEM':'FWR_TDEM\\OUT_DIR\\Grid65p5\\',
             'DC': '\\FWR_DC\\clean\\',
             'GTEM':'FWR_GTEM\\OUT_DIR\\',
             'Gxx':'FWR_Gz\\Ground\\FWR_Gxx\\',
             'Gxy':'FWR_Gz\\Ground\\FWR_Gxy\\',
             'Gxz':'FWR_Gz\\Ground\\FWR_Gxz\\',
             'Gyy':'FWR_Gz\\Ground\\FWR_Gyy\\',
             'Gyz':'FWR_Gz\\Ground\\FWR_Gyz\\',
             'Gzz':'FWR_Gz\\Ground\\FWR_Gzz\\'}
#work_dir = 'C:\\LC\\Private\\dominiquef\\Projects\\4414_Minsim\\Modeling\\VPmg\\FWR_Gz\\V2'

#work_dir = 'C:\\LC\\Private\\dominiquef\\Projects\\4414_Minsim\\Modeling\\VPmg\\
meshfile = 'C:\\Egnyte\\Private\\dominiquef\\Projects\\4414_Minsim\\Modeling\\Forward\\Mesh_20m.msh'
mesh = Mesh.TensorMesh.readUBC(meshfile)

ore = np.loadtxt(work_dir + 'Ore_Lenses.xyz')
dtype = 'DC'

outfile = 'Ground_Mag.dat'
#outfile = 'Ground_'+ dtype + '_Complex.dat'
label = 'TMI (nT)'#'mGal'#'$10^{-4}$ mGal/m'#'$\partial B_z/\partial z$'
vmin, vmax = -20, 20#-300, 1200#
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

    sfile.close()

    obs = np.asarray(obs)

    return np.asarray(obs)

#%% Body starts here
os.chdir(work_dir+local_dir[dtype])

if np.all([dtype != 'TDEM',dtype != 'GTEM', dtype != 'DC']):
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
#        r = ((gridx - np.median(tile[:, 0]))**2. +
#             (gridy - np.median(tile[:, 1]))**2.) + 1e-1
#
#        wght = 1./r**3.
        # Create dissolving tile
        nx = rowx.shape[0]
        ny = rowy.shape[0]
        wghtx, wghty = np.meshgrid(np.abs(np.asarray(range(nx))-np.ceil(nx)/2), np.abs(np.asarray(range(ny))-np.ceil(ny)/2))

        wght = np.max(np.c_[mkvc(wghtx),mkvc(wghty)], axis=1)

        # Normalize it and cosine it
        wght = 1-wght/wght.max() + 1e-6
        wght = -0.5 * np.cos(np.pi*wght) + 0.5
        wght = wght.reshape(gridx.shape, order='F')**2.
#        plt.figure()
#        plt.imshow(wght.reshape(gridx.shape))

        # Resample the grid
        data = griddata(tile[:, 0:2], tile[:, 3],
                        (gridx, gridy))

        # Check if level of current values matches the new one
        # otherwise sub-stract a DC
        temp_val = dataGrid[inDex]
        temp_wgt = wghtGrid[inDex]

        indTmp = temp_val != 0

        if indTmp.any():

            dc = np.mean(temp_val[indTmp]/temp_wgt[indTmp] - mkvc(data)[indTmp])
            data += dc

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
    plt.contourf(gCx,gCy,dataGrid.T,100, cmap = 'jet')
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    plt.colorbar(orientation='vertical', shrink=0.5, label=label)

    plt.contour(gCx,gCy,dataGrid.T,10,colors='k')
    plt.scatter(ore[:,0],ore[:,1],c='k', s = 1, linewidth=0)

    axs.set_aspect('equal')
    axs.set_title(outfile)
    plt.savefig(out_dir + outfile[:-4] + '.png')

    #%% Pplot true data

    # dataObs = np.loadtxt(work_dir+'\\Observed_Gridded50m.dat')
    # fig, axs = plt.figure(figsize=(5,9)), plt.subplot()

    # gCx = np.unique(dataObs[:,0])
    # gCy = np.unique(dataObs[:,1])

    # dataGrid = dataObs[:,-1].reshape((len(gCx), len(gCy)), order='F')
    # plt.contourf(gCx,gCy,dataGrid.T,100, cmap = 'jet', vmin=vmin, vmax=vmax)
    # plt.clim(vmin,vmax)
    # plt.xlim(xmin,xmax)
    # plt.ylim(ymin,ymax)
    # plt.colorbar(orientation='vertical', shrink=0.5, label=label)

    # plt.contour(gCx,gCy,dataGrid.T,10,colors='k')
    # plt.scatter(ore[:,0],ore[:,1],c='k', s = 1, linewidth=0)

    # axs.set_aspect('equal')
    # axs.set_title('Observed')
    # plt.savefig(out_dir + 'Observed.png')

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

elif dtype == 'GTEM':

    # Read the grid and extract survey lines
    # Re-organize data in to X,Y, t1, t2,... for each tile
    # Write to gocad pline
    os.chdir(work_dir + local_dir[dtype])

    files = glob.glob("*Tile*")
    nlines = 10
    nstn = 20
    times = []
    vrtxID  = 1
    dOut = []
    for file in files:
        print('Processing: ' + file)
        obj = read_h3d_pred(file)
        tID = int(re.search('\d+', file).group())

        if ~np.any(times):
            times = np.unique(obj[:, 3])
            ntimes = len(times)

        # Create survey line and interpolate data
        # ASSUME EW LINES FOR NOW
        locs = np.c_[obj[::ntimes, 0], obj[::ntimes, 1]]

        xlim = np.r_[locs[:, 0].min(), locs[:, 0].max()]
        ylim = np.r_[locs[:, 1].min(), locs[:, 1].max()]

        yline = np.linspace(ylim[0], ylim[1], nlines)
        xstn = np.linspace(xlim[0], xlim[1], nstn)

        [Y, X] = np.meshgrid(yline, xstn)
        X, Y = mkvc(X), mkvc(Y)
        nvrtx = X.shape[0]
        vIDs = np.linspace(vrtxID, vrtxID+nvrtx-1, nvrtx)
        loopID = np.kron(np.ones(nlines)*tID, np.ones(nstn))
        lIDs = tID*100 + np.kron(np.r_[range(nlines)]+1, np.ones(nstn))

        # Stack all the times, no topography for now
        d = np.c_[vIDs, X, Y, np.zeros(nvrtx), loopID, lIDs]

        for tt in range(ntimes):

            dtime = griddata(locs, -obj[tt::ntimes, -1], (X, Y), method='linear')
            d = np.c_[d, dtime]

        dOut += [d]
        vrtxID += nvrtx

#
#
#        for tt in range(ntimes):
#            dOut = np.c_[dOut,-obj[tt::ntimes,-1]]
#
#
#        fid = open(out_dir + '\\Tile_' + str(tID[-1]) + '_XYd.obs', 'w')
#        fid.write('X Y ')
#        for tt in range(ntimes): fid.write('t' + str(tt) +'_' + str(int(times[tt]*10**6.)) + 'us ')
#        fid.write('\n')
#        np.savetxt(fid, dOut, fmt='%e',delimiter=' ',newline='\n')
#        fid.close()

    # Write out a
    fid = open(work_dir + 'GTEM_FWRdata.pl', 'wb')

    fid.write(('GOCAD PLine 1\n').encode())
    fid.write(('HEADER {\n').encode())
    fid.write(('name:GTEM_FWRdata\n').encode())
    fid.write(('}\n').encode())
    fid.write(('PROPERTIES LoopID LineID dbdt\n').encode())
    fid.write(('PROP_LEGAL\n').encode())
    fid.write(('ESIZES 1 1 ' + str(ntimes)+ '\n').encode())
    dOut = np.vstack(dOut)
    lines = np.unique(dOut[:, 5])
    for ll in range(lines.shape[0]):

        fid.write(('ILINE\n').encode())

        dsub = dOut[dOut[:, 5] == lines[ll], :]
        for ii in range(dsub.shape[0]):
            fid.write(('PVRTX ').encode())
            np.savetxt(fid, dsub[ii,:].reshape((1,dsub.shape[1])), fmt=['%i', '%e', '%e', '%e', '%i', '%i']+['%e']*ntimes,delimiter=' ')

        for ii in range(dsub.shape[0]-1):
            fid.write(('SEG ').encode())
            np.savetxt(fid, np.c_[dsub[ii,0],dsub[ii,0]+1].reshape((1,2)), fmt='%i',delimiter=' ')

    fid.write(('END').encode())
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
        tID.append(int(re.search('\d+', file).group()))

    tID = np.asarray(tID)
    txLoc = np.asarray(txLoc)
    inds = np.argsort(tID)

#    fid = open(out_dir + '\\TxLocs.xyz', 'w')
#    np.savetxt(fid,np.c_[tID[inds],txLoc[inds]],fmt='%i',delimiter=' ',newline='\n')
#    fid.close()

    # Loop through tile pair and calculate potential + apparent resistivity
    ngridx = 9
    ngridy = 23
    nlines = 10
    nstn = 20

    count = -1
    for ii in range(ngridx):
        for jj in range(ngridy):

            count += 1
            
            if not np.any(tID == count):
                continue
            
            Aind = np.where(tID == count)[0][0]

            # loop through pairs NS and EW
            for pp in range(2):

                if np.all([pp == 0, jj+3 < ngridy]):
                    # Grab the north neighbour
                    
                    if not np.any(tID == (count + 3)):
                        continue
                    Bind = np.where(tID == count + 3)[0][0]
                    xExtent = np.r_[-1000, 1000]
                    yExtent = np.r_[475., 975]
                    offsets = np.r_[0, 25]

                elif np.all([pp == 1, ii+3 < ngridx]):

                    # Grab the north neighbour
                    if not np.any(tID == (count + 3*ngridy)):
                        continue
                    
                    Bind = np.where(tID == count + 3*ngridy)[0][0]
                    xExtent = np.r_[475., 975]
                    yExtent = np.r_[-1000, 1000]
                    offsets = np.r_[25, 0]

                else:
                    continue


                # # Grab the pole location for supplied grid ID
                # endl = np.c_[np.r_[X[tID]], np.r_[Y[tID]]]

                # # Create grid of survey electrodes
                # grid_x1, grid_y1 = np.meshgrid(np.linspace(X[tID[0]]-1000., X[tID[0]]+1000., 40.),
                #                                np.linspace(Y[tID[0]]+475., Y[tID[0]]+975., 10.))

                # grid_x2, grid_y2 = np.meshgrid(np.linspace(X[tID[0]]-1000., X[tID[0]]+1000., 40.),
                #                                np.linspace(Y[tID[0]]+525., Y[tID[0]]+1025., 10.))


                # Create survey line and interpolate data
                # ASSUME EW LINES FOR NOW
                loc_A = np.squeeze(txLoc[Aind, :])
                loc_B = np.squeeze(txLoc[Bind, :])

                xlim = np.r_[loc_A[0] + xExtent[0], loc_A[0] + xExtent[1]]
                ylim = np.r_[loc_A[1] + yExtent[0], loc_A[1] + yExtent[1]]

                if pp == 0:

                    ystn = np.linspace(ylim[0], ylim[1], nstn)
                    xline = np.linspace(xlim[0], xlim[1], nlines)

                    [X_1, Y_1] = np.meshgrid(xline, ystn)
                    [X_2, Y_2] = np.meshgrid(xline + offsets[0], ystn + offsets[1])

                else:
                    yline = np.linspace(ylim[0], ylim[1], nlines)
                    xstn = np.linspace(xlim[0], xlim[1], nstn)

                    [Y_1, X_1] = np.meshgrid(yline, xstn)
                    [Y_2, X_2] = np.meshgrid(yline + offsets[1], xstn + offsets[0])

#                X, Y = mkvc(X), mkvc(Y)

                lineID = np.kron(np.r_[range(nlines)]+1, np.ones(nstn))


                # Find the closest grids from source locations
                # tid = np.argmin(np.abs(endl[0, 0] - txLoc[:, 0])+np.abs(endl[0, 1] - txLoc[:, 1]))
                # Interpolate the potentials ar suvey points from pre-computed grids
                P1_phi1 = interpolate.griddata(DcData[Aind].srcList[0].rxList[0].locs[0][:, 0:2],
                                               DcData[Aind].dobs,
                                               (X_1, Y_1), method='linear')

                P2_phi1 = interpolate.griddata(DcData[Aind].srcList[0].rxList[0].locs[0][:, 0:2],
                                               DcData[Aind].dobs,
                                               (X_2, Y_2), method='linear')

                # tid = np.argmin(np.abs(endl[1, 0] - txLoc[:, 0])+np.abs(endl[1, 1] - txLoc[:, 1]))
                P1_phi2 = interpolate.griddata(DcData[Bind].srcList[0].rxList[0].locs[0][:, 0:2],
                                               DcData[Bind].dobs,
                                               (X_1, Y_1), method='linear')

                P2_phi2 = interpolate.griddata(DcData[Bind].srcList[0].rxList[0].locs[0][:, 0:2],
                                               DcData[Bind].dobs,
                                               (X_2, Y_2), method='linear')

                # Reshape the survey points into colums
                locRx1 = np.c_[mkvc(X_1), mkvc(Y_1)]
                locRx2 = np.c_[mkvc(X_2), mkvc(Y_2)]

                Pmid = (locRx1 + locRx2)/2.

                # Number of points
                nD = locRx1.shape[0]

                # Get distance between each electrodes
                rC1P1 = np.sqrt(np.sum((npm.repmat(loc_A[:2], nD, 1) - locRx1) ** 2, axis=1))
                rC2P1 = np.sqrt(np.sum((npm.repmat(loc_B[:2], nD, 1) - locRx1) ** 2, axis=1))
                rC1P2 = np.sqrt(np.sum((npm.repmat(loc_A[:2], nD, 1) - locRx2) ** 2, axis=1))
                rC2P2 = np.sqrt(np.sum((npm.repmat(loc_B[:2], nD, 1) - locRx2) ** 2, axis=1))

                # Compute potential difference
                volt = mkvc((P1_phi1 - P1_phi2) - (P2_phi1 - P2_phi2))

                # Compute apparent resistivity
                rho = np.abs(volt) * np.pi * 2. / (1/rC1P1 - 1/rC2P1 - 1/rC1P2 + 1/rC2P2)
                indx = np.isnan(rho) == False

                vIDs = np.linspace(0, np.sum(indx)-1, np.sum(indx))
                
                dOut = np.c_[vIDs, Pmid[indx, :], np.zeros(np.sum(indx)), lineID[indx], volt[indx], rho[indx]]

                # Write out a
                surveyName = 'Gradient_A_' + str(tID[Aind]) + '_B_' + str(tID[Bind])
                fid = open(out_dir + surveyName + '.pl', 'wb')

                fid.write(('GOCAD PLine 1\n').encode())
                fid.write(('HEADER {\n').encode())
                fid.write(('name:' + surveyName + '\n').encode())
                fid.write(('}\n').encode())
                fid.write(('PROPERTIES lineID volt AppRes\n').encode())
                fid.write(('PROP_LEGAL\n').encode())
                fid.write(('ESIZES 1 1 1\n').encode())
                dOut = np.vstack(dOut)
                lines = np.unique(dOut[:, 4])
                for ll in range(lines.shape[0]):

                    fid.write(('ILINE\n').encode())

                    dsub = dOut[dOut[:, 4] == lines[ll], :]
                    for rr in range(dsub.shape[0]):
                        fid.write(('PVRTX ').encode())
                        np.savetxt(fid, dsub[rr,:].reshape((1,dsub.shape[1])), fmt=['%i', '%e', '%e', '%e','%i', '%e', '%e'], delimiter=' ')

                    for rr in range(dsub.shape[0]-1):
                        fid.write(('SEG ').encode())
                        np.savetxt(fid, np.c_[dsub[rr,0],dsub[rr,0]+1].reshape((1,2)), fmt='%i',delimiter=' ')

                fid.write(('END').encode())
                fid.close()
