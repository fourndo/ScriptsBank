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
from SimPEG.EM.Static.Utils import drapeTopotoLoc, gen_DCIPsurvey, writeUBC_DCobs
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import time
import gc
import os
import glob

work_dir = 'C:\\LC\\Private\\dominiquef\\Projects\\4414_Minsim\\Modeling\\VPmg\\FWR_MAG'
gridFile = '..\\Mag_grid_50m.dat'
meshfile = '..\\Mesh_20m.msh'
mesh = Mesh.TensorMesh.readUBC(work_dir + '\\' + meshfile)
os.chdir(work_dir)


MagData = []
for file in glob.glob("*MAG_Tile*") :
    MagData.append(np.loadtxt(file))
        

# Grid data
gridLoc = np.loadtxt(work_dir + '\\' + gridFile)


 
# Get tile extents
Xlim = []
Ylim = []
for tile in MagData:
    
    Xlim.append(np.c_[np.min(tile[:,0]), np.mean(tile[:,0]), np.max(tile[:,0])]) 
    Ylim.append(np.c_[np.min(tile[:,1]), np.mean(tile[:,1]), np.max(tile[:,1])]) 


Xlim = np.squeeze(np.asarray(Xlim))
Ylim = np.squeeze(np.asarray(Ylim))

# Create data CC grid
offTile = 500
lenTile = 200
gCx = np.asarray(range(2))*lenTile + offTile + mesh.x0[0]
gCy = np.asarray(range(5))*lenTile + offTile + mesh.x0[1]

X,Y = np.meshgrid(gCx,gCy)
X, Y = mkvc(X), mkvc(Y)

fig, axs = plt.figure(figsize=(6,10)), plt.subplot(111)
axs.scatter(gridLoc[:,0],gridLoc[:,1],c=gridLoc[:,3],s=10,edgecolor=None)
axs.set_aspect('equal')

for ii in range(len(X)):
    axs.text(X[ii], Y[ii], str(ii),
        verticalalignment='center', horizontalalignment='center',
        color='r', fontsize=8)

#%% Progress function
def progress(iter, prog, final):
    """
    progress(iter,prog,final)

    Function measuring the progress of a process and print to screen the %.
    Useful to estimate the remaining runtime of a large problem.

    Created on Dec, 20th 2015

    @author: dominiquef
    """
    arg = np.floor(float(iter)/float(final)*100.)

    if arg > prog:

        print("Done " + str(arg) + " %")
        prog = arg

    return prog

#%% Create data grid
def grid_survey(spacing,width, x0 = (0,0,0), topo=None):
    """
        grid_survey(spacing,width)
        Generate concentric grid surveys centered at the origin
        :grid: List of grid spacing
        :width: Grid width
        
        return: rxLoc an n-by-3 receiver locations
        
    """
    rxLoc = []
    
    if len(spacing)!=len(width):
        raise 'Number of elements in spacing different than grid width'
    
    ii = -1
    for dx in spacing:

        ii += 1

        # Define survey grids centered at [0,0]
        
        # Number of grid points that fits in the width
        nC = int(width[ii]/dx)
        
        rxVec = -width[ii]/2. + dx/2. + np.asarray(range(nC))*dx
        
        rxGridx, rxGridy = np.meshgrid(rxVec,rxVec)
        
        rxGridx += x0[0]
        rxGridy += x0[1]
    
        if topo is not None:
            
            rxGridz = scipy.interpolate.griddata(topo[:,:2], topo[:,2], (rxGridx, rxGridy),
                                              method='linear') + x0[2]
                                             
        else:
            rxGridz = np.zeros_like(rxGridx) + x0[2]
            
        # Remove points if already inside inner grid
        if ii > 0:
            indx = np.logical_and(np.logical_and(rxGridy<rxLoc[:,1].max(), rxGridy>rxLoc[:,1].min()),np.logical_and(rxGridx>rxLoc[:,0].min(), rxGridx<rxLoc[:,0].max()))
        
            indx = indx.reshape(rxGridx.shape)
            rxGridx = rxGridx[indx==0]
            rxGridy = rxGridy[indx==0]
            rxGridz = rxGridz[indx==0]
            
            rxLoc = np.vstack([rxLoc, np.c_[mkvc(rxGridx),mkvc(rxGridy),mkvc(rxGridz)]])
            
        else:
            rxLoc = np.c_[mkvc(rxGridx),mkvc(rxGridy),mkvc(rxGridz)]
        

        
    return rxGridx, rxGridy, rxGridz
    
#%% Cycle through the grid and interpolate
ndat = len(X[:])
dOut = []
count = -1

print("Begin intepolation")
 
for ii in range(ndat):
    
    spacing=[20]
    width = [200]    
    rxGridx, rxGridy, rxGridz = grid_survey(spacing,width,x0 = (X[ii],Y[ii], 60.))
            
    # First figure out which tiles are touched by loc
    inDex = [False]*Xlim.shape[0]
    
    for jj in range(Xlim.shape[0]):
        
        indx = np.logical_and( np.any(Xlim[jj,0] < rxGridx), np.any(Xlim[jj,2] > rxGridx))
        indy = np.logical_and( np.any(Ylim[jj,0] < rxGridy), np.any(Ylim[jj,2] > rxGridy))
    
        if np.logical_and( indx, indy):
            inDex[jj] = True
            
    inDex =  np.where( inDex )
    inDex = np.array(inDex[0]).tolist()
    dat = 0.
    wgt = 0.
    
    for tID in inDex:
        
        # Compute distance to center of grid
        r = ((MagData[tID][:,0] - Xlim[tID,1])**2. + \
            (MagData[tID][:,1] - Ylim[tID,1])**2.) + 1e-8
        Idist = griddata( MagData[tID][:,0:2], 1./r**0.5 , (rxGridx,rxGridy), method= 'nearest' )
        
#        xid = np.abs(np.asarray(range(rxGridx.shape[0])) - rxGridx.shape[0]/2)
#        yid = np.abs(np.asarray(range(rxGridx.shape[1])) - rxGridx.shape[1]/2)
#        
#        Xid,Yid = np.meshgrid(-(xid-xid.max())+1,-(yid-yid.max())+1)
#        
#        Idist = np.minimum(Xid,Yid)
#        
#        # Normalize weights and cosine ramp off
#        Idist = Idist/float(np.max(Idist)+1e-6)
#        Idist = (-0.5*np.cos(-Idist*np.pi)+0.5) **0.5
        wgt += Idist
        
#        plt.figure(),plt.contourf(rxGridx,rxGridy,wgt)
        
        dat += griddata( MagData[tID][:,0:2], MagData[tID][:,3], (rxGridx,rxGridy), method= 'nearest' ) * Idist
#        plt.figure(),plt.contourf(rxGridx,rxGridy,dat)
    
    if len(inDex) > 0:
        dOut.append(np.c_[mkvc(rxGridx), mkvc(rxGridy), mkvc(dat / wgt)])
    
    count = progress(ii, count, ndat)
    
dOut = np.asarray(dOut).reshape((len(dOut)*dOut[0].shape[0],dOut[0].shape[1]))
#dOut = dOut.sort(axis=0,order=[0,1])

fid = open(work_dir + '\\MAGObsInterp.dat', 'w')
np.savetxt(fid, dOut, fmt='%e',delimiter=' ',newline='\n')
fid.close()

# Plot result grid
#xx = np.unique(dOut[:,0])
#yy = np.unique(dOut[:,1])
#
#X,Y = np.meshgrid(xx,yy)
#gridO = dOut[:,2].reshape((len(xx),len(yy)),order='F')
#
#plt.figure()
#plt.imshow(gridO)
