# -*- coding: utf-8 -*-
"""
Created on Sun Nov 13 08:31:12 2016

@author: dominiquef
"""

import scipy.io
from SimPEG import Maps, Mesh, mkvc, Utils
import numpy.matlib as npm
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os
import glob
import re
from matplotlib import animation
from JSAnimation import HTMLWriter


work_dir = 'C:\\LC\\Private\\dominiquef\\Projects\\4414_Minsim\\Modeling\\Forward\\Compilation'

ore = np.loadtxt(work_dir + '\\..\\Ore_Lenses.xyz')
inp_dir = {'DC':'\\FWR_DC',
           'TDEM':'\\TEM',
           'GTEM':'\\FWR_GTEM\\OUT_DIR'}

meshfile = '..\\Mesh_20m.msh'
mesh_global = Mesh.TensorMesh.readUBC(work_dir + '\\' + meshfile)

modelfile = {'MAG':'Mesh_20m_Susc\\Randn_Std_model.sus','DC':'Mesh_20m_Cond\\Randn_Std_model.con','Gz':'Mesh_20m_Dens\\Density_Randn_Std.den'}
#nullcell = 'Mesh_20m_Susc\\nullcell.dat'
#
#topo = np.genfromtxt(work_dir + '\\' + topofile,
#                                     skip_header=1)
#
#actvmod = mesh_global.readModelUBC(work_dir + '\\' + nullcell)
#actv = np.ones(mesh_global.nC,dtype=bool)
#actv[actvmod==0] = False

obsfile = 'Mesh_20m_Susc\\Mag_grid_50m.dat'
topofile = 'Mesh_20m_Susc\\ROT_DEM_30m.topo'

dType = 'TDEM'

ndv = {'MAG':0, 'Gz':0, 'DC':1e-8, 'FDEM':1e-8, 'TDEM':1e-8, 'GTEM': 1e-8}

npadxy = {'DC':10,'TDEM':7, }
npadz = {'DC':10, 'TDEM':10, }
expf = 1.3

offTile = 500
lenTile = 250

# Original grid from "DataSimulator"
#gCx = np.asarray(range(9))*lenTile + offTile + mesh_global.x0[0]
#gCy = np.asarray(range(23))*lenTile + offTile + mesh_global.x0[1]

gCx = np.asarray(range(18))*lenTile + offTile + mesh_global.x0[0]
gCy = np.asarray(range(46))*lenTile + offTile + mesh_global.x0[1]

X,Y = np.meshgrid(gCx,gCy)

indT = np.ones(X.shape, dtype=bool)
#indT[::2,::2] = True
#indT[::2,::2] = True

X, Y, indT = mkvc(X), mkvc(Y), mkvc(indT)

dx = [20., 40.]



core = {'DC':np.r_[np.ones(80)*dx[0]],
        'TDEM':np.r_[np.ones(5)*dx[1],[33,26],np.ones(10)*dx[0],[26,33],np.ones(5)*dx[1]]}

if dType == 'DC':

    # Fix bug in UBC export
    os.chdir(work_dir + inp_dir[dType])
    source_dir=os.getcwd()
    target_dir="clean"
    source_files = [fname for fname in glob.glob(os.path.join(source_dir,"*.obs"))]

    # check if target directory exists... if not, create it.
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    for source_file in source_files:
       target_file = os.path.join(target_dir,os.path.basename(source_file))
       with open(source_file,'r') as sfile:
          with open(target_file,'w') as tfile:
             lines = sfile.readlines()
             # do the replacement in the second line.
             # (remember that arrays are zero indexed)
             lines[2]=re.sub("E\+0210000",'E+02 10000',lines[2])
             tfile.writelines(lines)

    print "DONE"


    os.chdir(work_dir + inp_dir[dType] + '\\' + target_dir)
    DcData = []
    files = glob.glob("*Tile*")
    txLoc = []
    txId = []
    for file in files :


        obj = readUBC_DC3Dobs(file)

        DcData.append(obj['DCsurvey'])

        txLoc.append(obj['DCsurvey'].srcList[0].loc[0])

        txId.append(int(re.search('\d+',file).group()))

    txLoc = np.asarray(txLoc)

    # Sort tiles in increasing name and re-order list
    indx = np.argsort(np.asarray(txId))
    txLoc = txLoc[indx]
    DcData = [ DcData[i] for i in indx.tolist()]

    fid = open(work_dir + '\\FWR_DC\\Gradient_Tx_ALL.dat', 'w')
    np.savetxt(fid, txLoc, fmt='%e',delimiter=' ',newline='\n')
    fid.close()

    arrays = []
#    fid = open(work_dir + '\\FWR_DC\\Gradient_50m_ALL.dat', 'w')
    fig = plt.figure(figsize=(6,8))
    ax1 = plt.subplot()

    endl = []
    for ii in range(txLoc.shape[0]-3):

        endp = np.c_[txLoc[[ii,ii+3],:]]

        if endp[0,0]!=endp[1,0]:

            continue

        endl.append(endp)

    def animate(ii):

        removeFrame()

        grid_x1, grid_y1 = np.meshgrid(np.linspace(endl[ii][0,0]-250.,endl[ii][0,0]+250.,20.),np.linspace(endl[ii][0,1]+475.,endl[ii][0,1]+975.,10.))
        grid_x2, grid_y2 = np.meshgrid(np.linspace(endl[ii][0,0]-250.,endl[ii][0,0]+250.,20.),np.linspace(endl[ii][0,1]+525.,endl[ii][0,1]+1025.,10.))

        # Extract data from pre-computed grids
        # Grid points
    #    grid_x1, grid_y1 = np.meshgrid(np.unique(FWR_survey.srcList[0].rxList[0].locs[0][:,0]),
    #                                   np.unique(FWR_survey.srcList[0].rxList[0].locs[0][:,1]))
    #
    #    grid_x2, grid_y2 = np.meshgrid(np.unique(FWR_survey.srcList[0].rxList[0].locs[1][:,0]),
    #                                   np.unique(FWR_survey.srcList[0].rxList[0].locs[1][:,1]))

        tid = np.argmin(np.abs(endl[ii][0,0] - txLoc[:,0])+np.abs(endl[ii][0,1] - txLoc[:,1]))

        P1_phi1 = scipy.interpolate.griddata(DcData[tid].srcList[0].rxList[0].locs[0][:,0:2],
                                             DcData[tid].dobs,
                                             (grid_x1, grid_y1), method='linear')

        P2_phi1 = scipy.interpolate.griddata(DcData[tid].srcList[0].rxList[0].locs[0][:,0:2],
                                             DcData[tid].dobs,
                                             (grid_x2, grid_y2), method='linear')

        tid = np.argmin(np.abs(endl[ii][1,0] - txLoc[:,0])+np.abs(endl[ii][1,1] - txLoc[:,1]))
        P1_phi2 = scipy.interpolate.griddata(DcData[tid].srcList[0].rxList[0].locs[0][:,0:2],
                                             DcData[tid].dobs,
                                             (grid_x1, grid_y1), method='linear')

        P2_phi2 = scipy.interpolate.griddata(DcData[tid].srcList[0].rxList[0].locs[0][:,0:2],
                                             DcData[tid].dobs,
                                             (grid_x2, grid_y2), method='linear')


        locRx1 = np.c_[mkvc(grid_x1),mkvc(grid_y1)]
        locRx2 = np.c_[mkvc(grid_x2),mkvc(grid_y2)]

        Pmid = (locRx1 + locRx2)/2.

        nD = locRx1.shape[0]
        rC1P1 = np.sqrt( np.sum( (npm.repmat(endl[ii][0,:2],nD, 1) - locRx1)**2, axis=1 ))
        rC2P1 = np.sqrt( np.sum( (npm.repmat(endl[ii][1,:2],nD, 1) - locRx1)**2, axis=1 ))
        rC1P2 = np.sqrt( np.sum( (npm.repmat(endl[ii][0,:2],nD, 1) - locRx2)**2, axis=1 ))
        rC2P2 = np.sqrt( np.sum( (npm.repmat(endl[ii][1,:2],nD, 1) - locRx2)**2, axis=1 ))

        #rC1C2 = np.sqrt( np.sum( (npm.repmat(Tx[0][0:2]-Tx[0][3:5],Rx[0].shape[0], 1) )**2, axis=1 ))
        #rP1P2 = np.sqrt( np.sum( (Rx[0][:,0:2] - Rx[1][:,0:2])**2, axis=1 ))
        volt = ((P1_phi1 - P1_phi2) - (P2_phi1 - P2_phi2))


        rho = np.abs(mkvc(volt)) *np.pi *2./ ( 1/rC1P1 - 1/rC2P1 - 1/rC1P2 + 1/rC2P2 )#*((rC1P1)**2 / rP1P2)#

        global ax1, fig
        ax1 = plt.subplot()
        plt.gca().set_aspect('equal', adjustable='box')

        plt.contourf(grid_x1[0,:],
                     grid_y1[:,0],
                    rho.reshape(volt.shape,order='F'))
        plt.scatter(Pmid[:,0],Pmid[:,1],c=rho)
        plt.colorbar(ax = ax1)
        plt.title('Tiled')
        plt.show()

#        indx = np.isnan(rho) == False

#        np.savetxt(fid, np.c_[Pmid[indx,:],rho[indx]], fmt='%e',delimiter=' ',newline='\n')

#    fid.close()
    def removeFrame():

        global ax1, fig
        fig.delaxes(ax1)
        plt.draw()

    anim = animation.FuncAnimation(fig, animate,
                                   frames=len(endl) , interval=1000, repeat = False)
    #/home/dominiquef/3796_AGIC_Research/DCIP3D/MtISa
    anim.save('Data_slice.html', writer=HTMLWriter(embed_frames=True,fps=1))




if dType == 'TDEM':

    os.chdir(work_dir + inp_dir[dType])

    AllData = np.loadtxt('Airborne_TDEM.dat', dtype = float, skiprows = 1)

    xy = AllData[:,:2]
    data = AllData[:,2:]
    times = np.loadtxt('times.txt', dtype = float)
#    TDEMData = []
#    for file in glob.glob("*FWR*") :
#        TDEMData.append(np.loadtxt(file))
#
#    # Get the time channels from the first station
#    tchannels = TDEMData[0][:,3]
#
#    # Extract dBzdt for each time channel on all stations
#    dBzdt = []
#
#    for stn in TDEMData:
#
#        dBzdt.append(np.c_[stn[:,0:4], stn[:,-1]])
#
#    # Transform into a 3D array (stationID x times x (xyztd))
#    dBzdt = np.asarray(dBzdt)
#
#    # Special figures
#    plt.figure()
#
#    # Find all unique NS lines and grab the first one
#    stnX = np.unique(dBzdt[:,0,0])
#    indx = np.where(dBzdt[:,0,0] == stnX[1])[0]
#    tsort = dBzdt[indx,0,1].argsort()
#    for tt in range(len(tchannels)):
#        plt.semilogy(dBzdt[indx[tsort],tt,1],dBzdt[indx[tsort],tt,4])

    # Load in data matrix

    offTile = 500
    dx = 250


    gCx = np.linspace(xy[:,0].min(),xy[:,0].max(),100)
    gCy = np.linspace(xy[:,1].min(),xy[:,1].max(),100)

    X, Y = np.meshgrid(gCx, gCy)

    aa=np.asarray(range(18*46))
    aa = aa.reshape((18,46))
    ids = aa[::2,::2]
    
#    txOut = np.c_[mkvc(ids.T),mkvc(X[::2,::2]),mkvc(Y[::2,::2])]
#
#    fid = open(work_dir + '\\TxLocs.xyz', 'w')
#    np.savetxt(fid,txOut,fmt='%i',delimiter=' ',newline='\n')
#    fid.close()
#    
#    fig, axs = plt.figure(figsize=(6,10)), plt.subplot(111)
#    axs.scatter(txOut[:,1],txOut[:,2],c='r',s=20,edgecolor=None)
#    axs.set_aspect('equal')
#    
#    for ii in range(txLoc.shape[0]):
#
#        axs.text(txOut[ii,1],txOut[ii,2], str(int(txOut[ii,0])),
#            verticalalignment='bottom', horizontalalignment='center',
#            color='k', fontsize=10)
#    
#    axs.set_xlim(X.min()-500,X.max()+500)
#    axs.set_ylim(Y.min()-500,Y.max()+500)
    
    
    # Plot
    fig = plt.figure(figsize=(6,8))
    ax1 = plt.subplot()

    def animate(ii):

        removeFrame()

        dgrid = sp.interpolate.griddata(xy, data[:,ii],(X, Y))

        global ax1, fig

        ax1 = plt.subplot()
        plt.contourf(gCx,gCy,np.log10(dgrid),100)
        plt.scatter(xy[:,0], xy[:,1], c='k')
        ax1.set_xticks(np.linspace(xy[:,0].min(),xy[:,0].max(),3))
        ax1.set_yticks(np.linspace(xy[:,1].min(),xy[:,1].max(),3))
#        plt.colorbar(orientation='horizontal')
        ax1.set_aspect('equal')
        ax1.set_title('Time:' + str(int(times[ii]*1e+6)) + ' usec')


    def removeFrame():

        global ax1, fig
        fig.delaxes(ax1)
        plt.draw()

    anim = animation.FuncAnimation(fig, animate,
                                   frames=data.shape[1] , interval=1000, repeat = False)
    #/home/dominiquef/3796_AGIC_Research/DCIP3D/MtISa
    anim.save('Data_slice.html', writer=HTMLWriter(embed_frames=True,fps=1))

if dType == "GTEM":

    os.chdir(work_dir + inp_dir[dType])

    offTile = 500
    lenTile = 500


    gCx = np.asarray(range(9))*lenTile + offTile + mesh_global.x0[0]
    gCy = np.asarray(range(23))*lenTile + offTile + mesh_global.x0[1]

    X, Y = np.meshgrid(gCx, gCy)

    aa=np.asarray(range(9*23))
    ids = aa.reshape((9,23))

    txOut = np.c_[mkvc(ids.T),mkvc(X),mkvc(Y)]

    fid = open(work_dir + '\\TxLocs.xyz', 'w')
    np.savetxt(fid,txOut,fmt='%i',delimiter=' ',newline='\n')
    fid.close()
    
#    fig, axs = plt.figure(figsize=(6,10)), plt.subplot(111)
#    axs.scatter(txOut[:,1],txOut[:,2],c='r',s=20,edgecolor=None)
#    axs.set_aspect('equal')
#    
#    for ii in range(txLoc.shape[0]):
#
#        axs.text(txOut[ii,1],txOut[ii,2], str(int(txOut[ii,0])),
#            verticalalignment='bottom', horizontalalignment='center',
#            color='k', fontsize=10)
#    
#    axs.set_xlim(X.min()-500,X.max()+500)
#    axs.set_ylim(Y.min()-500,Y.max()+500)
    
    
    AllTiles = []
    for sourceFile in glob.glob("*Tile*"):

        with open(sourceFile,'r') as sfile:
            lines = sfile.readlines()
            stn = []
            for line in lines:

                if line[0]!='%':

                    temp = map(float,line.split())

                    if len(temp)>0:
                        stn.append(temp)


        # All ground loop data
        AllTiles.append(np.asarray(stn))

    times = np.unique(AllTiles[0][:,3])
    ntimes = len(times)

    # Re-order data over time


    # Plot
    fig = plt.figure(figsize=(6,8))
    ax1 = plt.subplot()

    def animate(ii):

        removeFrame()

        for tile in AllTiles:

            xloc = np.unique(tile[:,0])
            yloc = np.unique(tile[:,1])

            indx = tile[:,3] == times[ii]
            vmin, vmax = np.log10(np.min(np.abs(AllTiles[0][indx,-1])))-2., np.log10(np.max(np.abs(AllTiles[0][indx,-1])))+2.

            dgrid = np.abs(tile[indx,-1].reshape((len(xloc),len(yloc)), order='F'))

            global ax1, fig

            ax1 = plt.subplot()
            plt.contourf(xloc,yloc,np.log10(dgrid),100, clim=(vmin,vmax), vmin=vmin, vmax=vmax )
            plt.scatter(ore[:,0],ore[:,1],c=ore[:,2], s = 1, cmap='RdBu', linewidth=0)
#            ax1.set_xticks(np.linspace(xloc.min(),xloc.max(),3))
#            ax1.set_yticks(np.linspace(yloc.min(),yloc.max(),3))
    #        plt.colorbar(orientation='horizontal')
        ax1.set_aspect('equal')
        ax1.set_title('Time:' + str(int(times[ii]*1e+6)) + ' usec')


    def removeFrame():

        global ax1, fig
        fig.delaxes(ax1)
        plt.draw()

    anim = animation.FuncAnimation(fig, animate,
                                   frames=ntimes , interval=1000, repeat = False)
    #/home/dominiquef/3796_AGIC_Research/DCIP3D/MtISa
    anim.save('Data_slice.mp4', writer='ffmpeg')
#    anim.save('Data_slice.html', writer=HTMLWriter(embed_frames=True,fps=1))
