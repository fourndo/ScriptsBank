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


inp_dir = {'DC':'\\FWR_DC',
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


dType = 'DC'

ndv = {'MAG':0, 'Gz':0, 'DC':1e-8, 'FDEM':1e-8, 'TDEM':1e-8}

npadxy = {'DC':10,'TDEM':7, }
npadz = {'DC':10, 'TDEM':10, }
expf = 1.3

offTile = 500
lenTile = 250

# Original grid from "DataSimulator"
gCx = np.asarray(range(18))*lenTile + offTile + mesh_global.x0[0]
gCy = np.asarray(range(46))*lenTile + offTile + mesh_global.x0[1]

X,Y = np.meshgrid(gCx,gCy)

aa=np.asarray(range(18*46))
aa = aa.reshape((18,46))
ids = aa[::2,::2]

txOut = np.c_[mkvc(ids.T),mkvc(X[::2,::2]),mkvc(Y[::2,::2])]
    
    
X, Y = mkvc(X), mkvc(Y)

dx = [20., 40.]

tID = [118,124]

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
             lines[2]=re.sub("E\+0212100",'E+02 12100',lines[2])
             tfile.writelines(lines)

    print "DONE"


    os.chdir(work_dir + inp_dir[dType] + '\\' + target_dir)
    DcData = []
    files = glob.glob("*Tile*")
    txLoc = []
    for file in files :
        
        
        obj = readUBC_DC3Dobs(file)

        DcData.append(obj['DCsurvey'])
        
        txLoc.append(obj['DCsurvey'].srcList[0].loc[0])
        
    txLoc = np.asarray(txLoc)


    fid = open(work_dir + '\\TxLocs.xyz', 'w')
    np.savetxt(fid,txOut,fmt='%i',delimiter=' ',newline='\n')
    fid.close()
    
    fig, axs = plt.figure(figsize=(6,10)), plt.subplot(111)
    axs.scatter(txOut[:,1],txOut[:,2],c='r',s=20,edgecolor=None)
    axs.set_aspect('equal')
    
    for ii in range(txLoc.shape[0]):

        axs.text(txOut[ii,1],txOut[ii,2], str(int(txOut[ii,0])),
            verticalalignment='bottom', horizontalalignment='center',
            color='k', fontsize=10)
    
    axs.set_xlim(X.min()-500,X.max()+500)
    axs.set_ylim(Y.min()-500,Y.max()+500)

    #%%    ## TESTING CODE
    # Generate survey using two tiles, compute solution using UBC-DCIP3D
    # and compare by using the pre-computed grid
    padxy = np.r_[dx[1]*expf**(np.asarray(range(npadxy['DC']))+1)]
    padz = np.r_[dx[1]*expf**(np.asarray(range(npadz['DC']))+1)]
    
    hx = np.r_[padxy[::-1], core['DC'], padxy]
    hy = np.r_[padxy[::-1], core['DC'],core['DC'], padxy]
    hz = np.r_[padz[::-1],[33],np.ones(25)*30,[23,18,15,12], np.ones(15)*10,[12,15,18,23],30*expf**(np.asarray(range(npadz['DC'])))]

    mesh = Mesh.TensorMesh([hx,hy,hz], 'CC0')
        
    mesh._x0 = (mesh.x0[0] + np.mean(X[tID]), mesh.x0[1]+np.mean(Y[tID]), mesh.x0[2]-(np.sum(hz[:(npadz['DC']+20)])))
    
    # Extract model from global tolocal mesh
    actv2 = Utils.surface2ind_topo(mesh, topo, 'N')
        
    P = Maps.Mesh2MeshTopo([mesh_global,mesh],[actv,actv2],nIterpPts = 12)  

    model = mesh_global.readModelUBC(work_dir + '\\' + modelfile['DC'])
         
    model_Tile = P*(model[actv])

    m = np.ones(mesh.nC)*1e-8
    m[actv2] = model_Tile
#    mtemp = m.reshape(mesh.vnC, order='F')

    Mesh.TensorMesh.writeUBC(mesh,work_dir +'\\MeshTile.msh')
    Mesh.TensorMesh.writeModelUBC(mesh,work_dir +'\\MeshTile.mod',m)
            
    # Create survey
    a = 50
    b = 525
    n = None
    endl = np.c_[np.r_[X[tID]],np.r_[Y[tID]]]
    
    survey = gen_DCIPsurvey(endl, mesh, 'gradient', a, b, n)
            
    # Extract survey points
    RxLocs = survey.srcList[0].rxList[0].locs
    

    #
    locRx1 = drapeTopotoLoc(mesh, topo, np.c_[RxLocs[0][:,0],RxLocs[0][:,1]], airind=actv2)
    locRx2 = drapeTopotoLoc(mesh, topo, np.c_[RxLocs[1][:,0],RxLocs[1][:,1]], airind=actv2)

    locTx1 = drapeTopotoLoc(mesh, topo, np.c_[endl[0,0],endl[0,1]], airind=actv2)
    locTx2 = drapeTopotoLoc(mesh, topo, np.c_[endl[1,0],endl[1,1]], airind=actv2)

    fid = open(work_dir + '\\DCObslocs.dat', 'w')
    np.savetxt(fid, np.r_[locTx1[0],locTx2[0],locRx1.shape[0]].reshape((1,7)), fmt='%e',delimiter=' ',newline='\n')
    np.savetxt(fid, np.c_[locRx1,locRx2,np.ones_like(locRx1)], fmt='%e',delimiter=' ',newline='\n')
    fid.close()
    
    # Write ibs locs for plotting
    fid = open(work_dir + '\\XYZlocs.dat', 'w')
    np.savetxt(fid, locRx1, fmt='%e',delimiter=' ',newline='\n')
    fid.close()
    
    
    # Run forward model
    os.chdir(work_dir)
    os.system('dcipf3d dcipf3d.inp')
            
    FWR_dc = readUBC_DC3Dobs(work_dir + '\\dc3d.dat')
    FWR_survey = FWR_dc['DCsurvey']
    
    # Extract data from pre-computed grids
    # Grid points
    grid_x1, grid_y1 = np.meshgrid(np.unique(FWR_survey.srcList[0].rxList[0].locs[0][:,0]),
                                   np.unique(FWR_survey.srcList[0].rxList[0].locs[0][:,1]))

    grid_x2, grid_y2 = np.meshgrid(np.unique(FWR_survey.srcList[0].rxList[0].locs[1][:,0]),
                                   np.unique(FWR_survey.srcList[0].rxList[0].locs[1][:,1]))
    
    tid = np.argmin(np.abs(endl[0,0] - txLoc[:,0])+np.abs(endl[0,1] - txLoc[:,1]))
    
    P1_phi1 = scipy.interpolate.griddata(DcData[tid].srcList[0].rxList[0].locs[0][:,0:2], 
                                         DcData[tid].dobs,
                                         (grid_x1, grid_y1), method='linear')

    P2_phi1 = scipy.interpolate.griddata(DcData[tid].srcList[0].rxList[0].locs[0][:,0:2], 
                                         DcData[tid].dobs, 
                                         (grid_x2, grid_y2), method='linear')

    tid = np.argmin(np.abs(endl[1,0] - txLoc[:,0])+np.abs(endl[1,1] - txLoc[:,1]))
    P1_phi2 = scipy.interpolate.griddata(DcData[tid].srcList[0].rxList[0].locs[0][:,0:2], 
                                         DcData[tid].dobs,
                                         (grid_x1, grid_y1), method='linear')

    P2_phi2 = scipy.interpolate.griddata(DcData[tid].srcList[0].rxList[0].locs[0][:,0:2], 
                                         DcData[tid].dobs,
                                         (grid_x2, grid_y2), method='linear')


    locRx1 = np.c_[mkvc(grid_x1),mkvc(grid_y1)]
    locRx2 = np.c_[mkvc(grid_x2),mkvc(grid_y2)]
    
    Pmid = (locRx1 + locRx2)/2.
    
    rC1P1 = np.sqrt( np.sum( (npm.repmat(endl[0,:],FWR_survey.nD, 1) - locRx1)**2, axis=1 ))
    rC2P1 = np.sqrt( np.sum( (npm.repmat(endl[1,:],FWR_survey.nD, 1) - locRx1)**2, axis=1 ))
    rC1P2 = np.sqrt( np.sum( (npm.repmat(endl[0,:],FWR_survey.nD, 1) - locRx2)**2, axis=1 ))
    rC2P2 = np.sqrt( np.sum( (npm.repmat(endl[1,:],FWR_survey.nD, 1) - locRx2)**2, axis=1 ))

    #rC1C2 = np.sqrt( np.sum( (npm.repmat(Tx[0][0:2]-Tx[0][3:5],Rx[0].shape[0], 1) )**2, axis=1 ))
    #rP1P2 = np.sqrt( np.sum( (Rx[0][:,0:2] - Rx[1][:,0:2])**2, axis=1 ))
    volt = ((P1_phi1 - P1_phi2) - (P2_phi1 - P2_phi2))


    rho = np.abs(mkvc(volt)) *np.pi *2./ ( 1/rC1P1 - 1/rC2P1 - 1/rC1P2 + 1/rC2P2 )#*((rC1P1)**2 / rP1P2)#

    FWR_mid = np.c_[np.mean(np.c_[FWR_survey.srcList[0].rxList[0].locs[0][:,0],FWR_survey.srcList[0].rxList[0].locs[1][:,0]],axis=1),
                    np.mean(np.c_[FWR_survey.srcList[0].rxList[0].locs[0][:,1],FWR_survey.srcList[0].rxList[0].locs[1][:,1]],axis=1)]
    FWR_volt = scipy.interpolate.griddata(FWR_mid, 
                                          FWR_survey.dobs, 
                                          (Pmid[:,0].reshape(grid_x1.shape,order='F'), Pmid[:,1].reshape(grid_x1.shape,order='F')), method='linear')
 
    rho2 = np.abs(mkvc(FWR_volt)) *np.pi *2./ ( 1/rC1P1 - 1/rC2P1 - 1/rC1P2 + 1/rC2P2 )
#    survey.dobs = mkvc(rho)
#    survey.std = np.ones(survey.nD)
#
    


    plt.figure()
    ax_prim = plt.subplot(1,2,1)
    plt.gca().set_aspect('equal', adjustable='box')

    plt.contourf(np.unique(FWR_survey.srcList[0].rxList[0].locs[0][:,0]),
                 np.unique(FWR_survey.srcList[0].rxList[0].locs[0][:,1]),
                rho.reshape(volt.shape,order='F'), 20, vmin=rho2.min(), vmax=rho2.max(), cmap='magma_r')
    plt.clim(rho2.min(), rho2.max())
    plt.scatter(Pmid[:,0],Pmid[:,1],c=rho, vmin=rho2.min(), vmax=rho2.max(), cmap='magma_r')
    plt.colorbar(ax = ax_prim)
    plt.title('Tiled')
    plt.show()
    
    ax_prim = plt.subplot(1,2,2)
    plt.gca().set_aspect('equal', adjustable='box')

    plt.contourf(np.unique(FWR_survey.srcList[0].rxList[0].locs[0][:,0]),
                 np.unique(FWR_survey.srcList[0].rxList[0].locs[0][:,1]),
                rho2.reshape(volt.shape,order='F'), 20, vmin=rho2.min(), vmax=rho2.max(), cmap='magma_r' )
    plt.clim(rho2.min(), rho2.max())
    plt.scatter(Pmid[:,0],Pmid[:,1],c=rho2, vmin=rho2.min(), vmax=rho2.max(), cmap='magma_r' )
    plt.colorbar(ax = ax_prim)
    plt.title('Full')
    plt.show()





