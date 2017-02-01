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
from SimPEG.EM.Static.Utils import drapeTopotoLoc, gen_DCIPsurvey, writeUBC_DCobs
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import time
import gc
from test_TDEM_forward_Analytic import halfSpaceProblemAnaDiff
from shutil import copyfile
import os
#from pymatsolver import PardisoSolver as Solver
from pymatsolver import BicgJacobiSolver as Solver
#%% Input parameters

work_dir = 'C:\\LC\\Private\\dominiquef\\Projects\\4414_Minsim\\Modeling\\Forward'
meshfile = 'Mesh_20m.msh'
locfile = 'Mag_grid_50m.dat'
modelfile = {'MAG': 'Mesh_20m_Susc\\maginv3d_001.sus',
             'DC': 'Mesh_20m_Cond\\Cond_Randn_Std_model.con',
             'Gz': 'Mesh_20m_Dens\\Density_Randn_Std.den'}

nullcell = 'nullcell.dat'

obsfile = 'Mesh_20m_Susc\\Mag_grid_50m.dat'
topofile = 'ROT_DEM_30m.topo'

dtype = ['DC']
# dtype = ['DC']
gc.collect()

# %% Load files
mesh_global = Mesh.TensorMesh.readUBC(work_dir + '\\' + meshfile)
locs = np.loadtxt(work_dir + '\\' + locfile)
model = {'MAG': mesh_global.readModelUBC(work_dir + '\\' + modelfile['MAG']),
         'Gz': mesh_global.readModelUBC(work_dir + '\\' + modelfile['Gz']),
         'DC': mesh_global.readModelUBC(work_dir + '\\' + modelfile['DC'])}

topo = np.genfromtxt(work_dir + '\\' + topofile,
                     skip_header=1)

actvmod = mesh_global.readModelUBC(work_dir + '\\' + nullcell)
actv = np.ones(mesh_global.nC, dtype=bool)
actv[actvmod == 0] = False

srvy = PF.MagneticsDriver.MagneticsDriver_Inv()
srvy.basePath = work_dir
data = srvy.readMagneticsObservations('\\'+obsfile)

locs = data.srcField.rxList[0].locs

offTile = 500
lenTile = 500

dx = [20., 40.]
#
nCx = np.floor((500/dx[0]))-1

#core = np.r_[np.ones(12)*dx[1], [20,15], np.ones(10)*dx[0], [15,20], np.ones(12)*dx[1]]

core = {'DC': np.r_[np.ones(110)*dx[0]],
        'ATEM': np.r_[np.ones(5)*dx[1], [33, 26], np.ones(10)*dx[0],
                      [26, 33], np.ones(5)*dx[1]],
        'GTEM': np.r_[dx[1], [33, 26], np.ones(nCx)*dx[0],
                      [26, 33], dx[1]],
        'Gz': np.r_[np.ones(15)*dx[1], [20, 15], np.ones(nCx)*dx[0],
                    [15, 20], np.ones(15)*dx[1]],
        'MAG': np.r_[np.ones(15)*dx[1], [20, 15], np.ones(nCx)*dx[0],
                     [15, 20], np.ones(15)*dx[1]]}

ndv = {'MAG': 0, 'Gz': 0, 'DC': 1e-8, 'FDEM': 1e-8, 'ATEM': 1e-8, 'GTEM': 1e-8}

npadxy = {'DC': 10, 'ATEM': 7, 'GTEM': 7, 'Gz': 7, 'MAG': 7}
npadz = {'DC': 10, 'ATEM': 10, 'GTEM': 10, 'Gz': 10, 'MAG': 10}

gCx = np.asarray(range(9))*lenTile + offTile + mesh_global.x0[0]
gCy = np.asarray(range(23))*lenTile + offTile + mesh_global.x0[1]

X, Y = np.meshgrid(gCx, gCy)

indT = np.zeros(X.shape, dtype=bool)
indT[:, :] = True
indT[:, :] = True

X, Y, indT = mkvc(X), mkvc(Y), mkvc(indT)

# %% TEST script
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


#%% Loop through the tiles and create meshs
meshes = []
fig, axs = plt.figure(figsize=(6,10)), plt.subplot(111)
axs.scatter(locs[:,0],locs[:,1],c=data.dobs,s=10,edgecolor=None)

axs.scatter(mkvc(X),mkvc(Y),c='r',s=20,edgecolor=None)
axs.set_aspect('equal')

for ii in range(len(X)):

    if indT[ii]:
        axs.text(X[ii], Y[ii], str(ii),
            verticalalignment='center', horizontalalignment='center',
            color='k', fontsize=10)

#    axs.add_patch(
#    patches.Rectangle(
#        (mtemp.x0[0],mtemp.x0[1]),   # (x,y)
#        np.sum(mtemp.hx),          # width
#        np.sum(mtemp.hy), fill = False   ))

axs.set_xlim(X.min(),X.max())
axs.set_ylim(Y.min(),Y.max())
expf = 1.3

tID = 0

#%% Loop through the meshes and forward model

DcData = []
MagData = []
GzData = []
FdemData = []

os.chdir(work_dir)

tiles = range(len(X[:]))

if any(x in 'FDEM' for x in dtype):
    fid = open(work_dir + '\\FWR_FDEM\\FDEM_All' + str(tID) + '.obs', 'w')

# Initialize ckdTree and store after the first time it is called
ckdTree = None

TxLoop = []
for tID in tiles:
    # Create active cells for sub mesh

    # actv2 = np.asarray([inds for inds, elem in enumerate(actv2, 1)
    #                   if elem], dtype=int) - 1

    if not indT[tID]:
        continue

    # if not np.any(tID == np.r_[3218]):
    #     continue

    for dID in dtype:

        start_time = time.time()

        padxy = np.r_[dx[1]*expf**(np.asarray(range(npadxy[dID]))+1)]
        padz = np.r_[dx[1]*expf**(np.asarray(range(npadz[dID]))+1)]

        hx = np.r_[padxy[::-1], core[dID], padxy]
        hy = np.r_[padxy[::-1], core[dID], padxy]
#        hz = np.r_[padb*2.,[33,26],np.ones(25)*22,[18,15,12,10,8,7,6], np.ones(18)*5,5*expf**(np.asarray(range(2*npad)))]
        hz = np.r_[padz[::-1],[33],np.ones(25)*30,[23,18,15,12], np.ones(15)*10,[12,15,18,23],30*expf**(np.asarray(range(npadz[dID])))]

        mesh = Mesh.TensorMesh([hx,hy,hz], 'CC0')

    #    mtemp._x0 = [x0[ii]-np.sum(padb), y0[ii]-np.sum(padb), mesh.x0[2]]

        mesh._x0 = (mesh.x0[0] + X[tID], mesh.x0[1]+Y[tID], mesh.x0[2]-(np.sum(hz[:(npadz[dID]+20)])))
#        meshes.append(mtemp)


        # Extract model from global to local mesh
        actv2 = Utils.surface2ind_topo(mesh, topo, 'N')

        if ckdTree is None:
            ckdTree = Maps.Mesh2MeshTopo([mesh_global,mesh],[actv,actv2],nIterpPts = 12)

        P = Maps.Mesh2MeshTopo([mesh_global,mesh],[actv,actv2],nIterpPts = 12, tree=ckdTree.tree )

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
        print("Model Interpolated --- %s seconds ---" % (time.time() - start_time))


        if dID=='MAG':

            # OPTION TO WRITE OUT
            Mesh.TensorMesh.writeUBC(mesh,work_dir +'\\MeshTile1.msh')
            Mesh.TensorMesh.writeModelUBC(mesh,work_dir +'\\MeshTile1.mod',m)

            # Grab the data and forward ##NEED TO GENERATE DATA ON THE FLY
#            indx = (locs[:,0] > X[tID] - nCx/2*dx[0]) & (locs[:,0] < X[tID] + nCx/2*dx[0]) & (locs[:,1] > Y[tID] - nCx/2*dx[0]) & (locs[:,1] < Y[tID] + nCx/2*dx[0])

            # Create grid of points
            spacing=(20,40)
            width = (500,1500)
            rxLoc = grid_survey(spacing,width,x0 = (X[tID],Y[tID], 60.), topo = topo)

            # INTEGRAL FORMULATION
#            subrx = PF.BaseMag.RxObs(rxLoc)
#
#            srcParam = np.asarray((60000,79,31))
#            srcField = PF.BaseMag.SrcField([subrx],param=srcParam)
#            survey = PF.BaseMag.LinearSurvey(srcField)
#
#            fid = open(work_dir + '\\MAGObslocs.dat', 'w')
#            np.savetxt(fid, rxLoc, fmt='%e',delimiter=' ',newline='\n')
#            fid.close()
#
#            # Create PF problem
#            nC = int(np.sum(actv2))
#            idenMap = Maps.IdentityMap(nP = nC)
#            # Create active map to go from reduce set to full
#            actvMap = Maps.InjectActiveCells(mesh, actv2, -100)
#
#            prob = PF.Magnetics.MagneticIntegral(mesh, chiMap = idenMap, actInd = actv2, forwardOnly = True)
#
#            survey.pair(prob)
#            start_time = time.time()
#            d = prob.fields(model_Tile)
#            print("Solve Integral --- %s seconds ---" % (time.time() - start_time))

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

            fid = open(work_dir + '\\FWR_MAG\\MAG_Tile' + str(tID) + '.obs', 'w')
            np.savetxt(fid, np.c_[rxLoc,data], fmt='%e',delimiter=' ',newline='\n')
            fid.close()

#            PF.Magnetics.writeUBCobs(work_dir + '\\MAG_Tile' + str(tID) + '.obs', survey,  data)
#            p1 = PF.Magnetics.plot_obs_2D(rxLoc,d, vmin=d.min(), vmax=d.max())
#            p2 = PF.Magnetics.plot_obs_2D(rxLoc,dtmi_PDE, vmin=d.min(), vmax=d.max())
            del prob, u, survey, m_map, data
        elif dID == 'Gz':

            # Create grid of points
            spacing=(20,40)
            width = (500,1500)
            locXyz = grid_survey(spacing,width,x0 = (X[tID],Y[tID], 60.), topo = topo)

            Mesh.TensorMesh.writeUBC(mesh,work_dir +'\\MeshTile.msh')
            Mesh.TensorMesh.writeModelUBC(mesh,work_dir +'\\MeshTile.sus',m)
#
#            nC = int(np.sum(actv2))
#            idenMap = Maps.IdentityMap(nP = nC)
#
#            rxLoc = PF.BaseGrav.RxObs(locXyz)
#            srcField = PF.BaseGrav.SrcField([rxLoc])
#            survey = PF.BaseGrav.LinearSurvey(srcField)
#
#            prob_z = PF.Gravity.GravityIntegral(mesh, rhoMap=idenMap,
#                                                actInd = actv2,
#                                                forwardOnly=True,
#                                                rtype='z')
#
#
#            # Compute 3-component mag data
#            survey.pair(prob_z)
#            dInt = prob_z.fields(m[actv2])
#            print("Solve Integral --- %s seconds ---" % (time.time() - start_time))

            # Grab the data and forward ##NEED TO GENERATE DATA ON THE FLY
#            indx = (locs[:,0] > X[tID] - nCx/2*dx[0]) & (locs[:,0] < X[tID] + nCx/2*dx[0]) & (locs[:,1] > Y[tID] - nCx/2*dx[0]) & (locs[:,1] < Y[tID] + nCx/2*dx[0])



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

#            data = np.c_[locXyz,data]
#            GzData.append(data)

#            PF.Gravity.writeUBCobs(work_dir + '\\GRAV_Tile' + str(tID) + '.obs', survey,  data)
            fid = open(work_dir + '\\FWR_Gz\\GRAV_Tile_' + str(tID) + '.obs', 'w')
            np.savetxt(fid, np.c_[locXyz,data], fmt='%e',delimiter=' ',newline='\n')
            fid.close()

            del prob, u, survey, m_map, data

#            p1 = PF.Magnetics.plot_obs_2D(survey.srcField.rxList[0].locs,dInt, vmin=dInt.min(), vmax=dInt.max())
#            p2 = PF.Magnetics.plot_obs_2D(survey.srcField.rxList[0].locs,gz, vmin=gz.min(), vmax=gz.max())
#%%
        elif dID == 'DC':

            Mesh.TensorMesh.writeUBC(mesh,work_dir +'\\MeshTile.msh')
            Mesh.TensorMesh.writeModelUBC(mesh,work_dir +'\\MeshTile.mod',m)
            stype = 'gradient'

            # Survey parameters
#            a = 50
#            b = 500
#            n = 40
#
#            # Forward solver
#            slvr = 'BiCGStab' #'LU'
#
#            # Preconditioner
#            pcdr = 'Jacobi'#'Gauss-Seidel'#
#
#            # Inversion parameter
#            pct = 0.01
#            flr = 1e-4
#            chifact = 100
#            ref_mod = 1e-1
#
#            #%% Create system
#            #Set boundary conditions
#            mesh.setCellGradBC('neumann')
#
#            m = np.ones(mesh.nC)*1e-10
#            m[actv2] = model_Tile
#
#            Div = mesh.faceDiv
#            Grad = mesh.cellGrad
#            Msig = Utils.sdiag(1./(mesh.aveF2CC.T*(1./m)))
#
#            A = Div*Msig*Grad
#
#            # Change one corner to deal with nullspace
#            A[0,0] = 1./mesh.vol[0]
#            A = sp.csc_matrix(A)
#
#            if slvr=='BiCGStab':
#                # Create Jacobi Preconditioner
#                if pcdr=='Jacobi':
#                    dA = A.diagonal()
#                    PC = sp.spdiags(1/dA,0,A.shape[0],A.shape[0])
#
#                # Create Gauss-Seidel Preconditioner
#                elif pcdr=='Gauss-Seidel':
#                    LD = sp.tril(A,k=0)
#                    #LDinv = sp.linalg.splu(LD)
#
#            elif slvr=='LU':
#                # Factor A matrix
#                Ainv = sp.linalg.splu(A)
#                print("LU DECOMP--- %s seconds ---" % (time.time() - start_time))
#
#            elif slvr=='Pardiso':
#                Ainv = PardisoSolver(A)
#
##            Msig = mesh.getFaceInnerProduct(1/m)
##            iMsig= Utils.sdInv(Msig)
##            #Volume diagonal, Projection and inside DIV matrices
##            DIV = mesh.faceDiv
##            V = Utils.sdiag(mesh.vol)
##            Pbc, Pin, Pout = mesh.getBCProjWF('neumann')
##            Din = DIV*Pin.T*Pin
##            #Build the forward model matrix (fv for finite volume)
###            Afvdir = V*DIV*iMsig*DIV.T*V
##            Afvneu = V*Din*iMsig*Din.T*V
##            Afvneu[0,0] = 1 /mesh.vol[0]
##            #Afvneu[0,0] /= mesh.vol[0] it did not work... try it to check
##            #seems a too big pertubations does not allow for solving it (vol is big here)
##
##            #solving for Neumann
##            Ainvneu = Solver(Afvneu)
#
#            # Get the cell in center of mesh
#            midx, midy = mesh.vectorCCx[mesh.nCx/2], mesh.vectorCCy[mesh.nCy/2]
            tx = mkvc(drapeTopotoLoc(mesh, topo, np.c_[X[tID],Y[tID]], airind=actv2) )
#
##            rxLocx = mesh.vectorCCx[2*npad:-2*npad]
##            rxLocy = mesh.vectorCCy[2*npad:-2*npad]
#
#            # Create grid of points
            spacing=(20)
            width = (2200)
            rxLoc = grid_survey([spacing],[width],x0 = (X[tID]-5.,Y[tID]-5., 0.))
##            rxLocx = np.r_[mkvc(rxCorex),mkvc(rxMarginx)] + X[tID]
##            rxLocy = np.r_[mkvc(rxCorey),mkvc(rxMarginy)] + Y[tID]
#            rxLoc2 = grid_survey([spacing],[width],x0 = (X[tID],Y[tID]+20., 0.))
##            rxLocx,rxLocy = np.meshgrid(rxLocx,rxLocy)
#
            rxLocDrapped = drapeTopotoLoc(mesh, topo, rxLoc[:,:2], airind=actv2)
#            rxLocDrapped2 = drapeTopotoLoc(mesh, topo, rxLoc2[:,:2], airind=actv2)
#
#
#            tinf = np.array([mesh.vectorCCx[-1],mesh.vectorCCy[-1],mesh.vectorCCz[0]])
#            inds = Utils.closestPoints(mesh, np.c_[tx].T)
#
#            src = mesh.gridCC[inds]
#            RHS = mesh.getInterpolationMat(src, 'CC').T*( [-1] / mesh.vol[inds] )
#
#            # Solve for phi on pole locations
#            P1 = mesh.getInterpolationMat(rxLocDrapped, 'CC')
#            P2 = mesh.getInterpolationMat(rxLocDrapped2, 'CC')
#
#            if slvr=='BiCGStab':
#
#                if pcdr=='Jacobi':
#
#                    # Iterative Solve
#                    Ainvb = sp.linalg.bicgstab(PC*A,PC*RHS, tol=1e-5)
#
#                # Create Gauss-Seidel Preconditioner
#                elif pcdr=='Gauss-Seidel':
#                    LD = sp.tril(A,k=0)
#
#
#                phi = mkvc(Ainvb[0])
#
#            elif slvr=='LU':
#                #Direct Solve
#                phi = Ainv.solve(RHS)
#
#            elif slvr == 'Pardiso':
#                #Direct Solve
#                phi = Ainvneu*(RHS)
#
#            # Compute potential at each electrode
#            dtemp = np.c_[rxLocDrapped,(P2*phi-P1*phi)*np.pi]
#            DcData.append(dtemp)

            #%% USE UBC DCIPF3D
            # Write locations out
            fid = open(work_dir + '\\DCObslocs.dat', 'w')
            np.savetxt(fid, np.r_[tx,tx,rxLocDrapped.shape[0]].reshape((1,7)), fmt='%e',delimiter=' ',newline='\n')
            np.savetxt(fid, np.c_[rxLocDrapped,rxLocDrapped], fmt='%e',delimiter=' ',newline='\n')
            fid.close()

            # Run forward model
            os.system('dcipf3d dcipf3d.inp')

            copyfile(work_dir + '\\dc3d.dat', work_dir + '\\FWR_DC\\Tile_' + str(tID) + '.obs')
            # Write out model
#            Mesh.TensorMesh.writeUBC(mesh,work_dir +'\\MeshTile1.msh')
#            Mesh.TensorMesh.writeModelUBC(mesh,work_dir +'\\MeshTile1.con',m)

            print("Solve Tile %i --- %s seconds ---" % (tID, time.time() - start_time))
#%%
        elif dID == 'FDEM':

            # OPTION TO WRITE OUT
            Mesh.TensorMesh.writeUBC(mesh,work_dir +'\\MeshTile1.msh')
            Mesh.TensorMesh.writeModelUBC(mesh,work_dir +'\\MeshTile1.mod',m)


            start_time = time.time()

            # Create grid of points
#            spacing=(40,80)
#            width = (500,2000)
#            rxLoc = grid_survey(spacing,width,x0 = (X[tID],Y[tID], 60.), topo = topo)



            rxLoc = grid_survey([1],[1],x0 = (X[tID],Y[tID], 60.), topo=topo)
#            rxlocs = Utils.ndgrid([np.r_[50.], np.r_[0], np.r_[0.]])
            bzi = FDEM.Rx.Point_bSecondary(rxLoc, 'z', 'real')
            bzr = FDEM.Rx.Point_bSecondary(rxLoc, 'z', 'imag')

            freqs = np.asarray([400.])
#            freqs = np.asarray([400.,1500.,8000.,40000.,150000.])
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

            print("Solve Tile %i --- %s seconds ---" % (tID, time.time() - start_time))

#            dtemp =

            np.savetxt(fid, np.r_[srcLoc[0],surveyFD.dobs].reshape((1,13)), fmt='%e',delimiter=' ',newline='\n')

#            FdemData.append(dtemp)

        elif dID == 'ATEM':

            # OPTION TO WRITE OUT
            Mesh.TensorMesh.writeUBC(mesh,work_dir +'\\FWR_TDEM\\INP_DIR\\Tile_' + str(tID) + '.msh')
            Mesh.TensorMesh.writeModelUBC(mesh,work_dir +'\\FWR_TDEM\\INP_DIR\\Tile_' + str(tID) + '.con',m)

#            halfSpaceProblemAnaDiff
            start_time = time.time()

            rxLoc = grid_survey([1],[1],x0 = (X[tID],Y[tID], 30.), topo=topo)

            srcLoc = rxLoc
            # Create grid of points
#            spacing=(40,80)
#            width = (500,2000)
#            rxLoc = grid_survey(spacing,width,x0 = (X[tID],Y[tID], 60.), topo = topo)
#            times = np.logspace(-4, np.log10(2e-3), 10)
            times = np.loadtxt(work_dir + '\\VTEM_times.dat')
#            print('min diffusion distance ', 1.28*np.sqrt(times.min()/(sig_half*mu_0)),
#                  'max diffusion distance ', 1.28*np.sqrt(times.max()/(sig_half*mu_0)))
#            rx = TDEM.Rx(srcLoc, times, 'dbzdt')
#            src = TDEM.Src.MagDipole(
#                [rx],
#                waveform=TDEM.Src.StepOffWaveform(),
#                loc=rxLoc  # same src location as FDEM problem
#            )
#
#            mapping = Maps.IdentityMap(mesh)
#            surveyTD = TDEM.Survey([src])
#            prbTD = TDEM.Problem3D_e(mesh, sigmaMap=mapping, Solver=Solver)
#            prbTD.timeSteps = [(2.5e-5, 10), (1e-4, 10), (5e-4, 15)]
#            prbTD.pair(surveyTD)
#
#            std = 0.00
#            surveyTD.makeSyntheticData(m, std)

            print("Solve Tile %i --- %s seconds ---" % (tID, time.time() - start_time))

#            dtemp =
            fid = open(work_dir + '\\FWR_TDEM\\INP_DIR\\Trx_Tile_' + str(tID) + '.loc', 'w')
#            fid.write('IGNORE -9.9999 \n\n')
            fid.write('N_TRX 1\n\n')
            fid.write('TRX_LOOP\n')
            np.savetxt(fid, np.r_[srcLoc[0], 13., 0, 0].reshape((1,6)), fmt='%e',delimiter=' ',newline='\n\n')
            fid.write('N_RECV 1\n')
#            fid.write('N_TIME ' + str(len(times)) + '\n\n')

            np.savetxt(fid, srcLoc[0].reshape((1,3)), fmt='%e',delimiter=' ',newline='\n\n')
#            for tt in range(len(times)):
#
#                np.savetxt(fid, np.c_[np.r_[rxLoc[0],times[tt]].reshape((1,4)),np.ones((1,18))*1.], fmt='%e',delimiter=' ')
#
            fid.close()

        elif dID == 'GTEM':

            # OPTION TO WRITE OUT
#            Mesh.TensorMesh.writeUBC(mesh,work_dir +'\\FWR_GTEM\\INP_DIR\\Tile_' + str(tID) + '.msh')
#            Mesh.TensorMesh.writeModelUBC(mesh,work_dir +'\\FWR_GTEM\\INP_DIR\\Tile_' + str(tID) + '.con',m)

#            halfSpaceProblemAnaDiff
            start_time = time.time()

            rxLoc = grid_survey([20], [500],
                                x0=(X[tID], Y[tID], 1.), topo=topo)

            # Create grid of points
#            spacing=(40,80)
#            width = (500,2000)
#            rxLoc = grid_survey(spacing,width,x0 = (X[tID],Y[tID], 60.), topo = topo)
#            times = np.logspace(-4, np.log10(2e-3), 10)
#            times = np.loadtxt(work_dir + '\\VTEM_times.dat')
#            print('min diffusion distance ', 1.28*np.sqrt(times.min()/(sig_half*mu_0)),
#                  'max diffusion distance ', 1.28*np.sqrt(times.max()/(sig_half*mu_0)))
#            rx = TDEM.Rx(srcLoc, times, 'dbzdt')
#            src = TDEM.Src.MagDipole(
#                [rx],
#                waveform=TDEM.Src.StepOffWaveform(),
#                loc=rxLoc  # same src location as FDEM problem
#            )
#
#            mapping = Maps.IdentityMap(mesh)
#            surveyTD = TDEM.Survey([src])
#            prbTD = TDEM.Problem3D_e(mesh, sigmaMap=mapping, Solver=Solver)
#            prbTD.timeSteps = [(2.5e-5, 10), (1e-4, 10), (5e-4, 15)]
#            prbTD.pair(surveyTD)
#
#            std = 0.00
#            surveyTD.makeSyntheticData(m, std)

            print("Solve Tile %i --- %s seconds ---" % (tID, time.time() - start_time))

#            dtemp =
            fid = open(work_dir + '\\FWR_GTEM\\INP_DIR\\Trx_Tile_' + str(tID) + '.loc', 'w')
#            fid.write('IGNORE -9.9999 \n\n')
            fid.write('N_TRX 1\n\n')
            fid.write('TRX_ORIG\n')
            fid.write('5\n')

            # Define points for square loop centered on grid
            xOrig = X[tID] + lenTile * np.r_[-0.25, 0.25]
            yOrig = Y[tID] + lenTile * np.r_[-0.25, 0.25]

            xOrig, yOrig = np.meshgrid(xOrig, yOrig)

            # Flip the last two nodes for clockwise ordering
            yOrig = np.c_[yOrig[:,0],yOrig[::-1,1]]

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

            fid.write('N_RECV ' + str(rxLoc.shape[0]) + '\n')
#            fid.write('N_TIME ' + str(len(times)) + '\n\n')

            np.savetxt(fid, rxLoc, fmt='%e',delimiter=' ',newline='\n')
#            for tt in range(len(times)):
#
#                np.savetxt(fid, np.c_[np.r_[rxLoc[0],times[tt]].reshape((1,4)),np.ones((1,18))*1.], fmt='%e',delimiter=' ')
#
            fid.close()
#                np.savetxt(fid, , fmt='%e',delimiter=' ',newline='\n')
    # Clear the projection
#    del P
    gc.collect()

#fid.close()


#%% Write out GOcad polyline
fid = open(work_dir + '\\FWR_GTEM\\GroundLoops.pl','w')
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
## CREATE TEMPORARY MESH AND MODEL FOR DCIP3D
#hx = np.r_[padb, core, padf]
#hy = np.r_[padb, np.r_[np.ones(25)*dx[1], np.ones(4*nCx)*dx[0],np.ones(25)*dx[1]], padf]
#hz = np.r_[np.ones(25)*50,np.ones(20)*25,np.ones(10)*10, np.ones(30)*5,5*expf**(np.asarray(range(npad)))]
#
#mtemp = Mesh.TensorMesh([hx,hy,hz], 'CC0')
#
##    mtemp._x0 = [x0[ii]-np.sum(padb), y0[ii]-np.sum(padb), mesh.x0[2]]
#
#mtemp._x0 = (meshes[0].x0[0]-30,meshes[0].x0[1]-10,meshes[0].x0[2])
#
#actv2 = Utils.surface2ind_topo(mtemp, topo, 'N')
#    # actv2 = np.asarray([inds for inds, elem in enumerate(actv2, 1)
#    #                   if elem], dtype=int) - 1
#
#P = Maps.Mesh2MeshTopo([mesh,mtemp],[actv,actv2],nIterpPts = 12)
#model_DCIP3D = P*(model[dID][actv])
#m = np.ones(mtemp.nC)*1e-8
#m[actv2] = model_DCIP3D
#
#Mesh.TensorMesh.writeUBC(mtemp,work_dir +'\\MeshTile_TestDCIP3D.msh')
#Mesh.TensorMesh.writeModelUBC(mtemp,work_dir +'\\MeshTile_TestDCIP3D.sus',m)

##%% Create tiles
## Begin tiling
#max_mcell = 2e+5;
#min_Olap = 0e+2;
#
## In the x-direction
#ntile = 1
#Olap  = -1
#
## Cell size

#
#while Olap < min_Olap:
## for ii in range(5):
#
#    ntile += 1
#
#    # Set location of SW corners
#    x0 = np.asarray([mesh.vectorNx[0],mesh.vectorNx[-1] - nCx*dx[0]])
#
#    dx_t = np.round( ( x0[1] - x0[0] ) / ( (ntile-1) * dx[0]) )
#
#    if ntile>2:
#        x1 = np.r_[x0[0],
#                   x0[0] + np.cumsum(np.ones(ntile-2) * dx_t * dx[0]),
#                   x0[1]]
#    else:
#        x1 = np.asarray([x0[0],x0[1]])
#
#
#    x2 = x1 + nCx*dx[0];
#
##     y1 = np.ones(x1.shape[0])*np.min(locs[:,1]);
##     y2 = np.ones(x1.shape[0])*(np.min(locs[:,1]) + nCx*dx);
#
#    Olap = x1[0] + nCx*dx[0] - x1[1];
#
#
## Save x-corner location
#xtile = np.c_[x1,x2]
#
#
#
## In the Y-direction
#ntile = 1
#Olap  = -1
#
## Cell size
#
#while Olap < min_Olap:
## for ii in range(5):
#
#    ntile += 1
#
#    # Set location of SW corners
#    y0 = np.asarray([mesh.vectorNy[0],mesh.vectorNy[-1] - nCx*dx[0]])
#
#    dy_t = np.round( ( y0[1] - y0[0] ) / ( (ntile-1) * dx[0]) )
#
#    if ntile>2:
#        y1 = np.r_[y0[0],
#                   y0[0] + np.cumsum(np.ones(ntile-2) * dy_t * dx[0]),
#                   y0[1]]
#    else:
#        y1 = np.asarray([y0[0],y0[1]])
#
#
#    y2 = y1 + nCx*dx[0];
#
##     x1 = np.ones(y1.shape[0])*np.min(locs[:,0]);
##     x2 = np.ones(y1.shape[0])*(np.min(locs[:,0]) + nCx*dx);
#
#    Olap = y1[0] + nCx*dx[0] - y1[1];
#
#
## Save x-corner location
#ytile = np.c_[y1,y2]
#
#
#X,Y = np.meshgrid(x1,y1)
#
#x0 = mkvc(X)
#y0 = mkvc(Y)
