"""

This script runs an Magnetic Amplitude Inversion (MAI) from TMI data.
Magnetic amplitude data are weakly sensitive to the orientation of
magnetization, and can therefore better recover the location and geometry of
magnetic bodies in the presence of remanence. The algorithm is inspired from
Li & Shearer (2008), with an added iterative sensitivity weighting strategy to
counter the vertical streatchin that the old code had

This is done in three parts:

1- TMI data are inverted for an equivalent source layer.

2-The equivalent source layer is used to predict component data -> amplitude

3- Amplitude data are inverted in 3-D for an effective susceptibility model

Created on December 7th, 2016

@author: fourndo@gmail.com

"""
from SimPEG import Mesh, Directives, Maps, InvProblem, Optimization, DataMisfit, Inversion, Utils, Regularization
from SimPEG.Utils import mkvc
import SimPEG.PF as PF
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import os
import gc

#work_dir = "C:\\Users\\DominiqueFournier\\ownCloud\\Research\\Kevitsa\\Modeling\\MAG\\"
#work_dir = "C:\\Users\\DominiqueFournier\\ownCloud\\\Research\\Modelling\\Synthetic\\Block_Gaussian_topo\\"
# work_dir = "C:\\Users\\DominiqueFournier\\ownCloud\\Research\\Modelling\\Synthetic\\Nut_Cracker\\"
#work_dir = "C:\\Egnyte\\Private\\dominiquef\\Projects\\4559_CuMtn_ZTEM\\Modeling\\MAG\\A1_Fenton\\"
#work_dir = 'C:\\Users\\DominiqueFournier\\ownCloud\\Research\\Kevitsa\\Modeling\\MAG\\Aiborne\\'
work_dir = '/tera_raid/dfournier/Kevitsa/MAG/Aiborne/'
out_dir = "SimPEG_PF_Inv_MVI_S/"
tile_dirl2 = 'Tiles_l2/'
tile_dirlp = 'Tiles_lp/'
input_file = "SimPEG_MAG.inp"

max_mcell = 2.0e+5
min_Olap = 5e+2

# %%
# Read in the input file which included all parameters at once
# (mesh, topo, model, survey, inv param, etc.)
driver = PF.MagneticsDriver.MagneticsDriver_Inv(work_dir + input_file)

os.system('if not exist ' + work_dir + out_dir + ' mkdir ' + work_dir+out_dir)
os.system('if not exist ' + work_dir + out_dir + tile_dirlp + ' mkdir ' + work_dir+out_dir+ tile_dirlp)
os.system('if not exist ' + work_dir + out_dir +  tile_dirl2 +' mkdir ' + work_dir+out_dir+ tile_dirl2)

# Access the mesh and survey information
mesh = driver.mesh
survey = driver.survey
actv = driver.activeCells

# # TILE THE PROBLEM
lx = int(np.floor(np.sqrt(max_mcell/mesh.nCz)))
##%% Create tiles
## Begin tiling

# In the x-direction
ntile = 1
Olap  = -1
dx = [mesh.hx.min(), mesh.hy.min()]
# Cell size


while Olap < min_Olap:
# for ii in range(5):

    ntile += 1

    # Set location of SW corners
    x0 = np.asarray([mesh.vectorNx[0], mesh.vectorNx[-1] - lx*dx[0]])

    dx_t = np.round( ( x0[1] - x0[0] ) / ( (ntile-1) * dx[0]) )

    if ntile>2:
        x1 = np.r_[x0[0],
                   x0[0] + np.cumsum(np.ones(ntile-2) * dx_t * dx[0]),
                   x0[1]]
    else:
        x1 = np.asarray([x0[0],x0[1]])


    x2 = x1 + lx*dx[0];

#     y1 = np.ones(x1.shape[0])*np.min(locs[:,1]);
#     y2 = np.ones(x1.shape[0])*(np.min(locs[:,1]) + nCx*dx);

    Olap = x1[0] + lx*dx[0] - x1[1];


# Save x-corner location
xtile = np.c_[x1,x2]



# In the Y-direction
ntile = 1
Olap  = -1

# Cell size

while Olap < min_Olap:
# for ii in range(5):

    ntile += 1

    # Set location of SW corners
    y0 = np.asarray([mesh.vectorNy[0],mesh.vectorNy[-1] - lx*dx[0]])

    dy_t = np.round( ( y0[1] - y0[0] ) / ( (ntile-1) * dx[0]) )

    if ntile>2:
        y1 = np.r_[y0[0],
                   y0[0] + np.cumsum(np.ones(ntile-2) * dy_t * dx[0]),
                   y0[1]]
    else:
        y1 = np.asarray([y0[0],y0[1]])


    y2 = y1 + lx*dx[0];

#     x1 = np.ones(y1.shape[0])*np.min(locs[:,0]);
#     x2 = np.ones(y1.shape[0])*(np.min(locs[:,0]) + nCx*dx);

    Olap = y1[0] + lx*dx[0] - y1[1];


# Save x-corner location
ytile = np.c_[y1,y2]


X1,Y1 = np.meshgrid(x1,y1)
X1,Y1 = mkvc(X1), mkvc(Y1)
X2,Y2 = np.meshgrid(x2,y2)
X2,Y2 = mkvc(X2), mkvc(Y2)

# Plot data and tiles
rxLoc = survey.srcField.rxList[0].locs
# fig, ax1 = plt.figure(), plt.subplot()
# PF.Magnetics.plot_obs_2D(rxLoc,survey.dobs, ax=ax1)
# for ii in range(X1.shape[0]):
#     ax1.add_patch(Rectangle((X1[ii],Y1[ii]),X2[ii]-X1[ii],Y2[ii]-Y1[ii], facecolor = 'none', edgecolor='k'))
# ax1.set_aspect('equal')

# LOOP THROUGH TILES
npadxy = 8
core_x = np.ones(lx)*dx[0]
core_y = np.ones(lx)*dx[1]
expf = 1.3

os.system('if not exist ' + work_dir + out_dir + tile_dirl2 + ' mkdir ' + work_dir+out_dir+tile_dirl2)
os.system('if not exist ' + work_dir + out_dir + tile_dirlp + ' mkdir ' + work_dir+out_dir+tile_dirlp)

tiles = range(X1.shape[0])
tilestep = 14
for tt in tiles[:tilestep]:

    midX = np.mean([X1[tt], X2[tt]])
    midY = np.mean([Y1[tt], Y2[tt]])

    # Create new mesh
    padx = np.r_[dx[0]*expf**(np.asarray(range(npadxy))+1)]
    pady = np.r_[dx[1]*expf**(np.asarray(range(npadxy))+1)]

    hx = np.r_[padx[::-1], core_x, padx]
    hy = np.r_[pady[::-1], core_y, pady]
#        hz = np.r_[padb*2.,[33,26],np.ones(25)*22,[18,15,12,10,8,7,6], np.ones(18)*5,5*expf**(np.asarray(range(2*npad)))]
    hz = mesh.hz

    mesh_t = Mesh.TensorMesh([hx,hy,hz], 'CC0')

#    mtemp._x0 = [x0[ii]-np.sum(padb), y0[ii]-np.sum(padb), mesh.x0[2]]

    mesh_t._x0 = (mesh_t.x0[0] + midX, mesh_t.x0[1]+midY, mesh.x0[2])

    Mesh.TensorMesh.writeUBC(mesh_t,work_dir + out_dir + tile_dirl2 + "MAG_Tile"+ str(tt) +".msh")
    Mesh.TensorMesh.writeUBC(mesh_t,work_dir + out_dir + tile_dirlp + "MAG_Tile"+ str(tt) +".msh")

#        meshes.append(mtemp)
    # Grab the right data
    xlim = [mesh_t.vectorCCx[npadxy], mesh_t.vectorCCx[-npadxy]]
    ylim = [mesh_t.vectorCCy[npadxy], mesh_t.vectorCCy[-npadxy]]

    ind_t = np.all([rxLoc[:,0] > xlim[0], rxLoc[:,0] < xlim[1],
                    rxLoc[:,1] > ylim[0], rxLoc[:,1] < ylim[1]], axis=0)


    if np.sum(ind_t) < 20:
        continue

    rxLoc_t = PF.BaseMag.RxObs(rxLoc[ind_t,:])
    srcField = PF.BaseMag.SrcField([rxLoc_t], param=survey.srcField.param)
    survey_t = PF.BaseMag.LinearSurvey(srcField)
    survey_t.dobs = survey.dobs[ind_t]
    survey_t.std = survey.std[ind_t]

    # Extract model from global to local mesh
    if driver.topofile is not None:
        topo = np.genfromtxt(work_dir + driver.topofile,
                             skip_header=1)
        actv = Utils.surface2ind_topo(mesh_t, topo, 'N')
        actv = np.asarray(np.where(mkvc(actv))[0], dtype=int)
    else:
        actv = np.ones(mesh_t.nC, dtype='bool')

    surf = PF.MagneticsDriver.actIndFull2layer(mesh_t, actv)

    # Get the layer of cells directyl below topo
    #surf = Utils.actIndFull2layer(mesh_t, active)
    nC = len(surf)  # Number of active cells

    # Create active map to go from reduce set to full
    surfMap = Maps.InjectActiveCells(mesh_t, surf, -100)

    # Create identity map
    idenMap = Maps.IdentityMap(nP=nC)

    # Create static map
    prob = PF.Magnetics.MagneticIntegral(mesh_t, chiMap=idenMap, actInd=surf, equiSourceLayer=True)
    prob.solverOpts['accuracyTol'] = 1e-4

    # Pair the survey and problem
    survey.pair(prob)

    # Create a regularization function, in this case l2l2
    reg = Regularization.Simple(mesh_t, indActive=surf)
    reg.mref = np.zeros(nC)

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
    targetMisfit = Directives.TargetMisfit(chifact=0.1)

    # Put all the parts together
    inv = Inversion.BaseInversion(invProb,
                                  directiveList=[betaest, betaSchedule, targetMisfit])

    # Run the equivalent source inversion
    mstart = np.zeros(nC)
    mrec = inv.run(mstart)

    Mesh.TensorMesh.writeModelUBC(mesh_t, work_dir + out_dir + tile_dirl2 + "EQS_Tile" + str(tt) + ".sus", surfMap*mrec)

    # %% STEP 2: COMPUTE AMPLITUDE DATA
    # Now that we have an equialent source layer, we can forward model alh three
    # components of the field and add them up: |B| = ( Bx**2 + Bx**2 + Bx**2 )**0.5

    # Won't store the sensitivity and output 'xyz' data.
    prob.forwardOnly = True
    pred_x = prob.Intrgl_Fwr_Op(m=mrec, recType='x')
    pred_y = prob.Intrgl_Fwr_Op(m=mrec, recType='y')
    pred_z = prob.Intrgl_Fwr_Op(m=mrec, recType='z')

    ndata = survey.nD

    d_amp = np.sqrt(pred_x**2. +
                    pred_y**2. +
                    pred_z**2.)

    rxLoc = survey.srcField.rxList[0].locs
    # PF.Magnetics.plot_obs_2D(rxLoc,survey.dobs,varstr='TMI Data')
    # PF.Magnetics.plot_obs_2D(rxLoc,damp,varstr='Amplitude Data')

    # Write data out
    PF.Magnetics.writeUBCobs(work_dir + out_dir + 'Amplitude_data.obs', survey, d_amp)

    # %% STEP 3: RUN AMPLITUDE INVERSION
    # Now that we have |B| data, we can invert. This is a non-linear inversion,
    # which requires some special care for the sensitivity weights (see Directives)

    # Re-set the active cells to entire mesh
    # Create active map to go from reduce space to full
    actvMap = Maps.InjectActiveCells(mesh_t, active, -100)
    nC = len(active)

    # Create identity map
    idenMap = Maps.IdentityMap(nP=nC)

    mstart = np.ones(len(active))*1e-4

    # Create the forward model operator
    prob = PF.Magnetics.MagneticAmplitude(mesh_t, chiMap=idenMap,
                                          actInd=active)
    prob.chi = mstart


    # Pair the survey and problem
    survey.pair(prob)

    # Re-set the observations to |B|
    survey.dobs = d_amp

    # Create a sparse regularization
    reg = Regularization.Sparse(mesh_t, indActive=active, mapping=idenMap)
    reg.mref = driver.mref
    reg.norms = driver.lpnorms
    if driver.eps is not None:
        reg.eps_p = driver.eps[0]
        reg.eps_q = driver.eps[1]

    # Data misfit function
    dmis = DataMisfit.l2_DataMisfit(survey)
    dmis.W = 1/survey.std

    # Add directives to the inversion
    opt = Optimization.ProjectedGNCG(maxIter=100, lower=0., upper=1.,
                                     maxIterLS=20, maxIterCG=10,
                                     tolCG=1e-3)

    invProb = InvProblem.BaseInvProblem(dmis, reg, opt)

    # Here is the list of directives
    betaest = Directives.BetaEstimate_ByEig()

    # Specify the sparse norms
    IRLS = Directives.Update_IRLS(f_min_change=1e-3,
                                  minGNiter=3, coolingRate=1, chifact=0.25,
                                  maxIRLSiter=1)

    # Special directive specific to the mag amplitude problem. The sensitivity
    # weights are update between each iteration.
    update_Jacobi = Directives.Amplitude_Inv_Iter()

    saveModel = Directives.SaveUBCModelEveryIteration(mapping=actvMap)
    saveModel.fileName = work_dir + out_dir + 'ModelAmp'

    # Put all together
    inv = Inversion.BaseInversion(invProb,
                                  directiveList=[IRLS, update_Jacobi, betaest,
                                                 saveModel])

    # Invert
    mrec = inv.run(mstart)

    Mesh.TensorMesh.writeModelUBC(mesh_t, work_dir + out_dir + tile_dirl2 + "AMP_Tile" + str(tt) + ".sus", actvMap*reg.l2model)
    Mesh.TensorMesh.writeModelUBC(mesh_t, work_dir + out_dir + tile_dirlp + "AMP_Tile" + str(tt) + ".sus", actvMap*mrec)

    # Outputs
    PF.Magnetics.writeUBCobs(work_dir+out_dir + tile_dirlp + "AMP_Tile"+ str(tt) +"_Inv.pre", survey_t, invProb.dpred)


    del prob
    gc.collect()

#PF.Magnetics.plot_obs_2D(rxLoc,invProb.dpred,varstr='Amplitude Data')
