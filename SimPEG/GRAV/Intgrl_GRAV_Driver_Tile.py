#%%
from SimPEG import Mesh, Directives, Maps, InvProblem, Optimization, Utils
from SimPEG import DataMisfit, Inversion, Regularization, mkvc
import SimPEG.PF as PF
import pylab as plt
import os
import numpy as np
from matplotlib.patches import Rectangle
from SimPEG.ObjectiveFunction import ComboObjectiveFunction
work_dir = '.\\'
inpFile = 'Gabbs6_Unconstr_GravInv_input.inp'
out_dir = "Gabbs6_output\\"
dsep = os.path.sep

padLen = 2000
maxRAM = 0.5
n_cpu = 8

octreeObs = [5, 5, 3, 3, 3]  # Octree levels below observation points
octreeTopo = [0, 1]

meshType = 'TENSOR'
tileProblem = False
os.system('mkdir ' + work_dir+out_dir)


#%% User input
# Plotting parameter
#%%
# Read input file
#[mshfile, obsfile, topofile, mstart, mref, wgtfile, chi, alphas, bounds, lpnorms] = PF.Gravity.read_GRAVinv_inp(work_dir + dsep + inpFile)
driver = PF.GravityDriver.GravityDriver_Inv(work_dir + dsep + inpFile)
meshInput = driver.mesh
survey = driver.survey
topo = driver.topo
rxLoc = survey.srcField.rxList[0].locs

h = np.r_[meshInput.hx.min(), meshInput.hy.min(), meshInput.hz.min()]

# LOOP THROUGH TILES
surveyMask = np.ones(survey.nD, dtype='bool')
# Going through all problems:
# 1- Pair the survey and problem
# 2- Add up sensitivity weights
# 3- Add to the ComboMisfit

# Create first mesh outside the parallel process
padDist = np.r_[np.c_[padLen, padLen], np.c_[padLen, padLen], np.c_[padLen, 0]]

if meshType != meshInput._meshType:
    print("Creating Global Octree")
    mesh = Utils.modelutils.meshBuilder(
            rxLoc, h, padDist, meshType='TREE', meshGlobal=meshInput,
            verticalAlignment='center'
        )

    if topo is not None:
        mesh = Utils.modelutils.refineTree(
            mesh, topo, dtype='surface',
            nCpad=octreeTopo, finalize=False
        )

    mesh = Utils.modelutils.refineTree(
        mesh, rxLoc, dtype='surface',
        nCpad=octreeObs, finalize=True
    )

    if topo is not None:
        actv = Utils.surface2ind_topo(mesh, topo)
    else:
        actv = np.zeros(mesh.nC, dtype='bool')
        print(meshInput.vectorNz[-1])
        actv[mesh.gridCC[:, 2] < meshInput.vectorNz[-1]] = True

    m0 = np.ones(actv.sum()) * 1e-4

    print("Writing global Octree to file" + work_dir + out_dir + 'OctreeMeshGlobal.msh')
    Mesh.TreeMesh.writeUBC(
          mesh, work_dir + out_dir + 'OctreeMeshGlobal.msh',
          models={work_dir + out_dir + 'ActiveGlobal.act': actv}
        )


else:
    mesh = meshInput
    actv = np.zeros(mesh.nC, dtype='bool')
    actv[driver.activeCells] = True
    actvMap = Maps.InjectActiveCells(mesh, actv, 0)

    m0 = actvMap*driver.m0  # Starting model

wrGlobal = np.zeros(int(actv.sum()))
if tileProblem:


    # Loop over different tile size and break problem until memory usage is preserved
    usedRAM = np.inf
    count = 0
    while usedRAM > maxRAM:
        print("Tiling:" + str(count))

        tiles = Utils.modelutils.tileSurveyPoints(rxLoc, count)

        X1, Y1 = tiles[0][:, 0], tiles[0][:, 1]
        X2, Y2 = tiles[1][:, 0], tiles[1][:, 1]

        ind_t = np.all([rxLoc[:, 0] >= tiles[0][0, 0], rxLoc[:, 0] <= tiles[1][0, 0],
                        rxLoc[:, 1] >= tiles[0][0, 1], rxLoc[:, 1] <= tiles[1][0, 1],
                        surveyMask], axis=0)

        meshLocal = Utils.modelutils.meshBuilder(
            rxLoc, h, padDist, meshType='TREE', meshGlobal=meshInput,
            verticalAlignment='center'
        )

        if topo is not None:
            meshLocal = Utils.modelutils.refineTree(
                meshLocal, topo, dtype='surface',
                nCpad=octreeTopo, finalize=False
            )

        meshLocal = Utils.modelutils.refineTree(
            meshLocal, rxLoc[ind_t, :], dtype='surface',
            nCpad=octreeObs, finalize=True
        )

        nD, nC = ind_t.sum()*1., meshLocal.nC*1.

        nChunks = n_cpu # Number of chunks
        cSa, cSb = int(nD/nChunks), int(nC/nChunks) # Chunk sizes
        usedRAM = nD * nC * 8. * 1e-9
        count
        print(nD, nC, usedRAM)

    # Plot data and tiles
    fig, ax1 = plt.figure(), plt.subplot()
    Utils.PlotUtils.plot2Ddata(rxLoc, survey.dobs, ax=ax1)
    for ii in range(X1.shape[0]):
        ax1.add_patch(Rectangle((X1[ii], Y1[ii]),
                                X2[ii]-X1[ii],
                                Y2[ii]-Y1[ii],
                                facecolor='none', edgecolor='k'))
    ax1.set_xlim([X1.min()-20, X2.max()+20])
    ax1.set_ylim([Y1.min()-20, Y2.max()+20])
    ax1.set_aspect('equal')
    plt.show()

    def createLocalProb(rxLoc, wrGlobal, lims, ind):

        # Grab the data for current tile
        ind_t = np.all([rxLoc[:, 0] >= lims[0], rxLoc[:, 0] <= lims[1],
                        rxLoc[:, 1] >= lims[2], rxLoc[:, 1] <= lims[3],
                        surveyMask], axis=0)

        # Remember selected data in case of tile overlap
        surveyMask[ind_t] = False

        # Create new survey
        rxLoc_t = PF.BaseGrav.RxObs(rxLoc[ind_t, :])
        srcField = PF.BaseGrav.SrcField([rxLoc_t])
        survey_t = PF.BaseGrav.LinearSurvey(srcField)
        survey_t.dobs = survey.dobs[ind_t]
        survey_t.std = survey.std[ind_t]
        survey_t.ind = ind_t

        meshLocal = Utils.modelutils.meshBuilder(
            rxLoc, h, padDist, meshType='TREE', meshGlobal=meshInput,
            verticalAlignment='center'
        )

        if topo is not None:
            meshLocal = Utils.modelutils.refineTree(
                meshLocal, topo, dtype='surface',
                nCpad=octreeTopo, finalize=False
            )

        # Refine the mesh around loc
        meshLocal = Utils.modelutils.refineTree(
            meshLocal, rxLoc[ind_t, :], dtype='surface',
            nCpad=octreeObs, finalize=True
        )

        actv_t = np.ones(meshLocal.nC, dtype='bool')

        Mesh.TreeMesh.writeUBC(
              meshLocal, work_dir + out_dir + 'OctreeMesh' + str(tt) + '.msh',
              models={work_dir + out_dir + 'Active' + str(tt) + '.act': actv_t}
            )

        print(meshLocal.nC)
        # Create reduced identity map
        tileMap = Maps.Tile((mesh, actv), (meshLocal, actv_t))

        # Create the forward model operator
        prob =  PF.Gravity.GravityIntegral(
            meshLocal, rhoMap=tileMap, actInd=actv_t, parallelized=True,
            Jpath=work_dir + out_dir + "Tile" + str(ind) + ".zarr")

        survey_t.pair(prob)

        # Data misfit function
        dmis = DataMisfit.l2_DataMisfit(survey_t)
        dmis.W = 1./survey_t.std

        wrGlobal += prob.getJtJdiag(np.ones(tileMap.P.shape[1]), W=dmis.W)

        del meshLocal

        # Create combo misfit function
        return dmis, wrGlobal


    for tt in range(X1.shape[0]):

        print("Tile " + str(tt+1) + " of " + str(X1.shape[0]))

        dmis, wrGlobal = createLocalProb(rxLoc, wrGlobal, np.r_[X1[tt], X2[tt], Y1[tt], Y2[tt]], tt)

        if tt == 0:
            ComboMisfit = dmis

        else:
            ComboMisfit += dmis
else:

    # Create the forward model operator
    ## Create identity map
    nC = int(actv.sum())
    idenMap = Maps.IdentityMap(nP=nC)

    prob =  PF.Gravity.GravityIntegral(
        mesh, rhoMap=idenMap, actInd=actv, parallelized=True,
        Jpath=work_dir + out_dir + "Sensitivity.zarr")

    survey.pair(prob)

    # Data misfit function
    ComboMisfit = DataMisfit.l2_DataMisfit(survey)
    ComboMisfit.W = 1./survey.std

    wrGlobal += prob.getJtJdiag(np.ones(nC), W=ComboMisfit.W)
    actvGlobal = actv
# print('Sum of all problems:' + str(probSize*1e-6) + ' Mb')
# Scale global weights for regularization

# Check if global mesh has regions untouched by local problem
actvGlobal = wrGlobal != 0
if actvGlobal.sum() < actv.sum():

    if isinstance(ComboMisfit, ComboObjectiveFunction):
        for ind, dmis in enumerate(ComboMisfit.objfcts):
            dmis.prob.rhoMap.index = actvGlobal
            dmis.prob.rhoMap._P = None
            dmis.prob._model = None
            dmis.prob.gtgdiag = None
    else:
        ComboMisfit.prob.rhoMap.index = actvGlobal
        ComboMisfit.prob.rhoMap._P = None
        ComboMisfit.prob.model = np.zeros(actvGlobal.sum())
        ComboMisfit.prob.gtgdiag = None

wrGlobal = wrGlobal[actvGlobal]**0.5
wrGlobal = (wrGlobal/np.max(wrGlobal))

#%% Create a regularization
actvMap = Maps.InjectActiveCells(mesh, actv, 0)

actv = np.all([actv, actvMap*actvGlobal], axis=0)
actvMap = Maps.InjectActiveCells(mesh, actv, -100)

m0 = m0[actv]
idenMap = Maps.IdentityMap(nP=int(np.sum(actv)))
reg = Regularization.Sparse(mesh, indActive=actv, mapping=idenMap)
reg.norms = np.c_[driver.lpnorms].T

if driver.eps is not None:
    reg.eps_p = driver.eps[0]
    reg.eps_q = driver.eps[1]

reg.cell_weights = wrGlobal
reg.mref = np.zeros(mesh.nC)[actv]

# Add directives to the inversion
opt = Optimization.ProjectedGNCG(maxIter=20, lower=-10, upper=10.,
                                 maxIterLS=20, maxIterCG=20, tolCG=1e-4)
invProb = InvProblem.BaseInvProblem(ComboMisfit, reg, opt)
betaest = Directives.BetaEstimate_ByEig(beta0_ratio=1e0)

# Here is where the norms are applied
# Use pick a treshold parameter empirically based on the distribution of
#  model parameters
IRLS = Directives.Update_IRLS(f_min_change=1e-3, minGNiter=1,
                              maxIRLSiter=10)

IRLS.target = driver.survey.nD
update_Jacobi = Directives.UpdatePreconditioner()
saveModel = Directives.SaveUBCModelEveryIteration(mapping=actvMap, fileName=work_dir + out_dir + dsep + "GRAV_Tile")
inv = Inversion.BaseInversion(invProb,
                              directiveList=[betaest, saveModel,
                                             IRLS, update_Jacobi])

# Run the inversion
mrec = inv.run(m0)


# Outputs
if mesh._meshType == 'TENSOR':
    Mesh.TensorMesh.writeUBC(mesh, work_dir + out_dir + "GRAV_Tile.msh")
    Mesh.TensorMesh.writeModelUBC(mesh, work_dir + out_dir + "GRAV_Tile_lp.sus",
                                  actvMap*invProb.model)
    Mesh.TensorMesh.writeModelUBC(mesh, work_dir + out_dir + "GRAV_Tile_l2.sus",
                                  actvMap*invProb.l2model)
else:

    if invProb.l2model is not None:
        Mesh.TreeMesh.writeUBC(mesh, work_dir + out_dir + 'GRAV_Tile.msh',
                               models={work_dir + out_dir + 'Model_l2.sus': actvMap*invProb.l2model})
    Mesh.TreeMesh.writeUBC(mesh, work_dir + out_dir + 'GRAV_Tile.msh',
                           models={work_dir + out_dir + 'Model_lp.sus': actvMap*mrec})


# Get predicted data for each tile and write full predicted to file
if getattr(ComboMisfit, 'objfcts', None) is not None:
    dpred = np.zeros(survey.nD)
    for ind, dmis in enumerate(ComboMisfit.objfcts):
        dpred[dmis.survey.ind] += dmis.survey.dpred(mrec)
else:
    dpred = ComboMisfit.survey.dpred(mrec)
    # PF.Magnetics.writeUBCobs(work_dir+out_dir + "Tile" + str(ind) + ".pre",
    #                          survey, survey.dpred(mrec))

Utils.io_utils.writeUBCgravityObservations(
    work_dir+out_dir + "Predicted_lp.pre", survey, dpred
    )
#Utils.io_utils.writeUBCgravityObservations(
#    work_dir + out_dir + 'Predicted_l2.pre', survey,
#    survey.dpred(invProb.l2model)
#    )

