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
from SimPEG import Mesh, Directives, Maps, InvProblem, Optimization
from SimPEG import DataMisfit, Inversion, Utils, Regularization
from SimPEG.Utils import mkvc
import SimPEG.PF as PF
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.spatial import cKDTree
from discretize.utils import closestPoints
import os
import multiprocessing

if __name__ == '__main__':

    #work_dir = "C:\\Users\\DominiqueFournier\\ownCloud\\Research\\Kevitsa\\Modeling\\MAG\\"
    #work_dir = "C:\\Users\\DominiqueFournier\\ownCloud\\\Research\\Synthetic\\Block_Gaussian_topo\\"
    #work_dir = "C:\\Users\\DominiqueFournier\\ownCloud\\Research\\Synthetic\\Nut_Cracker\\"
    #work_dir = 'C:\\Users\\DominiqueFournier\\ownCloud\\Research\\Kevitsa\\Modeling\\MAG\\Airborne\\'
    # work_dir = "C:\\Users\\DominiqueFournier\\ownCloud\\Research\\CraigModel\\MAG\\"
    #work_dir = "C:\\Users\\DominiqueFournier\\ownCloud\\Research\\Yukon\\Modeling\\MAG\\"
    # work_dir = 'C:\\Users\\DominiqueFournier\\ownCloud\\Research\\Kevitsa\\Modeling\\MAG\\Airborne\\'
    #work_dir = "C:\\Users\\DominiqueFournier\\ownCloud\\Research\\Synthetic\\Triple_Block_lined\\"
    work_dir = "C:\\Users\\DominiqueFournier\\ownCloud\\Research\\Synthetic\\Nut_Cracker\\"
    out_dir = "SimPEG_MVI_S_TileInv\\"
    input_file = "SimPEG_MAG.inp"
    meshType = 'TreeMesh'
    padLen = 3000
    dwnFact = .25
    maxNpoints = 100
    numProcessors = 8

    # %%
    # Read in the input file which included all parameters at once
    # (mesh, topo, model, survey, inv param, etc.)
    driver = PF.MagneticsDriver.MagneticsDriver_Inv(work_dir + input_file)

    os.system('if not exist ' + work_dir + out_dir + ' mkdir ' + work_dir+out_dir)

    # Access the mesh and survey information
    meshInput = driver.mesh
    survey = driver.survey
    xyzLocs = survey.srcField.rxList[0].locs.copy()

    topo = None
    if driver.topofile is not None:
        topo = np.genfromtxt(driver.basePath + driver.topofile,
                             skip_header=1)
    else:
        # Grab the top coordinate and make a flat topo
        indTop = meshInput.gridCC[:, 2] == meshInput.vectorCCz[-1]
        topo = meshInput.gridCC[indTop, :]
        topo[:, 2] += meshInput.hz.min()/2. + 1e-8

    if meshType == 'TreeMesh':
        if isinstance(meshInput, Mesh.TensorMesh):
            # Define an octree mesh based on the provided tensor
            h = np.r_[meshInput.hx.min(), meshInput.hy.min(), meshInput.hz.min()]
            coreX, coreY, coreZ = meshInput.hx == h[0], meshInput.hy == h[1], meshInput.hz == h[2]
            padx, pady, padz = meshInput.hx[~coreX].sum(), meshInput.hy[~coreY].sum(), meshInput.hz[~coreZ].sum()

            padDist = np.r_[np.c_[padx, padx], np.c_[pady, pady], np.c_[padz, padz]]

            print("Creating TreeMesh. Please standby...")
            mesh = Utils.modelutils.meshBuilder(topo, h, padDist,
                                                meshGlobal=meshInput,
                                                meshType='TREE',
                                                gridLoc='CC')

            mesh = Utils.modelutils.refineTree(mesh, topo, dtype='surface',
                                               nCpad=[0, 10, 5], finalize=False)

            mesh = Utils.modelutils.refineTree(mesh, xyzLocs, dtype='surface',
                                               nCpad=[10, 0, 0], finalize=True)

        else:
            mesh = Mesh.TreeMesh.readUBC(driver.basePath + driver.mshfile)
    else:
        mesh = meshInput
    actv = Utils.surface2ind_topo(mesh, topo)

    if isinstance(mesh, Mesh.TreeMesh):
        Mesh.TreeMesh.writeUBC(mesh, work_dir + out_dir + 'OctreeMesh.msh',
                               models={work_dir + out_dir + 'ActiveOctree.dat': actv})
    else:
        mesh.writeModelUBC(mesh, work_dir + out_dir + 'ActiveOctree.dat', actv)

    actvMap = Maps.InjectActiveCells(mesh, actv, 0)

    tiles = Utils.modelutils.tileSurveyPoints(xyzLocs, maxNpoints)

    X1, Y1 = tiles[0][:, 0], tiles[0][:, 1]
    X2, Y2 = tiles[1][:, 0], tiles[1][:, 1]


    # Plot data and tiles
    fig, ax1 = plt.figure(), plt.subplot()
    PF.Magnetics.plot_obs_2D(xyzLocs, survey.dobs, ax=ax1)
    for ii in range(X1.shape[0]):
        ax1.add_patch(Rectangle((X1[ii], Y1[ii]),
                                X2[ii]-X1[ii],
                                Y2[ii]-Y1[ii],
                                facecolor='none', edgecolor='k'))
    ax1.set_xlim([X1.min()-20, X2.max()+20])
    ax1.set_ylim([Y1.min()-20, Y2.max()+20])
    ax1.set_aspect('equal')
    plt.show()

    # LOOP THROUGH TILES
    # expf = 1.3
    # dx = [mesh.hx.min(), mesh.hy.min()]
    surveyMask = np.ones(survey.nD, dtype='bool')
    # Going through all problems:
    # 1- Pair the survey and problem
    # 2- Add up sensitivity weights
    # 3- Add to the ComboMisfit

    nC = int(actv.sum())
    wrGlobal = np.zeros(3*nC)
    probSize = 0

    padDist = np.r_[np.c_[padLen, padLen], np.c_[padLen, padLen], np.c_[padLen, 0]]
    def createLocalProb(rxLoc, wrGlobal, lims):

        # Grab the data for current tile
        ind_t = np.all([rxLoc[:, 0] >= lims[0], rxLoc[:, 0] <= lims[1],
                        rxLoc[:, 1] >= lims[2], rxLoc[:, 1] <= lims[3],
                        surveyMask], axis=0)

        # Remember selected data in case of tile overlap
        surveyMask[ind_t] = False

        # Create new survey
        rxLoc_t = PF.BaseMag.RxObs(rxLoc[ind_t, :])
        srcField = PF.BaseMag.SrcField([rxLoc_t], param=survey.srcField.param)
        survey_t = PF.BaseMag.LinearSurvey(srcField)
        survey_t.dobs = survey.dobs[ind_t]
        survey_t.std = survey.std[ind_t]
        survey_t.ind = ind_t

        # mesh_t = meshTree.copy()
        mesh_t = Utils.modelutils.meshBuilder(rxLoc[ind_t, :], h, padDist,
                                            meshGlobal=meshInput,
                                            meshType='TREE',
                                            gridLoc='CC')

        mesh_t = Utils.modelutils.refineTree(mesh_t, topo, dtype='surface',
                                           nCpad=[0, 10, 2], finalize=False)

        mesh_t = Utils.modelutils.refineTree(mesh_t, rxLoc[ind_t, :], dtype='surface',
                                           nCpad=[10, 0, 0], finalize=False)

        center = np.mean(rxLoc[ind_t, :], axis=0)
        tileCenter = np.r_[np.mean(lims[0:2]), np.mean(lims[2:]), center[2]]

        ind = closestPoints(mesh, tileCenter, gridLoc='CC')

        shift = np.squeeze(mesh.gridCC[ind, :]) - center

        mesh_t.x0 += shift
        mesh_t.finalize()

        actv_t = Utils.surface2ind_topo(mesh_t, topo)

        # Create reduced identity map
        tileMap = Maps.Tile((mesh, actv), (mesh_t, actv_t))
        tileMap.nCell = 40
        tileMap.nBlock = 3

        # Create the forward model operator
        prob = PF.Magnetics.MagneticVector(mesh_t, chiMap=tileMap, actInd=actv_t,
                                           memory_saving_mode=True, parallelized=True)
        survey_t.pair(prob)

        # Data misfit function
        dmis = DataMisfit.l2_DataMisfit(survey_t)
        dmis.W = 1./survey_t.std

        wrGlobal += prob.getJtJdiag(np.ones(tileMap.P.shape[1]))

        # Create combo misfit function
        return dmis, wrGlobal

    for tt in range(X1.shape[0]):

        print("Tile " + str(tt+1) + " of " + str(X1.shape[0]))

        dmis, wrGlobal = createLocalProb(xyzLocs, wrGlobal, np.r_[X1[tt], X2[tt], Y1[tt], Y2[tt]])

        # Create combo misfit function

        if tt == 0:
            ComboMisfit = dmis

        else:
            ComboMisfit += dmis

        # Add problem size
    #    probSize += prob.F.shape[0] * prob.F.shape[1] * 32 / 4

    #ComboMisfit = ComboMisfit*1
    #print('Sum of all problems:' + str(probSize*1e-6) + ' Mb')
    # Scale global weights for regularization
    # Check if global mesh has regions untouched by local problem
    actvGlobal = wrGlobal != 0
    actvMeshGlobal = (wrGlobal[:nC]) != 0
    if actvMeshGlobal.sum() < actv.sum():

        for ind, dmis in enumerate(ComboMisfit.objfcts):
            dmis.prob.chiMap.index = actvMeshGlobal
            dmis.prob.gtgdiag = None

    wrGlobal = wrGlobal[actvGlobal]**0.5
    wrGlobal = (wrGlobal/np.max(wrGlobal))

    #%% Create a regularization
    actv = np.all([actv, actvMap*actvMeshGlobal], axis=0)
    actvMap = Maps.InjectActiveCells(mesh, actv, 0)
    actvMapAmp = Maps.InjectActiveCells(mesh, actv, -100)

    nC = int(np.sum(actv))

    mstart = np.ones(3*nC)*1e-4
    mref = np.zeros(3*nC)

    # Create a block diagonal regularization
    wires = Maps.Wires(('p', nC), ('s', nC), ('t', nC))

    # Create a regularization
    reg_p = Regularization.Sparse(mesh, indActive=actv, mapping=wires.p)
    reg_p.cell_weights = (wires.p * wrGlobal)
    reg_p.norms = np.c_[2, 2, 2, 2]
    reg_p.mref = mref

    reg_s = Regularization.Sparse(mesh, indActive=actv, mapping=wires.s)
    reg_s.cell_weights = (wires.s * wrGlobal)
    reg_s.norms = np.c_[2, 2, 2, 2]
    reg_s.mref = mref

    reg_t = Regularization.Sparse(mesh, indActive=actv, mapping=wires.t)
    reg_t.cell_weights = (wires.t * wrGlobal)
    reg_t.norms = np.c_[2, 2, 2, 2]
    reg_t.mref = mref

    reg = reg_p + reg_s + reg_t
    reg.mref = mref

    # Add directives to the inversion
    opt = Optimization.ProjectedGNCG(maxIter=7, lower=-10., upper=10.,
                                     maxIterCG=20, tolCG=1e-3)

    invProb = InvProblem.BaseInvProblem(ComboMisfit, reg, opt)
    betaest = Directives.BetaEstimate_ByEig()

    # Here is where the norms are applied
    IRLS = Directives.Update_IRLS(f_min_change=1e-3,
                                  minGNiter=1)

    update_Jacobi = Directives.UpdateJacobiPrecond()
    targetMisfit = Directives.TargetMisfit()

    saveModel = Directives.SaveUBCModelEveryIteration(mapping=actvMap)
    saveModel.fileName = work_dir + out_dir + 'MVI_C'
    inv = Inversion.BaseInversion(invProb,
                                  directiveList=[betaest, IRLS, update_Jacobi,
                                                 saveModel])

    mrec_MVI = inv.run(mstart)

    x = actvMap * (wires.p * mrec_MVI)
    y = actvMap * (wires.s * mrec_MVI)
    z = actvMap * (wires.t * mrec_MVI)

    amp = (np.sum(np.c_[x, y, z]**2., axis=1))**0.5

    if isinstance(mesh, Mesh.TreeMesh):
        Mesh.TreeMesh.writeUBC(
          mesh, work_dir + out_dir + 'OctreeMesh.msh',
          models={work_dir + out_dir + 'MVI_C_amp.sus': amp}
        )
    else:
        mesh.writeModelUBC(mesh, work_dir+out_dir + 'MVI_C_amp.sus', amp)

    # Get predicted data for each tile and write full predicted to file
    if getattr(ComboMisfit, 'objfcts', None) is not None:
        dpred = np.zeros(survey.nD)
        for ind, dmis in enumerate(ComboMisfit.objfcts):
            dpred[dmis.survey.ind] += dmis.survey.dpred(mrec_MVI)
    else:
        dpred = ComboMisfit.survey.dpred(mrec_MVI)

    PF.Magnetics.writeUBCobs(
      work_dir+out_dir + 'MVI_C_pred.pre', survey, dpred
    )

    beta = invProb.beta

    # %% RUN MVI-S WITH SPARSITY

    # # STEP 3: Finish inversion with spherical formulation
    mstart = Utils.matutils.xyz2atp(mrec_MVI.reshape((nC,3),order='F'))

    for misfit in ComboMisfit.objfcts:
        misfit.prob.coordinate_system = 'spherical'
        misfit.prob.model = mstart

    # Create a block diagonal regularization
    wires = Maps.Wires(('amp', nC), ('theta', nC), ('phi', nC))

    # Create a regularization
    reg_a = Regularization.Sparse(mesh, indActive=actv,
                                  mapping=wires.amp, gradientType='total')
    reg_a.norms = np.c_[driver.lpnorms[:4]].T
    reg_a.mref = mref
    if driver.eps is not None:
        reg_a.eps_p = driver.eps[0]
        reg_a.eps_q = driver.eps[1]

    reg_t = Regularization.Sparse(mesh, indActive=actv,
                                  mapping=wires.theta, gradientType='total')
    reg_t.alpha_s = 0.  # No reference angle
    reg_t.space = 'spherical'
    reg_t.norms = np.c_[driver.lpnorms[4:8]].T
    reg_t.eps_q = 5e-2
    reg_t.mref = mref
    # reg_t.alpha_x, reg_t.alpha_y, reg_t.alpha_z = 0.25, 0.25, 0.25

    reg_p = Regularization.Sparse(mesh, indActive=actv,
                                  mapping=wires.phi, gradientType='total')
    reg_p.alpha_s = 0.  # No reference angle
    reg_p.space = 'spherical'
    reg_p.norms = np.c_[driver.lpnorms[8:]].T
    reg_p.eps_q = 5e-2
    reg_p.mref = mref

    reg = reg_a + reg_t + reg_p
    reg.mref = mref

    Lbound = np.kron(np.asarray([0, -np.inf, -np.inf]), np.ones(nC))
    Ubound = np.kron(np.asarray([10, np.inf, np.inf]), np.ones(nC))


    # Add directives to the inversion
    opt = Optimization.ProjectedGNCG(maxIter=40,
                                     lower=Lbound,
                                     upper=Ubound,
                                     maxIterLS=10,
                                     maxIterCG=20, tolCG=1e-3,
                                     stepOffBoundsFact=1e-8)

    invProb = InvProblem.BaseInvProblem(ComboMisfit, reg, opt, beta=beta)
    #  betaest = Directives.BetaEstimate_ByEig()

    # Here is where the norms are applied
    IRLS = Directives.Update_IRLS(f_min_change=1e-4, maxIRLSiter=20,
                                  minGNiter=1, beta_tol = 0.5, prctile=100,
                                  coolingRate=1, coolEps_q=True,
                                  betaSearch=True)

    # Special directive specific to the mag amplitude problem. The sensitivity
    # weights are update between each iteration.
    ProjSpherical = Directives.ProjSpherical()
    update_SensWeight = Directives.UpdateSensWeighting()
    update_Jacobi = Directives.UpdateJacobiPrecond()
    saveModel = Directives.SaveUBCModelEveryIteration(mapping=actvMap)
    saveModel.fileName = work_dir+out_dir + 'MVI_S'

    inv = Inversion.BaseInversion(invProb,
                                  directiveList=[ProjSpherical, IRLS, update_SensWeight,
                                                 update_Jacobi, saveModel])

    mrec_MVI_S = inv.run(mstart)


    # Get predicted data for each tile and write full predicted to file
    if getattr(ComboMisfit, 'objfcts', None) is not None:
        dpred = np.zeros(survey.nD)
        for ind, dmis in enumerate(ComboMisfit.objfcts):
            dpred[dmis.survey.ind] += dmis.survey.dpred(mrec_MVI_S)
    else:
        dpred = ComboMisfit.survey.dpred(mrec_MVI_S)

    PF.Magnetics.writeUBCobs(work_dir+out_dir + 'MVI_S_pred.pre', survey, dpred)

    ##%%
    #JtJdiag = np.zeros_like(invProb.model)
    #for obj in ComboMisfit.objfcts:
    #    f = obj.prob.fields(mstart)
    #    JtJdiag += np.sum((obj.prob.getJ(mstart, f))**2., axis=0)
