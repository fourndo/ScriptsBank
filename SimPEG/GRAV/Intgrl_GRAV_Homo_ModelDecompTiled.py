#%%
from SimPEG import Mesh, Directives, Maps, InvProblem, Optimization, Utils
from SimPEG import DataMisfit, Inversion, Regularization
import SimPEG.PF as PF
import pylab as plt
import os
import numpy as np
from scipy.spatial import cKDTree
from matplotlib.patches import Rectangle

if __name__ == '__main__':
    work_dir = "C:\\Users\\DominiqueFournier\\Dropbox\\Projects\\Synthetic\\Block_Gaussian_topo\\GRAV\\"
    inpfile = 'SimPEG_GRAV.inp'
    out_dir = "SimPEG_GRAV_ModelDecomp\\"
    dsep = os.path.sep

    padLen = 200
    maxRAM = 0.05
    n_cpu = 4

    octreeObs = [5, 5, 3, 3, 3]  # Octree levels below observation points
    octreeTopo = [0, 1]
    ndv = -100
    meshType = 'TREE'
    tileProblem = True
    parallization = "dask" # "dask" ||  "multiprocessing"

    os.system('if not exist ' + work_dir + out_dir + ' mkdir ' + work_dir+out_dir)

    # Choice for the homogeneous model
    useMrefValues = True

    # Read input file
    driver = PF.GravityDriver.GravityDriver_Inv(work_dir + dsep + inpfile)
    meshInput = driver.mesh
    survey = driver.survey
    actv = driver.activeCells
    topo = driver.topo
    rxLoc = survey.srcField.rxList[0].locs

    # Tile the forward problem
    tree = cKDTree(meshInput.gridCC)

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

        if isinstance(driver.mstart, float):
            m0 = np.ones(mesh.nC) * driver.mstart

        else:
            print("Interpolating the starting model")

            _, ind = tree.query(mesh.gridCC)

            m0 = driver.m0
            m0[m0 == ndv] = 0
            m0 = m0[ind]

        if isinstance(driver._mrefInput, float):
            mref = np.ones(mesh.nC) * driver._mrefInput

        else:
            print("Interpolating the reference model")
            _, ind = tree.query(mesh.gridCC)

            mref = driver.mref
            mref[mref == ndv] = 0
            mref = mref[ind]

        print("Writing global Octree to file" + work_dir + out_dir + 'OctreeMeshGlobal.msh')
        Mesh.TreeMesh.writeUBC(
              mesh, work_dir + out_dir + 'OctreeMeshGlobal.msh',
              models={work_dir + out_dir + 'ActiveGlobal.act': actv}
            )

        # Create new inteprolation tree for the tiles
        tree = cKDTree(mesh.gridCC)

    else:
        mesh = meshInput
        actv = np.zeros(mesh.nC, dtype='bool')
        actv[driver.activeCells] = True
        actvMap = Maps.InjectActiveCells(mesh, actv, 0)

        m0 = driver.m0  # Starting model
        mref = driver.mref  # Starting model

    # Create active map to go from reduce set to full
    fullMap = Maps.InjectActiveCells(mesh, actv, 0)

    # Set the global homogeneous mapping
    m0 = m0[actv]
    mgeo = mref[actv]

    # Get unique geo units
    geoUnits = np.unique(mgeo).tolist()

    # Compute an a median value for each homogeneous units
    mUnit = np.asarray([np.median(mgeo[mgeo==unit]) for unit in geoUnits])



    # Build list of indecies for the geounits
    index = []
    for unit in geoUnits:
    #    if unit!=0:
        index += [mgeo==unit]

    nC = len(index)

    # Creat reduced identity map
    homogMap = Maps.SurjectUnits(index)

    # Create a wire map to link the homogeneous and heterogeneous spaces
    wires = Maps.Wires(('homo', nC), ('hetero', int(actv.sum())))

    # Tile the problem
    wrGlobal = np.zeros(int(nC + actv.sum()))
    if tileProblem:

        # Loop over different tile size and break
        # problem fit in memory
        usedRAM = np.inf
        count = 0
        while usedRAM > maxRAM:
            print("Tiling:" + str(count))

            tiles, binCount = Utils.modelutils.tileSurveyPoints(rxLoc, count)

            # Grab the smallest bin and generate a temporary mesh
            indMin = np.argmin(binCount)

            X1, Y1 = tiles[0][:, 0], tiles[0][:, 1]
            X2, Y2 = tiles[1][:, 0], tiles[1][:, 1]

            ind_t = np.all([rxLoc[:, 0] >= tiles[0][indMin, 0], rxLoc[:, 0] <= tiles[1][indMin, 0],
                            rxLoc[:, 1] >= tiles[0][indMin, 1], rxLoc[:, 1] <= tiles[1][indMin, 1],
                            surveyMask], axis=0)

            # Create the mesh and refine the same as the global mesh
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

            # Calculate approximate problem size
            nD_t, nC_t = ind_t.sum()*1., meshLocal.nC*1.

            nChunks = n_cpu # Number of chunks
            # cSa, cSb = int(nD_t/nChunks), int(nC_t/nChunks) # Chunk sizes
            usedRAM = nD_t * nC_t * 8. * 1e-9
            count += 1
            print(nD_t, nC_t, usedRAM)

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

        def createLocalProb(rxLoc, wrGlobal, lims, tileID):

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
                  models={work_dir + out_dir + 'Active' + str(tileID) + '.act': actv_t}
                )

            # Create reduced identity map
            tileMap = Maps.Tile((mesh, actv), (meshLocal, actv_t))

            actv_t = tileMap.activeLocal

            # Interpolate the geo model
            _, ind = tree.query(meshLocal.gridCC)

            mgeo_t = ((fullMap * mgeo)[ind])[actv_t]

            # # Get unique geo units
            # geoUnits = np.unique(mgeo_t).tolist()

            # Build list of indecies for the geounits
            index = []
            for unit in geoUnits:
                index += [mgeo_t == unit]

            # Creat reduced identity map
            homogMap_t = Maps.SurjectUnits(index)

            # Create Sum map
            sumMap = Maps.SumMap([homogMap_t*wires.homo, tileMap*wires.hetero])

            # Create the forward model operator
            prob = PF.Gravity.GravityIntegral(
                meshLocal, rhoMap=sumMap, actInd=actv_t, parallelized=parallization,
                Jpath=work_dir + out_dir + "Tile" + str(tileID) + ".zarr", n_cpu=n_cpu)

            survey_t.pair(prob)

            # Data misfit function
            dmis = DataMisfit.l2_DataMisfit(survey_t)
            dmis.W = 1./survey_t.std

            wrGlobal += prob.getJtJdiag(np.ones(wires.homo.shape[1]), W=dmis.W)

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

        sumMap = Maps.SumMap([homogMap*wires.homo, wires.hetero])

        prob = PF.Gravity.GravityIntegral(
            mesh, rhoMap=sumMap, actInd=actv, parallelized=parallization,
            Jpath=work_dir + out_dir + "Sensitivity.zarr", n_cpu=n_cpu)

        survey.pair(prob)

        # Data misfit function
        ComboMisfit = DataMisfit.l2_DataMisfit(survey)
        ComboMisfit.W = 1./survey.std

        wrGlobal += prob.getJtJdiag(np.ones(nC), W=ComboMisfit.W)
        actvGlobal = actv


    # Create augmented mstart and mref
    mref = np.r_[mUnit, mgeo]
    mstart = np.r_[mUnit*0., m0]

    # Create Sum map
    sumMap = Maps.SumMap([homogMap*wires.homo, wires.hetero])


    #%% Run inversion
    # prob = PF.Gravity.GravityIntegral(mesh, rhoMap=sumMap, actInd=actv)

    # survey.pair(prob)

    #%%
    # Load weighting  file

    # dmis = DataMisfit.l2_DataMisfit(survey)
    # dmis.W = 1./survey.std

    # wr = prob.getJtJdiag(np.ones(int(nC + len(actv))), W=dmis.W)
    wrGlobal[wires.homo.index] /= (np.max((wires.homo*wrGlobal)))
    wrGlobal[wires.hetero.index] /= (np.max(wires.hetero*wrGlobal))
    wrGlobal = wrGlobal**0.5


    ## Create a regularization
    # For the homogeneous model
    regMesh = Mesh.TensorMesh([nC])

    reg_m1 = Regularization.Sparse(regMesh, mapping=wires.homo)
    reg_m1.cell_weights = wires.homo*wrGlobal
    reg_m1.mref = mref

    # Regularization for the voxel model
    reg_m2 = Regularization.Sparse(mesh, indActive=actv, mapping=wires.hetero)
    reg_m2.cell_weights = wires.hetero*wrGlobal
    reg_m2.norms = np.c_[driver.lpnorms].T
    reg_m2.mref =  mref

    reg = reg_m1 + reg_m2

    opt = Optimization.ProjectedGNCG(maxIter=10, lower=driver.bounds[0],
                                     upper=driver.bounds[1],
                                     maxIterLS = 20, maxIterCG= 30,
                                     tolCG = 1e-4)
    invProb = InvProblem.BaseInvProblem(dmis, reg, opt)

    betaest = Directives.BetaEstimate_ByEig(beta0_ratio = 1.)
    IRLS = Directives.Update_IRLS(f_min_change=1e-4, minGNiter=2)
    update_Jacobi = Directives.UpdatePreconditioner()
    #saveModel = Directives.SaveUBCModelEveryIteration(mapping=actvMap*sumMap)
    #saveModel.fileName = work_dir + dsep + out_dir + 'GRAV'

    saveDict = Directives.SaveOutputDictEveryIteration()
    inv = Inversion.BaseInversion(invProb, directiveList=[betaest, IRLS, saveDict,
                                                          update_Jacobi])
    # Run inversion
    mrec = inv.run(mstart)

    # Plot predicted
#    pred = prob.fields(mrec)

    # PF.Gravity.plot_obs_2D(survey, 'Observed Data')
#    print("Final misfit:" + str(np.sum(((survey.dobs-pred)/survey.std)**2.)))
        # Create active map to go from reduce set to full
    outMap = Maps.InjectActiveCells(mesh, actv, ndv)

    #%% Write result
    # if getattr(invProb, 'l2model', None) is not None:

    #     m_l2 = actvMap*(sumMap*invProb.l2model)
    #     Mesh.TensorMesh.writeModelUBC(mesh, work_dir + dsep + out_dir + 'Total_inv_l2l2.den', m_l2)

    #     m_l2 = actvMap*(homogMap*wires.homo*invProb.l2model)
    #     Mesh.TensorMesh.writeModelUBC(mesh, work_dir + dsep + out_dir + 'Homoge_inv_l2l2.den', m_l2)

    #     m_l2 = actvMap*(wires.hetero*invProb.l2model)
    #     Mesh.TensorMesh.writeModelUBC(mesh, work_dir + dsep + out_dir + 'Hetero_inv_l2l2.den', m_l2)

    #     PF.Gravity.writeUBCobs(work_dir + out_dir + dsep + 'Predicted_l2.pre',
    #                          survey, d=survey.dpred(invProb.l2model))

if meshType == "TENSOR":
    m_lp = outMap*(sumMap*invProb.model)
    Mesh.TensorMesh.writeModelUBC(mesh, work_dir + dsep + out_dir + 'Total_inv_lp.den', m_lp)

    m_lp = outMap*(homogMap*wires.homo*invProb.model)
    Mesh.TensorMesh.writeModelUBC(mesh, work_dir + dsep + out_dir + 'Homoge_inv_lp.den', m_lp)

    m_lp = outMap*(wires.hetero*invProb.model)
    Mesh.TensorMesh.writeModelUBC(mesh, work_dir + dsep + out_dir + 'Hetero_inv_lp.den', m_lp)

else:
    m_lp = outMap*(sumMap*invProb.model)
    Mesh.TreeMesh.writeUBC(mesh, work_dir + dsep + out_dir + "TreeMesh.msh", models={work_dir + dsep + out_dir + 'Total_inv_lp.den': m_lp})

    m_lp = outMap*(homogMap*wires.homo*invProb.model)
    Mesh.TreeMesh.writeUBC(mesh, work_dir + dsep + out_dir + "TreeMesh.msh", models={work_dir + dsep + out_dir + 'Homoge_inv_lp.den': m_lp})

    m_lp = outMap*(wires.hetero*invProb.model)
    Mesh.TreeMesh.writeUBC(mesh, work_dir + dsep + out_dir + "TreeMesh.msh", models={work_dir + dsep + out_dir + 'Hetero_inv_lp.den': m_lp})

    PF.Gravity.writeUBCobs(work_dir + out_dir + dsep + 'Predicted_lp.pre',
                             survey, d=invProb.dpred)
