"""

This script runs an Magnetic Vector Inversion (MVI) tiled with
octree meshes. The scipt uses Dask to parallelize and store individual
sensitivities to disk, reducing the RAM footprint.

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
from SimPEG.ObjectiveFunction import ComboObjectiveFunction
import os


if __name__ == '__main__':
#    work_dir = "C:\\Users\\DominiqueFournier\\Dropbox\\Projects\\Synthetic\\Triple_Block_lined\\"
    work_dir = "C:\\Users\\DominiqueFournier\\Dropbox\\Projects\\Synthetic\\Nut_Cracker\\"
    #    work_dir = "C:\\Users\\DominiqueFournier\\Downloads\\Ruapehu\\"

#    work_dir = "C:\\Users\\DominiqueFournier\\ownCloud\\Research\\Osborne\\Inversion\\UTM\\"
    out_dir = "SimPEG_MVI_S_TileInv\\"
    input_file = "SimPEG_MAG.inp"

    dsep = os.path.sep


    padLen = 300  # Padding distance around the convex haul of the data
    maxRAM = .07  # Maximum memory size allowed for each tiles
    n_cpu = 4     # Number of processors allocated

    octreeObs = [5, 5, 10]  # Octree levels below data points [n1*dz, n2*dz**2, ...]
    octreeTopo = [0, 0 ,1]   # Octree levels below topography [n1*dz, n2*dz**2, ...]

    meshType = 'TREE'   # Inversion mesh type
    tileProblem = True  # Tile the forward True | False
    parallelized = "dask"  #"multiprocessing" # "dask" ||  "multiprocessing"
    ###############################################################################
    # Inversion script starts here
    # ----------------------------
    # Read in the input file which included all parameters at once
    # (mesh, topo, model, survey, inv param, etc.)
    #

    os.system('if not exist ' + work_dir + out_dir + ' mkdir ' + work_dir+out_dir)
    driver = PF.MagneticsDriver.MagneticsDriver_Inv(work_dir + input_file)

    # Access the mesh and survey information
    meshInput = driver.mesh
    survey = driver.survey
    rxLoc = survey.srcField.rxList[0].locs

    topo = None
    if driver.topofile is not None:
        topo = np.genfromtxt(driver.basePath + driver.topofile,
                             skip_header=1)
    else:
        # Grab the top coordinate of mesh and make a flat topo
        indTop = meshInput.gridCC[:, 2] == meshInput.vectorCCz[-1]
        topo = meshInput.gridCC[indTop, :]
        topo[:, 2] += meshInput.hz.min()/2. + 1e-8

    # Remember selected data in case of tile overlap
    surveyMask = np.ones(survey.nD, dtype='bool')

    # Define the octree mesh based on the provided tensor
    h = np.r_[meshInput.hx.min(), meshInput.hy.min(), meshInput.hz.min()]
    padDist = np.r_[np.c_[padLen, padLen], np.c_[padLen, padLen], np.c_[padLen, 0]]

    if meshType != meshInput._meshType:
        print("Creating Global Octree")
        mesh = Utils.modelutils.meshBuilder(
                rxLoc, h, padDist, meshType='TREE', meshGlobal=meshInput,
                verticalAlignment='center'
            )

        # Refine the octree mesh below topography
        if topo is not None:
            mesh = Utils.modelutils.refineTree(
                mesh, topo, dtype='surface',
                nCpad=octreeTopo, finalize=False
            )

        # Refine the octree mesh around the obs and finalize
        mesh = Utils.modelutils.refineTree(
            mesh, rxLoc, dtype='surface',
            nCpad=octreeObs, finalize=True
        )

        # Compute active cells
        if topo is not None:
            actv = Utils.surface2ind_topo(mesh, topo)
        else:
            actv = np.zeros(mesh.nC, dtype='bool')
            print(meshInput.vectorNz[-1])
            actv[mesh.gridCC[:, 2] < meshInput.vectorNz[-1]] = True

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
    wrGlobal = np.zeros(int(actv.sum()*3))


    if tileProblem:

        # Loop over different tile size and break problem until
        # memory footprint false below maxRAM
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
            nD, nC = ind_t.sum()*1., meshLocal.nC*1.

            nChunks = n_cpu # Number of chunks
            cSa, cSb = int(nD/nChunks), int(nC/nChunks) # Chunk sizes
            usedRAM = nD * nC * 8. * 1e-9 * 3
            count += 1
            print(nD, nC, usedRAM, binCount.min())

        # After tiling:
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
            # createLocalProb(rxLoc, wrGlobal, lims, ind)
            # Generate a problem, calculate/store sensitivities for
            # given data points

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



            # Create reduced identity map
            tileMap = Maps.Tile((mesh, actv), (meshLocal, actv_t))
            tileMap.nBlock = 3

            # Create the forward model operator
            prob = PF.Magnetics.MagneticIntegral(
                        meshLocal, chiMap=tileMap, actInd=actv_t, parallelized=parallelized,
                        Jpath=work_dir + out_dir + "Tile" + str(ind) + ".zarr",
                        modelType='vector', n_cpu=n_cpu
            )

            survey_t.pair(prob)

            # Data misfit function
            dmis = DataMisfit.l2_DataMisfit(survey_t)
            dmis.W = 1./survey_t.std

            wr = prob.getJtJdiag(np.ones(int(tileMap.P.shape[1])), W=dmis.W)
            Mesh.TreeMesh.writeUBC(
                  meshLocal, work_dir + out_dir + 'OctreeMesh' + str(tt) + '.msh',
                  models={work_dir + out_dir + 'Wr_' + str(tt) + '.act': (tileMap*wr)[:meshLocal.nC]+1e-8}
                )
            wrGlobal += wr

            del meshLocal
            # Create combo misfit function
            return dmis, wrGlobal

        # Loop through the tiles and generate all sensitivities
        for tt in range(X1.shape[0]):

            print("Tile " + str(tt+1) + " of " + str(X1.shape[0]))

            dmis, wrGlobal = createLocalProb(rxLoc, wrGlobal, np.r_[X1[tt], X2[tt], Y1[tt], Y2[tt]], tt)

            # Add the problems to a Combo Objective function
            if tt == 0:
                ComboMisfit = dmis

            else:
                ComboMisfit += dmis

    # If not tiled, just a single problem
    else:

        # Create the forward model operator
        # Create identity map
        nC = int(actv.sum())
        idenMap = Maps.IdentityMap(nP=int(3*nC))

        prob = PF.Magnetics.MagneticIntegral(
            mesh, chiMap=idenMap, actInd=actv, parallelized=parallelized,
            Jpath=work_dir + out_dir + "Sensitivity.zarr",
            modelType='vector', n_cpu=n_cpu)

        survey.pair(prob)

        # Data misfit function
        ComboMisfit = DataMisfit.l2_DataMisfit(survey)
        ComboMisfit.W = 1./survey.std

        wrGlobal += prob.getJtJdiag(np.ones(int(3*nC)), W=ComboMisfit.W)
        actvGlobal = actv

    # Scale global weights for regularization
    # Check if global mesh has regions untouched by local problem
    # NEED REVIEW FOR THIS CASE
    # (ONLY AN ISSUE FOR TENSOR INVERSION MESH)
    nC = int(actv.sum())
    actvGlobal = wrGlobal != 0
    actvMeshGlobal = (wrGlobal[:nC]) != 0
    if actvMeshGlobal.sum() < actv.sum():

        if isinstance(ComboMisfit, ComboObjectiveFunction):
            for ind, dmis in enumerate(ComboMisfit.objfcts):
                dmis.prob.chiMap.index = actvMeshGlobal
                dmis.prob.gtgdiag = None

        else:
            ComboMisfit.prob.chiMap.index = actvGlobal
            ComboMisfit.prob.rhoMap._P = None
            ComboMisfit.prob.model = np.zeros(actvGlobal.sum())
            ComboMisfit.prob.gtgdiag = None

    # Global sensitivity weights (linear)
    wrGlobal = wrGlobal[actvGlobal]**0.5
    wrGlobal = (wrGlobal/np.max(wrGlobal))

    # Create global active set
    actv = np.all([actv, actvMap*actvMeshGlobal], axis=0)
    actvMap = Maps.InjectActiveCells(mesh, actv, 0)  # For re-projection
    actvMapAmp = Maps.InjectActiveCells(mesh, actv, -100)  # For final output

    nC = int(np.sum(actv))

    mstart = np.ones(3*nC) * 1e-4

    # Assumes amplitude reference, distributed on 3 components
    mref = np.ones(3*nC) * (np.mean(driver.mref)**2./3)**0.5

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

    # Assemble the 3-component regularizations
    reg = reg_p + reg_s + reg_t
    reg.mref = mref


    opt = Optimization.ProjectedGNCG(maxIter=5, lower=-10., upper=10.,
                                     maxIterCG=20, tolCG=1e-3)

    invProb = InvProblem.BaseInvProblem(ComboMisfit, reg, opt)
    betaest = Directives.BetaEstimate_ByEig()

    # Add directives to the inversion
    # Here is where the norms are applied
    IRLS = Directives.Update_IRLS(f_min_change=1e-3,
                                  minGNiter=1)

    # Pre-conditioner
    update_Jacobi = Directives.UpdatePreconditioner()

    # Output models between each iteration
    saveModel = Directives.SaveUBCModelEveryIteration(mapping=actvMap, vector=True)
    saveModel.fileName = work_dir + out_dir + 'MVI_C'

    # Put it all together
    inv = Inversion.BaseInversion(invProb,
                                  directiveList=[betaest, IRLS, update_Jacobi,
                                                 saveModel])

    # Invert
    mrec_MVI = inv.run(mstart)


    ###############################################################################
    # MVI-Spherical with sparsity
    # ---------------------------
    #
    # Finish inversion with spherical formulation
    #

    # Extract the vector components for the MVI-S
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
        mesh.writeModelUBC(work_dir+out_dir + 'MVI_C_amp.sus', amp)

    # Get predicted data for each tile and form/write full predicted to file
    if getattr(ComboMisfit, 'objfcts', None) is not None:
        dpred = np.zeros(survey.nD)
        for ind, dmis in enumerate(ComboMisfit.objfcts):
            dpred[dmis.survey.ind] += dmis.survey.dpred(mrec_MVI)
    else:
        dpred = ComboMisfit.survey.dpred(mrec_MVI)

    Utils.io_utils.writeUBCmagneticsObservations(
      work_dir+out_dir + 'MVI_C_pred.pre', survey, dpred
    )

    beta = invProb.beta

    # Change the starting model from Cartesian to Spherical
    mstart = Utils.matutils.xyz2atp(mrec_MVI.reshape((nC, 3), order='F'))
    mref = np.kron(np.r_[np.mean(driver.mref), -np.deg2rad(survey.srcField.param[1]), np.deg2rad(survey.srcField.param[2])], np.ones(nC))

    # Flip the problem from Cartesian to Spherical
    if getattr(ComboMisfit, 'objfcts', None) is not None:
        for misfit in ComboMisfit.objfcts:
            misfit.prob.coordinate_system = 'spherical'
            misfit.prob.model = mstart
    else:
        ComboMisfit.prob.coordinate_system = 'spherical'
        ComboMisfit.prob.model = mstart

    # Create a block diagonal regularization
    wires = Maps.Wires(('amp', nC), ('theta', nC), ('phi', nC))

    # Create a regularization
    reg_a = Regularization.Sparse(mesh, indActive=actv,
                                  mapping=wires.amp, gradientType='component')
    reg_a.norms = np.c_[driver.lpnorms[:4]].T
    reg_a.mref = mref


    reg_t = Regularization.Sparse(mesh, indActive=actv,
                                  mapping=wires.theta, gradientType='component')
    reg_t.alpha_s = 1e-8  # No reference angle
    reg_t.space = 'spherical'
    reg_t.norms = np.c_[driver.lpnorms[4:8]].T
    reg_t.eps_q = np.deg2rad(15)
    reg_t.mref = mref

    reg_p = Regularization.Sparse(mesh, indActive=actv,
                                  mapping=wires.phi, gradientType='component')
    reg_p.alpha_s = 1e-8  # No reference angle
    reg_p.space = 'spherical'
    reg_p.norms = np.c_[driver.lpnorms[8:]].T
    reg_p.eps_q = np.deg2rad(45)
    reg_p.mref = mref

    # Assemble the three regularization
    reg = reg_a + reg_t + reg_p
    reg.mref = mref

    Lbound = np.kron(np.asarray([0, -np.inf, -np.inf]), np.ones(nC))
    Ubound = np.kron(np.asarray([10, np.inf, np.inf]), np.ones(nC))


    # Add directives to the inversion
    opt = Optimization.ProjectedGNCG(maxIter=40,
                                     lower=Lbound,
                                     upper=Ubound,
                                     maxIterLS=20,
                                     maxIterCG=20, tolCG=1e-3,
                                     stepOffBoundsFact=1e-8,
                                     LSshorten=0.25)

    invProb = InvProblem.BaseInvProblem(ComboMisfit, reg, opt, beta=beta*3)
    #  betaest = Directives.BetaEstimate_ByEig()

    # Here is where the norms are applied
    IRLS = Directives.Update_IRLS(f_min_change=1e-4, maxIRLSiter=40,
                                  minGNiter=1, beta_tol=0.5, prctile=100,
                                  coolingRate=1, coolEps_q=True,
                                  betaSearch=False)

    # Special directive specific to the mag amplitude problem. The sensitivity
    # weights are update between each iteration.
    ProjSpherical = Directives.ProjSpherical()
    update_SensWeight = Directives.UpdateSensitivityWeights()
    update_Jacobi = Directives.UpdatePreconditioner()
    saveModel = Directives.SaveUBCModelEveryIteration(mapping=actvMap, vector=True)
    saveModel.fileName = work_dir+out_dir + 'MVI_S'

    inv = Inversion.BaseInversion(invProb,
                                  directiveList=[
                                    ProjSpherical, IRLS, update_SensWeight,
                                    update_Jacobi, saveModel
                                    ])

    # Run the inversion
    mrec_MVI_S = inv.run(mstart)


    # Get predicted data for each tile and write full predicted to file
    if getattr(ComboMisfit, 'objfcts', None) is not None:
        dpred = np.zeros(survey.nD)
        for ind, dmis in enumerate(ComboMisfit.objfcts):
            dpred[dmis.survey.ind] += dmis.survey.dpred(mrec_MVI_S)
    else:
        dpred = ComboMisfit.survey.dpred(mrec_MVI_S)

    Utils.io_utils.writeUBCmagneticsObservations(work_dir+out_dir + 'MVI_S_pred.pre', survey, dpred)

    # Export tensor model if inputMesh is tensor
    if isinstance(meshInput, Mesh.TensorMesh):

        tree = cKDTree(mesh.gridCC)

        dd, ind = tree.query(meshInput.gridCC)

        # Model lp out
        vec_xyz = Utils.matutils.atp2xyz(
            invProb.model.reshape((nC, 3), order='F')).reshape((nC, 3), order='F')
        vec_x = actvMap * vec_xyz[:, 0]
        vec_y = actvMap * vec_xyz[:, 1]
        vec_z = actvMap * vec_xyz[:, 2]

        vec_xyzTensor = np.zeros((meshInput.nC, 3))
        vec_xyzTensor = np.c_[vec_x[ind], vec_y[ind], vec_z[ind]]

        Utils.io_utils.writeVectorUBC(
            meshInput, work_dir+out_dir + 'MVI_S_TensorLp.fld', vec_xyzTensor)
        amp = np.sum(vec_xyzTensor**2., axis=1)**0.5
        meshInput.writeModelUBC(work_dir+out_dir + 'MVI_S_TensorLp.amp', amp)

        # Model l2 out
        vec_xyz = Utils.matutils.atp2xyz(
            invProb.l2model.reshape((nC, 3), order='F')).reshape((nC, 3), order='F')
        vec_x = actvMap * vec_xyz[:, 0]
        vec_y = actvMap * vec_xyz[:, 1]
        vec_z = actvMap * vec_xyz[:, 2]

        vec_xyzTensor = np.zeros((meshInput.nC, 3))
        vec_xyzTensor = np.c_[vec_x[ind], vec_y[ind], vec_z[ind]]

        Utils.io_utils.writeVectorUBC(
            meshInput, work_dir+out_dir + 'MVI_S_TensorL2.fld', vec_xyzTensor)

        amp = np.sum(vec_xyzTensor**2., axis=1)**0.5
        meshInput.writeModelUBC(work_dir+out_dir + 'MVI_S_TensorL2.amp', amp)
