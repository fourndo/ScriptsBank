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
import os

#work_dir = "C:\\Users\\DominiqueFournier\\ownCloud\\Research\\Kevitsa\\Modeling\\MAG\\"
#work_dir = "C:\\Users\\DominiqueFournier\\ownCloud\\\Research\\Synthetic\\Block_Gaussian_topo\\"
#work_dir = "C:\\Users\\DominiqueFournier\\ownCloud\\Research\\Synthetic\\Nut_Cracker\\"
#work_dir = "C:\\Egnyte\\Private\\dominiquef\\Projects\\4559_CuMtn_ZTEM\\Modeling\\MAG\\A1_Fenton\\"
#work_dir = 'C:\\Users\\DominiqueFournier\\ownCloud\\Research\\Kevitsa\\Modeling\\MAG\\Aiborne\\'
work_dir = "C:\\Users\\DominiqueFournier\\ownCloud\\Research\\CraigModel\\MAG\\"
out_dir = "SimPEG_Susc_TileInv\\"
input_file = "SimPEG_MAG.inp"

# %%
# Read in the input file which included all parameters at once
# (mesh, topo, model, survey, inv param, etc.)
driver = PF.MagneticsDriver.MagneticsDriver_Inv(work_dir + input_file)

os.system('if not exist ' + work_dir + out_dir + ' mkdir ' + work_dir+out_dir)

# Access the mesh and survey information
mesh = driver.mesh
survey = driver.survey
actv = np.zeros(mesh.nC, dtype='bool')
actv[driver.activeCells] = True

rxLoc = survey.srcField.rxList[0].locs
#tree = cKDTree(np.c_[mesh.gridCC[actv, 0],
#                     mesh.gridCC[actv, 1],
#                     mesh.gridCC[actv, 2]])


# # TILE THE PROBLEM
# Define core mesh properties
h = np.r_[[np.min(np.r_[mesh.hx.min(), mesh.hy.min(), mesh.hz.min()])]*3]

maxNpoints = 100

tiles = Utils.modelutils.tileSurveyPoints(rxLoc, maxNpoints)

X1, Y1 = tiles[0][:,0], tiles[0][:, 1]
X2, Y2 = tiles[1][:,0], tiles[1][:, 1]


# Plot data and tiles
fig, ax1 = plt.figure(), plt.subplot()
PF.Magnetics.plot_obs_2D(rxLoc, survey.dobs, ax=ax1)
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
expf = 1.3
dx = [mesh.hx.min(), mesh.hy.min()]
surveyMask = np.ones(driver.survey.nD, dtype='bool')
# Going through all problems:
# 1- Pair the survey and problem
# 2- Add up sensitivity weights
# 3- Add to the ComboMisfit
wrGlobal = np.zeros(int(actv.sum()))
probSize = 0
for tt in range(X1.shape[0]):
    
    print("Tile " + str(tt+1) + " of " + str(X1.shape[0]))
    if tt == 0:
        tree = cKDTree(np.c_[mesh.gridCC[actv, 0],
                             mesh.gridCC[actv, 1],
                             mesh.gridCC[actv, 2]])

    # Grab the data for current tile
    ind_t = np.all([rxLoc[:, 0] >= X1[tt], rxLoc[:, 0] <= X2[tt],
                    rxLoc[:, 1] >= Y1[tt], rxLoc[:, 1] <= Y2[tt],
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

    padDist = np.r_[np.c_[1500, 1500], np.c_[1500, 1500], np.c_[1500, 0]]
    mesh_t = Utils.modelutils.meshBuilder(rxLoc[ind_t, :], h,
                                          padDist, meshGlobal=mesh,
                                          meshType='TREE',
                                          padCore=np.r_[3, 3, 3])

    # Extract model from global to local mesh
#    if driver.topofile is not None:
#        topo = np.genfromtxt(work_dir + driver.topofile,
#                             skip_header=1)
#        actv_t = Utils.surface2ind_topo(mesh_t, topo, 'N')
#        # actv_t = np.asarray(np.where(mkvc(actv))[0], dtype=int)
#    else:
    actv_t = np.ones(mesh_t.nC, dtype='bool')

    # Create reduced identity map
    tileMap = Maps.Tile((mesh, actv), (mesh_t, actv_t))
    tileMap.nCell = 27
    # Create the forward model operator
    prob = PF.Magnetics.MagneticIntegral(mesh_t, chiMap=tileMap, actInd=actv_t)
    survey_t.pair(prob)

    # Data misfit function
    dmis = DataMisfit.l2_DataMisfit(survey_t)
    dmis.W = 1./survey_t.std
#    ncell = Utils.mkvc(np.sum(prob.chiMap.deriv(0) > 0, axis=1))
#    ncell[ncell>0] = (1./ncell[ncell>0]).copy()
#    SCALE = Utils.sdiag(mesh_t.vol**(0))
    for ii in range(prob.G.shape[0]):
        wrGlobal += ((prob.G[ii, :])*(prob.chiMap.deriv(0)))**2.
    
#    wrGlobal += np.abs(prob.Jtvec(0, prob.Jvec(0, np.ones(mesh.nC)*1e-4)))
#    wrGlobal += prob.chiMap.deriv(0).T*wr

    # Create combo misfit function
    if tt == 0:
        ComboMisfit = dmis

    else:
        ComboMisfit += dmis

    # Add problem size
    probSize += prob.G.shape[0] * prob.G.shape[1] * 32 / 4

print('Sum of all problems:' + str(probSize*1e-6) + ' Mb')
# Scale global weights for regularization

# Check if global mesh has regions untouched by local problem
actvGlobal = wrGlobal != 0
if actvGlobal.sum() < actv.sum():

    for ind, dmis in enumerate(ComboMisfit.objfcts):
        dmis.prob.chiMap.index = actvGlobal

wrGlobal = wrGlobal[actvGlobal]**0.5
wrGlobal = (wrGlobal/np.max(wrGlobal))

# Create a regularization
actvMap = Maps.InjectActiveCells(mesh, actv, 0)
actv = np.all([actv, actvMap*actvGlobal], axis=0)
actvMap = Maps.InjectActiveCells(mesh, actv, -100)
idenMap = Maps.IdentityMap(nP=int(np.sum(actv)))
reg = Regularization.Sparse(mesh, indActive=actv, mapping=idenMap)
reg.norms = driver.lpnorms

if driver.eps is not None:
    reg.eps_p = driver.eps[0]
    reg.eps_q = driver.eps[1]

reg.cell_weights = wrGlobal
reg.mref = np.zeros(mesh.nC)[actv]

# Add directives to the inversion
opt = Optimization.ProjectedGNCG(maxIter=40, lower=0., upper=10.,
                                 maxIterLS=20, maxIterCG=10, tolCG=1e-4)
invProb = InvProblem.BaseInvProblem(ComboMisfit, reg, opt)
betaest = Directives.BetaEstimate_ByEig()

# Here is where the norms are applied
# Use pick a treshold parameter empirically based on the distribution of
#  model parameters
IRLS = Directives.Update_IRLS(f_min_change=1e-3, minGNiter=3,
                              maxIRLSiter=10)

IRLS.target = driver.survey.nD
update_Jacobi = Directives.UpdateJacobiPrecond()

inv = Inversion.BaseInversion(invProb,
                              directiveList=[betaest, IRLS, update_Jacobi])

# Run the inversion
m0 = np.ones(mesh.nC)[actv]*1e-4  # Starting model
mrec = inv.run(m0)

# Outputs
Mesh.TensorMesh.writeUBC(mesh, work_dir + out_dir + "MAG_Tile.msh")
Mesh.TensorMesh.writeModelUBC(mesh, work_dir + out_dir + "MAG_Tile_lp.sus",
                              actvMap*invProb.model)
Mesh.TensorMesh.writeModelUBC(mesh, work_dir + out_dir + "MAG_Tile_l2.sus",
                              actvMap*IRLS.l2model)


# Get predicted data for each tile and write full predicted to file
if getattr(ComboMisfit, 'objfcts', None) is not None:
    dpred = np.zeros(driver.survey.nD)
    for ind, dmis in enumerate(ComboMisfit.objfcts):
        dpred[dmis.survey.ind] += dmis.survey.dpred(mrec)
else:
    dpred = ComboMisfit.survey.dpred(mrec)
    # PF.Magnetics.writeUBCobs(work_dir+out_dir + "Tile" + str(ind) + ".pre",
    #                          survey, survey.dpred(mrec))

PF.Magnetics.writeUBCobs(work_dir+out_dir + "MAG_Tile_Inv.pre",
                         driver.survey, dpred)

# # %% OUTPUT models for each tile
# model = Mesh.TensorMesh.readModelUBC(driver.mesh,
#                                      work_dir + '\Effec_sus_20mGrid.sus')
# for ind, prob in enumerate(problems):
#     Mesh.TensorMesh.writeUBC(prob.mesh,
#                              work_dir + out_dir + "Tile" + str(ind) + ".msh")
#     Mesh.TensorMesh.writeModelUBC(prob.mesh,
#                                   work_dir + out_dir + "Tile" + str(ind) + ".sus", prob.mapping() * model)

