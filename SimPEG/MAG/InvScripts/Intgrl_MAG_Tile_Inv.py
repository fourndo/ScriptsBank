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
work_dir = "C:\\Users\\DominiqueFournier\\ownCloud\\\Research\\Synthetic\\Block_Gaussian_topo\\"
#work_dir = "C:\\Users\\DominiqueFournier\\ownCloud\\Research\\Synthetic\\Nut_Cracker\\"
#work_dir = "C:\\Egnyte\\Private\\dominiquef\\Projects\\4559_CuMtn_ZTEM\\Modeling\\MAG\\A1_Fenton\\"
#work_dir = 'C:\\Users\\DominiqueFournier\\ownCloud\\Research\\Kevitsa\\Modeling\\MAG\\Aiborne\\'
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

actvMap = Maps.InjectActiveCells(mesh, actv, -100)

rxLoc = survey.srcField.rxList[0].locs
tree = cKDTree(np.c_[mesh.gridCC[actv, 0],
                     mesh.gridCC[actv, 1],
                     mesh.gridCC[actv, 2]])


# # TILE THE PROBLEM
# #%% Create tiles
# # Begin tiling

# In the x-direction
nNx = 2
nNy = 1
maxObs = 30
dx = [mesh.hx.min(), mesh.hy.min()]
# Cell size

xlim = [rxLoc[:, 0].min(), rxLoc[:, 0].max()]
ylim = [rxLoc[:, 1].min(), rxLoc[:, 1].max()]

nObs = 1e+8
countx = 0
county = 0

while nObs > maxObs:

    nObs = 0

    if countx > county:
        nNx += 1
    else:
        nNy += 1

    countx = 0
    county = 0
    xtiles = np.linspace(xlim[0], xlim[1], nNx)
    ytiles = np.linspace(ylim[0], ylim[1], nNy)

    inTile = []
    # Number of points in each tiles
    for ii in range(xtiles.shape[0]-1):
        for jj in range(ytiles.shape[0]-1):
            maskx = np.all([rxLoc[:, 0] >= xtiles[ii],
                           rxLoc[:, 0] <= xtiles[int(ii+1)]], axis=0)
            masky = np.all([rxLoc[:, 1] >= ytiles[jj],
                           rxLoc[:, 1] <= ytiles[int(jj+1)]], axis=0)

            countx = np.max([np.sum(maskx), countx])
            county = np.max([np.sum(masky), county])

            mask = np.all([maskx, masky], axis=0)
            nObs = np.max([nObs, np.sum(mask)])
            inTile += [mask]

x1, x2 = xtiles[:-1], xtiles[1:]
y1, y2 = ytiles[:-1], ytiles[1:]

X1, Y1 = np.meshgrid(x1, y1)
X1, Y1 = mkvc(X1), mkvc(Y1)
X2, Y2 = np.meshgrid(x2, y2)
X2, Y2 = mkvc(X2), mkvc(Y2)

# Plot data and tiles
fig, ax1 = plt.figure(), plt.subplot()
PF.Magnetics.plot_obs_2D(rxLoc, survey.dobs, ax=ax1)
for ii in range(X1.shape[0]):
    ax1.add_patch(Rectangle((X1[ii], Y1[ii]),
                            X2[ii]-X1[ii],
                            Y2[ii]-Y1[ii],
                            facecolor='none', edgecolor='k'))
ax1.set_xlim([x1[0]-20, x2[-1]+20])
ax1.set_ylim([y1[0]-20, y2[-1]+20])
ax1.set_aspect('equal')
plt.show()

# LOOP THROUGH TILES
npadxy = 8
expf = 1.3
problems = []
surveys = []
surveyMask = np.ones(driver.survey.nD, dtype='bool')
for tt in range(X1.shape[0]):

    if tt == 0:
        tree = cKDTree(np.c_[mesh.gridCC[actv, 0],
                             mesh.gridCC[actv, 1],
                             mesh.gridCC[actv, 2]])

    midX = np.mean([X1[tt], X2[tt]])
    midY = np.mean([Y1[tt], Y2[tt]])

    # Make sure the core has odd number of cells + 3 buffer cells
    nCx = int((X2[tt] - X1[tt]) / dx[0])
    nCx += 7 - nCx % 2
    nCy = int((Y2[tt] - Y1[tt]) / dx[1])
    nCy += 7 - nCy % 2

    core_x = np.ones(nCx)*dx[0]
    core_y = np.ones(nCy)*dx[1]

    # Create new mesh
    padx = np.r_[dx[0]*expf**(np.asarray(range(npadxy))+1)]
    pady = np.r_[dx[1]*expf**(np.asarray(range(npadxy))+1)]

    # Add paddings
    hx = np.r_[padx[::-1], core_x, padx]
    hy = np.r_[pady[::-1], core_y, pady]
    hz = mesh.hz

    mesh_t = Mesh.TensorMesh([hx, hy, hz], 'CC0')

    # Shift tile center to closest cell
    shiftx = mesh.vectorCCx - midX
    shifty = mesh.vectorCCy - midY

    shiftx = shiftx[np.argmin(np.abs(shiftx))]
    shifty = shifty[np.argmin(np.abs(shifty))]

    mesh_t._x0 = (mesh_t.x0[0] + midX + shiftx,
                  mesh_t.x0[1] + midY + shifty,
                  mesh.x0[2])

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

    # Extract model from global to local mesh
    if driver.topofile is not None:
        topo = np.genfromtxt(work_dir + driver.topofile,
                             skip_header=1)
        actv_t = Utils.surface2ind_topo(mesh_t, topo, 'N')
        # actv_t = np.asarray(np.where(mkvc(actv))[0], dtype=int)
    else:
        actv_t = np.ones(mesh_t.nC, dtype='bool')

    nC = len(actv_t)

    # Create active map to go from reduce set to full
    # actvMap = Maps.InjectActiveCells(mesh_t, actv, -100)

    # Creat reduced identity map
    tileMap = Maps.Tile((mesh, actv), (mesh_t, actv_t), tree=tree)

    # Create the forward model operator
    prob = PF.Magnetics.MagneticIntegral(mesh_t, chiMap=tileMap, actInd=actv_t)

    problems += [prob]
    surveys += [survey_t]


# Going through all problems:
# 1- Pair the survey and problem
# 2- Add up sensitivity weights
# 3- Create a ComboProblem
wrGlobal = np.zeros(int(actv.sum()))
dmisfits = []
for prob, survey in zip(problems, surveys):

    survey.pair(prob)

    # Data misfit function
    dmis = DataMisfit.l2_DataMisfit(survey)
    dmis.W = 1./survey.std

    dmisfits += [dmis]

    wr = np.sum(prob.G**2., axis=0)

    wrGlobal += prob.chiMap.deriv(0).T*wr

# Scale global weights for regularization
wrGlobal = wrGlobal**0.5
wrGlobal = (wrGlobal/np.max(wrGlobal))

# Create combo misfit function
if len(dmisfits) > 1:
    ComboMisfit = dmisfits[0] + dmisfits[1]
    for dmisfit in dmisfits[2:]:
        ComboMisfit += dmisfit
else:
    ComboMisfit = dmisfits[0]

# Create a regularization
idenMap = Maps.IdentityMap(nP=int(np.sum(actv)))
reg = Regularization.Sparse(mesh, indActive=actv, mapping=idenMap)
reg.norms = driver.lpnorms

if driver.eps is not None:
    reg.eps_p = driver.eps[0]
    reg.eps_q = driver.eps[1]

reg.cell_weights = wrGlobal
reg.mref = np.zeros(mesh.nC)[actv]

# Add directives to the inversion
opt = Optimization.ProjectedGNCG(maxIter=20, lower=0., upper=10.,
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
Mesh.TensorMesh.writeModelUBC(mesh, work_dir + out_dir + "MAG_Tile_l2.sus",
                              actvMap*IRLS.l2model)
Mesh.TensorMesh.writeModelUBC(mesh, work_dir + out_dir + "MAG_Tile_lp.sus",
                              actvMap*invProb.model)

# Get predicted data for each tile and write full predicted to file
dpred = np.zeros(driver.survey.nD)
for ind, survey in enumerate(surveys):
    dpred[survey.ind] += survey.dpred(mrec)
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

