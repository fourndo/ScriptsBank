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
import SimPEG.PF as PF
import numpy as np
import matplotlib.pyplot as plt
import os

# # STEP 1: Setup and data simulation # #

# Magnetic inducing field parameter (A,I,D)
B = [50000, 90, 0]

# Create a mesh
dx = 10.

hxind = [(dx, 5, -1.3), (dx, 15), (dx, 5, 1.3)]
hyind = [(dx, 5, -1.3), (dx, 15), (dx, 5, 1.3)]
hzind = [(dx, 5, -1.3), (dx, 10)]

mesh = Mesh.TensorMesh([hxind, hyind, hzind], 'CC0')
mesh.x0[2] -= mesh.vectorNz[-1]

susc = 0.1
nX = 1

# Get index of the center of block
locx = [int(mesh.nCx/2)-3, int(mesh.nCx/2)+3]
midy = int(mesh.nCy/2)
midz = -4

# Lets create a simple flat topo and set the active cells
[xx, yy] = np.meshgrid(mesh.vectorNx, mesh.vectorNy)
zz = np.ones_like(xx)*mesh.vectorNz[-1]
topo = np.c_[Utils.mkvc(xx), Utils.mkvc(yy), Utils.mkvc(zz)]

# Go from topo to actv cells
actv = Utils.surface2ind_topo(mesh, topo, 'N')
actv = np.asarray([inds for inds, elem in enumerate(actv, 1) if elem],
                  dtype=int) - 1

# Create active map to go from reduce space to full
actvMap = Maps.InjectActiveCells(mesh, actv, -100)
nC = int(len(actv))

# Create and array of observation points
xr = np.linspace(-80., 80., 10)
yr = np.linspace(-80., 80., 10)
X, Y = np.meshgrid(xr, yr)

# Move the observation points 5m above the topo
Z = np.ones_like(X) * mesh.vectorNz[-1] + dx

# Create a MAGsurvey
rxLoc = np.c_[Utils.mkvc(X.T), Utils.mkvc(Y.T), Utils.mkvc(Z.T)]
rxLoc = PF.BaseMag.RxObs(rxLoc)
srcField = PF.BaseMag.SrcField([rxLoc], param=(B[0], B[1], B[2]))
survey = PF.BaseMag.LinearSurvey(srcField)

# We can now create a susceptibility model and generate data
# Here a simple block in half-space
model = np.zeros((mesh.nCx, mesh.nCy, mesh.nCz))
for midx in locx:
    model[(midx-nX):(midx+nX+1), (midy-nX):(midy+nX+1), (midz-nX):(midz+nX+1)] = susc
model = Utils.mkvc(model)
model = model[actv]

# We create a magnetization model different than the inducing field
# to simulate remanent magnetization. Let's do something simple [45,90]
M = PF.Magnetics.dipazm_2_xyz(np.ones(nC) * 0., np.ones(nC) * 90.)
#M = PF.Magnetics.dipazm_2_xyz(np.ones(nC) * B[1], np.ones(nC) * B[2])


# Create active map to go from reduce set to full
actvMap = Maps.InjectActiveCells(mesh, actv, -100)

# Create reduced identity map
idenMap = Maps.IdentityMap(nP=nC)

# Create the forward problem (forwardOnly)
prob = PF.Magnetics.MagneticIntegral(mesh, chiMap=idenMap, actInd=actv,
                                     M=M)

# Pair the survey and problem
survey.pair(prob)

# Compute forward model some data
d = prob.fields(model)

# Add noise and uncertainties
# We add some random Gaussian noise (1nT)
d_TMI = d + np.random.randn(len(d))*0.
wd = np.ones(len(d_TMI))  # Assign flat uncertainties
survey.dobs = d_TMI
survey.std = wd

# For comparison, let's run the inversion assuming an induced response
M = PF.Magnetics.dipazm_2_xyz(np.ones(nC) * B[1], np.ones(nC) * B[2])

# Reset the magnetization
prob.M = M
prob._G = None

# Create a regularization function, in this case l2l2
wr = np.sum(prob.G**2., axis=0)**0.5
wr = (wr/np.max(wr))


# Create a regularization
reg_Susc = Regularization.Sparse(mesh, indActive=actv, mapping=idenMap)
reg_Susc.cell_weights = wr

# Data misfit function
dmis = DataMisfit.l2_DataMisfit(survey)
dmis.Wd = 1/wd

# Add directives to the inversion
opt = Optimization.ProjectedGNCG(maxIter=100, lower=0., upper=1.,
                                 maxIterLS=20, maxIterCG=10, tolCG=1e-3)
invProb = InvProblem.BaseInvProblem(dmis, reg_Susc, opt)
betaest = Directives.BetaEstimate_ByEig()

# Here is where the norms are applied
# Use pick a treshold parameter empirically based on the distribution of
#  model parameters
IRLS = Directives.Update_IRLS(norms=([2, 2, 2, 2]),
                              f_min_change=1e-3, minGNiter=3)
update_Jacobi = Directives.Update_lin_PreCond()
inv = Inversion.BaseInversion(invProb,
                              directiveList=[IRLS, betaest, update_Jacobi])

# Run the inversion
m0 = np.ones(nC)*1e-4  # Starting model
mrec_sus = inv.run(m0)


# # STEP 2: Equivalent source inversion and amplitude data # #

# Get the layer of cells directly below topo
surf = Utils.surface2ind_topo(mesh, topo, 'N', layer=True)
nC = int(np.sum(surf))  # Number of active cells

# Create active map to go from reduce set to full
surfMap = Maps.InjectActiveCells(mesh, surf, -100)

# Create identity map
idenMap = Maps.IdentityMap(nP=nC)

# Create MAG equivalent layer problem
prob = PF.Magnetics.MagneticIntegral(mesh, chiMap=idenMap, actInd=surf,
                                     equiSourceLayer=True)
prob.solverOpts['accuracyTol'] = 1e-4

# Pair the survey and problem
survey.pair(prob)

# Create a regularization function, in this case l2l2
reg = Regularization.Simple(mesh, indActive=surf)
reg.mref = np.zeros(nC)

# Specify how the optimization will proceed
opt = Optimization.ProjectedGNCG(maxIter=150, lower=-np.inf,
                                 upper=np.inf, maxIterLS=20,
                                 maxIterCG=20, tolCG=1e-3)

# Define misfit function (obs-calc)
dmis = DataMisfit.l2_DataMisfit(survey)
dmis.Wd = 1./survey.std

# Create the default L2 inverse problem from the above objects
invProb = InvProblem.BaseInvProblem(dmis, reg, opt)

# Specify how the initial beta is found
betaest = Directives.BetaEstimate_ByEig()

# Beta schedule for inversion
betaSchedule = Directives.BetaSchedule(coolingFactor=2., coolingRate=1)

# Target misfit to stop the inversion
targetMisfit = Directives.TargetMisfit()

# Put all the parts together
inv = Inversion.BaseInversion(invProb,
                              directiveList=[betaest, betaSchedule,
                                             targetMisfit])

# Run the equivalent source inversion
mstart = np.zeros(nC)
mrec = inv.run(mstart)

# COMPUTE AMPLITUDE DATA
# Now that we have an equialent source layer, we can forward model all
# three components of the field and add them up:
# |B| = ( Bx**2 + Bx**2 + Bx**2 )**0.5

# Won't store the sensitivity and output 'xyz' data.
prob.forwardOnly = True
prob.rtype = 'xyz'
pred = prob.Intrgl_Fwr_Op(m=mrec)

ndata = survey.nD

d_amp = np.sqrt(pred[:ndata]**2. +
                pred[ndata:2*ndata]**2. +
                pred[2*ndata:]**2.)

rxLoc = survey.srcField.rxList[0].locs

# # STEP 3: RUN AMPLITUDE INVERSION ##

# Now that we have |B| data, we can invert. This is a non-linear inversion,
# which requires some special care for the sensitivity weighting
# (see Directives)

# Create active map to go from reduce space to full
actvMap = Maps.InjectActiveCells(mesh, actv, -100)
nC = int(len(actv))

# Create identity map
idenMap = Maps.IdentityMap(nP=nC)

# Create the forward model operator
prob = PF.Magnetics.MagneticAmplitude(mesh, chiMap=idenMap,
                                      actInd=actv)

# Define starting model
mstart = np.ones(len(actv))*1e-4
prob.chi = mstart

# Change the survey to xyz components
survey.srcField.rxList[0].rxType = 'xyz'

# Pair the survey and problem
survey.pair(prob)

# Re-set the observations to |B|
survey.dobs = d_amp

# Create a sparse regularization
reg = Regularization.Sparse(mesh, indActive=actv, mapping=idenMap)
reg.mref = mstart*0.
reg.cell_weights=wr
# Data misfit function
dmis = DataMisfit.l2_DataMisfit(survey)
dmis.Wd = 2./(d_amp.min())

# Add directives to the inversion
opt = Optimization.ProjectedGNCG(maxIter=100, lower=0., upper=1.,
                                 maxIterLS=20, maxIterCG=10,
                                 tolCG=1e-3)

invProb = InvProblem.BaseInvProblem(dmis, reg, opt)

# Here is the list of directives
betaest = Directives.BetaEstimate_ByEig()

# Specify the sparse norms
IRLS = Directives.Update_IRLS(norms=([2, 2, 2, 2]),
                              eps=(1e-3, 1e-3), f_min_change=1e-3,
                              minGNiter=3, chifact=0.25)

# Special directive specific to the mag amplitude problem. The sensitivity
# weights are update between each iteration.
update_Jacobi = Directives.Amplitude_Inv_Iter()
update_Jacobi.test = True

# Put all together
inv = Inversion.BaseInversion(invProb,
                              directiveList=[IRLS, update_Jacobi, betaest])

# Invert
mrec_MAI = inv.run(mstart)

# # REPEAT WITH SENSITIVITY REWEIGHTING
# Create a sparse regularization
reg = Regularization.Sparse(mesh, indActive=actv, mapping=idenMap)
reg.mref = mstart*0.
#reg.cell_weights=wr
# Data misfit function
dmis = DataMisfit.l2_DataMisfit(survey)
dmis.Wd = 2./(d_amp.min())

# Add directives to the inversion
opt = Optimization.ProjectedGNCG(maxIter=100, lower=0., upper=1.,
                                 maxIterLS=20, maxIterCG=10,
                                 tolCG=1e-3)

invProb = InvProblem.BaseInvProblem(dmis, reg, opt)

# Here is the list of directives
betaest = Directives.BetaEstimate_ByEig()

# Specify the sparse norms
IRLS = Directives.Update_IRLS(norms=([2, 2, 2, 2]),
                              eps=(1e-3, 1e-3), f_min_change=1e-3,
                              minGNiter=3, chifact=0.25)

# Special directive specific to the mag amplitude problem. The sensitivity
# weights are update between each iteration.
update_Jacobi = Directives.Amplitude_Inv_Iter()

# Put all together
inv = Inversion.BaseInversion(invProb,
                              directiveList=[IRLS, update_Jacobi, betaest])

# Invert
mrec_MAIS = inv.run(mstart)


#%% Plot models
from matplotlib.patches import Rectangle
vmax = 0.05
contours = [0.02]

plt.figure(figsize=(15, 5))

ax3 = plt.subplot(231)
PF.Magnetics.plotModelSections(mesh, mrec_sus, normal='z', ind=-3, subFact=2, scale=0.25, xlim=[-75, 75], ylim=[-30, 30],
                      title="Esus Model", axs=ax3, vmin=0, vmax=vmax, contours = contours)
ax3.xaxis.set_visible(False)

ax2 = plt.subplot(234)
PF.Magnetics.plotModelSections(mesh, mrec_sus, normal='y', ind=midy, subFact=2, scale=0.25, xlim=[-75, 75], ylim=[-85, 5],
                      axs=ax2, vmin=0, vmax=vmax, contours = contours)
for midx in locx:
    ax2.add_patch(Rectangle((mesh.vectorCCx[midx-nX]-dx/2.,mesh.vectorCCz[midz-nX]-dx/2.),3*dx,3*dx, facecolor = 'none', edgecolor='w'))

# Draw a box

ax3 = plt.subplot(232)
PF.Magnetics.plotModelSections(mesh, mrec_MAI, normal='z', ind=-3, subFact=2, scale=0.25, xlim=[-75, 75], ylim=[-30, 30],
                      title="Esus Model", axs=ax3, vmin=0, vmax=vmax, contours = contours)
ax3.xaxis.set_visible(False)

ax2 = plt.subplot(235)
PF.Magnetics.plotModelSections(mesh, mrec_MAI, normal='y', ind=midy, subFact=2, scale=0.25, xlim=[-75, 75], ylim=[-85, 5],
                      axs=ax2, vmin=0, vmax=vmax, contours = contours)
for midx in locx:
    ax2.add_patch(Rectangle((mesh.vectorCCx[midx-nX]-dx/2.,mesh.vectorCCz[midz-nX]-dx/2.),3*dx,3*dx, facecolor = 'none', edgecolor='w'))                

ax3 = plt.subplot(233)
PF.Magnetics.plotModelSections(mesh, mrec_MAIS, normal='z', ind=-3, subFact=2, scale=0.25, xlim=[-75, 75], ylim=[-30, 30],
                      title="Esus Model", axs=ax3, vmin=0, vmax=vmax, contours = contours)
ax3.xaxis.set_visible(False)

ax2 = plt.subplot(236)
PF.Magnetics.plotModelSections(mesh, mrec_MAIS, normal='y', ind=midy, subFact=2, scale=0.25, xlim=[-75, 75], ylim=[-85, 5],
                      axs=ax2, vmin=0, vmax=vmax, contours = contours)
for midx in locx:
    ax2.add_patch(Rectangle((mesh.vectorCCx[midx-nX]-dx/2.,mesh.vectorCCz[midz-nX]-dx/2.),3*dx,3*dx, facecolor = 'none', edgecolor='w'))