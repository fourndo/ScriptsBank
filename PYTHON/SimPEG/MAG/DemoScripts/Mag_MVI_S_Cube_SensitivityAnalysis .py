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
from SimPEG import mkvc
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import os

# # STEP 1: Setup and data simulation # #
out_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Two_Blocks_test'

# Magnetic inducing field parameter (A,I,D)
B = [50000, 90, 0]

# Create a mesh
dx = 5.

hxind = [(dx, 5, -1.3), (dx, 15), (dx, 5, 1.3)]
hyind = [(dx, 5, -1.3), (dx, 15), (dx, 5, 1.3)]
hzind = [(dx, 5, -1.3), (dx, 10)]

mesh = Mesh.TensorMesh([hxind, hyind, hzind], 'CC0')
mesh.x0[2] -= mesh.vectorNz[-1]

susc = 0.05
inc = [0, 0]
dec = [90, 90]
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
xr = np.linspace(-40., 40., 20)
yr = np.linspace(-30., 30., 20)
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

I = np.zeros((mesh.nCx, mesh.nCy, mesh.nCz))
count = -1
for midx in locx:
    count += 1
    I[(midx-nX):(midx+nX+1), (midy-nX):(midy+nX+1), (midz-nX):(midz+nX+1)] = inc[count]
I = Utils.mkvc(I)
I = I[actv]

D = np.zeros((mesh.nCx, mesh.nCy, mesh.nCz))
count = -1
for midx in locx:
    count += 1
    D[(midx-nX):(midx+nX+1), (midy-nX):(midy+nX+1), (midz-nX):(midz+nX+1)] = dec[count]
D = Utils.mkvc(D)
D = D[actv]

M = PF.Magnetics.dipazm_2_xyz(I, D)
#M = PF.Magnetics.dipazm_2_xyz(np.ones(nC) * B[1], np.ones(nC) * B[2])
m = mkvc(sp.diags(model, 0) * M)

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

#%% # RUN MVI INVERSION
# Create active map to go from reduce set to full
# Creat reduced identity map
idenMap = Maps.IdentityMap(nP=3*nC)

# Create the forward model operator
prob = PF.Magnetics.MagneticVector(mesh, chiMap=idenMap,
                                     actInd=actv)


survey.pair(prob)

# Create sensitivity weights from our linear forward operator
wr = np.sum(prob.G**2., axis=0)**0.5
wr = (wr/np.max(wr))

# Create a regularization
reg = Regularization.Sparse(mesh, indActive=actv, mapping=idenMap, nSpace=3)
reg.cell_weights = wr
reg.mref = np.zeros(3*nC)

# Data misfit function
dmis = DataMisfit.l2_DataMisfit(survey)
dmis.Wd = 1/wd

# Add directives to the inversion
opt = Optimization.ProjectedGNCG(maxIter=30,lower=-10.,upper=10., maxIterCG= 20, tolCG = 1e-3)


invProb = InvProblem.BaseInvProblem(dmis, reg, opt, beta=3.5e+7)
betaest = Directives.BetaEstimate_ByEig()

# Here is where the norms are applied
IRLS = Directives.Update_IRLS(norms=([2, 2, 2, 2]),
                              eps=None, f_min_change=1e-2,
                              minGNiter=1, beta_tol=1e-2,
                              chifact = 10.)

update_Jacobi = Directives.Update_lin_PreCond()

inv = Inversion.BaseInversion(invProb,
                              directiveList=[IRLS, update_Jacobi])
   
mrec_MVI = inv.run(np.ones(3*len(actv))*1e-4)

beta = invProb.beta

#%% # RUN MVI-S 

# # STEP 3: Finish inversion with spherical formulation
mstart = PF.Magnetics.xyz2atp(mrec_MVI)
prob.ptype = 'Spherical'
prob.chi = mstart

# Create a block diagonal regularization
reg = Regularization.Sparse(mesh, indActive=actv, mapping=idenMap,
                            nSpace=3)
reg.mref = np.zeros(3*nC)
reg.cell_weights = np.ones(3*nC)
reg.alpha_s = [1., 0., 0.]
reg.mspace = ['lin', 'sph', 'sph']

# Data misfit function
dmis = DataMisfit.l2_DataMisfit(survey)
dmis.Wd = 1./survey.std

# Add directives to the inversion
opt = Optimization.ProjectedGNCG_nSpace(maxIter=30, lower=[0., -np.pi/2.,-np.pi],
                                        upper=[10., np.pi/2., np.pi], maxIterLS=2,
                                        maxIterCG=40, tolCG=1e-3, LSreduction=1e-1,
                                        ptype=['lin', 'sph', 'sph'], nSpace=3)

invProb = InvProblem.BaseInvProblem(dmis, reg, opt, beta=beta)
#betaest = Directives.BetaEstimate_ByEig()

# Here is where the norms are applied
IRLS = Directives.Update_IRLS(norms=([0, 2, 2, 2]),
                              f_min_change=1e-4,
                              minGNiter=3, beta_tol=1e-2,
                              coolingRate=3)
IRLS.eps = [[2e-3,5e-3],[2e-3,1e-2],[1e-4, 1e-2]]

# Special directive specific to the mag amplitude problem. The sensitivity
# weights are update between each iteration.
update_Jacobi = Directives.Amplitude_Inv_Iter()
#update_Jacobi.test = True
update_Jacobi.ptype = 'MVI-S'

inv = Inversion.BaseInversion(invProb,
                              directiveList=[IRLS, update_Jacobi, ])

mrec_MVI_S = inv.run(mstart)
mrec_MVI_S = PF.Magnetics.atp2xyz(mrec_MVI_S)

Mesh.TensorMesh.writeVectorUBC(mesh,out_dir + '\Vector_CAR.sus',mrec_MVI.reshape(mesh.nC,3,order='F'))

Mesh.TensorMesh.writeVectorUBC(mesh,out_dir + '\Vector_SPH.sus',mrec_MVI_S.reshape(mesh.nC,3,order='F'))
#%% Plot models
from matplotlib.patches import Rectangle

contours = [0.02]

plt.figure(figsize=(16, 8))
ax1 = plt.subplot(231)
ax2 = plt.subplot(234)
ax3 = plt.subplot(232)
ax4 = plt.subplot(235)
ax5 = plt.subplot(233)
ax6 = plt.subplot(236)

ypanel = int(mesh.nCy/2)
zpanel = -4

vmin = 0.
xlim = [-60, 60]
cmap = 'magma_r'

# # FIRST MODEL # #
vmax = model.max()
ax2, im2, cbar = PF.Magnetics.plotModelSections(mesh, m, normal='y',
                               ind=ypanel, axs=ax2,
                               xlim=xlim, scale = 0.75, vec ='w',
                               ylim=(mesh.vectorNz[3], mesh.vectorNz[-1]+dx),
                               vmin=vmin, vmax=vmax, cmap = cmap)

for midx in locx:
    ax2.add_patch(Rectangle((mesh.vectorCCx[midx-nX-1]+dx/2.,mesh.vectorCCz[midz-nX-1]+dx/2.),3*dx,3*dx, lw=2, facecolor = 'none', edgecolor='r'))
ax2.set_title('(a) True')

ax1, im2, cbar = PF.Magnetics.plotModelSections(mesh, m, normal='z',
                               ind=zpanel, axs=ax1,
                               xlim=xlim, scale = 0.75, vec ='w',
                               ylim=xlim,
                               vmin=vmin, vmax=vmax, cmap = cmap)

for midx in locx:
    ax1.add_patch(Rectangle((mesh.vectorCCx[midx-nX-1]+dx/2.,mesh.vectorCCy[midy-nX-1]+dx/2.),3*dx,3*dx, lw=2, facecolor = 'none', edgecolor='r'))
ax1.set_title('(a) True')

ax1.xaxis.set_visible(False)

# # SECOND MODEL # #
vmax = mrec_MVI.max()
scale = mrec_MVI.max()/m.max()*0.75

ax3, im2, cbar = PF.Magnetics.plotModelSections(mesh, mrec_MVI, normal='z',
                               ind=zpanel, axs=ax3,
                               xlim=xlim, scale = scale, vec ='w',
                               ylim=xlim,
                               vmin=vmin, vmax=vmax, cmap = cmap)

for midx in locx:
    ax3.add_patch(Rectangle((mesh.vectorCCx[midx-nX-1]+dx/2.,mesh.vectorCCy[midy-nX-1]+dx/2.),3*dx,3*dx, lw=2, facecolor = 'none', edgecolor='r'))
ax3.set_title('(a) True')
ax3.xaxis.set_visible(False)


ax4, im2, cbar = PF.Magnetics.plotModelSections(mesh, mrec_MVI, normal='y',
                               ind=ypanel, axs=ax4,
                               xlim=xlim, scale = scale, vec ='w',
                               ylim=(mesh.vectorNz[3], mesh.vectorNz[-1]+dx),
                               vmin=vmin, vmax=vmax, cmap = cmap)

for midx in locx:
    ax4.add_patch(Rectangle((mesh.vectorCCx[midx-nX-1]+dx/2.,mesh.vectorCCz[midz-nX-1]+dx/2.),3*dx,3*dx, lw=2, facecolor = 'none', edgecolor='r'))
ax4.set_title('(a) True')
ax4.xaxis.set_visible(False)

# # THIRD MODEL # #
vmax = mrec_MVI_S.max()
scale = mrec_MVI_S.max()/m.max()*0.75

ax5, im2, cbar = PF.Magnetics.plotModelSections(mesh, mrec_MVI_S, normal='z',
                               ind=zpanel, axs=ax5,
                               xlim=xlim, scale = scale, vec ='w',
                               ylim=xlim,
                               vmin=vmin, vmax=vmax, cmap = cmap)

for midx in locx:
    ax5.add_patch(Rectangle((mesh.vectorCCx[midx-nX-1]+dx/2.,mesh.vectorCCy[midy-nX-1]+dx/2.),3*dx,3*dx, lw=2, facecolor = 'none', edgecolor='r'))
ax5.set_title('(a) True')
ax5.xaxis.set_visible(False)

ax6, im2, cbar = PF.Magnetics.plotModelSections(mesh, mrec_MVI_S, normal='y',
                               ind=ypanel, axs=ax6,
                               xlim=xlim, scale = scale, vec ='w',
                               ylim=(mesh.vectorNz[3], mesh.vectorNz[-1]+dx),
                               vmin=vmin, vmax=vmax, cmap = cmap)

for midx in locx:
    ax6.add_patch(Rectangle((mesh.vectorCCx[midx-nX-1]+dx/2.,mesh.vectorCCz[midz-nX-1]+dx/2.),3*dx,3*dx, lw=2, facecolor = 'none', edgecolor='r'))
ax6.set_title('(a) True')
ax6.xaxis.set_visible(False)