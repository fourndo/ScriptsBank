"""
This script runs a Magnetic Vector Inversion - Cartesian formulation. The code
is used to get the magnetization orientation and makes no induced assumption.
The inverse problem is three times larger than the usual susceptibility
inversion, but won't suffer the same issues in the presence of remanence. The
drawback of having 3x the number of model parameters is in non-uniqueness of
the solution. Models have a tendency to be overlly complicated.

To counter this issue, the inversion can be done cooperatively with the
Magnetic Amplitude Inverion (MAI) that uses sparse norms to recover compact
bodies. The information about the location and shape of magnetic bodies is
used by the MVI code through a cell-based wieght.

This script will run both inversions in sequence in order to compare the
results if the flag CMI is activated

-------------------------------------------------------------------------------
 !! IMPORTANT: PLease run Intgrl_MAG_AMP_Inv.py before running this script !!
-------------------------------------------------------------------------------


Created on Thu Sep 29 10:11:11 2016

@author: dominiquef
"""
from SimPEG import Mesh, Directives, Maps, InvProblem, Optimization
from SimPEG import DataMisfit, Inversion, Regularization
import SimPEG.PF as PF
import numpy as np
import os

# Define the inducing field parameter
# work_dir = "C:\\Users\\DominiqueFournier\\ownCloud\\\Research\\Modelling\\Synthetic\\Block_Gaussian_topo\\"
work_dir = "C:\\Users\\DominiqueFournier\\Dropbox\\Two_Blocks_test\\SingleBlock\\"
out_dir = "SimPEG_PF_Inv\\"
input_file = "SimPEG_MAG.inp"


# %% INPUTS
# Read in the input file which included all parameters at once
# (mesh, topo, model, survey, inv param, etc.)
driver = PF.MagneticsDriver.MagneticsDriver_Inv(work_dir + input_file)
os.system('if not exist ' +work_dir+out_dir + ' mkdir ' + work_dir+out_dir)
#%% Access the mesh and survey information
mesh = driver.mesh
survey = driver.survey

# Extract active region
actv = driver.activeCells

# Set starting mdoel
mstart = np.ones(3*len(actv))*1e-4

# Create active map to go from reduce space to full
actvMap = Maps.InjectActiveCells(mesh, actv, 0)
nC = int(len(actv))

# Create identity map
idenMap = Maps.IdentityMap(nP=3*nC)

mstart = np.ones(3*len(actv))*1e-4

# Create the forward model operator
prob = PF.Magnetics.MagneticVector(mesh, chiMap=idenMap,
                                     actInd=actv)

# Explicitely set starting model
prob.chi = mstart

# Pair the survey and problem
survey.pair(prob)


# Create sensitivity weights from our linear forward operator
wr = np.sum(prob.F**2., axis=0)**0.5
wr = (wr/np.max(wr))

# Create a block diagonal regularization
wires = Maps.Wires(('p', nC), ('s', nC), ('t', nC))

# Create a regularization
reg_p = Regularization.Sparse(mesh, indActive=actv, mapping=wires.p)
reg_p.cell_weights = wires.p * wr
reg_p.norms = [2, 2, 2, 2]

reg_s = Regularization.Sparse(mesh, indActive=actv, mapping=wires.s)
reg_s.cell_weights = wires.s * wr
reg_s.norms = [2, 2, 2, 2]

reg_t = Regularization.Sparse(mesh, indActive=actv, mapping=wires.t)
reg_t.cell_weights = wires.t * wr
reg_t.norms = [2, 2, 2, 2]

reg = reg_p + reg_s + reg_t
reg.mref = np.zeros(3*nC)

# Data misfit function
dmis = DataMisfit.l2_DataMisfit(survey)
dmis.W = 1./survey.std

# Add directives to the inversion
opt = Optimization.ProjectedGNCG(maxIter=7, lower=-10., upper=10.,
                                 maxIterCG=20, tolCG=1e-3)

invProb = InvProblem.BaseInvProblem(dmis, reg, opt)
betaest = Directives.BetaEstimate_ByEig()

# Here is where the norms are applied
IRLS = Directives.Update_IRLS(f_min_change=1e-4,
                              minGNiter=3, beta_tol=1e-2)

update_Jacobi = Directives.Update_lin_PreCond()
targetMisfit = Directives.TargetMisfit()

saveModel = Directives.SaveUBCVectorsEveryIteration(mapping=actvMap)
saveModel.fileName = work_dir + out_dir + 'MVI_C'
inv = Inversion.BaseInversion(invProb,
                              directiveList=[betaest, IRLS, update_Jacobi,
                                              saveModel])

mrec_MVI = inv.run(mstart)

beta = invProb.beta

# %% RUN MVI-S

# # STEP 3: Finish inversion with spherical formulation
mstart = PF.Magnetics.xyz2atp(mrec_MVI)
prob.ptype = 'Spherical'
prob.chi = mstart

# Create a block diagonal regularization
wires = Maps.Wires(('amp', nC), ('theta', nC), ('phi', nC))

# Create a regularization
reg_a = Regularization.Sparse(mesh, indActive=actv, mapping=wires.amp)
reg_a.norms = driver.lpnorms[:4]
if driver.eps is not None:
    reg_a.eps_p = driver.eps[0]
    reg_a.eps_q = driver.eps[1]
else:
    reg_a.eps_p = np.percentile(np.abs(mstart[:nC]), 95)


reg_t = Regularization.Sparse(mesh, indActive=actv, mapping=wires.theta)
reg_t.alpha_s = 0.
reg_t.space = 'spherical'
reg_t.norms = driver.lpnorms[4:8]
reg_t.eps_q = 2e-2
# reg_t.alpha_x, reg_t.alpha_y, reg_t.alpha_z = 0.25, 0.25, 0.25

reg_p = Regularization.Sparse(mesh, indActive=actv, mapping=wires.phi)
reg_p.alpha_s = 0.
reg_p.space = 'spherical'
reg_p.norms = driver.lpnorms[8:]
reg_p.eps_q = 1e-2

reg = reg_a + reg_t + reg_p
reg.mref = np.zeros(3*nC)

# Data misfit function
dmis = DataMisfit.l2_DataMisfit(survey)
dmis.W = 1./survey.std

# Add directives to the inversion
opt = Optimization.ProjectedGNCG(maxIter=40,
                                 lower=[0., -np.inf, -np.inf],
                                 upper=[10., np.inf, np.inf],
                                 maxIterLS=10,
                                 maxIterCG=20, tolCG=1e-3,
                                 alwaysPass=True,
                                 stepOffBoundsFact=1e-8)

invProb = InvProblem.BaseInvProblem(dmis, reg, opt, beta=beta)
#  betaest = Directives.BetaEstimate_ByEig()

# Here is where the norms are applied
IRLS = Directives.Update_IRLS(f_min_change=1e-4,
                              minGNiter=5, beta_tol=1e-2,
                              coolingRate=5)

invProb = InvProblem.BaseInvProblem(dmis, reg, opt, beta=beta)
# betaest = Directives.BetaEstimate_ByEig()


# Special directive specific to the mag amplitude problem. The sensitivity
# weights are update between each iteration.
update_Jacobi = Directives.Amplitude_Inv_Iter()
# update_Jacobi.test = True
update_Jacobi.ptype = 'MVI-S'

ProjSpherical = Directives.ProjSpherical()
saveModel = Directives.SaveUBCVectorsEveryIteration(mapping=actvMap,
                                                    saveComp=True,
                                                    spherical=True)

saveModel.fileName = work_dir+out_dir + 'MVI_S'

inv = Inversion.BaseInversion(invProb,
                              directiveList=[ProjSpherical, IRLS,
                                             update_Jacobi, saveModel])

mrec_MVI_S = inv.run(mstart)

# %%# Output the vector model and amplitude

# m_l2 = actvMap * reg.l2model[0:nC]
# m_l2[m_l2==-100] = np.nan

m_lpx = actvMap * mrec_MVI[0:nC]
m_lpy = actvMap * mrec_MVI[nC:2*nC]
m_lpz = actvMap * mrec_MVI[2*nC:]

mvec = np.c_[m_lpx, m_lpy, m_lpz]

m_lpx[m_lpx == -100] = 0

m_lpy[m_lpy == -100] = 0
m_lpz[m_lpz == -100] = 0

amp = np.sqrt(m_lpx**2. + m_lpy**2. + m_lpz**2.)

# Mesh.TensorMesh.writeModelUBC(mesh,work_dir + out_dir + "MVI_C_amp.amp",amp)
# PF.Magnetics.writeUBCobs(work_dir+out_dir + 'MVI.pre',survey,invProb.dpred)

# Mesh.TensorMesh.writeModelUBC(mesh,work_dir+out_dir + '\MVI_S_amp.sus',actvMap*mrec_MVI_S[:nC])
# Mesh.TensorMesh.writeModelUBC(mesh,work_dir+out_dir + '\MVI_S_phi.sus',actvMap*mrec_MVI_S[2*nC:])
# Mesh.TensorMesh.writeModelUBC(mesh,work_dir+out_dir + '\MVI_S_theta.sus',actvMap*mrec_MVI_S[nC:2*nC])

mrec_MVI_S_xyz = PF.Magnetics.atp2xyz(mrec_MVI_S)

# PF.MagneticsDriver.writeVectorUBC(mesh,work_dir+out_dir + '\Vector_MVI_S.fld',actvMap*(mrec_MVI_S_xyz).reshape(mesh.nC,3,order='F'))

# obs_loc = survey.srcField.rxList[0].locs

# PF.Magnetics.plot_obs_2D(obs_loc, invProb.dpred, varstr='MVI Predicted data')
