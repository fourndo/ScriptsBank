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

#work_dir = "C:\\Users\\DominiqueFournier\\ownCloud\\Research\\Kevitsa\\Modeling\\MAG\\"
#work_dir = "C:\\Users\\DominiqueFournier\\ownCloud\\\Research\\Modelling\\Synthetic\\Block_Gaussian_topo\\"
work_dir = "C:\\Users\\DominiqueFournier\\ownCloud\\Research\\Modelling\\Synthetic\\SingleBlock\Simpeg\\"
#work_dir = "C:\\Egnyte\\Private\\dominiquef\\Projects\\4559_CuMtn_ZTEM\\Modeling\\MAG\\A1_Fenton\\"
out_dir = "SimPEG_PF_Inv\\"
input_file = "SimPEG_MAG.inp"
# %%
# Read in the input file which included all parameters at once
# (mesh, topo, model, survey, inv param, etc.)
driver = PF.MagneticsDriver.MagneticsDriver_Inv(work_dir + input_file)

os.system('if not exist ' + work_dir + out_dir + ' mkdir ' + work_dir+out_dir)

# Access the mesh and survey information
mesh = driver.mesh
survey = driver.survey
actv = driver.activeCells


nC = len(actv)


# Create active map to go from reduce set to full
actvMap = Maps.InjectActiveCells(mesh, actv, -100)

# Creat reduced identity map
idenMap = Maps.IdentityMap(nP=nC)

# Create the forward model operator
prob = PF.Magnetics.MagneticIntegral(mesh, chiMap=idenMap, actInd=actv)

# Pair the survey and problem
survey.pair(prob)

# Create sensitivity weights from our linear forward operator
rxLoc = survey.srcField.rxList[0].locs
wr = np.zeros(prob.F.shape[1])
for ii in range(survey.nD):
    wr += (prob.F[ii, :]/survey.std[ii])**2.
wr = wr**0.5
wr = (wr/np.max(wr))


# Create a regularization
reg = Regularization.Sparse(mesh, indActive=actv, mapping=idenMap)
reg.norms = driver.lpnorms

if driver.eps is not None:
    reg.eps_p = driver.eps[0]
    reg.eps_q = driver.eps[1]

reg.cell_weights = wr
reg.mref = driver.mref
# Data misfit function
dmis = DataMisfit.l2_DataMisfit(survey)
dmis.W = 1./survey.std

# Add directives to the inversion
opt = Optimization.ProjectedGNCG(maxIter=20, lower=0., upper=10.,
                                 maxIterLS=20, maxIterCG=10, tolCG=1e-4)
invProb = InvProblem.BaseInvProblem(dmis, reg, opt)
betaest = Directives.BetaEstimate_ByEig()

# Here is where the norms are applied
# Use pick a treshold parameter empirically based on the distribution of
#  model parameters
IRLS = Directives.Update_IRLS(f_min_change=1e-3, minGNiter=3, maxIRLSiter=10)
update_Jacobi = Directives.Update_lin_PreCond()

saveModel = Directives.SaveUBCModelEveryIteration(mapping=actvMap)
saveModel.fileName = work_dir + out_dir + 'ModelSus'

inv = Inversion.BaseInversion(invProb,
                              directiveList=[betaest, IRLS, update_Jacobi,  saveModel])

# Run the inversion
m0 = driver.m0  # Starting model
mrec = inv.run(m0)

# Outputs
Mesh.TensorMesh.writeModelUBC(mesh,work_dir + out_dir + "SensWeight.sus", actvMap*(wr**0.5))
if reg.l2model is not None:
	Mesh.TensorMesh.writeModelUBC(mesh,work_dir + out_dir + "MAG_SuS_l2l2.sus", actvMap*reg.l2model)
Mesh.TensorMesh.writeModelUBC(mesh,work_dir + out_dir + "MAG_SuS_lplq.sus", actvMap*invProb.model)
PF.Magnetics.writeUBCobs(work_dir+out_dir + 'MAG_SuS_Inv.pre', survey, invProb.dpred)

#PF.Magnetics.plot_obs_2D(rxLoc,invProb.dpred,varstr='Amplitude Data')
