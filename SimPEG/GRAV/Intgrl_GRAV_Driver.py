#%%
from SimPEG import Mesh, Directives, Maps, InvProblem, Optimization, Utils
from SimPEG import DataMisfit, Inversion, Regularization
import SimPEG.PF as PF
import pylab as plt
import os
import numpy as np

#work_dir = 'C:\\Users\\DominiqueFournier\\ownCloud\\Research\\Modelling\\Synthetic\\SingleBlock\\GRAV'
work_dir = 'C:\\Users\\DominiqueFournier\\Downloads'
inpfile = 'SimPEG_GRAV.inp'
out_dir = "SimPEG_PF_Inv\\"
dsep = '\\'
dsep = os.path.sep
plt.close('all')

os.system('mkdir ' + work_dir + dsep + out_dir)


# Read input file
driver = PF.GravityDriver.GravityDriver_Inv(work_dir + dsep + inpfile)
mesh = driver.mesh
survey = driver.survey

rxLoc = survey.srcField.rxList[0].locs
d = survey.dobs
wd = survey.std

ndata = survey.srcField.rxList[0].locs.shape[0]

actv = driver.activeCells
nC = len(actv)

# Create active map to go from reduce set to full
actvMap = Maps.InjectActiveCells(mesh, actv, -100)

# Create static map
static = driver.staticCells
dynamic = driver.dynamicCells

staticCells = Maps.InjectActiveCells(None, dynamic, driver.m0[static], nC=nC)
mstart = driver.m0[dynamic]


#%% Plot obs data
# PF.Gravity.plot_obs_2D(survey.srcField.rxList[0].locs, survey.dobs,'Observed Data')

#%% Run inversion
prob = PF.Gravity.GravityIntegral(mesh, rhoMap=staticCells, actInd=actv)
prob.solverOpts['accuracyTol'] = 1e-4

survey.pair(prob)

# Write out the predicted file and generate the forward operator
pred = prob.fields(mstart)

PF.Gravity.writeUBCobs(work_dir + dsep + 'Pred0.dat',survey,pred)

# Load weighting  file
if driver.wgtfile is None:
    # wr = PF.Magnetics.get_dist_wgt(mesh, rxLoc, actv, 3., np.min(mesh.hx)/4.)
    # wr = wr**2.

    # Make depth weighting
    wr = np.sum(prob.F**2., axis=0)**0.5
    wr = (wr/np.max(wr))
    # wr_out = actvMap * wr

else:
    wr = Mesh.TensorMesh.readModelUBC(mesh, work_dir + dsep + wgtfile)
    wr = wr[actv]
    wr = wr**2.

# % Create inversion objects
reg = Regularization.Sparse(mesh, indActive=actv, mapping=staticCells)
reg.mref = driver.mref[dynamic]
reg.cell_weights = wr
reg.norms = driver.lpnorms
if driver.eps is not None:
    reg_a.eps_p = driver.eps[0]
    reg_a.eps_q = driver.eps[1]


opt = Optimization.ProjectedGNCG(maxIter=100, lower=driver.bounds[0],upper=driver.bounds[1], maxIterLS = 20, maxIterCG= 10, tolCG = 1e-3)
dmis = DataMisfit.l2_DataMisfit(survey)
dmis.W = 1./wd
invProb = InvProblem.BaseInvProblem(dmis, reg, opt)

betaest = Directives.BetaEstimate_ByEig()
IRLS = Directives.Update_IRLS(f_min_change=1e-4, minGNiter=3)
update_Jacobi = Directives.UpdatePreCond()
saveModel = Directives.SaveUBCModelEveryIteration(mapping=actvMap)
saveModel.fileName = work_dir + dsep + out_dir + 'GRAV'
inv = Inversion.BaseInversion(invProb, directiveList=[betaest, IRLS,
                                                      update_Jacobi, saveModel])
# Run inversion
mrec = inv.run(mstart)

# Plot predicted
pred = prob.fields(mrec)

survey.dobs = pred
# PF.Gravity.plot_obs_2D(survey, 'Observed Data')
print("Final misfit:" + str(np.sum(((d-pred)/wd)**2.)))

m_out = actvMap*staticCells*invProb.l2model

# Write result
Mesh.TensorMesh.writeModelUBC(mesh, work_dir + dsep + out_dir + 'SimPEG_inv_l2l2.den',m_out)

m_out = actvMap*staticCells*mrec
# Write result
Mesh.TensorMesh.writeModelUBC(mesh, work_dir + dsep + out_dir + 'SimPEG_inv_lplq.den',m_out)
