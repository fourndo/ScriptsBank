#%%
from SimPEG import *
import simpegPF as PF
import pylab as plt
import os
home_dir = 'C:\\LC\\Private\\dominiquef\\Projects\\4414_Minsim\\Modeling\\GRAV\\SimPEG_Compact'

#inpfile = 'PYGRAV3D_inv.inp'
dsep = '\\'
dsep = os.path.sep
plt.close('all')

#%% User input
# Plotting parameter
vmin = -0.3
vmax = 0.3
#%%
# Read input file
#[mshfile, obsfile, topofile, mstart, mref, wgtfile, chi, alphas, bounds, lpnorms] = PF.Gravity.read_GRAVinv_inp(home_dir + dsep + inpfile)
driver = PF.GravityDriver.GravityDriver_Inv(home_dir + dsep + 'PYGRAV3D_inv.inp')
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


# Get index of the center
midx = int(mesh.nCx/2)
midy = int(mesh.nCy/2)


#%% Plot obs data
PF.Gravity.plot_obs_2D(survey,'Observed Data')

#%% Run inversion
prob = PF.Gravity.GravityIntegral(mesh, mapping = staticCells, actInd = actv)
prob.solverOpts['accuracyTol'] = 1e-4

survey.pair(prob)

# Write out the predicted file and generate the forward operator
pred = prob.fields(mstart)

PF.Gravity.writeUBCobs(home_dir + dsep + 'Pred0.dat',survey,pred)

# Load weighting  file
if driver.wgtfile == 'DEFAULT':  
    wr = PF.Magnetics.get_dist_wgt(mesh, rxLoc, actv, 3., np.min(mesh.hx)/4.)
    wr = wr**2.
    
    # Make depth weighting
    #wr = np.sum(prob.G**2.,axis=0)**0.5 / mesh.vol[actv]
    #wr = ( wr/np.max(wr) )
    #wr_out = actvMap * wr
    
else:
    wr = Mesh.TensorMesh.readModelUBC(mesh, home_dir + dsep + wgtfile)   
    wr = wr[actv]
    wr = wr**2.

#%% Plot depth weighting
#plt.figure()
#ax = plt.subplot()
#mesh.plotSlice(actvMap*wr, ax = ax, normal = 'Y', ind=midx+1 ,clim = (0, wr.max()))
#plt.title('Distance weighting')
#plt.xlabel('x');plt.ylabel('z')
#plt.gca().set_aspect('equal', adjustable='box')

#%% Create inversion objects

reg = Regularization.Sparse(mesh, indActive=actv, mapping=staticCells)
reg.mref = driver.mref[dynamic]
reg.cell_weights = wr*mesh.vol[actv]
    

opt = Optimization.ProjectedGNCG(maxIter=100 ,lower=driver.bounds[0],upper=driver.bounds[1], maxIterLS = 20, maxIterCG= 10, tolCG = 1e-3)
dmis = DataMisfit.l2_DataMisfit(survey)
dmis.Wd = 1./wd
invProb = InvProblem.BaseInvProblem(dmis, reg, opt)
#beta = Directives.BetaSchedule(coolingFactor=1, coolingRate=1)
#update_beta = Directives.Scale_Beta(tol = 0.05, coolingRate=5)
betaest = Directives.BetaEstimate_ByEig()
IRLS = Directives.Update_IRLS( norms=driver.lpnorms,  eps=driver.eps, f_min_change = 1e-4, minGNiter=3)
update_Jacobi = Directives.Update_lin_PreCond()
inv = Inversion.BaseInversion(invProb, directiveList=[IRLS,betaest,update_Jacobi])

# Run inversion
mrec = inv.run(mstart)

# Plot predicted
pred = prob.fields(mrec)
#PF.Magnetics.plot_obs_2D(rxLoc,pred,wd,'Predicted Data')
#PF.Magnetics.plot_obs_2D(rxLoc,(d-pred),wd,'Residual Data')
survey.dobs = pred
PF.Gravity.plot_obs_2D(survey,'Observed Data')
print "Final misfit:" + str(np.sum( ((d-pred)/wd)**2. ) ) 

#%% Plot out a section of the model

yslice = midx

m_out = actvMap*staticCells*reg.l2model

# Write result
Mesh.TensorMesh.writeModelUBC(mesh,home_dir + dsep + 'SimPEG_inv_l2l2.den',m_out)

# Nan aircells for plotting
m_out[m_out==-100] = np.nan

plt.figure()
ax = plt.subplot(221)
mesh.plotSlice(m_out, ax = ax, normal = 'Z', ind=-10, clim = (mrec.min(), mrec.max()))
plt.plot(np.array([mesh.vectorCCx[0],mesh.vectorCCx[-1]]), np.array([mesh.vectorCCy[yslice],mesh.vectorCCy[yslice]]),c='w',linestyle = '--')
plt.title('Z: ' + str(mesh.vectorCCz[-5]) + ' m')
plt.xlabel('x');plt.ylabel('z')
plt.gca().set_aspect('equal', adjustable='box')

ax = plt.subplot(222)
mesh.plotSlice(m_out, ax = ax, normal = 'Z', ind=-15, clim = (mrec.min(), mrec.max()))
plt.plot(np.array([mesh.vectorCCx[0],mesh.vectorCCx[-1]]), np.array([mesh.vectorCCy[yslice],mesh.vectorCCy[yslice]]),c='w',linestyle = '--')
plt.title('Z: ' + str(mesh.vectorCCz[-15]) + ' m')
plt.xlabel('x');plt.ylabel('z')
plt.gca().set_aspect('equal', adjustable='box')


ax = plt.subplot(212)
mesh.plotSlice(m_out, ax = ax, normal = 'Y', ind=yslice, clim = (mrec.min(), mrec.max()))
plt.title('Cross Section')
plt.xlabel('x');plt.ylabel('z')
plt.gca().set_aspect('equal', adjustable='box')

plt.figure()
ax = plt.subplot(121)
plt.hist(reg.l2model,100)
plt.yscale('log', nonposy='clip')
plt.title('Histogram of model values - Smooth')
ax = plt.subplot(122)
plt.hist(reg.regmesh.cellDiffxStencil*(staticCells*reg.l2model),100)
plt.yscale('log', nonposy='clip')
plt.title('Histogram of model gradient values - Smooth')

#%% Plot out a section of the model

yslice = midx

m_out = actvMap*staticCells*mrec
# Write result
Mesh.TensorMesh.writeModelUBC(mesh,home_dir + dsep + 'SimPEG_inv_lplq.den',m_out)

m_out[m_out==-100] = np.nan

plt.figure()
ax = plt.subplot(221)
mesh.plotSlice(m_out, ax = ax, normal = 'Z', ind=-10, clim = (mrec.min(), mrec.max()))
plt.plot(np.array([mesh.vectorCCx[0],mesh.vectorCCx[-1]]), np.array([mesh.vectorCCy[yslice],mesh.vectorCCy[yslice]]),c='w',linestyle = '--')
plt.title('Z: ' + str(mesh.vectorCCz[-5]) + ' m')
plt.xlabel('x');plt.ylabel('z')
plt.gca().set_aspect('equal', adjustable='box')

ax = plt.subplot(222)
mesh.plotSlice(m_out, ax = ax, normal = 'Z', ind=-15, clim = (mrec.min(), mrec.max()))
plt.plot(np.array([mesh.vectorCCx[0],mesh.vectorCCx[-1]]), np.array([mesh.vectorCCy[yslice],mesh.vectorCCy[yslice]]),c='w',linestyle = '--')
plt.title('Z: ' + str(mesh.vectorCCz[-8]) + ' m')
plt.xlabel('x');plt.ylabel('z')
plt.gca().set_aspect('equal', adjustable='box')

ax = plt.subplot(212)
mesh.plotSlice(m_out, ax = ax, normal = 'Y', ind=yslice, clim = (mrec.min(), mrec.max()))
plt.title('Cross Section')
plt.xlabel('x');plt.ylabel('z')
plt.gca().set_aspect('equal', adjustable='box')

plt.figure()
ax = plt.subplot(121)
plt.hist(mrec,100)
plt.yscale('log', nonposy='clip')
plt.title('Histogram of model values - Compact')
ax = plt.subplot(122)
plt.hist(reg.regmesh.cellDiffxStencil*(staticCells*mrec),100)
plt.yscale('log', nonposy='clip')
plt.title('Histogram of model gradient values - Compact')

plt.show()