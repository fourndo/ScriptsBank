#%%
from SimPEG import Mesh, Directives, Maps, InvProblem, Optimization, Utils
from SimPEG import DataMisfit, Inversion, Regularization, mkvc
import SimPEG.PF as PF
import pylab as plt
import os
import numpy as np

# work_dir = 'C:\\Users\\DominiqueFournier\\ownCloud\\Research\\Modelling\\Synthetic\\SingleBlock\\GRAV\\'
work_dir = 'C:\\Users\\DominiqueFournier\\ownCloud\\Research\\Kevitsa\\Modeling\\GRAV\\'
tile_dirl2 = 'Tiles_l2\\'
tile_dirlp = 'Tiles_lp\\'
inpfile = 'SimPEG_GRAV.inp'

dsep = '\\'
dsep = os.path.sep
plt.close('all')

out_dir = "SimPEG_PF_Inv\\"

max_mcell = 5.0e+4
min_Olap = 1e+3

os.system('mkdir ' + work_dir+out_dir)
os.system('mkdir ' + work_dir+out_dir + tile_dirlp)
os.system('mkdir ' + work_dir+out_dir + tile_dirl2)


#%% User input
# Plotting parameter
#%%
# Read input file
#[mshfile, obsfile, topofile, mstart, mref, wgtfile, chi, alphas, bounds, lpnorms] = PF.Gravity.read_GRAVinv_inp(work_dir + dsep + inpfile)
driver = PF.GravityDriver.GravityDriver_Inv(work_dir + dsep + inpfile)
mesh = driver.mesh
survey = driver.survey

actv = driver.activeCells

# # TILE THE PROBLEM
lx = int(np.floor(np.sqrt(max_mcell/mesh.nCz)))
# #%% Create tiles
# # Begin tiling

# In the x-direction
ntile = 1
Olap = -1
dx = [mesh.hx.min(), mesh.hy.min()]
# Cell size


while Olap < min_Olap:

    ntile += 1

    # Set location of SW corners
    x0 = np.asarray([mesh.vectorNx[0], mesh.vectorNx[-1] - lx*dx[0]])

    dx_t = np.round((x0[1] - x0[0]) / ((ntile-1) * dx[0]))

    if ntile > 2:
        x1 = np.r_[x0[0],
                   x0[0] + np.cumsum(np.ones(ntile-2) * dx_t * dx[0]),
                   x0[1]]
    else:
        x1 = np.asarray([x0[0], x0[1]])

    x2 = x1 + lx*dx[0]

    Olap = x1[0] + lx*dx[0] - x1[1]

# Save x-corner location
xtile = np.c_[x1, x2]

# In the Y-direction
ntile = 1
Olap = -1

# Cell size

while Olap < min_Olap:

    ntile += 1

    # Set location of SW corners
    y0 = np.asarray([mesh.vectorNy[0], mesh.vectorNy[-1] - lx*dx[0]])

    dy_t = np.round((y0[1] - y0[0]) / ((ntile-1) * dx[0]))

    if ntile > 2:
        y1 = np.r_[y0[0],
                   y0[0] + np.cumsum(np.ones(ntile-2) * dy_t * dx[0]),
                   y0[1]]
    else:
        y1 = np.asarray([y0[0], y0[1]])

    y2 = y1 + lx*dx[0]

    Olap = y1[0] + lx*dx[0] - y1[1]

# Save x-corner location
ytile = np.c_[y1, y2]


X1, Y1 = np.meshgrid(x1, y1)
X1, Y1 = mkvc(X1), mkvc(Y1)
X2, Y2 = np.meshgrid(x2, y2)
X2, Y2 = mkvc(X2), mkvc(Y2)

# Plot data and tiles
rxLoc = survey.srcField.rxList[0].locs
# fig, ax1 = plt.figure(), plt.subplot()
# PF.Magnetics.plot_obs_2D(rxLoc,survey.dobs, ax=ax1)
# for ii in range(X1.shape[0]):
#     ax1.add_patch(Rectangle((X1[ii],Y1[ii]),X2[ii]-X1[ii],Y2[ii]-Y1[ii], facecolor = 'none', edgecolor='k'))
# ax1.set_aspect('equal')

# LOOP THROUGH TILES
npadxy = 8
core_x = np.ones(lx)*dx[0]
core_y = np.ones(lx)*dx[1]
expf = 1.3

os.system('if not exist ' + work_dir + out_dir + tile_dirl2 + ' mkdir ' + work_dir+out_dir+tile_dirl2)
os.system('if not exist ' + work_dir + out_dir + tile_dirlp + ' mkdir ' + work_dir+out_dir+tile_dirlp)

tiles = range(X1.shape[0])
tilestep = 3
for tt in tiles[:tilestep]:

    midX = np.mean([X1[tt], X2[tt]])
    midY = np.mean([Y1[tt], Y2[tt]])

    # Create new mesh
    padx = np.r_[dx[0]*expf**(np.asarray(range(npadxy))+1)]
    pady = np.r_[dx[1]*expf**(np.asarray(range(npadxy))+1)]

    hx = np.r_[padx[::-1], core_x, padx]
    hy = np.r_[pady[::-1], core_y, pady]
#        hz = np.r_[padb*2.,[33,26],np.ones(25)*22,[18,15,12,10,8,7,6], np.ones(18)*5,5*expf**(np.asarray(range(2*npad)))]
    hz = mesh.hz

    mesh_t = Mesh.TensorMesh([hx, hy, hz], 'CC0')

#    mtemp._x0 = [x0[ii]-np.sum(padb), y0[ii]-np.sum(padb), mesh.x0[2]]

    mesh_t._x0 = (mesh_t.x0[0] + midX, mesh_t.x0[1]+midY, mesh.x0[2])

    Mesh.TensorMesh.writeUBC(mesh_t, work_dir + out_dir + tile_dirl2 + "MVI_S_Tile" + str(tt) + ".msh")
    Mesh.TensorMesh.writeUBC(mesh_t, work_dir + out_dir + tile_dirlp + "MVI_S_Tile" + str(tt) + ".msh")

#        meshes.append(mtemp)
    # Grab the right data
    xlim = [mesh_t.vectorCCx[npadxy], mesh_t.vectorCCx[-npadxy]]
    ylim = [mesh_t.vectorCCy[npadxy], mesh_t.vectorCCy[-npadxy]]

    ind_t = np.all([rxLoc[:, 0] > xlim[0], rxLoc[:, 0] < xlim[1],
                    rxLoc[:, 1] > ylim[0], rxLoc[:, 1] < ylim[1]], axis=0)

    if np.sum(ind_t) < 20:
        continue

    rxLoc_t = PF.BaseGrav.RxObs(rxLoc[ind_t, :])
    srcField = PF.BaseGrav.SrcField([rxLoc_t])
    survey_t = PF.BaseGrav.LinearSurvey(srcField)
    survey_t.dobs = survey.dobs[ind_t]
    survey_t.std = survey.std[ind_t]

    # Extract model from global to local mesh
    if driver.topofile is not None:
        topo = np.genfromtxt(work_dir + driver.topofile,
                             skip_header=1)
        actv = Utils.surface2ind_topo(mesh_t, topo, 'N')
        actv = np.asarray(np.where(mkvc(actv))[0], dtype=int)
    else:
        actv = np.ones(mesh_t.nC, dtype='bool')

    nC = len(actv)
    print("Tile "+str(tt))
    print(nC, np.sum(ind_t))
    # Create active map to go from reduce space to full
    actvMap = Maps.InjectActiveCells(mesh_t, actv, 0)

    # Create identity map
    idenMap = Maps.IdentityMap(nP=3*nC)

    mstart = np.ones(3*len(actv))*1e-4
    #%% Run inversion
    prob = PF.Gravity.GravityIntegral(mesh_t, rhoMap=idenMap, actInd=actv)
    prob.solverOpts['accuracyTol'] = 1e-4

    survey.pair(prob)

    # Write out the predicted file and generate the forward operator
    pred = prob.fields(mstart)

    PF.Gravity.writeUBCobs(work_dir + dsep + 'Pred0.dat', survey, pred)

    wr = np.sum(prob.F**2., axis=0)**0.5
    wr = (wr/np.max(wr))

    #%% Create inversion objects
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
    saveModel.fileName = work_dir + out_dir + 'SimPEG_GRAV'
    inv = Inversion.BaseInversion(invProb, directiveList=[betaest, IRLS,
                                                          update_Jacobi,
                                                          saveModel])

    # Run inversion
    mrec = inv.run(mstart)

    # Plot predicted
    pred = prob.fields(mrec)
    #PF.Magnetics.plot_obs_2D(rxLoc,pred,wd,'Predicted Data')
    #PF.Magnetics.plot_obs_2D(rxLoc,(d-pred),wd,'Residual Data')
    survey.dobs = pred
    # PF.Gravity.plot_obs_2D(survey, 'Observed Data')
    print("Final misfit:" + str(np.sum(((d-pred)/wd)**2.)))

    #%% Plot out a section of the model

    yslice = midx

    m_out = actvMap*staticCells*invProb.l2model

    # Write result
    Mesh.TensorMesh.writeModelUBC(mesh, work_dir + dsep + 'SimPEG_inv_l2l2.den', m_out)

    m_out = actvMap*staticCells*mrec
    # Write result
    Mesh.TensorMesh.writeModelUBC(mesh, work_dir + dsep + 'SimPEG_inv_lplq.den', m_out)


# Nan aircells for plotting
# m_out[m_out==-100] = np.nan

# plt.figure()
# ax = plt.subplot(221)
# mesh.plotSlice(m_out, ax = ax, normal = 'Z', ind=-10, clim = (mrec.min(), mrec.max()))
# plt.plot(np.array([mesh.vectorCCx[0],mesh.vectorCCx[-1]]), np.array([mesh.vectorCCy[yslice],mesh.vectorCCy[yslice]]),c='w',linestyle = '--')
# plt.title('Z: ' + str(mesh.vectorCCz[-5]) + ' m')
# plt.xlabel('x');plt.ylabel('z')
# plt.gca().set_aspect('equal', adjustable='box')

# ax = plt.subplot(222)
# mesh.plotSlice(m_out, ax = ax, normal = 'Z', ind=-15, clim = (mrec.min(), mrec.max()))
# plt.plot(np.array([mesh.vectorCCx[0],mesh.vectorCCx[-1]]), np.array([mesh.vectorCCy[yslice],mesh.vectorCCy[yslice]]),c='w',linestyle = '--')
# plt.title('Z: ' + str(mesh.vectorCCz[-15]) + ' m')
# plt.xlabel('x');plt.ylabel('z')
# plt.gca().set_aspect('equal', adjustable='box')


# ax = plt.subplot(212)
# mesh.plotSlice(m_out, ax = ax, normal = 'Y', ind=yslice, clim = (mrec.min(), mrec.max()))
# plt.title('Cross Section')
# plt.xlabel('x');plt.ylabel('z')
# plt.gca().set_aspect('equal', adjustable='box')

# plt.figure()
# ax = plt.subplot(121)
# plt.hist(reg.l2model,100)
# plt.yscale('log', nonposy='clip')
# plt.title('Histogram of model values - Smooth')
# ax = plt.subplot(122)
# plt.hist(reg.regmesh.cellDiffxStencil*(staticCells*reg.l2model),100)
# plt.yscale('log', nonposy='clip')
# plt.title('Histogram of model gradient values - Smooth')

# #%% Plot out a section of the model

# yslice = midx

