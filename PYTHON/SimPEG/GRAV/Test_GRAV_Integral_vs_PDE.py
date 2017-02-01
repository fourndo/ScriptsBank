# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 09:47:38 2016


@author: dominiquef
"""

from SimPEG import Mesh, Utils, np, PF, Maps, Problem, Survey, mkvc
import matplotlib.pyplot as plt
from pymatsolver import PardisoSolver
import time

# Define inducing field and sphere parameters
H0 = (50000., 60., 270.)
rad = 2.
rho = 0.1

ntrials = 2

# Initialize metrics
clock = np.zeros((2,ntrials))
psize = np.zeros(ntrials)
mresi = np.zeros(ntrials)

for ii in range(ntrials):
    
    # Define a mesh
    cs = 0.2
    nC = 11 + ii*5
    
    hxind = [(cs,10,-1.3),(cs, nC),(cs,10,1.3)]
    hyind = [(cs,10,-1.3),(cs, nC),(cs,10,1.3)]
    hzind = [(cs,10,-1.3),(cs, nC),(cs,10,1.3)]
    mesh = Mesh.TensorMesh([hxind, hyind, hzind], 'CCC')
    
    # Get cells inside the sphere
    sph_ind = PF.MagAnalytics.spheremodel(mesh, 0., 0., 0., rad)
    
    # Adjust susceptibility for volume difference
    Vratio = (4./3.*np.pi*rad**3.) / (np.sum(sph_ind)*cs**3.)
    model = np.zeros(mesh.nC)
    model[sph_ind] = rho*Vratio
    m = model[sph_ind]
    
    # Creat reduced identity map for Linear Pproblem
    idenMap = Maps.IdentityMap(nP=int(sum(sph_ind)))
    
    # Create plane of observations
    xr = np.linspace(-cs*nC/2., cs*nC/2., nC)
    yr = np.linspace(-cs*nC/2., cs*nC/2., nC)
    X, Y = np.meshgrid(xr, yr)
    
    # Move obs plane 2 radius away from sphere
    Z = np.ones((xr.size, yr.size))*2.*rad
    locXyz = np.c_[Utils.mkvc(X), Utils.mkvc(Y), Utils.mkvc(Z)]
    rxLoc = PF.BaseGrav.RxObs(locXyz)
    srcField = PF.BaseGrav.SrcField([rxLoc])
    survey = PF.BaseGrav.LinearSurvey(srcField)
    
    prob_z = PF.Gravity.GravityIntegral(mesh, rhoMap=idenMap,
                                                  actInd=sph_ind,
                                                  forwardOnly=True,
                                                  rtype='z')
    
    
    # Compute 3-component mag data
    survey.pair(prob_z)
    
    start_time = time.time()
    
    d = prob_z.fields(m)
    
    print("Solve Integral --- %s seconds ---" % (time.time() - start_time))
    
    clock[0,ii] = time.time() - start_time
    
    ndata = locXyz.shape[0]
    dbx = d[0:ndata]
    dby = d[ndata:2*ndata]
    dbz = d[2*ndata:]
    
    # Compute analytical response from a magnetized sphere
    bxa, bya, bza = PF.GravAnalytics.GravSphereFreeSpace(locXyz[:, 0],
                                                         locXyz[:, 1],
                                                         locXyz[:, 2],
                                                         rad, 0, 0, 0,
                                                         rho)
    
    # Projection matrix
    err_xyz = (np.linalg.norm(d-bza) /
               np.linalg.norm(bza))
    
    #fig = plt.figure()
    #axs = plt.subplot(111)
    #mesh.plotSlice(model, normal='Z', ind = mesh.nCz/2, ax= axs)
    #axs.set_aspect('equal')
    
    #PF.Magnetics.plot_obs_2D(locXyz,dtmi)
    
    #%% Repeat using PDE solve
    m_rho = PF.BaseGrav.BaseGravMap(mesh)
    prob = PF.Gravity.Problem3D_Diff(mesh, rhoMap=m_rho)
    
    rxLoc = PF.BaseGrav.RxObs(locXyz)
    srcField = PF.BaseGrav.SrcField([rxLoc])
    survey = PF.BaseGrav.LinearSurvey(srcField)
    
    start_time = time.time()
    
    prob.pair(survey)
    u = prob.fields(model)
    
    print("Solve PDE --- %s seconds ---" % (time.time() - start_time))
    
    gg = survey.projectFields(u)
    gz = gg['gz']
    
    clock[1,ii] = time.time() - start_time
    psize[ii] = mesh.nC * survey.nD
    mresi[ii] = np.linalg.norm(gz-d)

fig, ax1, ax2, ax3 = plt.figure(), plt.subplot(1,3,1), plt.subplot(1,3,2), plt.subplot(1,3,3)

PF.Gravity.plot_obs_2D(locXyz,d, title='Integral Solution',axs=ax1)
PF.Gravity.plot_obs_2D(locXyz,gz, title='PDE Solution',axs=ax2)
PF.Gravity.plot_obs_2D(locXyz,gz-d,title='Residual',axs=ax3)

fig, ax1, ax2, ax3 = plt.figure(), plt.subplot(3,3,1), plt.subplot(3,3,2), plt.subplot(3,3,3)
ax5, ax6, ax9 =  plt.subplot(3,3,5), plt.subplot(3,3,6), plt.subplot(3,3,9)

PF.Gravity.plot_obs_2D(locXyz,gg['gxx'], title='PDE gxx',axs=ax1)
PF.Gravity.plot_obs_2D(locXyz,gg['gxy'], title='PDE gxy',axs=ax2)
PF.Gravity.plot_obs_2D(locXyz,gg['gxz'], title='PDE gxz',axs=ax3)
PF.Gravity.plot_obs_2D(locXyz,gg['gyy'], title='PDE gyy',axs=ax5)
PF.Gravity.plot_obs_2D(locXyz,gg['gyz'], title='PDE gyz',axs=ax6)
PF.Gravity.plot_obs_2D(locXyz,gg['gzz'], title='PDE gzz',axs=ax9)

fig, axs = plt.figure(), plt.subplot(1,1,1)
plt.semilogx(psize,clock[0,:])
plt.semilogx(psize,clock[1,:])
axs.set_ylabel('Solve time (sec)')
axs.set_xlabel('Problem size (nC * nD)')
plt.legend(['Integral','PDE'])
plt.grid()
#Mesh.TensorMesh.writeUBC(mesh,'MeshGrav.msh')
#Mesh.TensorMesh.writeModelUBC(mesh,'MeshGrav.den',model)
#PF.Gravity.writeUBCobs("Obs.grv",survey,d)