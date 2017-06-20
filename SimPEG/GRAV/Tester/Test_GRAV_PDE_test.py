# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 09:47:38 2016


@author: dominiquef
"""

from SimPEG import Mesh, Utils, np, PF, Maps, Problem, Survey, mkvc
import matplotlib.pyplot as plt
from pymatsolver import PardisoSolver
import time

# Define sphere parameters
rad = 1.
rho = 0.1

ntrials = 2

# Initialize metrics
clock = np.zeros((2,ntrials))
psize = np.zeros(ntrials)
mresi = np.zeros(ntrials)

    
# Define a mesh
cs = 0.1
nC = 20

hxind = [(cs,15,-1.3),(cs, nC),(cs,15,1.3)]
hyind = [(cs,15,-1.3),(cs, nC),(cs,15,1.3)]
hzind = [(cs,15,-1.3),(cs, 2*nC),(cs,15,1.3)]
mesh = Mesh.TensorMesh([hxind, hyind, hzind], 'CCC')
mesh.x0[2] += rad

# Get cells inside the sphere
sph_ind = PF.MagAnalytics.spheremodel(mesh, 0., 0., 0., rad)

# Adjust susceptibility for volume difference
Vratio = (4./3.*np.pi*rad**3.) / (np.sum(sph_ind)*cs**3.)
model = np.zeros(mesh.nC)
model[sph_ind] = rho*Vratio
m = model[sph_ind]

# Create plane of observations
xr = np.linspace(-cs*nC/2., cs*nC/2., 11)
yr = np.linspace(-cs*nC/2., cs*nC/2., 11)
X, Y = np.meshgrid(xr, yr)

# Move obs plane 2 radius away from sphere
Z = np.ones((xr.size, yr.size))*2.*rad
locXyz = np.c_[Utils.mkvc(X), Utils.mkvc(Y), Utils.mkvc(Z)]

ndata = locXyz.shape[0]

# Compute analytical response from a magnetized sphere
AnaSphere = PF.GravAnalytics.GravSphereFreeSpace(locXyz[:, 0],
                                                 locXyz[:, 1],
                                                 locXyz[:, 2],
                                                 rad, 0, 0, 0,
                                                 rho)


# Save analytic to file
fid = open('GRAV_gg_Sphere.obs', 'w')
np.savetxt(fid, np.c_[locXyz,AnaSphere['gxx'],AnaSphere['gxy'],
                      AnaSphere['gxz'],AnaSphere['gyy'],AnaSphere['gyz'],
                      AnaSphere['gzz']], fmt='%e',delimiter=' ',newline='\n')
fid.close()
            

#%% Compute data using PDE solve
m_rho = PF.BaseGrav.BaseGravMap(mesh)
prob = PF.Gravity.Problem3D_PDE(mesh, rhoMap=m_rho)

rxLoc = PF.BaseGrav.RxObs(locXyz)
srcField = PF.BaseGrav.SrcField([rxLoc])
survey = PF.BaseGrav.LinearSurvey(srcField)

start_time = time.time()

prob.pair(survey)
u = prob.fields(model)

print("Solve PDE --- %s seconds ---" % (time.time() - start_time))

gg = survey.projectFields(u)

gx = gg['gx']
gy = gg['gy']
gz = gg['gz']

fig, ax1, ax2, ax3 = plt.figure(figsize=(15,12)), plt.subplot(3,3,1), plt.subplot(3,3,2), plt.subplot(3,3,3)
ax4, ax5, ax6 = plt.subplot(3,3,4), plt.subplot(3,3,5), plt.subplot(3,3,6)
ax7, ax8, ax9 = plt.subplot(3,3,7), plt.subplot(3,3,8), plt.subplot(3,3,9)

PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gx'], title='Analytic gx',axs=ax1)
PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gy'], title='Analytic gy',axs=ax2)
PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gz'], title='Analytic gz',axs=ax3)

PF.Gravity.plot_obs_2D(locXyz,gx, title='PDE gx',axs=ax4)
PF.Gravity.plot_obs_2D(locXyz,gy, title='PDE gy',axs=ax5)
PF.Gravity.plot_obs_2D(locXyz,gz, title='PDE gz',axs=ax6)

PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gx']-gx, title='Residual gx',axs=ax7)
PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gy']-gy, title='Residual gy',axs=ax8)
PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gz']-gz, title='Residual gz',axs=ax9)

fig, ax1, ax2, ax3 = plt.figure(figsize=(15,12)), plt.subplot(3,3,1), plt.subplot(3,3,2), plt.subplot(3,3,3)
ax5, ax6, ax9 =  plt.subplot(3,3,5), plt.subplot(3,3,6), plt.subplot(3,3,9)

PF.Gravity.plot_obs_2D(locXyz,gg['gxx'], title='PDE gxx',axs=ax1)
PF.Gravity.plot_obs_2D(locXyz,gg['gxy'], title='PDE gxy',axs=ax2)
PF.Gravity.plot_obs_2D(locXyz,gg['gxz'], title='PDE gxz',axs=ax3)
PF.Gravity.plot_obs_2D(locXyz,gg['gyy'], title='PDE gyy',axs=ax5)
PF.Gravity.plot_obs_2D(locXyz,gg['gyz'], title='PDE gyz',axs=ax6)
PF.Gravity.plot_obs_2D(locXyz,gg['gzz'], title='PDE gzz',axs=ax9)

fig, ax1, ax2, ax3 = plt.figure(figsize=(15,12)), plt.subplot(3,3,1), plt.subplot(3,3,2), plt.subplot(3,3,3)
ax5, ax6, ax9 =  plt.subplot(3,3,5), plt.subplot(3,3,6), plt.subplot(3,3,9)

PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gxx'], title='Analytic gxx',axs=ax1)
PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gxy'], title='Analytic gxy',axs=ax2)
PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gxz'], title='Analytic gxz',axs=ax3)
PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gyy'], title='Analytic gyy',axs=ax5)
PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gyz'], title='Analytic gyz',axs=ax6)
PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gzz'], title='Analytic gzz',axs=ax9)

fig, ax1, ax2, ax3 = plt.figure(figsize=(15,12)), plt.subplot(3,3,1), plt.subplot(3,3,2), plt.subplot(3,3,3)
ax5, ax6, ax9 =  plt.subplot(3,3,5), plt.subplot(3,3,6), plt.subplot(3,3,9)

PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gxx']-gg['gxx'], title='Residual PDE gxx',axs=ax1)
PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gxy']-gg['gxy'], title='Residual PDE gxy',axs=ax2)
PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gxz']-gg['gxz'], title='Residual PDE gxz',axs=ax3)
PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gyy']-gg['gyy'], title='Residual PDE gyy',axs=ax5)
PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gyz']-gg['gyz'], title='Residual PDE gyz',axs=ax6)
PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gzz']-gg['gzz'], title='Residual PDE gzz',axs=ax9)


Mesh.TensorMesh.writeUBC(mesh,'Mesh.msh')
Mesh.TensorMesh.writeModelUBC(mesh,'Sphere.den',model)
survey.std = np.ones(survey.nD)
PF.Gravity.writeUBCobs('Obs.loc',survey,gg['gxy'])

#%% Load UBC file
surveyUBC = PF.Gravity.readUBCgravObs('ggfor3d.gg',gravGrad=True)

fig, ax1, ax2, ax3 = plt.figure(figsize=(15,12)), plt.subplot(3,3,1), plt.subplot(3,3,2), plt.subplot(3,3,3)
ax5, ax6, ax9 =  plt.subplot(3,3,5), plt.subplot(3,3,6), plt.subplot(3,3,9)

PF.Gravity.plot_obs_2D(locXyz,surveyUBC.dobs[:,0], title='UBC gxx',axs=ax1)
PF.Gravity.plot_obs_2D(locXyz,surveyUBC.dobs[:,1], title='UBC gxy',axs=ax2)
PF.Gravity.plot_obs_2D(locXyz,surveyUBC.dobs[:,2], title='UBC gxz',axs=ax3)
PF.Gravity.plot_obs_2D(locXyz,surveyUBC.dobs[:,3], title='UBC gyy',axs=ax5)
PF.Gravity.plot_obs_2D(locXyz,surveyUBC.dobs[:,4], title='UBC gyz',axs=ax6)
PF.Gravity.plot_obs_2D(locXyz,surveyUBC.dobs[:,5], title='UBC gzz',axs=ax9)

fig, ax1, ax2, ax3 = plt.figure(figsize=(15,12)), plt.subplot(3,3,1), plt.subplot(3,3,2), plt.subplot(3,3,3)
ax5, ax6, ax9 =  plt.subplot(3,3,5), plt.subplot(3,3,6), plt.subplot(3,3,9)

PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gxx']-surveyUBC.dobs[:,0], title='Residual UBC gxx',axs=ax1)
PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gxy']-surveyUBC.dobs[:,1], title='Residual UBC gxy',axs=ax2)
PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gxz']-surveyUBC.dobs[:,2], title='Residual UBC gxz',axs=ax3)
PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gyy']-surveyUBC.dobs[:,3], title='Residual UBC gyy',axs=ax5)
PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gyz']-surveyUBC.dobs[:,4], title='Residual UBC gyz',axs=ax6)
PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gzz']-surveyUBC.dobs[:,5], title='Residual UBC gzz',axs=ax9)

#%% Load VPmg file
surveyVPmg = np.loadtxt('vp.txt',skiprows=1)

fig, ax1, ax2, ax3 = plt.figure(figsize=(15,12)), plt.subplot(3,3,1), plt.subplot(3,3,2), plt.subplot(3,3,3)
ax5, ax6, ax9 =  plt.subplot(3,3,5), plt.subplot(3,3,6), plt.subplot(3,3,9)

PF.Gravity.plot_obs_2D(locXyz,surveyVPmg[:,3], title='VPmg gxx',axs=ax1)
PF.Gravity.plot_obs_2D(locXyz,surveyVPmg[:,4], title='VPmg gxy',axs=ax2)
PF.Gravity.plot_obs_2D(locXyz,surveyVPmg[:,5], title='VPmg gxz',axs=ax3)
PF.Gravity.plot_obs_2D(locXyz,surveyVPmg[:,6], title='VPmg gyy',axs=ax5)
PF.Gravity.plot_obs_2D(locXyz,surveyVPmg[:,7], title='VPmg gyz',axs=ax6)
PF.Gravity.plot_obs_2D(locXyz,surveyVPmg[:,8], title='VPmg gzz',axs=ax9)

fig, ax1, ax2, ax3 = plt.figure(figsize=(15,12)), plt.subplot(3,3,1), plt.subplot(3,3,2), plt.subplot(3,3,3)
ax5, ax6, ax9 =  plt.subplot(3,3,5), plt.subplot(3,3,6), plt.subplot(3,3,9)

PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gxx']-surveyVPmg[:,3], title='Residual VPmg gxx',axs=ax1)
PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gxy']-surveyVPmg[:,4], title='Residual VPmg gxy',axs=ax2)
PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gxz']-surveyVPmg[:,5], title='Residual VPmg gxz',axs=ax3)
PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gyy']-surveyVPmg[:,6], title='Residual VPmg gyy',axs=ax5)
PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gyz']-surveyVPmg[:,7], title='Residual VPmg gyz',axs=ax6)
PF.Gravity.plot_obs_2D(locXyz,AnaSphere['gzz']-surveyVPmg[:,8], title='Residual VPmg gzz',axs=ax9)