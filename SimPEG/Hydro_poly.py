from SimPEG import *
import simpegEM as EM

#from pymatsolver import MumpsSolver
from numpy.polynomial import polynomial

cs, ncx, ncy, ncz, npad = 40., 15, 10, 15, 2

hx = [(cs,npad,-1.4), (cs,ncx), (cs,npad,1.4)]
hy = [(cs,npad,-1.4), (cs,ncy), (cs,npad,1.4)]
hz = [(cs,npad,-1.4), (cs,ncz), (cs,npad,1.4)]
mesh = Mesh.TensorMesh([hx,hy,hz], 'CCC')
print ("Padding distance x: %10.5f m") % (np.sum(mesh.hx[:npad]))
print ("Padding distance z: %10.5f m") % (np.sum(mesh.hz[:npad]))
print ("Min dx: %10.5f m") % (mesh.hx.min())
print ("Min dz: %10.5f m") % (mesh.hz.min())

# sigma3D = Utils.meshutils.readUBCTensorModel('../sigma_realistic.con', mesh)
# sigma3D = Utils.meshutils.readUBCTensorModel('../examples/sigma_realistic.con', mesh)

sigma3D = np.ones(mesh.nC)*1e-8
sigma3D[(mesh.gridCC[:,2] >-60) & (mesh.gridCC[:,2] <0.)] = 3e-3
sigma3D[(mesh.gridCC[:,2] >-160) & (mesh.gridCC[:,2] <-60)] = 2e-3
sigma3D[mesh.gridCC[:,2] <-160 ] = 1e-3

x1 = np.arange(14)*20 - 310.
y1 = np.arange(14)*20 - 150.+20
xyz1 = Utils.ndgrid(x1, y1, np.r_[0.])

x2 = np.arange(14)*20 - 30.
y2 = np.arange(14)*20 - 150.+20
xyz2 = Utils.ndgrid(x2, y2, np.r_[0.])

ntx = 2
time = np.logspace(-4, -2, 31)

rx1 = EM.TDEM.RxTDEM(xyz1, time, 'bz')
tx1 = EM.TDEM.SrcTDEM_CircularLoop_MVP([rx1], np.array([-180., 0.,  0.]), radius=250.)
rx2 = EM.TDEM.RxTDEM(xyz2, time, 'bz')
tx2 = EM.TDEM.SrcTDEM_CircularLoop_MVP([rx2], np.array([100., 0.,  0.]), radius=250.)


actind = np.logical_and((mesh.gridCC[:,2] >-160) & (mesh.gridCC[:,2] <-60), 
                       (mesh.gridCC[:,1] >-400) & (mesh.gridCC[:,1] <400))
actMap = Maps.ActiveCells(mesh, actind, sigma3D)

ind1z = (mesh.vectorCCz >-60) & (mesh.vectorCCz <0.)
ind2z = (mesh.vectorCCz >-160) & (mesh.vectorCCz <-60)
indy = (mesh.vectorCCy >-400) & (mesh.vectorCCy <400)
meshact = Mesh.TensorMesh([mesh.hx, mesh.hy[indy], mesh.hz[ind2z]], x0 = np.r_[mesh.x0[0],mesh.vectorNy[indy][0], mesh.vectorNz[ind2z][0]])
order = [1, 1]
YZ = Utils.ndgrid(meshact.vectorCCy, meshact.vectorCCz)
V = polynomial.polyvander2d(YZ[:,0], YZ[:,1], order)

polymap = Maps.PolyMap(meshact, order, normal='X', logSigma=True)
m1D = Mesh.TensorMesh([(order[0]+1)*(order[1]+1)+2])
weight = (1./(V**2).sum(axis=0))**0.5
weight = weight / weight.max()
weightmap = Maps.Weighting(m1D, weights=np.r_[1., 1., weight])

actindpoly = np.ones(6,dtype=bool)
actindpoly[:2] = False
m0 = np.r_[np.log(2e-3), np.log(2e-3), np.zeros(4) / weight ]

mapping = actMap * polymap * weightmap
survey = EM.TDEM.SurveyTDEM([tx1, tx2])
prb = EM.TDEM.ProblemTDEM_b(mesh, mapping=mapping, verbose=False)
prb.Solver = MumpsSolver
prb.solverOpts = {"symmetric":True}
prb.timeSteps = [(1e-4/10, 15), (1e-3/10, 15), (1e-2/10, 15)]
if prb.ispaired:
    prb.unpair()
if survey.ispaired:
    survey.unpair()
prb.pair(survey)


dobs = np.load('bzobs_realistic.npy')

std = 0.05
survey.dobs = Utils.mkvc(dobs)
survey.std = survey.dobs*0 + std
dmis = DataMisfit.l2_DataMisfit(survey)
dmis.Wd = 1/(abs(survey.dobs)*std)
opt = Optimization.InexactGaussNewton(maxIter = 40, maxIterLS=20)
reg = Regularization.BaseRegularization(m1D)
beta = Directives.BetaSchedule(coolingFactor = 1, coolingRate = 1)
invprb = InvProblem.BaseInvProblem(dmis, reg, opt)
invprb.beta = 0.
savemodel = Directives.SaveModelEveryIteration()
saveinvres = Directives.SaveOutputEveryIteration()
target = Directives.TargetMisfit()
inv = Inversion.BaseInversion(invprb, directiveList = [beta, target, savemodel, saveinvres])
reg.mref = m0
C =  Utils.Counter()
prb.counter = C
opt.counter = C
mopt = inv.run(m0)
opt.counter.summary()
