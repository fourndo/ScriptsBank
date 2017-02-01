import Mag as MAG
import simpegPF as PF
import simpegCoordUtils as Utils
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import animation
from JSAnimation import HTMLWriter
from scipy import interpolate

#%% Load the data
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Processing'
dbase = np.loadtxt(work_dir + '\Gridded_DIGHEM\Grid_25m_Filters.dat',delimiter=',')

lake = ["lakes\\lakes1.txt","lakes\\lakes1b.txt","lakes\\lakes2.txt","lakes\\lakes3.txt","lakes\\lakes4.txt"]

header = ["X","Y","Z","Analytical Signal","1th Vertical Derivative","Total Horizontal Derivative","Upward Continuation 100m","Residual Total Field (Detrended)","$1^{th}$ Order Polynomial","magnetic_1st_order_trend","Residual Total Field (Raw)"]
# X,Y,Z,Analytic_Signal,HD_filtered,Upward_100m,VD_filtered,magnetic,magnetic_1st_order_trend,magnetic_trended

clims = np.asarray([[0,0],[-2,6],[-75,300],[-5,10],[-5,75],[-10,75],[-75,300]])
         
inc, dec = 83.8, 25.4
obsloc = np.loadtxt(work_dir +'\DIGHEM_Loc.dat');

dmat = np.asarray(dbase)

xloc = np.unique(np.round(dmat[:,0]))
yloc = np.unique(np.round(dmat[:,1]))
zloc = np.unique(np.round(dmat[:,2]))


Xloc = dmat[:,0].reshape((yloc.shape[0],xloc.shape[0]))
Yloc = dmat[:,1].reshape((yloc.shape[0],xloc.shape[0]))

# Interpolation line
xin = np.asarray([558400,559275])-1500
yin = np.asarray([7134250,7134640])-750

xlin = np.linspace(xin[0],xin[1],30)
ylin = np.linspace(yin[0],yin[1],30)

#%% Simiulation
fig = plt.figure(figsize=(10,5))
ax1 = plt.subplot(121, projection='3d')
ax2 = plt.subplot(224)
ax3 = plt.subplot(224)
# Define the dimensions of the prism (m)
dx, dy, dz = 300.,200., 50.
# Set the depth of burial (m)
depth = 110
pinc, pdec = 20., 225.
npts2D, xylim = 40., 500.
rx_h, View_elev, View_azim = 30., 20, -115
Einc, Edec, Bigrf = 45., 45., 50000.
x1, x2, y1, y2 = x1, x2, y1, y2 = -xylim, xylim, 0., 0.
comp = 'tf'
irt = 'total'
Q, rinc,rdec = 0., 0., 90.
susc = 0.05


ax1.axis('equal')
ax1.set_title(irt+' '+comp)
# Define the problem interactively
p = MAG.definePrism()
p.dx, p.dy, p.dz = dx, dy, dz
p.z0 = -depth
p.pinc, p.pdec = pinc, pdec

srvy = PF.survey()
srvy.rx_h, srvy.npts2D, srvy.xylim = rx_h, npts2D, xylim
#srvy._xr, srvy._yr = xloc, yloc

# Create problem
prob = PF.problem()
prob.prism = p
prob.survey = srvy

X, Y = np.meshgrid(prob.survey.xr, prob.survey.yr)

x, y = MAG.linefun(x1, x2, y1, y2, prob.survey.npts2D)
xyz_line = np.c_[x, y, np.ones_like(x)*prob.survey.rx_h]



ax1.plot(x,y,xyz_line[:,2], 'w.', ms=3,lw=2)
ax1.text(-1,0., -3.25, 'B-field', fontsize=14, color='k', horizontalalignment='left')
#im1 = ax1.contourf(X,Y,X)
#im2 = ax2.plot(x,y)
#im3 = ax2.plot(x,y)
#im4 = ax2.plot(x,y)
#im5 = ax1.text(0,0,0,'')
clim = np.asarray([-75,300])
#def animate(ii):
#
#removePlt()

#if ii<19:
#    dec = 0
#    inc = -90. + ii*10.
#    
#else:
#
#    dec = 90.
#    inc = 90. - (ii-19)*10.
    
MAG.plotObj3D(p, rx_h, View_elev, View_azim, npts2D, xylim, profile="X", fig= fig, axs = ax1, plotSurvey=False)
plt.show()
block_xyz = np.asarray([[-.2, -.2, .2, .2, 0],
                       [-.25, -.25, -.25, -.25, 0.5],
                       [-.2, .2, .2, -.2, 0]])

# rot = Utils.mkvc(Utils.dipazm_2_xyz(pinc, pdec))

# xyz = Utils.rotatePointsFromNormals(block_xyz.T, np.r_[0., 1., 0.], rot,
#                                     np.r_[p.xc, p.yc, p.zc])

R = Utils.rotationMatrix(inc, dec)

xyz = R.dot(block_xyz).T

#print xyz
# Face 1
ax1.add_collection3d(Poly3DCollection([zip(xyz[:4, 0]-1,
                                           xyz[:4, 1],
                                           xyz[:4, 2]-4)], facecolors='w'))

ax1.add_collection3d(Poly3DCollection([zip(xyz[[1, 2, 4], 0]-1,
                                           xyz[[1, 2, 4], 1],
                                           xyz[[1, 2, 4], 2]-4)], facecolors='b'))

ax1.add_collection3d(Poly3DCollection([zip(xyz[[0, 1, 4], 0]-1,
                                           xyz[[0, 1, 4], 1],
                                           xyz[[0, 1, 4], 2]-4)], facecolors='r'))

ax1.add_collection3d(Poly3DCollection([zip(xyz[[2, 3, 4], 0]-1,
                                           xyz[[2, 3, 4], 1],
                                           xyz[[2, 3, 4], 2]-4)], facecolors='k'))

ax1.add_collection3d(Poly3DCollection([zip(xyz[[0, 3, 4], 0]-1,
                                       xyz[[0, 3, 4], 1],
                                       xyz[[0, 3, 4], 2]-4)], facecolors='y'))




# Create problem
prob = PF.problem()
prob.prism = p
prob.survey = srvy

prob.Bdec, prob.Binc, prob.Bigrf = dec, inc, Bigrf
prob.Q, prob.rinc, prob.rdec = Q, rinc, rdec
prob.uType, prob.mType = comp, 'total'
prob.susc = susc

# Compute fields from prism
b_ind, b_rem = prob.fields()

if irt == 'total':
    out = b_ind + b_rem

elif irt == 'induced':
    out = b_ind

else:
    out = b_rem



#out = plogMagSurvey2D(prob, susc, Einc, Edec, Bigrf, x1, y1, x2, y2, comp, irt,  Q, rinc, rdec, fig=fig, axs1=ax2, axs2=ax3)


data = dmat[:,7]
#dat = axs1.contourf(X,Y, np.reshape(out, (X.shape)).T
#global im1
im1 = ax1.contourf(X,Y,np.reshape(out, (X.shape)).T,zdir='z',offset=rx_h-0.1, clim=clim, vmin=clim[0],vmax=clim[1])
ax1.set_zlim(-xylim,0)

pos = ax1.get_position()
ax3 = fig.add_axes([pos.x0 + 0.5, pos.y0+0.45,  pos.width*0.6, pos.height*0.5]) 
    
im7 = ax3.contourf(xloc,yloc,np.reshape(data, (xloc.shape[0],yloc.shape[0])).T,zdir='z',offset=rx_h-0.1, alpha=0.75, clim=clim, vmin=clim[0],vmax=clim[1])
ax3.scatter(xin,yin,color='k')
ax3.plot(xin,yin,'k--')
ax3.axis('equal')
ax3.set_xlim(xloc.min(),xloc.max())
ax3.set_ylim(yloc.min(),yloc.max())
xpos = np.asarray([557000, 559000])
ax3.set_xticks(map(int, xpos))
ax3.set_xticklabels(map(str, map(int, xpos)),size=12)

#global im5
im5 = ax1.text(-.5,0,-4,'I: ' + str(inc) + ' D: ' + str(dec), horizontalalignment='left')

# Create problem
prob1D = PF.problem()
srvy1D = PF.survey()
srvy1D._rxLoc = xyz_line

prob1D.prism = prob.prism
prob1D.survey = srvy1D

prob1D.Bdec, prob1D.Binc, prob1D.Bigrf = dec, inc, Bigrf
prob1D.Q, prob1D.rinc, prob1D.rdec = Q, rinc, rdec
prob1D.uType, prob1D.mType = comp, 'total'
prob1D.susc = susc

# Compute fields from prism
out_linei, out_liner = prob1D.fields()

#out_linei, out_liner = getField(p, xyz_line, comp, 'total')
#out_linei = getField(p, xyz_line, comp,'induced')
#out_liner = getField(p, xyz_line, comp,'remanent')

out_linet = out_linei+out_liner

distance = np.sqrt((x-x1)**2.+(y-y1)**2.)


# Interpolate true data
F = interpolate.interp2d(xloc,yloc,np.reshape(data, (xloc.shape[0],yloc.shape[0])).T)
interp = np.zeros(xlin.shape)
for ii in range(xlin.shape[0]):
    interp[ii] = F(xlin[ii],ylin[ii])

#global im2
im2 =ax2.plot(distance, out_linei, 'b.-')

#global im3
dist = np.sqrt((xin[0]-xlin)**2.+(yin[0]-ylin)**2.)
im3 =ax2.plot(dist, interp, 'r.-')

#global im4
#im4 =ax2.plot(distance, out_linei-interp, 'k.-')

ax2.set_xlim(distance.min(), distance.max())
ax2.set_xlabel("Distance (m)")
ax2.set_ylabel("Magnetic field (nT)")
ax2.set_ylim(clim)

ax2.legend(("parametric", "Observed data"),fontsize=10)
ax2.grid(True)
plt.show()




#def removePlt():
#    #global im1
#    #im1.remove()
#
#    global im1
#    for coll in im1.collections:
#        coll.remove()
#
#    for cc in range(6):
#        for coll in ax1.collections:
#            ax1.collections.remove(coll)
#
#
#    global im2
#    im2.pop(0).remove()
#
#
#    global im3
#    im3.pop(0).remove()
#
#    global im4
#    im4.pop(0).remove()
#
#    global im5
#    im5.remove()

#anim = animation.FuncAnimation(fig, animate,
#                               frames=37 , interval=200,repeat=False)
#
#anim.save('animation.html', writer=HTMLWriter(embed_frames=True,fps=10,default_mode = 'reflect'))
