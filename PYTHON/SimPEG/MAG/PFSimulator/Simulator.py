import Mag as MAG
import Utils
import SimPEG.PF as PF
from SimPEG import mkvc
from scipy.constants import mu_0
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import notebook
import ipywidgets as widgets
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.interpolate import griddata

# monFile = "data2015/StudentData2015_Monday.csv"
# monData = pd.DataFrame(pd.read_csv(filename, header = 0))

# filename = "data2014/HZrebarProfile.csv"
# data = pd.DataFrame(pd.read_csv(filename, header = 0))
# loc = data["Distance"].values

diameter = 1.4e-2
length = 3.
xlim = np.r_[5., 25.]
rx_h = 1.9

Bincd = 70.205
Bdecd = 16.63
Bigrfd = 54399

# Depth of burial: Monday was 35cm. I believe that Wednesday was ~45cm


class definePrism(object):
    """
        Define a prism and its attributes

        Prism geometry:
            - dx, dy, dz: width, length and height of prism
            - depth : depth to top of prism
            - susc : susceptibility of prism
            - x0, y0 : center of prism in horizontal plane
            - pinc, pdec : inclination and declination of prism
    """

    x0, y0, z0, dx, dy, dz = 0., 0., 0., 1., 1., 1.
    pinc, pdec = 0., 0.


    # Define the nodes of the prism
    @property
    def xn(self):
        xn = np.asarray([-self.dx/2. + self.x0, self.dx/2. + self.x0])

        return xn

    @property
    def yn(self):
        yn = np.asarray([-self.dy/2. + self.y0, self.dy/2. + self.y0])

        return yn

    @property
    def zn(self):
        zn = np.asarray([-self.dz + self.z0, self.z0])

        return zn

    @property
    def xc(self):
        xc = (self.xn[0] + self.xn[1]) / 2.

        return xc

    @property
    def yc(self):
        yc = (self.yn[0] + self.yn[1]) / 2.

        return yc

    @property
    def zc(self):
        zc = (self.zn[0] + self.zn[1]) / 2.

        return zc


def ForwardMagnetic(prism, param):


    nx = npts2D
    xr = np.linspace(-xylim, xylim, nx)

    X, Y = np.meshgrid(xr, xr)
    Z = np.ones(np.shape(X))*rx_h

    rxLoc = np.c_[Utils.mkvc(X), Utils.mkvc(Y), Utils.mkvc(Z)]

    rxLoc = PF.BaseMag.RxObs(rxLoc)
    srcField = PF.BaseMag.SrcField([rxLoc], param=survey.srcField.param)
    survey = PF.BaseMag.LinearSurvey(srcField)

    # Create problem
    prob = PF.Magnetics.problem()

    survey.pair(prob)

    nx, ny = 100, 1
    shape = (nx, ny)
    xLoc = np.linspace(xlim[0], xlim[1], nx)

    zLoc = np.ones(np.shape(xLoc))*rx_h
    yLoc = np.zeros(np.shape(xLoc))

    #xpl, ypl, zpl = fatiandoGridMesh.regular(surveyArea,shape, z=z)
    rxLoc = np.c_[Utils.mkvc(xLoc), Utils.mkvc(yLoc), Utils.mkvc(zLoc)]

    prob1D = MAG.problem()
    srvy1D = MAG.survey()
    srvy1D._rxLoc = rxLoc

    prob1D.prism = p
    prob1D.survey = srvy1D

    prob1D.Bdec, prob1D.Binc, prob1D.Bigrf = Bdec, Binc, Bigrf
    prob1D.Q, prob1D.rinc, prob1D.rdec = Q, rinc, rdec
    prob1D.uType, prob1D.mType = 'tf', 'total'
    prob1D.susc = susc

    # Compute fields from prism
    magi, magr = prob1D.fields()

    #out_linei, out_liner = getField(p, xyz_line, comp, 'total')
    #out_linei = getField(p, xyz_line, comp,'induced')
    #out_liner = getField(p, xyz_line, comp,'remanent')

    # distance = np.sqrt((x-x1)**2.+(y-y1)**2.)

    f = plt.figure(figsize = (10, 5))
    gs = gridspec.GridSpec(2, 1,height_ratios=[2,1])

    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])

    ax1.plot(p.x0, p.z0, 'ko')
    ax1.text(p.x0+0.5, p.z0, 'Rebar', color='k')
    ax1.text(xlim[0]+1.,-1.2, 'Magnetometer height (1.9 m)', color='b')
    ax1.plot(xlim, np.r_[-rx_h, -rx_h], 'b--')

    # magi,magr = getField(p, rxLoc, 'bz', 'total')

    ax1.plot(xlim, np.r_[0., 0.], 'k--')
    ax1.set_xlim(xlim)
    ax1.set_ylim(-2.5, 2.5)

    ax0.scatter(loc,tfa,c=teams)
    ax0.errorbar(loc,tfa,yerr=std,linestyle = "None",color="k")
    ax0.set_xlim(xlim)
    ax0.grid(which="both")

    ax0.plot(xLoc, magi, 'b', label='induced')
    ax0.plot(xLoc, magr, 'r', label='remnant')
    ax0.plot(xLoc, magi+magr, 'k', label='total')
    ax0.legend(loc=2)
    # ax[1].plot(loc-8, magnT[::-1], )

    ax1.set_xlabel("Northing (m)")
    ax1.set_ylabel("Depth (m)")

    ax0.set_ylabel("Total field anomaly (nT)")

    ax0.grid(True)
    ax0.set_xlabel("Northing (m)")

    ax1.grid(True)
    ax1.set_xlabel("Northing (m)")

    ax1.invert_yaxis()

    plt.tight_layout()
    plt.show()

    return True





def linefun(x1, x2, y1, y2, nx, tol=1e-3):
    dx = x2-x1
    dy = y2-y1

    if np.abs(dx) < tol:
        y = np.linspace(y1, y2, nx)
        x = np.ones_like(y)*x1
    elif np.abs(dy) < tol:
        x = np.linspace(x1, x2, nx)
        y = np.ones_like(x)*y1
    else:
        x = np.linspace(x1, x2, nx)
        slope = (y2-y1)/(x2-x1)
        y = slope*(x-x1)+y1
    return x, y


def plotMagSurvey2D(survey, xlim, ylim, a, b, npts, fig=None, axs1=None, axs2=None):
    """
    Plot the data and line profile inside the spcified limits
    """

    if fig is None:
        fig = plt.figure(figsize=(9, 12))

        plt.rcParams.update({'font.size': 14})


    if axs1 is None:
        axs1 = plt.subplot(2, 1, 1)


    if axs2 is None:
        axs2 = plt.subplot(2, 1, 2)



    rxLoc = survey.srcField.rxList[0].locs

    # Impose bound limits
    ind = np.all([rxLoc[:, 0] > xlim[0], rxLoc[:, 0] < xlim[1],
                  rxLoc[:, 1] > ylim[0], rxLoc[:, 1] < ylim[1]], axis=0)

    # Subsample the full survey
    rxLoc = rxLoc[ind, :]

    data = survey.dobs[ind]

    # Use SimPEG.PF ploting function
    PF.Magnetics.plot_obs_2D(rxLoc, data, fig=fig, ax=axs1)
    # pos = axs1.get_position()
    # axs1.set_position([pos.x0, pos.y0,  pos.height*1.1, pos.width*1.1])
    # Interpolate the data
    # F = LinearNDInterpolator(zip(rxLoc[:, :2]), data)

    x, y = linefun(a[0], b[0], a[1], b[1], npts)
    # print(rxLoc[:,:2],x,y)
    # F(np._r[x, y])

    dline = griddata(rxLoc[:, :2], data, (x, y), method='linear')

    distance = np.sqrt((x-a[0])**2.+(y-a[1])**2.)

    axs1.plot(x, y, 'w.', ms=10)

    axs1.text(x[0], y[0], 'A', fontsize=16, color='w', horizontalalignment='left')
    axs1.text(x[-1], y[-1], 'B', fontsize=16,
             color='w', horizontalalignment='right')


    axs2.plot(distance, dline, 'b.-')
    # axs2.plot(distance, out_liner, 'r.-')
    # axs2.plot(distance, out_linet, 'k.-')
    axs2.set_xlim(distance.min(), distance.max())

    axs2.set_xlabel("Distance (m)")
    axs2.set_ylabel("Magnetic field (nT)")

    #axs2.text(distance.min(), dline.max()*0.8, 'A', fontsize = 16)
    # axs2.text(distance.max()*0.97, out_linei.max()*0.8, 'B', fontsize = 16)
    axs2.legend(("induced", "remanent", "total"), bbox_to_anchor=(0.5, -0.3))
    axs2.grid(True)
    plt.show()

    out = {'survey': survey, 'xlim': xlim, 'ylim': ylim,
           'a': a, 'b': b, 'npts': npts, }
    return out


# def fitline(Box):

#     def profiledata(data, Binc, Bdec, Bigrf, x0, depth, susc, Q, rinc, rdec):

#         prob = Box.result[1]
#         prism = prob.prism
#         prism.x0, prism.z0 = x0-prism.dx/2, -depth

#         return plotProfile(prism, data, Binc, Bdec, Bigrf, susc, Q, rinc, rdec)

#     Q = widgets.interactive(profiledata, data=widgets.ToggleButtons(options=['MonSt','WedTA','WedSt']),\
#              Binc=widgets.FloatSlider(min=-90.,max=90,step=5,value=90,continuous_update=False),\
#              Bdec=widgets.FloatSlider(min=-90.,max=90,step=5,value=0,continuous_update=False),\
#              Bigrf=widgets.FloatSlider(min=54000.,max=55000,step=10,value=54500,continuous_update=False),\
#              x0=widgets.FloatSlider(min=5., max=25., step=0.1, value=15.), \
#              depth=widgets.FloatSlider(min=0.,max=2.,step=0.05,value=0.5), \
#              susc=widgets.FloatSlider(min=0., max=800.,step=5., value=1.),\
#              Q=widgets.FloatSlider(min=0., max=10.,step=0.1, value=0.),\
#              rinc=widgets.FloatSlider(min=-180., max=180.,step=1., value=0.),\
#              rdec=widgets.FloatSlider(min=-180., max=180.,step=1., value=0.),
#              )
#     return Q


def ViewMagSurvey2D(survey):

    def MagSurvey2D(East, North, dx, dy, azm, length, npts):

        # Get the line extent from the 2D survey for now
        azm /= 180./np.pi
        length /= 2.*0.98
        a = [East - np.cos(-azm)*length, North - np.sin(-azm)*length]
        b = [East + np.cos(-azm)*length, North + np.sin(-azm)*length]

        xlim = East + np.asarray([-dx/2., dx/2.])
        ylim = North + np.asarray([-dy/2., dy/2.])

        return plotMagSurvey2D(survey, xlim, ylim, a, b, npts)

    locs = survey.srcField.rxList[0].locs
    xlim = np.asarray([locs[:, 0].min(), locs[:, 0].max()])
    ylim = np.asarray([locs[:, 1].min(), locs[:, 1].max()])

    Lx = xlim[1] - xlim[0]
    Ly = ylim[1] - ylim[0]
    diag = (Lx**2. + Ly**2.)**0.5

    East = np.mean(xlim)
    North = np.mean(ylim)
    cntr = [East, North]


    out = widgets.interactive(MagSurvey2D,
                    East=widgets.FloatSlider(min=xlim[0], max=xlim[1], step=10, value=cntr[0],continuous_update=False),
                    North=widgets.FloatSlider(min=ylim[0], max=ylim[1], step=10, value=cntr[1],continuous_update=False),
                    dx=widgets.FloatSlider(min=10, max=Lx, step=10, value=Lx, continuous_update=False),
                    dy=widgets.FloatSlider(min=10, max=Ly, step=10, value=Ly, continuous_update=False),
                    azm=widgets.FloatSlider(min=-90, max=90, step=5, value=0, continuous_update=False),
                    length=widgets.FloatSlider(min=10, max=diag, step=10, value= Ly, continuous_update=False),
                    npts=widgets.BoundedFloatText(min=10, max=100, step=1, value=20, continuous_update=False))
    return out


def Prism(dx, dy, dz, depth, pinc, pdec, xylim, View_elev, View_azim):
    #p = definePrism(dx, dy, dz, depth,pinc=pinc, pdec=pdec, susc = 1., Einc=90., Edec=0., Bigrf=1e-6)
    prism = definePrism()
    prism.dx, prism.dy, prism.dz, prism.z0 = dx, dy, dz, -depth
    prism.pinc, prism.pdec = pinc, pdec

    return plotObj3D(prism, xylim, View_elev, View_azim)

def plotObj3D(prism, xylim, elev,  azim,
              profile=None, x0=15., y0=0., fig=None, axs=None, plotSurvey=True, title=None):

    """
    Plot the prism in 3D
    """
    surveyArea = (-xylim, xylim)

    depth = prism.z0
    x1, x2 = prism.xn[0]-prism.xc, prism.xn[1]-prism.xc
    y1, y2 = prism.yn[0]-prism.yc, prism.yn[1]-prism.yc
    z1, z2 = prism.zn[0]-prism.zc, prism.zn[1]-prism.zc
    pinc, pdec = prism.pinc, prism.pdec

    if fig is None:
        fig = plt.figure(figsize=(7, 7))

    if axs is None:
        axs = fig.add_subplot(111, projection='3d')

    if title is not None:
        axs.set_title(title)

    plt.rcParams.update({'font.size': 13})

    axs.set_xlim3d(surveyArea)
    axs.set_ylim3d(surveyArea)
#     axs.set_zlim3d(depth+np.array(surveyArea[:2]))
    axs.set_zlim3d(-surveyArea[-1]*2, 3)

    # Create a rectangular prism, rotate and plot
    block_xyz = np.asarray([[x1, x1, x2, x2, x1, x1, x2, x2],
                           [y1, y2, y2, y1, y1, y2, y2, y1],
                           [z1, z1, z1, z1, z2, z2, z2, z2]])

    # rot = Utils.mkvc(Utils.dipazm_2_xyz(pinc, pdec))

    # xyz = Utils.rotatePointsFromNormals(block_xyz.T, np.r_[0., 1., 0.], rot,
    #                                     np.r_[p.xc, p.yc, p.zc])

    R = Utils.rotationMatrix(pinc, pdec)

    xyz = R.dot(block_xyz).T

    #print xyz
    # Face 1
    axs.add_collection3d(Poly3DCollection([zip(xyz[:4, 0]+prism.xc,
                                               xyz[:4, 1]+prism.yc,
                                               xyz[:4, 2]+prism.zc)], facecolors='w'))

    # Face 2
    axs.add_collection3d(Poly3DCollection([zip(xyz[4:, 0]+prism.xc,
                                               xyz[4:, 1]+prism.yc,
                                               xyz[4:, 2]+prism.zc)], facecolors='w'))

    # Face 3
    axs.add_collection3d(Poly3DCollection([zip(xyz[[0, 1, 5, 4], 0]+prism.xc,
                                               xyz[[0, 1, 5, 4], 1]+prism.yc,
                                               xyz[[0, 1, 5, 4], 2]+prism.zc)], facecolors='w'))

    # Face 4
    axs.add_collection3d(Poly3DCollection([zip(xyz[[3, 2, 6, 7], 0]+prism.xc,
                                               xyz[[3, 2, 6, 7], 1]+prism.yc,
                                               xyz[[3, 2, 6, 7], 2]+prism.zc)], facecolors='w'))

   # Face 5
    axs.add_collection3d(Poly3DCollection([zip(xyz[[0, 4, 7, 3], 0]+prism.xc,
                                               xyz[[0, 4, 7, 3], 1]+prism.yc,
                                               xyz[[0, 4, 7, 3], 2]+prism.zc)], facecolors='w'))

   # Face 6
    axs.add_collection3d(Poly3DCollection([zip(xyz[[1, 5, 6, 2], 0]+prism.xc,
                                               xyz[[1, 5, 6, 2], 1]+prism.yc,
                                               xyz[[1, 5, 6, 2], 2]+prism.zc)], facecolors='w'))

    plt.show()
    axs.set_xlabel('Easting (X; m)')
    axs.set_ylabel('Northing (Y; m)')
    axs.set_zlabel('Depth (Z; m)')
    # axs.invert_zaxis()
    # axs.invert_yaxis()

    # if plotSurvey:
    #     axs.plot(rxLoc[:, 0], rxLoc[:, 1], rxLoc[:, 2], '.g', alpha=0.5)

    # if profile == "X":
    #     axs.plot(np.r_[surveyArea[:2]], np.r_[0., 0.], np.r_[rx_h, rx_h], 'r-')
    # elif profile == "Y":
    #     axs.plot(np.r_[0., 0.], np.r_[surveyArea[2:]], np.r_[rx_h, rx_h], 'r-')
    # elif profile == "XY":
    #     axs.plot(np.r_[0., 0.], np.r_[surveyArea[:2]], np.r_[rx_h, rx_h], 'r-')
    #     axs.plot(np.r_[surveyArea[2:]], np.r_[0., 0.], np.r_[rx_h, rx_h], 'r-')

    axs.view_init(elev, azim)
    plt.show()

    return True


def ViewPrism(param):
    elev, azim = 20, 250
    npts2D = param.result['npts']
    survey = param.result['survey']
    a, b = np.asarray(param.result['a']), np.asarray(param.result['b'])
    xylim = np.sum((a-b)**2.)**0.5

    Q = widgets.interactive(Prism \
                            , dx=widgets.FloatSlider(min=1., max=1000., step=1., value=100, continuous_update=False) \
                            , dy=widgets.FloatSlider(min=1., max=1000., step=1., value=100, continuous_update=False) \
                            , dz=widgets.FloatSlider(min=1., max=1000., step=1., value=100, continuous_update=False) \
                            , depth=widgets.FloatSlider(min=0., max=1000., step=1., value=-50, continuous_update=False)\
                            , pinc=(-90., 90., 5.)\
                            , pdec=(-90., 90., 5.)\
                            , xylim=widgets.FloatSlider(min=1, max=1000, step=1, value=xylim, continuous_update=False) \
                            , View_elev=widgets.FloatSlider(min=-90, max=90, step=5, value=elev, continuous_update=False) \
                            , View_azim=widgets.FloatSlider(min=0, max=360, step=5, value=azim, continuous_update=False)
                            )

    return Q
