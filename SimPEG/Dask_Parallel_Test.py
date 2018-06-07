import numpy as np
import SimPEG.PF as PF
import dask as dd
import dask.array as da
import scipy.constants as constants
from SimPEG import Problem
from SimPEG import Utils
from SimPEG import Props
from SimPEG.Utils import mkvc
import scipy.sparse as sp
import os

import time

# Create functions for Dask
def calcTrow(Xn, Yn, Zn, xyzLoc, rxType='z'):
    """
    Load in the active nodes of a tensor mesh and computes the gravity tensor
    for a given observation location xyzLoc[obsx, obsy, obsz]

    INPUT:
    Xn, Yn, Zn: Node location matrix for the lower and upper most corners of
                all cells in the mesh shape[nC,2]
    M
    OUTPUT:
    Tx = [Txx Txy Txz]
    Ty = [Tyx Tyy Tyz]
    Tz = [Tzx Tzy Tzz]

    where each elements have dimension 1-by-nC.
    Only the upper half 5 elements have to be computed since symetric.
    Currently done as for-loops but will eventually be changed to vector
    indexing, once the topography has been figured out.

    """

    NewtG = constants.G*1e+8  # Convertion from mGal (1e-5) and g/cc (1e-3)
    eps = 1e-8  # add a small value to the locations to avoid

    # Pre-allocate space for 1D array
    row = np.zeros((1, Xn.shape[0]))

    dz = xyzLoc[2] - Zn

    dy = Yn - xyzLoc[1]

    dx = Xn - xyzLoc[0]

    # Compute contribution from each corners
    for aa in range(2):
        for bb in range(2):
            for cc in range(2):

                r = (
                        mkvc(dx[:, aa]) ** 2 +
                        mkvc(dy[:, bb]) ** 2 +
                        mkvc(dz[:, cc]) ** 2
                    ) ** (0.50)

                if rxType == 'x':
                    row = row - NewtG * (-1) ** aa * (-1) ** bb * (-1) ** cc * (
                        dy[:, bb] * np.log(dz[:, cc] + r + eps) +
                        dz[:, cc] * np.log(dy[:, bb] + r + eps) -
                        dx[:, aa] * np.arctan(dy[:, bb] * dz[:, cc] /
                                              (dx[:, aa] * r + eps)))

                elif rxType == 'y':
                    row = row - NewtG * (-1) ** aa * (-1) ** bb * (-1) ** cc * (
                        dx[:, aa] * np.log(dz[:, cc] + r + eps) +
                        dz[:, cc] * np.log(dx[:, aa] + r + eps) -
                        dy[:, bb] * np.arctan(dx[:, aa] * dz[:, cc] /
                                              (dy[:, bb] * r + eps)))

                else:
                    row -= NewtG * (-1) ** aa * (-1) ** bb * (-1) ** cc * (
                        dx[:, aa] * np.log(dy[:, bb] + r + eps) +
                        dy[:, bb] * np.log(dx[:, aa] + r + eps) -
                        dz[:, cc] * np.arctan(dx[:, aa] * dy[:, bb] /
                                              (dz[:, cc] * r + eps)))

    return row

if __name__ == '__main__':

    work_dir = "C:\\Users\\DominiqueFournier\\ownCloud\\Research\\Synthetic\\Block_Gaussian_topo\\GRAV\\"

    # work_dir = "C:\\Users\\DominiqueFournier\\Documents\\GIT\\InnovationGeothermal\\FORGE\\SyntheticModel\\"

    inpfile = 'SimPEG_GRAV.inp'
    out_dir = "SimPEG_GRAV_Inv\\"
    dsep = '\\'
    dsep = os.path.sep

    # Read input file
    driver = PF.GravityDriver.GravityDriver_Inv(work_dir + dsep + inpfile)
    mesh = driver.mesh
    survey = driver.survey

    prob = PF.Gravity.GravityIntegral(mesh=mesh,
                                      n_cpu=None, parallelized=True)
    prob.solverOpts['accuracyTol'] = 1e-4

    survey.pair(prob)

    tc = time.time()

    F = prob.F

    print("Multirpocessing CPU Time: " + str(time.time()-tc) + " s")

    # Repeat with dask
    rxLoc = survey.srcField.rxList[0].locs

    # Setup the forward calc
    bsw = (mesh.gridCC -
           np.kron(mesh.vol.T**(1/3)/2,
                   np.ones(3)).reshape((mesh.nC, 3)))
    tne = (mesh.gridCC +
           np.kron(mesh.vol.T**(1/3)/2,
                   np.ones(3)).reshape((mesh.nC, 3)))

    xn1, xn2 = bsw[:, 0], tne[:, 0]
    yn1, yn2 = bsw[:, 1], tne[:, 1]
    zn1, zn2 = bsw[:, 2], tne[:, 2]

    Yn = np.c_[Utils.mkvc(yn1), Utils.mkvc(yn2)]
    Xn = np.c_[Utils.mkvc(xn1), Utils.mkvc(xn2)]
    Zn = np.c_[Utils.mkvc(zn1), Utils.mkvc(zn2)]

    tc = time.time()

    def appendTrows(Xn, Yn, Zn, row):
        row = calcTrow(Xn, Yn, Zn, row)
        return row

    buildMat = [dd.delayed(appendTrows)(Xn, Yn, Zn, rxLoc[ii, :]) for ii in range(rxLoc.shape[0])]

    T = np.squeeze(np.vstack(dd.compute(buildMat)))

    print("Dask CPU Time: " + str(time.time()-tc) + " s")

    x = np.random.randn(T.shape[1])

    tc = time.time()
    # Try dot product
    np.dot(T, x)
    print("Numpy Dot CPU Time: " + str(time.time()-tc) + " s")

    tc = time.time()
    # Try dot product
    da.dot(T, x)
    print("Dask Dot CPU Time: " + str(time.time()-tc) + " s")



