'''
Created on Jul 17, 2013

@author: dominiquef

MAG_lp_Driver_PyDev_v1.0
Author: D Fournier
Last Update: September 27th, 2015

This program sets up the magnetic forward operator for the integral formulation. 


def get_T_mat(xn,yn,zn):
'''
from SimPEG import Maps, Survey, Utils, np, sp
from scipy.constants import mu_0

'''   
Old scripts     
from get_UBC_mesh import get_UBC_mesh
from read_input import read_input
from read_UBC_log_NEWCODE import read_UBC_log_NEWCODE
from numpy import *
from numpy.linalg import *
from comp_dmdx import comp_dmdx
from comp_dmdy import comp_dmdy
from comp_dmdz import comp_dmdz
from make_depthW import make_depthW
from get_input_files import get_input_files
from get_Ws3D import get_Ws3D
'''

nx = np.size(y);
ny = np.size(y);
nz = np.size(y);

Tx = zeros((1,3*mcell), dtype=float);