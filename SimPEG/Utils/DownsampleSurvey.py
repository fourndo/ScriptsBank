# -*- coding: utf-8 -*-
"""
Script to subsample a survey with distance

Created on Wed Sep  6 09:57:20 2017

@author: DominiqueFournier
"""

from SimPEG import Utils
from SimPEG.Utils import mkvc
import SimPEG.PF as PF
import numpy as np
import matplotlib.pyplot as plt

# # USER INPUTS # #
workDir = "C:\\Users\\DominiqueFournier\\Desktop\\Craig"
dFile = 'mesh_trimmed_data.obs'
dType = 'MAG'
method = ('radius', 50)

dFileOut = 'mesh_trimmed_data_Filter50m.obs'

# # SCRIPT STARTS HERE # #
if dType == 'MAG':
    survey = PF.Magnetics.readMagneticsObservations(workDir + '\\' + dFile)
elif dType == 'GRAV':
    survey = PF.Gravity.readUBCgravObs(workDir + '\\' + dFile)

else:
    assert dType in ['MAG', 'GRAV'], "dType must be 'MAG' or 'GRAV'"

# Downsample the survey using specified method
assert method[0] in ['radius', 'random'], "Downsample method should be 'radius' or 'random' "


locXYZ = survey.srcField.rxList[0].locs

if method[0] == 'radius':

    nstn = locXYZ.shape[0]
    # Initialize the filter
    indx = np.ones(nstn, dtype='bool')
    for ii in range(nstn):

        if indx[ii]:

            rad = ((locXYZ[ii, 0] - locXYZ[:, 0])**2 +
                   (locXYZ[ii, 1] - locXYZ[:, 1])**2)**0.5

            indx[rad < method[1]] = False
            indx[ii] = True


elif method[0] == 'random':

    nD = int(survey.nD*0.5)
    print("nD ratio:" + str(nD) + '\\' + str(survey.nD))
    indx = np.random.randint(0, high=survey.nD, size=nD)


# Create a new downsampled survey
if dType == 'MAG':

    rxLoc = PF.BaseGrav.RxObs(locXYZ[indx, :])
    srcField = PF.BaseMag.SrcField([rxLoc], param=survey.srcField.param)
    survey_dwnS = PF.BaseMag.LinearSurvey(srcField)
    survey_dwnS.dobs = survey.dobs[indx]
    survey_dwnS.std = survey.std[indx]

    PF.Magnetics.writeUBCobs(workDir + '\\' + dFileOut, survey_dwnS)

elif dType == 'GRAV':

    rxLoc = BaseGrav.RxObs(locXYZ[indx, :])
    srcField = BaseGrav.SrcField([rxLoc])
    survey_dwnS = BaseGrav.LinearSurvey_dwnS(srcField)
    survey_dwnS.dobs = survey.dobs[indx]
    survey_dwnS.std = survey.std[indx]

    PF.Gravity.writeUBCobs(workDir + '\\' + dFileOut, survey_dwnS)
