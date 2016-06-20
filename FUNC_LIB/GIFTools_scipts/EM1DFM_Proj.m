% Test EM1DFM proj creation
clear all

inp_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\DIGHEM\1D';
dsep = '\';

% Create project
proj = GIFproject;

proj.setWorkDir(inp_dir)

% Load EM1D data
EM1DFMinversion.readDobs(proj,[inp_dir dsep 'DIGHEM_TKC_DO27_FLT_25m.obs'])