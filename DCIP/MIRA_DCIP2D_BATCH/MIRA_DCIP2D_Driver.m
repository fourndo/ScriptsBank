% MIRA - DCIP2D Batch Inversion 
% This program uses a series of sub-routines to extract data from the 
% TQIOdb output files and invert using the UBC code.
% Input arguments:
% 'all': Process all lines
% or
% '## ##' : Process only lines ## (sequential order from the top)
clear all
close all

addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB
addpath functions

%% Extract files and create directories
work_dir = 'C:\LC\Private\dominiquef\Projects\4329_Goldcorp_Wabamisk_DCIP3D\Modelling\Inversion\2D\Geosig_East';
input_dir = 'C:\LC\Private\dominiquef\Projects\4329_Goldcorp_Wabamisk_DCIP3D\Data\Processing\TMC_West';
output_dir = 'C:\LC\Private\dominiquef\Projects\4329_Goldcorp_Wabamisk_DCIP3D\Modelling\Inversion\2D';

lines = 'all';
% lines = [16 29];
% lines = [5 11];
% lines = num2str(lines);
% Different pre-formating routines are available

% For Titan type
% TQIP_2_folder(input_dir,work_dir,0.02,1e-2,0.02,1e-4,lines);
% UBC_2_folder_FWR_BCK(input_dir,work_dir,'all');
% CSV_2_UBC(input_dir,work_dir)

% For output from *.gdb files in format <line, T1X, T2X, R1X, R2X, data>
% gdb_2_UBC(input_folder,raw_folder,'TECK_Data_gdb_DC.dat','TECK_Data_gdb_IP.dat')
% UBC_2_folder(edt_folder,input_dir,lines);
%% Run DC inversions
% DC_Batch_Inv(work_dir,lines);

%% Run IP inversions
% IP_Batch_Inv(work_dir,lines);

%% Re-compute standard deviations from the log files

% DCIP_Batch_STD(work_dir,lines);



%% Georeference data
% 1- Extract survey information (format might need to change)

% survey_full = get_survey_v2('C:\LC\Private\dominiquef\Projects\4329_Goldcorp_Wabamisk_DCIP3D\Modelling\Inversion\2D\Geosig_UTM.xyz');
% save('C:\LC\Private\dominiquef\Projects\4329_Goldcorp_Wabamisk_DCIP3D\Modelling\Inversion\2D\GeoSig_East_UTM.XYZ','-ascii','survey_full');
survey_full =load('C:\LC\Private\dominiquef\Projects\4329_Goldcorp_Wabamisk_DCIP3D\Modelling\Inversion\2D\GeoSig_East_UTM.XYZ');
Georef_data(survey_full,work_dir,lines)


%% Create depth weighting
% work_dir = 'C:\Projects\4180_Wallbridge_Wisner_IP\Modelling\Inversion\Titan\Inv2\Tile2';
% meshfile = [work_dir '\ROT40_Titan_25m_15m_Tile2.msh'];
% mfile = [work_dir '\dcinv3d_01.con'];
% 
% 
% layers = [16 8 4 2];
% 
% make_depthw(work_dir,meshfile,mfile,layers);


%% Invert three new models for computing the Depth of Investigation
% Cond_model(1) = 3e-4;
% Cond_model(2) = 6e-4;
% Charg_model(1) = 5e-3;
% 
% DCIP_Batch_DOI(Cond_model,Charg_model,output_folder,lines);
