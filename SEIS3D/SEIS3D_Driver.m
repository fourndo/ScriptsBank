% Seismic refraction program for given 3D velocity model.

clear all
close all

addpath ..\FUNC_LIB\;
% Project folders
work_dir    = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\UBC_backup\EOSC350\Modeling';
meshfile    = 'Mesh_0p5m.msh';
modelfile   = 'Tumuli.dat';
recfile     = 'Receivers.dat';

% Load mesh
[xn,yn,zn] = read_UBC_mesh([work_dir '\' meshfile]);
dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);

% Load velocity model
m = load([work_dir '\' modelfile]); 

% Get active cells from air
nullcell = m ~= -100;

% Load receiver location
rx = load([work_dir '\' recfile]);

% Define source location
q = [mean(xn) mean(yn)];

% Need to define snell's law


% Need kernel operator
