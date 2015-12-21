function [Q] = make_EM1D_Q_3D(work_dir,meshfile,nullcell,xyz)
% Function make_1D_P_3D
% Create querry matrix (Q) to transfer 1D models information to 3D UBC-mesh
%
% INPUT
% work_dir  : directory of files
% meshfile  : Mesh file in UBC format
% nullcell  : Binary topo model in UBC format
% obsfile   : 1D Data used for EM1DFM in UBC format
%
% SUB-FUNCTION CALLS
% get_UBC_mesh.m    : Load mesh file and output [5xn] matrix
% load_EM1DFM_obs   : Load observation file in UBC-EM1DFM format

%% FOR DEV ONLY

% clear all 
% close all
% 
% work_dir    = 'C:\Users\dominiquef\Dropbox\DIGHEM\Processed_data\1DEMInversions\ALL_DF';
% meshfile    = 'UBC_mesh.msh';
% load([work_dir '\nullcell.dat']);
% obsfile = 'DIGHEM_TKC_ALL.obs';
% nnodes=4;

%% SCRIPT STARTS HERE
%% Load observation data file
% data = load_EM1DFM_obs(work_dir,obsfile);

% Load nullcell file
topo = nullcell;

% Create observation vectors
obsx = xyz(:,1);
obsy = xyz(:,2);
% obsz = xyz(:,3);

% Extract data location
nsnds  = length(obsx);

%% Load mesh
[meshfile]=get_UBC_mesh([work_dir '\' meshfile]);
nx = meshfile(1,1); %size(X,1);    %number of cell in X
ny = meshfile(1,2); %size(X,2);    %number of cell in Y
nz = meshfile(1,3); %size(X,3);    %number of cell in Z

mcell = nx*ny*nz;

% Cell size array
dx = meshfile(3,1:nx);
dy = meshfile(4,1:ny);
dz = meshfile(5,1:nz);

% Corner grid
x0 = meshfile(2,1);
y0 = meshfile(2,2);
z0 = meshfile(2,3);

% Create 3D cell center array
xc = x0 + cumsum(dx) - dx/2;
yc = y0 + cumsum(dy) - dy/2;
zc = z0 - cumsum(dz) + dz/2;

[Xc,Yc] = ndgrid(xc,yc);

topo = reshape(topo,nz,nx,ny);

%% Create Projector and Querry matrices
% Pre-allocate memory for Querry matrix
% Store the index of the first cell in the 3D mesh for each station
% location and index of topography in column i,j. Information will be used
% later for running the 1D inversions and pass model values from 3D to 1D
% Q = [K,I,J,ktopo]

Q = zeros(nsnds,3);
for ss = 1 : nsnds
    
    r = ( ( Xc - obsx(ss) ).^2 +...
                ( Yc - obsy(ss) ).^2 ) .^ 0.5;

    % Sort points by distance and keep the smallest
    [I,J] = find(min(min(r))==r);
    I = I(1);
    J = J(1);
    
    Q(ss,1:2) = [I J];
    
    % Look at column of cells in topo model and find first air-ground cell
    K = 0;    
    for ii = 1 : nz
        
        if topo(ii,I,J)==1
            
            break;
            
        end
        
    end
%     K = K + 1;
    
    Q(ss,3) = ii;
    
%     % Move K-index up until reaches the station location vertically 
%     while zc(K)- zc(Q(ss,4)) < obsz(ss)
%         
%         K = K-1;
%         
%     end
%     K = K + 1;
%     
%     Q(ss,1) = K;
    
end

