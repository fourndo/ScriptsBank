function [nullcell,topocell,ztopo_n] = topocheck_surf(Xn,Yn,Zn,topo)
% [nullcell,topocell,toponode] = topocheck(Xn,Yn,Zn,topo)
% Create active cell matrix from discretized topography and topocell for
% all cells intersected by the toposurface
%
% TOPOCHECK SURF
% nullcell only takes the top layer below topography to be active
% Inputs//
% Xn, Yn, Zn: 3D arrays for the X, Y and Z location of all nodes in mesh
% topo: Topography array 3-by-points [x(:) y(:) z(:)]
%
% Output//
% nullcell: 1D vector, in UBC format, of binary values 1 (active) or
% 0(inactive) depending if below or above topography
% 
% topocell: 1D vector, in UBC format, of binary values 1 (yes) or 0 (no)
% if cell is intersected by topo surface
% 

%% FOR DEV
% clear all
% close all
% 
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Documents\Projects\Research\Modelling\Topo_adjust';
% 
% meshfile = 'Mesh_5m.msh';
% topofile = 'Gaussian.topo';

%% SCRIPT START HERE
nxn = size(Xn,2);
nyn = size(Xn,3);
nzn = size(Xn,1);

nx = nxn-1;
ny = nyn-1;
nz = nzn-1;

Xn = Xn(:);
Yn = Yn(:);
  
% Pre-allocate space
nullcell = ones(nz,nx,ny);
topocell = zeros(nz,nx,ny);

if sum(topo(:,3))==0
    
    fprintf('Topo file name is empty\n')
    fprintf('Nullcell model of all ones created\n')
    ztopo_n = ones(nzn,nxn,nyn)*z0;
    nullcell = reshape(nullcell,nx*ny*nz,1);
       return 
end
    
% Start looking for elevation of topography at each cell
%     fprintf('Start computing topography over the mesh ... \n')

% Create topo surface
T = scatteredInterpolant( topo(:,1) , topo(:,2) , topo(:,3) );

% Querry topo on horizontal vertices
ztopo_n = T(Xn,Yn);
ztopo_n = reshape(ztopo_n,nzn,nxn,nyn);    

% Look at 8 corner of each cell and form a logical matrices 
% depending on their location with respect to topography
% below=1 , above = 0;
N = (Zn(1:end-1,1:end-1,1:end-1) <= ztopo_n(1:end-1,1:end-1,1:end-1))*1;
N = N + (Zn(2:end,1:end-1,1:end-1) <= ztopo_n(2:end,1:end-1,1:end-1))*1;
N = N + (Zn(1:end-1,2:end,1:end-1) <= ztopo_n(1:end-1,2:end,1:end-1))*1;
N = N + (Zn(2:end,2:end,1:end-1) <= ztopo_n(2:end,2:end,1:end-1))*1;    
N = N + (Zn(1:end-1,1:end-1,2:end) <= ztopo_n(1:end-1,1:end-1,2:end))*1;
N = N + (Zn(2:end,1:end-1,2:end) <= ztopo_n(2:end,1:end-1,2:end))*1;
N = N + (Zn(1:end-1,2:end,2:end) <= ztopo_n(1:end-1,2:end,2:end))*1;
N = N + (Zn(2:end,2:end,2:end) <= ztopo_n(2:end,2:end,2:end))*1;

% Create nullcell and topocell
% If N==8, then all nodes of cell are below topo
% If N~=8 && N~=0, then topo intersect cell
nullcell = (N==8)*1; 

for ii = 1: nx
    
    for jj = 1:ny
        
        index = find(nullcell(:,ii,jj)==1);
        
        if length(index)~=1
            nullcell(index(2):end,ii,jj)=0;
        end
        
    end
end
topocell = find((N~=8)&(N~=0));
    
nullcell = reshape(nullcell,nx*ny*nz,1);
