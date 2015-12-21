function indx = null_2_idx(nullcell)
% topo_2_idx(topofile,meshfile,out_dir)
% Function to convert topocheck file to UBC-topo.idx format.
%
% INPUT
% topofile: ...\topo_model.txt (topocheck file)
% meshfile: Mesh file in UBC format associated to the topofile
% out_dir : Output directory where the topo.idx will be created

%% FOR DEV ONLY
% clear all
% close all
% out_dir = 'C:\Projects\4180_Wallbridge_Wisner_IP\Data\DataProcessing\DHIP\FWR';
% topofile = [out_dir '\topo_model.txt'];
% meshfile = [out_dir '\DHIP_BrokenHammer.msh'];
% model = load('C:\Projects\4180_Wallbridge_Wisner_IP\Data\DataProcessing\DHIP\FWR\Inv\dcinv3d_03.con');

%% SCRIPT BEGGING HERE
% addpath C:\Users\dominiquef\Dropbox\Master\Miscellaneous

% Load meshfile
% mesh = get_UBC_mesh(meshfile);
nx = size(nullcell,2);
ny = size(nullcell,3);
nz = size(nullcell,1);

% Reshape topofile and find the number cells above topography, then write
% to Topo.idx format

% topomodel = load(topofile);
idz = zeros(nx,ny);

for ii = 1 : ny
    
    for jj = 1 : nx
        
        count = 1;
        
        while nullcell(count,jj,ii)==0 && count < nz
            
            count = count + 1;
            
        end
        
        idz(jj,ii) = count;
    end
    
end

[idx,idy] = ndgrid(1:nx,1:ny);

indx = [idz(:) idx(:) idy(:)];
