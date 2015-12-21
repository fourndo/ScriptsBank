% function topo_2_idx(topofile,meshfile,out_dir)
% topo_2_idx(topofile,meshfile,out_dir)
% Function to convert topocheck file to UBC-topo.idx format.
%
% INPUT
% topofile: ...\topo_model.txt (topocheck file)
% meshfile: Mesh file in UBC format associated to the topofile
% out_dir : Output directory where the topo.idx will be created

%% FOR DEV ONLY
clear all
close all
out_dir = 'C:\Projects\4180_Wallbridge_Wisner_IP\Data\DataProcessing\DHIP\FWR';
topofile = [out_dir '\topo_model.txt'];
meshfile = [out_dir '\DHIP_BrokenHammer.msh'];
model = load('C:\Projects\4180_Wallbridge_Wisner_IP\Data\DataProcessing\DHIP\FWR\Inv\dcinv3d_03.con');

%% SCRIPT BEGGING HERE
addpath C:\Users\dominiquef\Dropbox\Master\Miscellaneous

% Load meshfile
mesh = get_UBC_mesh(meshfile);
nx = mesh(1,1);
ny = mesh(1,2);
nz = mesh(1,3);

% Reshape topofile and find the number cells above topography, then write
% to Topo.idx format

% topomodel = load(topofile);
topomodel = model;
topomodel(model==1e-8) = 0;
topomodel(model~=1e-8) = 1;
topomodel = reshape(topomodel,nz,nx,ny); 

fid = fopen([out_dir '\topo.dat'],'w');
fprintf(fid,'%i %i\n',nx,ny);

for ii = 1 : ny
    
    for jj = 1 : nx
        
        count = 0;
        
        for kk = 1 : nz
            
            if topomodel(kk,jj,ii)==0
                
                count = count+1;
                
            end
            
        end
        
        fprintf(fid,'%i %i %i\n',jj,ii,count);
    end
    
end

fclose(fid);

topomodel = reshape(topomodel,nx*ny*nz,1);
save([out_dir '\topomodel.dat'],'-ascii','topomodel');