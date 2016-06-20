% Function copy_cells_up
%
% Reads in a mesh and a model and propagate the top cells all the way to
% the top of the mesh.


clear all
close all

%% USER INPUT
modelfile   = 'MAG3D_TMI_lplq.sus';
work_dir    = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Block_Gaussian_topo';
meshfile    = 'Mesh_10m.msh';

% Specify the no-data value in the original model (-100:Mag , 1e-8:Cond)
ndv = -100;

out_file    = 'Model_FullSpace.dat';
%% SCRIPT STARTS HERE
mesh=get_UBC_mesh([work_dir '\' meshfile]);


[xn,yn,zn] = read_UBC_mesh([work_dir '\' meshfile]);
dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);

% Load model
m = load([work_dir '\' modelfile]);
m = reshape(m,nz,nx,ny);

% Loop the XY plane and grab one column        
for ii = 1:nx;
    
    for jj = 1:ny
        
       airc = m(:,ii,jj)~=ndv;
       
       % If there is no ground cells, then leave it empty
       if isempty(airc)
           
           m(:,ii,jj) = ndv;
           
       % Otherwise copy the last cell up      
       else 
           
           ground = m(airc,ii,jj);
           m(airc==0,ii,jj) = ground(1);
           
       end
       
    end
end

%Output model
m = m(:);

save([work_dir '\' out_file],'-ascii','m');