clear all
close all

modelfile = 'TMI_p2q2l1_iter_32.sus';
work_dir     = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Research\Modelling\Topo_adjust\Gaussian_topo_plane_model';
meshfile    = 'Mesh_5m.msh';

mesh=get_UBC_mesh([work_dir '\' meshfile]);
ndv = -100;

nx = mesh(1,1); %size(X,1);    %number of cell in X
ny = mesh(1,2); %size(X,2);    %number of cell in Y
nz = mesh(1,3); %size(X,3);    %number of cell in Z

dx = mesh(3,1:nx);
dy = mesh(4,1:ny);
dz = mesh(5,1:nz);

x0 = mesh(2,1);
y0 = mesh(2,2);
z0 = mesh(2,3);

% Load model
m = load([work_dir '\' modelfile]);
m = reshape(m,nz,nx,ny);

for ii = 1:nx;
    for jj = 1:ny
        
       airc = sum(m(:,ii,jj)==ndv);
       
       if airc~=nz
           m(1:airc,ii,jj)=m(airc+1,ii,jj);
       end
       
    end
end

% model(nullcell==0) = 0;
% 
% model = reshape(model,ny*nx*nz,1);
m = m(:);

save([work_dir '\' modelfile 'noair'],'-ascii','m');