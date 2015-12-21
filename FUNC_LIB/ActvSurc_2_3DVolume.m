% Active surface to 3D volume

clear all
close all

work_dir = 'C:\LC\Private\dominiquef\Projects\3796_AGIC_Research\HyperCube\Processing\MtDore_files_for_Hypercube\Targeting_2015';

meshfile = 'CEM_to_2p5km.msh';
modelfile = 'Pseudo_litho_3D.dat';
actvfile = 'Mesh100m_activecell.dat';

ndv = -99999;

% lOad mesh
[xn,yn,zn] = read_UBC_mesh([work_dir '\' meshfile]);
dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);


% Load model and active cells
m = load([work_dir '\' modelfile]);
m = reshape(m,nz,nx,ny);

% actv = load([work_dir '\' actvfile]);
% actv = reshape(actv,nz,nx,ny);

for ii = 1 : nx
    
    for jj = 1 : ny
        
        idx = m(:,ii,jj)~=ndv;

        if sum(idx)~=0
            m(idx==0,ii,jj) = m(idx==1,ii,jj);
        end
        
    end

end

% m = m(actv==1);

save([work_dir '\' modelfile],'-ascii','m');