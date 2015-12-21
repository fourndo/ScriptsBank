% MAG3D_Tiles_Driver
%
% Inversion code for MAG3D style inversion. The large problem is divided in
% small inversion tiles with overlap.
% Sub-function calls:
% MAGSEN_Tiled: Compute sensitivity matrix for sub-problem
% MAGINV_Tiled: Run the inversion code MAG3D style and save model
% MERGE_Tiles: Merge all models to the underlying mesh
%
% Inversion code with lp,lq norm for 3D magnetostatic problem
% Read input file and proceede with inversion
% 
% DEVELOPMENT CODE
% Other sub-functions are required to run the code
% Written by: Dominique Fournier 
% Last update: June 6th, 2015
clear all
close all

addpath '.\func_lib';

dsep = '\';
% Project folders
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Nut_Cracker\Tiled_AMI\All_Tiles';
% out_dir = [work_dir dsep 'ALL_Tiles'];
% mkdir(out_dir);


% meshfile= '..\..\Mesh_50m_v2.msh';
% obsfile = '..\..\Obs_Paul_Lake_SUB_1pc_5nT.dat';
% magfile = 'Merged_M_model.fld';
topofile = [];

meshfile= '..\Mesh_20m.msh';
obsfile = '..\Obs_REM_GRID_TMI.obs';
magfile = 'Merged_M_model.fld';

% mtrue = load([work_dir '\..\Effec_sus.sus']);
% Load mesh file and convert to vectors (UBC format)
[xn,yn,zn] = read_UBC_mesh([work_dir dsep meshfile]);
dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);

%% 3D vertical nodal location
[Zn,Xn,Yn] = ndgrid(zn,xn,yn);

if isempty(topofile)==1
    
    nxn = length(xn);
    nyn = length(yn);
    topo = [reshape(Xn(1,:,:),nxn*nyn,1) reshape(Yn(1,:,:),nxn*nyn,1) reshape(Zn(1,:,:),nxn*nyn,1)+1e-8];
    
else
    % Load topo
    topo = read_UBC_topo([work_dir '\' topofile]);
end

% Create nullcell
[nullcell,tcellID,ztopo_n] = topocheck(Xn,Yn,Zn,topo);
% save([work_dir '\nullcell.dat'],'-ascii','nullcell');

% Get index of active cells
cellID = find(nullcell==1);

% Load susceptibility model
mcell = nx*ny*nz;

% Get nodal discretization for octree levels
% celln = MAG3C_RTC_OctNodes(Xn,Yn,Zn,cellID,ztopo_n,0);
celln = zeros(mcell,1,6);

celln(:,1,1) = kron(ones(1,ny),kron(ones(1,nx),[zn(1:end-1)]));
celln(:,1,2) = kron(ones(1,ny),kron([xn(1:end-1)],ones(1,nz)));
celln(:,1,3) = kron([yn(1:end-1)],kron(ones(1,nx),ones(1,nz)));

celln(:,1,4) = kron(ones(1,ny),kron(ones(1,nx),[zn(2:end)]));
celln(:,1,5) = kron(ones(1,ny),kron([xn(2:end)],ones(1,nz)));
celln(:,1,6) = kron([yn(2:end)],kron(ones(1,nx),ones(1,nz)));




% Create selector matrix for active cells
X = spdiags(nullcell,0,mcell,mcell);
X = X(nullcell==1,:);

%%
% Load observation loc  
[H, I, Dazm, D, obsx, obsy, obsz, d, wd] = read_MAG3D_obs([work_dir '\' obsfile]);
% plot_mag3C(obsx,obsy,oldB,I,D,'Observed 3C-Data data')
% obsx = -222.5:5:-50;
% obsy = -162.5:5:-50;
% [obsx,obsy] = ndgrid(obsx,obsy);
% 
% obsx = reshape(obsx,size(obsx,1) * size(obsx,2),1);
% obsy = reshape(obsy,size(obsx,1) * size(obsx,2),1);
% 
% obsz = obsz - 25;

ndata = length(obsx);

% Load magnetization model
if isempty(magfile) == 1
    
    mag_model = [ones(mcell,1)*Dazm ones(mcell,1)*I];
    
else
%     mag_model = [ones(mcell,1)*Dazm ones(mcell,1)*I];
    mag_model = load([work_dir '\' magfile]);
    
end

mag_model = X*mag_model;
mcell = sum(nullcell);

M = H*[mag_model(:,1);mag_model(:,2);mag_model(:,3)];

%% Forward model data
TMI = [(cosd(I) * cosd(D)) (cosd(I) * sind(D)) sind(I)];
d_pre = [];


% matlabpool OPEN 12
% Create sub-problems
nworkers = 1;
nsub = round(ndata/nworkers);

for jj = 1 
    
    dsub = zeros(nsub,1);
    
    if jj == nworkers
    
        indx = (1+(jj-1)*nsub):length(d);
        
    else
        
        indx = (1+(jj-1)*nsub):(jj*nsub);
        
    end
    %%
    progress = -1;
    tic 
    
    for ii = 1:length(indx)

        
        % compute kernel for active cells
        [Tx,Ty,Tz] = MAG3C_T(obsx(indx(ii)),obsy(indx(ii)),obsz(indx(ii)),celln);

        dsub(ii) = TMI * [Tx;Ty;Tz] * M;

        d_iter = floor(ii/nsub*100);
        if  d_iter > progress

            fprintf('Computed %i pct of data in %8.5f sec\n',d_iter,toc)
            progress = d_iter;

        end

    end
    
    write_MAG3D_TMI([work_dir '\PRED_TMI' num2str(jj) '.pre'],H,I,Dazm,...
    obsx(indx),obsy(indx),obsz(indx),dsub',wd(indx));
    
    write_MAG3D_TMI([work_dir '\OBS_TMI' num2str(jj) '.pre'],H,I,Dazm,...
    obsx(indx),obsy(indx),obsz(indx),d(indx),wd(indx));
end
% matlabpool CLOSE


% plot_TMI(obsx,obsy,d,d_pre,wd,'Observed vs Predicted Magnitude');