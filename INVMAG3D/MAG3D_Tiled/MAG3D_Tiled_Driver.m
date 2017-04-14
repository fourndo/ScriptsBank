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

addpath 'func_lib';

dsep = '\';
% Project folders
% work_dir = '/tera_raid/dfournier/TKC/Mag/DIGHEM_lp';
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Dipping_Prism_Li';
% work_dir = 'C:\LC\Private\dominiquef\Projects\Gervais_Solado\Processing\Inversion';
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Nut_Cracker\Tiled_MAG3D';
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Triple_Block_lined';
work_dir = 'C:\Users\DominiqueFournier\ownCloud\Research\Modelling\Synthetic\Dipping_Plane';

out_dir = [work_dir dsep 'ALL_Tiles'];
mkdir(out_dir);

inpfile   = 'MAG3D_Tile.inp'; 

[meshfile,obsfile,topofile,mstart,mref,magfile,weightfile,chi_target,alphas,beta,bounds,norm_vec,eps_p,eps_q,FLAG1,FLAG2,Rot,tilefile] = MAG3D_Tile_read_inp([work_dir dsep inpfile]);

% mtrue = load([work_dir '\..\Effec_sus.sus']);
% Load mesh file and convert to vectors (UBC format)
[xn,yn,zn] = read_UBC_mesh([work_dir dsep meshfile]);
dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);

xc = (xn(2:end) + xn(1:end-1))/2;
yc = (yn(2:end) + yn(1:end-1))/2;
zc = (zn(2:end) + zn(1:end-1))/2;

[Zc,Xc,Yc] = ndgrid(zc,xc,yc);

mcell = (length(xn)-1) * (length(yn)-1) * (length(zn)-1);

expanfac = 1.4;
npadx = 8;
npady = 8;

R = 3;
R0 = 0;

azm = 0;
dip = 0;

%% Load observation file (3C UBC-MAG format)
[H, BI, BD, MI, MD, Obsx, Obsy, Obsz, data, wd_full] = read_MAG3D_obs([work_dir dsep obsfile]);
% plot_mag3C(obsx,obsy,d,I,D,'Observed 3C-data')
% plot_TMI(obsx,obsy,d,d,wd,'Observed vs Predicted Magnitude');

% % Filter data
% flt_r = 10;
% [indx] = Filter_xy(Obsx,Obsy,flt_r);
% 
% Obsx = Obsx(indx);
% Obsy = Obsy(indx);
% Obsz = Obsz(indx);
% data = data(indx);
% wd_full = wd_full(indx);
% 
ndata = length(data);

% write_MAG3D_TMI([work_dir dsep obsfile(1:end-4) '_FLT' num2str(flt_r) '.dat'],H, BI, BD, MI, MD,Obsx,Obsy,Obsz,data,wd_full);
%% Load input models and weights
% Load magnetization model
if isempty(magfile) == 1
    
    mag_azmdip = [ones(mcell,1)*MD ones(mcell,1)*MI];
    mag_xyz = azmdip_2_xyz( mag_azmdip(:,1) , mag_azmdip(:,2) );
    
else
    
    mag_xyz = load([work_dir dsep magfile]);
    
end

% Load weights model
if isempty(weightfile) == 1
    
    w = ones(4*mcell,1);
    
else
    
    w = load([work_dir dsep weightfile]);
    
end

% Create or load reference model
if ischar(mref)==1
    
    mref = load([work_dir dsep mref]);
    
else
    
    mref = ones(mcell,1)*mref;
    
end

% Create or load reference model
if ischar(mstart)==1
    
    mstart = load([work_dir dsep mstart]);
    
else
    
    mstart = ones(mcell,1)*mstart;
    
end

%% Load topography
if isempty(topofile)==1
    [Zn,Xn,Yn] = ndgrid(zn,xn,yn);
    nxn = length(xn);
    nyn = length(yn);
    Topo = [reshape(Xn(1,:,:),nxn*nyn,1) reshape(Yn(1,:,:),nxn*nyn,1) reshape(Zn(1,:,:),nxn*nyn,1)+1e-8];
    clear Zn Yn Xn
else
    % Load topo
    Topo = read_UBC_topo([work_dir dsep topofile]);
end

%% Generate s and t vectors for the different lp zones
if ischar(norm_vec)==1
    
    lpmat = load([work_dir dsep norm_vec]);
    [s,LP] = find_zones(lpmat);

    % Smooth out the regions with 8-point averager
    % Power determine the transition length
    A = get_AVG_8pt(dx,dy,dz);
    A = A*(A*A);
    A = spdiags(1./sum(A,2),0,mcell,mcell) *A;
    
    trans = A*s;
    t = trans;

    
else
    
    LP = norm_vec;
    t= ones(mcell,1);

end

%% Load in tile files and loop through the inversions
% Format of tile file is [xo, yo, x_max, y_max]
% Each row defines a seperate tile
if isempty(tilefile)
    ntiles = 1;%T = [xn(1) yn(1) xn(end) yn(end)];
else
    T = load([work_dir dsep tilefile]);
    ntiles = size(T,1);
end


%parpool(3)

for ii = 1 : ntiles

    % Define working directory for tile
    tile_dir = [work_dir dsep 'Tile' num2str(ii)];
    mkdir(tile_dir);
    
    if ~isempty(tilefile)
        % Create core mesh from tile extent
        xn_t = T(ii,1):min(dx):T(ii,3);
        yn_t = T(ii,2):min(dy):T(ii,4);
        zn_t = zn;

        % Grab data within the mesh boundary
        indx = Obsx >= xn_t(1) & Obsx <= xn_t(end) & Obsy >= yn_t(1) & Obsy <= yn_t(end);

        % Add padding cells
        padx = min(dx)*expanfac.^(0:npadx);
        pady = min(dy)*expanfac.^(0:npady);

        xn_t = [(xn_t(1) - fliplr(cumsum(padx))) xn_t (xn_t(end) + cumsum(padx))];
        yn_t = [(yn_t(1) - fliplr(cumsum(pady))) yn_t (yn_t(end) + cumsum(pady))];
        
        mcell = (length(xn_t)-1) * (length(yn_t)-1) * (length(zn_t)-1);
        
        w_t = ones(4*mcell,1);
        mref_t = zeros(mcell,1);
        mstart_t = ones(mcell,1)*1e-4;
        t = ones(mcell,1);
    else
        
        xn_t = xn;
        yn_t = yn;
        zn_t = zn;
        
        w_t = w;
    
        mref_t = mref;

        mstart_t = mstart;

        %t = ones(mcell_t,1);
    
        indx = ones(ndata,1) == 1;
        
    end
    
    obsx = Obsx(indx);
    obsy = Obsy(indx);
    obsz = Obsz(indx);
    d = data(indx);
    wd = wd_full(indx); 
    
    % Create cell center location for interp
    xc_t = (xn_t(2:end) + xn_t(1:end-1))/2;
    yc_t = (yn_t(2:end) + yn_t(1:end-1))/2;
    zc_t = (zn_t(2:end) + zn_t(1:end-1))/2;
        
    dx_t = xn_t(2:end) - xn_t(1:end-1); 
    dy_t = yn_t(2:end) - yn_t(1:end-1); 
    dz_t = zn_t(1:end-1) - zn_t(2:end); 

    write_MAG3D_TMI([tile_dir dsep 'Tile_data.dat'],H, BI, BD, MI, MD,obsx,obsy,obsz,d,wd)
    write_UBC_mesh(out_dir,['Tile' num2str(ii) '.msh'],dsep,xn_t(1),yn_t(1),zn_t(1),dx_t,dy_t,dz_t)
    write_UBC_mesh(tile_dir,['Tile' num2str(ii) '.msh'],dsep,xn_t(1),yn_t(1),zn_t(1),dx_t,dy_t,dz_t)
    
    mcell_t = length(xc_t) * length(yc_t) * length(zc_t);
    
%     [Zc_t,Xc_t,Yc_t] = ndgrid(zc_t,xc_t,yc_t);

    % Take a subset of topography over the mesh
    indx = Topo(:,1) >= xn_t(1)-100 & Topo(:,1) <= xn_t(end)+100 & Topo(:,2) >= yn_t(1)-100 & Topo(:,2) <= yn_t(end)+100;
    topo = Topo(indx,:);
    
    % Interpolate input models to mesh (cell-center)
%     indx = xc >= xc_t(1) & xc <= xc_t(end);
%     indy = yc >= yc_t(1) & yc <= yc_t(end);
%     indz = zc <= zc_t(1) & zc >= zc_t(end);
    
    %% NEED TO FIX THE TILED INTERPOLATION   
    [nullcell,tcell,~] = topocheck(xn_t,yn_t,zn_t,topo+1e-5);
    save([tile_dir dsep 'nullcell.dat'],'-ascii','nullcell');

    mactv = sum(nullcell);
    
    mag_azmdip = [ones(mactv,1)*MD ones(mactv,1)*MI];
    mag_xyz_t = azmdip_2_xyz( mag_azmdip(:,1) , mag_azmdip(:,2) );
    
    
    M = [spdiags(H * mag_xyz_t(:,1),0,mactv,mactv);spdiags(H * mag_xyz_t(:,2),0,mactv,mactv);spdiags(H * mag_xyz_t(:,3),0,mactv,mactv)];


    %% CALL MAGSEN_Tile and create sensitivity
    G = MAGSEN_Tiled(tile_dir,dsep,xn_t,yn_t,zn_t,H, BI, BD, MI, MD, obsx, obsy, obsz,nullcell, 'NONE', 'G', R, R0, M);
    
    %% CALL MAGINV_Tile and run inversion
    MAGINV_sus_Tiled(tile_dir,out_dir,dsep,ii,xn_t,yn_t,zn_t,H, BI, BD, MI, MD, obsx, obsy, obsz, G, d, wd,w_t,mref_t,mstart_t,bounds,chi_target,alphas,beta,LP,eps_p,eps_q,t,FLAG1,FLAG2,Rot)
    
end
% 
% % matlabpool CLOSE
% %% Merging results back to core grid
% TileMerge(out_dir,['..' dsep meshfile],500,-100,dsep,'topofile',['..' dsep topofile],'rem_pad',[7 7])

%% Forward model the merged model
% MAGFWR_Tiled(work_dir,meshfile,obsfile,['ALL_Tiles' dsep 'VOXMERGE_model.dat'],topofile,dsep)