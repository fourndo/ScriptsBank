% MAG3D_Tiles_Driver
%
% 3-D Magnetic Vector Inversion in spherical coordinates with sparsity
% constraints applied on the amplitude and angles.
% 
% Main Functions:
% MAGSEN_Tiled      : Compute sensitivity matrix for sub-problem
% MVI3D_PST_Tiled   : Run the MVI inversion in cartesian coordinates
% MVI3D_ATP_Tiled   : Run the MVI inversion in spherical coordinates
% MERGE_Tiles       : Merge all models to the underlying mesh
%
%
% DEVELOPMENT CODE
% Other sub-functions are required to run the code
% Written by: Dominique Fournier 
% Last update: June 6th, 2016

clear all
close all

addpath 'func_lib';

dsep = '\';
% Project folders
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Nut_Cracker\Tiled_CMI';
% work_dir ='C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\CraigModel\MAG';
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\SingleBlock\CMI';
% work_dir ='C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Triple_Block_lined';
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\MAG';
% work_dir = 'C:\LC\Private\dominiquef\Projects\4414_Minsim\Modeling\MAG\Composite\North_Main';
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\PLexample';
% work_dir ='C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Dipping_Prism_Li';
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Shared\KrisDom\BlockModel\SingleBlock';
% work_dir = 'C:\LC\Private\dominiquef\Projects\Gervais_Solado\Processing\Inversion';
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Shared\KrisDom\Nutcracker';
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Osborne\Inversion\ROT40\ATP';
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Desktop\Paul_Lake\Modeling\Inversion\Misery';
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\SingleBlock\Simpeg';
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Two_Blocks_test\SingleBlock';

out_dir = [work_dir dsep 'ALL_Tiles'];
mkdir(out_dir);

inpfile   = 'MAG3D_ATP_Tile.inp'; 

% Read input file
[meshfile,obsfile,topofile,mstart,mref,m_vec,chi_target,alphas,beta,bounds,lp_MVIvec,lp_tresh,weightfile,FLAG1,FLAG2,tilefile,ROT] = MVI_ATP_Tile_read_inp([work_dir dsep inpfile]);

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

% Define mesh parameters for local tile meshes
expanfac = 1.4;

npadx = 8;
npady = 8;

R = 3;
R0 = min(dx)/4;

%% Load observation file (3C UBC-MAG format)
[H, HI, HD, MI, MD, Obsx, Obsy, Obsz, data, wd_full] = read_MAG3D_obs([work_dir dsep obsfile]);
ndata = length(data);


%% Load input models and weights
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

% Create or load magnetization vector model
if ischar(m_vec)==1
    
    m_vec = load([work_dir dsep m_vec]);
    
    if size(m_vec,2) == 2
        
       m_vec = azmdip_2_xyz( m_vec(:,1) , m_vec(:,2) );
        
    end
    % Make sure the vectors are unity
%     temp = sqrt( sum( m_vec.^2 , 2 ) );
%     m_vec = spdiags( 1./ temp(:) ,0,mcell,mcell) * m_vec;
    
else
    
    m_vec = azmdip_2_xyz( ones(mcell,1)*MD , ones(mcell,1)*MI );
    
end

M_xyz = m_vec;
beta_in = 5e+6;
%% Load topography
if isempty(topofile)==1
    [Zn,Xn,Yn] = ndgrid(zn,xn,yn);
    nxn = length(xn);
    nyn = length(yn);
    Topo = [reshape(Xn(1,:,:),nxn*nyn,1) reshape(Yn(1,:,:),nxn*nyn,1) reshape(Zn(1,:,:),nxn*nyn,1)];
    
else
    % Load topo
    Topo = read_UBC_topo([work_dir dsep topofile]);
end

%% Generate t vectors for the different lp zones if provided a model of lp-norms
% Otherwise, all the same everywhere.
if ischar(lp_MVIvec)==1
    
    lpmat = load([work_dir dsep lp_MVIvec]);
    [s,LP] = find_zones(lpmat);

    % Smooth out the regions with 8-point averager
    % Power determine the transition length
    A = get_AVG_8pt(dx,dy,dz);
    A = A*(A*A);
    A = spdiags(1./sum(A,2),0,mcell,mcell) *A;
    
    t = A*s;

    
else
    
    LP = lp_MVIvec;
    t= ones(mcell,1);

end

%% Load in tile files and loop through the inversions
% Format of tile file is [xo, yo, x_max, y_max]
% Each row defines a seperate tile
if isempty(tilefile)
    T = [xn(1) yn(1) xn(end) yn(end)];
else
    T = load([work_dir dsep tilefile]);
end
ntiles = size(T,1);

% Can be runned in a parfor (in parallel)
% matlabpool OPEN 2

for ii = 1:ntiles
    
    addpath 'func_lib';
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
        
        obsx = Obsx(indx);
        obsy = Obsy(indx);
        obsz = Obsz(indx);
        
        
%         indx = kron(ones(3,1),indx);
        d = data(indx);
        wd = wd_full(indx); 
        
    else
        
        % If no tiles, than use the whole input mesh
        xn_t = xn;
        yn_t = yn;
        zn_t = zn;
        
        obsx = Obsx;
        obsy = Obsy;
        obsz = Obsz;
        d = data;
        wd = wd_full; 
        
    end
    

    
    % Create cell center location for interp
    xc_t = (xn_t(2:end) + xn_t(1:end-1))/2;
    yc_t = (yn_t(2:end) + yn_t(1:end-1))/2;
    zc_t = (zn_t(2:end) + zn_t(1:end-1))/2;
        
    dx_t = xn_t(2:end) - xn_t(1:end-1); 
    dy_t = yn_t(2:end) - yn_t(1:end-1); 
    dz_t = zn_t(1:end-1) - zn_t(2:end); 

    write_MAG3D_TMI([tile_dir dsep 'Tile_data.dat'],H, HI, HD, MI, MD,obsx,obsy,obsz,d,wd)
    write_UBC_mesh(out_dir,['Tile' num2str(ii) '.msh'],dsep,xn_t(1),yn_t(1),zn_t(1),dx_t,dy_t,dz_t)
    write_UBC_mesh(tile_dir,['Tile' num2str(ii) '.msh'],dsep,xn_t(1),yn_t(1),zn_t(1),dx_t,dy_t,dz_t)
    
    mcell_t = length(xc_t) * length(yc_t) * length(zc_t);
    
    [Zc_t,Xc_t,Yc_t] = ndgrid(zc_t,xc_t,yc_t);

    % Take a subset of topography over the mesh
    indx = Topo(:,1) >= xn_t(1)-100 & Topo(:,1) <= xn_t(end)+100 & Topo(:,2) >= yn_t(1)-100 & Topo(:,2) <= yn_t(end)+100;
    topo = Topo(indx,:);
    
    %% NEED TO FIX THE TILED INTERPOLATION   
    [Zn,Xn,Yn] = ndgrid(zn_t,xn_t,yn_t);
    [nullcell,tcell,~] = topocheck(xn_t,yn_t,zn_t,topo+1e-5);
    save([tile_dir dsep 'nullcell.dat'],'-ascii','nullcell');
    
%     nullcell = load([work_dir '\nullcell.dat']);
    
    mactv = sum(nullcell);
    
    mag_azmdip = [ones(mcell_t,1)*MD ones(mcell_t,1)*MI];
    m_vec_t = azmdip_2_xyz( mag_azmdip(:,1) , mag_azmdip(:,2) );
     
    M = [spdiags(H * m_vec_t(:,1),0,mactv,mactv);spdiags(H * m_vec_t(:,2),0,mactv,mactv);spdiags(H * m_vec_t(:,3),0,mactv,mactv)];

    w_t = ones(4*mcell_t,1);
    
    t = ones(mcell_t,1);


%% CALCULATE THE SENSITIVITIES

    G = MAGSEN_Tiled(tile_dir,dsep,xn_t,yn_t,zn_t,H, HI, HD, MI, MD, obsx, obsy, obsz,nullcell, 'NONE', 'Guvw', R, R0, M);

%% RUN MAGNETIC VECTOR INVERSION (MVI)
    % Run one round of cartesian formulation as a primer
    
    mref_t = zeros(mcell_t,1);
    mstart_t = ones(mcell_t,1)*1e-4;
    esus = ones(mcell_t,1); % Effective Susceptibility model as constraint
    bounds_PST = ones(3,2)*10;
    bounds_PST(:,1) = -10;
    
%     [M_xyz,beta_in] = MVI3D_PST_Tiled(tile_dir,out_dir,dsep,xn_t,yn_t,zn_t,H, HI,...
%         HD, obsx, obsy, obsz, G, d, wd, mstart_t, mref_t, esus,...
%         1,alphas,[],bounds_PST,kron([1 1 1],[2 2 2 2 1]),t,lp_tresh{1},...
%         lp_tresh{2} ,FLAG1,FLAG2,5 );
%     
%     
    
    % Convert from xyz to atp
    aa = sum(M_xyz.^2,2).^0.5;
    
    tt = zeros(mactv,1);
    pp = zeros(mactv,1);
    
    tt(aa>0) = asin(M_xyz(aa>0,3)./(aa(aa>0)));
    
    
    dydx = M_xyz(aa>0,2)./M_xyz(aa>0,1);
    pp(aa>0) = atan(dydx);
    pp(M_xyz(:,2)>0 & M_xyz(:,1)<0) = pp(M_xyz(:,2)>0 & M_xyz(:,1)<0) + pi;
    pp(M_xyz(:,2)<0 & M_xyz(:,1)<0) = pp(M_xyz(:,2)<0 & M_xyz(:,1)<0) - pi;

    amp = aa; amp(nullcell==0) = -1;

    save([tile_dir dsep 'Tile' num2str(ii) '_MVI_PST.fld'],'-ascii','M_xyz')
    save([tile_dir dsep 'Tile' num2str(ii) '_MVI_PST.amp'],'-ascii','amp')
    %% Second run the spherical
    mref_t = zeros(3*mcell_t,1);
    
    temp = randn(mactv,1);
    temp = temp/max(abs(temp));
    mstart_t = [aa;tt;pp];
    
    bounds(2,:) = [-pi/2 pi/2];
    bounds(3,:) = [-pi pi];
    
    esus = ones(mcell_t,1);
% 	beta_in = [];
    max_iterMVI = 45;
    [M] = MVI3D_ATP_Tiled_v6(tile_dir,...
        out_dir,dsep,ii,xn_t,yn_t,zn_t,H, HI, HD, obsx, obsy, obsz, G,...
        d,wd,mstart_t,mref_t,chi_target,alphas,beta_in,bounds,LP,t,...
        lp_tresh{1},lp_tresh{2} ,FLAG1,FLAG2,max_iterMVI,ROT);

    % Normalize magnetization vector
%     esus_MVI = sqrt( sum( M.^2 , 2 ) );
%     m_vec = spdiags( 1./ esus_MVI ,0,mcell_t,mcell_t) * M;

    % All zeros vector magnetization are put back to induced
%     indx = esus_MVI==0;

%     m_vec(indx,:) = kron(ones(sum(indx),1),mp');
    fprintf('END CMI ITERATION # %i\n',ii);
%     end
    
end

fclose all;
