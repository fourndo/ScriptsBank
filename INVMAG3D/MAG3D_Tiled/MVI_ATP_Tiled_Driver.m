% MAG3D_Tiles_Driver
%
% Magnetic Vector Inversion code. The large problem is divided in
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
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Nut_Cracker\Tiled_CMI';
work_dir ='C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\CraigModel\MAG';
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\SingleBlock\CMI';
% work_dir ='C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Triple_Block_lined';
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\MAG';
% work_dir = 'C:\LC\Private\dominiquef\Projects\4414_Minsim\Modeling\MAG\Composite\North_Main';
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Dipping_Prism_Li';
out_dir = [work_dir dsep 'ALL_Tiles'];
mkdir(out_dir);

inpfile   = 'MAG3D_ATP_Tile.inp'; 

% Read input file
[meshfile,obsfile,topofile,mstart,mref,m_vec,chi_target,alphas,beta,bounds,lp_MAIvec,lp_MVIvec,lp_tresh,weightfile,FLAG1,FLAG2,tilefile,ROT] = CMI_Tile_read_inp([work_dir dsep inpfile]);

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
R0 = min(dx)/4;

%% Load observation file (3C UBC-MAG format)
[H, HI, HD, MI, MD, Obsx, Obsy, Obsz, data, wd_full] = read_MAG3D_obs([work_dir dsep obsfile]);
% plot_mag3C(obsx,obsy,d,I,D,'Observed 3C-data')
% plot_TMI(obsx,obsy,d,d,wd,'Observed vs Predicted Magnitude');

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
    temp = sqrt( sum( m_vec.^2 , 2 ) );
    m_vec = spdiags( 1./ temp(:) ,0,mcell,mcell) * m_vec;
    
else
    
    m_vec = azmdip_2_xyz( ones(mcell,1)*MD , ones(mcell,1)*MI );
    
end

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

%% Generate s and t vectors for the different lp zones
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
    
    mactv = sum(nullcell);
    
    mag_azmdip = [ones(mcell_t,1)*MD ones(mcell_t,1)*MI];
    m_vec_t = azmdip_2_xyz( mag_azmdip(:,1) , mag_azmdip(:,2) );
     
    M = [spdiags(H * m_vec_t(:,1),0,mactv,mactv);spdiags(H * m_vec_t(:,2),0,mactv,mactv);spdiags(H * m_vec_t(:,3),0,mactv,mactv)];

    w_t = ones(4*mcell_t,1);
    
    
    
            
    t = ones(mcell_t,1);


%% RUN MAGNETIC VECTOR INVERSION (MVI)

    G = MAGSEN_Tiled(tile_dir,dsep,xn_t,yn_t,zn_t,H, HI, HD, MI, MD, obsx, obsy, obsz,nullcell, 'NONE', 'Guvw', R, R0, M);

    %% Run Inversion
    % First run the cartesian as primer
    
    mref_t = zeros(mcell_t,1);
    mstart_t = ones(mcell_t,1)*1e-4;
    esus_MAI = ones(mcell_t,1)*1e-4;
    bounds_PST = ones(3,2);
    bounds_PST(:,1) = -1;
    
    [m_xyz,~,~,~,~,~] = MVI3D_Tiled_v2(tile_dir,out_dir,dsep,ii,xn_t,yn_t,zn_t,H, HI, HD, obsx, obsy, obsz, G, d, wd,mstart_t,mref_t*0,esus_MAI,chi_target,[],[],alphas,[],bounds_PST,LP,t,lp_tresh{1},lp_tresh{2} ,0,FLAG1,FLAG2,1 );
    
    % Convert from xyz to atp
    aa = sum(m_xyz.^2,2).^0.5;
    tt = asin(m_xyz(:,3)./(aa));
    pp = acos((m_xyz(:,1)./(aa)./cos(tt)));
%     m_xyz2 = reshape(m_uvw(aa,tt,real(pp)),mactv,3);
%     save([work_dir dsep 'Tile_MVIrot.fld'],'-ascii','m_xyz2')
    
%     m_uvw = @(a,t,p) [a.*cos(t).*cos(p);...
%     a.*cos(t).*sin(p);...
%     a.*sin(t)];

    %% Second run the sperical
    mref_a = zeros(mcell_t,1);
    mref_t = ones(mcell_t,1)*(HI/180); 
    mref_p = ones(mcell_t,1)*(mod(450-HD,360)/180);

    mref_t = [mref_a;mref_t;mref_p];
    
    temp = randn(mactv,1);
    temp = temp/max(abs(temp));
    mstart_t = [aa;tt/pi;pp/pi];
    
    
    esus_MAI = ones(mcell_t,1)*1e-4;
	beta_in = [];
    max_iterMVI = 20;
    [M,pred_TMI,beta_out,phid_out,switcher] = MVI3D_ATP_Tiled_v2(tile_dir,...
        out_dir,dsep,ii,xn_t,yn_t,zn_t,H, HI, HD, obsx, obsy, obsz, G, d,...
        wd,mstart_t,mref_t,esus_MAI,chi_target,alphas,beta_in,bounds,LP,t,...
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
