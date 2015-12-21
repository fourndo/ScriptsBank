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
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Nut_Cracker\Tiled_CMI';
out_dir = [work_dir dsep 'ALL_Tiles'];
mkdir(out_dir);

inpfile   = 'MAG3D_CMI_Tile.inp'; 

[meshfile,obsfile,topofile,mstart,mref,m_vec,chi_target,alphas,beta,bounds,lp_MAIvec,lp_MVIvec,lp_tresh,weightfile,FLAG1,FLAG2,tilefile] = AMI_Tile_read_inp([work_dir dsep inpfile]);

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

npadx = 7;
npady = 7;

R = 3;
R0 = min(dx)/4;
%% Load observation file (3C UBC-MAG format)
[H, I, Dazm, D, Obsx, Obsy, Obsz, data, wd_full] = read_MAG3D_obs([work_dir dsep obsfile]);
% plot_mag3C(obsx,obsy,d,I,D,'Observed 3C-data')
% plot_TMI(obsx,obsy,d,d,wd,'Observed vs Predicted Magnitude');

ndata = length(data);


%% Load input models and weights
% Load magnetization model

    
% mag_azmdip = [ones(mcell,1)*Dazm ones(mcell,1)*I];
% mag_xyz = azmdip_2_xyz( mag_azmdip(:,1) , mag_azmdip(:,2) );
    


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
    
    m_vec = azmdip_2_xyz( ones(mcell,1)*Dazm , ones(mcell,1)*I );
    
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
    
    if ntiles ~=1
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
        
    else
        
        xn_t = xn;
        yn_t = yn;
        zn_t = zn;
        
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

    write_MAG3D_TMI([tile_dir dsep 'Tile_data.dat'],H,I,Dazm,obsx,obsy,obsz,d,wd)
    write_UBC_mesh(out_dir,['Tile' num2str(ii) '.msh'],dsep,xn_t(1),yn_t(1),zn_t(1),dx_t,dy_t,dz_t)
    write_UBC_mesh(tile_dir,['Tile' num2str(ii) '.msh'],dsep,xn_t(1),yn_t(1),zn_t(1),dx_t,dy_t,dz_t)
    
    mcell_t = length(xc_t) * length(yc_t) * length(zc_t);
    
    [Zc_t,Xc_t,Yc_t] = ndgrid(zc_t,xc_t,yc_t);

    % Take a subset of topography over the mesh
    indx = Topo(:,1) >= xn_t(1)-100 & Topo(:,1) <= xn_t(end)+100 & Topo(:,2) >= yn_t(1)-100 & Topo(:,2) <= yn_t(end)+100;
    topo = Topo(indx,:);
    
    
    % Interpolate input models to mesh (cell-center)
%     indx = xc >= xc_t(1) & xc <= xc_t(end);
%     indy = yc >= yc_t(1) & yc <= yc_t(end);
%     indz = zc <= zc_t(1) & zc >= zc_t(end);
    
    %% NEED TO FIX THE TILED INTERPOLATION   
%     F = interp3(zc, xc, yc,reshape(mag_xyz(:,1),nz,nx,ny),zc_t,xc_t,yc_t);
%     mag_xyz_t = F(Zc_t,Xc_t,Yc_t);
%     
%     F = scatteredInterpolant(Zc(indz,indx,indy), Xc(indz,indx,indy),Yc(indz,indx,indy),mstart);
%     mstart_t = F(Zc_t,Xc_t,Yc_t);
%     
%     F = scatteredInterpolant(Zc(indz,indx,indy), Xc(indz,indx,indy),Yc(indz,indx,indy),reshape(mref,nz,nx,ny));
%     mref_t = F(Zc_t,Xc_t,Yc_t);
%     
%     F = scatteredInterpolant(Zc(indz,indx,indy), Xc(indz,indx,indy),Yc(indz,indx,indy),w);
%     w_t = F(Zc_t,Xc_t,Yc_t);
    [Zn,Xn,Yn] = ndgrid(zn_t,xn_t,yn_t);
    [nullcell,tcell,~] = topocheck(xn_t,yn_t,zn_t,topo+1e-5);
    save([tile_dir dsep 'nullcell.dat'],'-ascii','nullcell');
    mag_azmdip = [ones(mcell_t,1)*Dazm ones(mcell_t,1)*I];
    m_vec_t = azmdip_2_xyz( mag_azmdip(:,1) , mag_azmdip(:,2) );
    
    w_t = ones(4*mcell_t,1);
    
    mref_t = ones(mcell_t,1)*0;
    
    mstart_t = ones(mcell_t,1)*1e-4;
    
    m_esus_t = ones(mcell_t,1)*1e-4;
    
    t = ones(mcell_t,1);
%     wr = get_wrAmp(obsx, obsy, obsz, d, D, I, xn, yn, zn, nullcell, 'DISTANCE', 3.0,0);
%     s = ones(nz,nx,ny);
%     s(1:nz-npadx,npadx:nx-npadx+1,npadx:ny-npadx+1)= 0;
%     
%     t = [s(:)==0 s(:)~=0]*1;
%     
%     lp_MAIvec = [lp_MAIvec;0 2 2 2 1];
    %% CALL MAGSEN_Tile and create sensitivity
    [Tx,Ty,Tz] = MAGSEN_Tiled(tile_dir,dsep,xn_t,yn_t,zn_t,I, D, obsx, obsy, obsz,nullcell, 'DISTANCE', 'TxTyTz', R, R0);

    %% Remove regional field
%     d = RegRem_Tiled(tile_dir,dsep,ii,xn_t,yn_t,zn_t,H, I, D,Dazm, obsx, obsy,obsz, Tx, Ty, Tz, d, wd,mstart_t,mref_t,w_t,1,alphas,beta,[0 100],FLAG1,FLAG2,[npadx npady]);
   
    %% RUN EQUIVALENT SOURCE
    [d_amp,~] = EMS3D_Tiled(tile_dir,dsep,ii,xn_t,yn_t,zn_t,H, I, Dazm, D, obsx, obsy, obsz, d, wd,nullcell,mstart_t,mref_t,1,alphas,beta,1,FLAG1);

%% RUN AMPLITUDE INVERSION (MAI)
    betaMAI = [];
    betaMVI = [];
    
    % Pre-initiate phid
    phidMAI = ndata*2;
    phidMVI = ndata*2;
    MOF_MAI = speye(size(Tx,2)/3);
    MOF_MVI = speye(size(Tx,2));
    esus_MAI = m_esus_t;
    switchMAI = 0;
    switchMVI = 0;
    count = 0;
    beta_tol = 0.25;
    
    wr = get_wrAmp(obsx, obsy, obsz, d_amp/2, Dazm, I, xn, yn, zn, nullcell, 'DISTANCE', 3.0,0);
    wr = wr(:);
    save([tile_dir dsep 'wr_amp.dat'],'-ascii','wr');
    
    count = count + 1;
    fprintf('TILE %i: BEGIN AMI ITERATION # %i\n',ii,count);

    if switchMAI == 1

        max_iterMAI = 30;
        fprintf('BEGIN IRLS FOR MAI\n') 
    else

        max_iterMAI = 30;

    end

    if switchMAI~=3
        
        fprintf('\n### MAI STEP\n')
        [esus_MAI,pred_ampB,betaMAI,phidMAI,MOF_MAI,switchMAI] = MAI3D_Tiled(tile_dir,out_dir,dsep,ii,xn_t,yn_t,zn_t,H, I, Dazm, obsx, obsy,obsz, Tx, Ty, Tz,d_amp,wd/2,m_esus_t,mref_t,m_vec_t,w_t,chi_target,phidMAI,MOF_MAI,alphas,betaMAI,[0 100],lp_MAIvec,t,lp_tresh{1},lp_tresh{2},switchMAI,FLAG1,FLAG2,max_iterMAI,beta_tol );
    
    end

    %% PASS INFORMATION FROM MAI TO MVI
    % MAI assumes all in the direction of magnetization.
%     [mp,ms,mt] = azmdip_2_pst(Dazm,I,1);
%     Esus = spdiags(m_esus_t,0,mcell_t,mcell_t);
% 
%     % PST formulation
%     Mp = Esus * (mp'*m')' ;
%     Ms = Esus * (ms'*m')' ;
%     Mt = Esus * (mt'*m')' ;
% % 
% %     mstart = [Mp;Ms;Mt];
% 
%     % Sperical formulation
% %     Mp = ones(mcell_t,1) * 1e-4 ;
% %     Ms = ones(mcell_t,1) * D * pi /180 ;
% %     Mt = ones(mcell_t,1) * I * pi /180 ;
% 
%     mstart = [Mp;Ms;Mt];
    % mref = [Mp;Ms;Mt];

    %% RUN MAGNETIC VECTOR INVERSION (MVI)

    if switchMVI == 1 && switchMAI == 3

        max_iterMVI = 30;
        fprintf('BEGIN IRLS FOR MVI\n')

    else

        max_iterMVI = 30;
        fprintf('\nMVI STEP\n')
    end

    [M,pred_TMI,betaMVI,phidMVI,MOF_MVI,switchMVI] = MVI3D_Tiled(tile_dir,out_dir,dsep,ii,xn_t,yn_t,zn_t,H, I, Dazm, D, obsx, obsy, obsz, Tx, Ty, Tz, d, wd,mstart_t,mref_t*0,esus_MAI,chi_target,phidMVI,MOF_MVI,alphas,betaMVI,bounds,LP,t,lp_tresh{1},lp_tresh{2} ,switchMVI,FLAG1,FLAG2,max_iterMVI );

    % Normalize magnetization vector
    esus_MVI = sqrt( sum( M.^2 , 2 ) );
    m_vec = spdiags( 1./ esus_MVI ,0,mcell_t,mcell_t) * M;

    % All zeros vector magnetization are put back to induced
%     indx = esus_MVI==0;

%     m_vec(indx,:) = kron(ones(sum(indx),1),mp');
    fprintf('END AMI ITERATION # %i\n',ii);
%     end
    
end

% matlabpool CLOSE
%% Merging results back to core grid
% TileMerge(out_dir,[meshfile],500,-100,dsep,'topofile',['..' dsep topofile],'rem_pad',7)
% TileMerge_AMI(out_dir,['..' dsep meshfile],['..' dsep tilefile],50,-100,dsep,'none',['..' dsep topofile],'rem_pad',[10 10])
