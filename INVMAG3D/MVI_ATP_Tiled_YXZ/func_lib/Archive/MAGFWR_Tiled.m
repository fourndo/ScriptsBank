function MAGFWR_Tiled(work_dir,meshfile,obsfile,suscfile,topofile,dsep)
% Forward model magnetic data from susceptibility model
% Dominique Fournier 2013/01/23


% Load mesh file and convert to vectors (UBC format)
[xn,yn,zn]=read_UBC_mesh([work_dir dsep meshfile]);

%% 3D vertical nodal location
[Zn,Xn,Yn] = ndgrid(zn,xn,yn);

if isempty(topofile)==1
    
    nxn = length(xn);
    nyn = length(yn);
    topo = [reshape(Xn(1,:,:),nxn*nyn,1) reshape(Yn(1,:,:),nxn*nyn,1) reshape(Zn(1,:,:),nxn*nyn,1)+1e-8];
    
else
    % Load topo
    topo = read_UBC_topo([work_dir dsep topofile]);
end

% Create nullcell
[nullcell,tcellID,ztopo_n] = topocheck(Xn,Yn,Zn,topo);
% save([work_dir '\nullcell.dat'],'-ascii','nullcell');

% Get index of active cells
cellID = find(nullcell==1);


% Get nodal discretization for octree levels
celln = MAG3C_RTC_OctNodes(Xn,Yn,Zn,cellID,ztopo_n,0);

% Load susceptibility model
m_sus = load([work_dir dsep suscfile]); %m_sus = m_sus;
m_sus(m_sus==-100) = 0;
mcell = length(m_sus);


% Create selector matrix for active cells
X = spdiags(nullcell,0,mcell,mcell);
X = X(nullcell==1,:);

m_sus = X*m_sus;
% Load observation loc  
[H, I, Dazm, D, obsx, obsy, obsz, ~, wd] = read_MAG3D_obs([work_dir dsep obsfile]);
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


mag_model = [ones(mcell,1)*Dazm ones(mcell,1)*I];
    
mag_model = X*mag_model;
mcell = sum(nullcell);

%% Create model magnetization vectors
% Azimuth and dip of magnitization
    
mag_xyz = azmdip_2_xyz( mag_model(:,1) , mag_model(:,2) );

M = [spdiags(H * mag_xyz(:,1),0,mcell,mcell);spdiags(H * mag_xyz(:,2),0,mcell,mcell);spdiags(H * mag_xyz(:,3),0,mcell,mcell)];

% [p,s,t] = azmdip_2_pst(Dazm,I,mcell);

% Create projection matrix for TMI
P = [cosd(I) * cosd(D) cosd(I) * sind(D) sind(I)];

% Compute sensitivities
% Pre-allocate
Tx = zeros(1,3*mcell);
Ty = zeros(1,3*mcell);
Tz = zeros(1,3*mcell);
data_TMI = zeros(ndata,1);


% count = 0;

% matlabpool OPEN 8
for ii = 1:ndata
   
    % compute kernel for active cells
    [Tx,Ty,Tz] = MAG3C_T(obsx(ii),obsy(ii),obsz(ii),celln);

    
    bx = Tx*(M*m_sus);%data_3C(1:ndata); 
    by = Ty*(M*m_sus);%data_3C(ndata+1:2*ndata); 
    bz = Tz*(M*m_sus);%data_3C(2*ndata+1:3*ndata); 

    % Write TMI data
    data_TMI(ii) = P * [bx;by;bz];

%     d_iter = floor(count/ndata*100);
%     if  d_iter > progress
% 
        fprintf('Computed data %i\n',ii)
%         progress = d_iter;
% 
%     end
         
end
% matlabpool CLOSE

%% Write fwr file with noise
% pct_noise   = 0.00;
% rand_noise = pct_noise*randn(ndata,1);
% floor = 1.0;
% 
% 
% 
% % Amplitude data
% lBl = sqrt( bx.^2 + by.^2 + bz.^2 );
% 
% % Write 3-components file
% % floor       = (pct_noise.*max(abs(bx)));
% noise       = floor .*randn(ndata,1);
% bx = bx + noise;
% wdx = floor*ones(ndata,1);%abs(bx)*pct_noise + pct_noise*std(bx);
% 
% % floor       = (pct_noise.*max(abs(by)));
% noise       = floor .*randn(ndata,1);
% by = by + noise;
% wdy = floor*ones(ndata,1);%abs(by)*pct_noise + pct_noise*std(by);
% 
% % floor       = (pct_noise.*max(abs(bz)));
% noise       = floor .*randn(ndata,1);
% bz = bz + noise;
% wdz = floor*ones(ndata,1);%abs(bz)*pct_noise + pct_noise*std(bz);
% 
% d_3C = [bx;by;bz];
% wd_3C = [wdx;wdy;wdz];
% 
% % Write TMI data
% d_TMI = P * [bx;by;bz];
% % floor       = (pct_noise.*max(abs(data_TMI)));
% % noise       = floor .*rand(ndata,1);
% % d_TMI = data_TMI + noise;
% wd_TMI = floor*ones(ndata,1);%abs(d_TMI)*pct_noise + pct_noise*std(d_TMI); %./ floor;


% write_MAG3D_3C([work_dir dsep obsfile(1:end-4) '_3C.obs'],H,I,Dazm,...
%     obsx,obsy,obsz,bx,by,bz,wdx,wdy,wdz);

% plot_mag3C(obsx,obsy,d_3C,I,D,'Observed 3C-Data data')

write_MAG3D_TMI([work_dir dsep 'PRED_' obsfile(1:end-4) '_TMI.pre'],H,I,Dazm,...
    obsx,obsy,obsz,data_TMI,wd);

% plot_TMI(obsx,obsy,data_TMI,d_TMI,wd_TMI,'Observed vs Predicted Magnitude');

% noise       = floor .*randn(ndata,1);
% lBl = lBl + noise;
% wd = floor*ones(ndata,1);
% lBl = sqrt( bx.^2 + by.^2 + bz.^2 );
% wd = sqrt( wdx.^2 + wdy.^2 + wdz.^2 );

% write_MAG3D_TMI([work_dir '\' obsfile(1:end-4) '_lBl.obs'],H,I,Dazm,obsx,obsy,obsz,lBl,wd)

