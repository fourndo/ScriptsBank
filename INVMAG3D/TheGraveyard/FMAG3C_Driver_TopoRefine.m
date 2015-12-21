% Generate model and observations for 3D gravity
% Dominique Fournier 2013/01/23
close all
clear all

% addpath C:\Users\dominiquef\Dropbox\Master\INVMAG3D\

addpath ..\FUNC_LIB\;
% Project folders
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Documents\Projects\Research\Modelling\Topo_adjust';

obsfile = 'Obs_loc.dat';
meshfile = 'Mesh_5m.msh';
model_sus = 'Model_0p01_SI.sus';
% model_azm = 'Dual_azm.sus';
% model_dip = 'Dual_dip.sus';
topofile = 'Gaussian.topo';

% Set octree levels
OctLev = [1 2 3];

% Load mesh file and convert to vectors (UBC format)
[xn,yn,zn]=read_UBC_mesh([work_dir '\' meshfile]);
% nx = mesh(1,1); %size(X,1);    %number of cell in X
% ny = mesh(1,2); %size(X,2);    %number of cell in Y
% nz = mesh(1,3); %size(X,3);    %number of cell in Z
% 
% dx = mesh(3,1:nx);
% dy = mesh(4,1:ny);
% dz = mesh(5,1:nz);
% 
% x0 = mesh(2,1);
% y0 = mesh(2,2);
% z0 = mesh(2,3);
% 
% % Create 3D cell node array
% xn = [x0 x0 + cumsum(dx)]; 
% yn = [y0 y0 + cumsum(dy)]; 
% zn = [z0 (z0 - cumsum(dz))];

% 3D vertical nodal location
[Zn,Xn,Yn] = ndgrid(zn,xn,yn);
mcell = (length(zn)-1) * (length(xn)-1) * (length(yn)-1);
% Load topo
topo = read_topo(work_dir,topofile);

% Create nullcell
[nullcell,topocell,ztopo_n] = topocheck(Xn,Yn,Zn,topo);

% Get nodal discretization for octree levels
OctNodes = MAG3C_RTC_OctNodes(Xn,Yn,Zn,topocell,ztopo_n,OctLev);

% Compute center of mass


% Find the index of the closest active cell for each topocell
% sprcell = MAG3C_RTC_link(dx,dy,dz,nullcell,topocell);

% test = ones(nx*ny*nz,1);
% for ii = 1:length(sprcell)
% test(sprcell(ii)) = 2;
% 
% end
% test(nullcell==0) = -1;
% 
% save([work_dir '\supercell.dat'],'-ascii','test');
% 
% save([work_dir '\nullcell.dat'],'-ascii','nullcell')
% save([work_dir '\topocell.dat'],'-ascii','topocell')
% load([work_dir '\nullcell.dat']);
% load([work_dir '\topocell.dat']);
% nullcell = ones(length(m),1);



% Load observation loc  
[H, I, Dazm, D, obsx, obsy, obsz, ~, ~] = read_MAG3D_obs([work_dir '\' obsfile]);
% plot_mag3C(obsx,obsy,oldB,I,D,'Observed 3C-Data data')
ndata = length(obsx);

% % Shift for test
% obsz = obsz - 4.5;

%% Get depth weighting
dwz = compdwz(obsx, obsy, obsz, xn,yn,zn, nullcell, 'distance');
% wr = (dwz).^(1/2);
% wr = wr./(max(wr));

save([work_dir '\dwz'],'dwz');

%% Create model magnetization vectors
m_azm_IND = ones(mcell,1)*Dazm;
m_dip_IND = ones(mcell,1)*I;
mv_IND = azmdip_2_xyz(m_azm_IND,m_dip_IND,mcell);

M = [spdiags(H * mv_IND(:,1),0,mcell,mcell);spdiags(H * mv_IND(:,2),0,mcell,mcell);spdiags(H * mv_IND(:,3),0,mcell,mcell)];


% m_azm_REM = load([work_dir '\' model_azm]);
% m_dip_REM = load([work_dir '\' model_dip]);
% mv_REM = azmdip_2_xyz(m_azm_REM,m_dip_REM,mcell);

% Create diagonal matrix of remenant mag
% M_REM =  [spdiags(H * mv_REM(:,1),0,mcell,mcell);spdiags(H * mv_REM(:,2),0,mcell,mcell);spdiags(H * mv_REM(:,3),0,mcell,mcell)];

% Pre-allocate to store fields
Tx = zeros(ndata,3*mcell);
Ty = zeros(ndata,3*mcell);
Tz = zeros(ndata,3*mcell);

progress = -1;
tic


for ii = 1:ndata

    dobsx = xn - obsx(ii);
    
    dobsy = yn - obsy(ii);
    
    dobsz = obsz(ii) - zn;
    
    % compute kernel for active cells
    [Tx(ii,:),Ty(ii,:),Tz(ii,:)] = MAG3C_T_row(dobsx,dobsy,dobsz);

    % refine kernel for topocell
%     [Tx(ii,:),Ty(ii,:),Tz(ii,:)] = MAG3C_T_octree(Tx(ii,:),Ty(ii,:),Tz(ii,:), Xn,Yn,Zn,nullcell,topocell,toponode,sprcell,obsx(ii),obsy(ii),obsz(ii));
    
    d_iter = floor(ii/ndata*100);
    if  d_iter > progress

        fprintf('Computed %i pct of data in %8.5f sec\n',d_iter,toc)
        progress = d_iter;

    end
            
end
% save([work_dir '\Tx'],'Tx');
% save([work_dir '\Ty'],'Ty');
% save([work_dir '\Tz'],'Tz');
% load([work_dir '\Tx']);
% load([work_dir '\Ty']);
% load([work_dir '\Tz']);

% Compute depth weighting


% Form full tensor matrix
% T = [Tx;Ty;Tz];

% Forward operator
% G       = T * M;


% Create projection matrix for TMI
P = [spdiags(ones(ndata,1)* (cosd(I) * cosd(D)),0,ndata,ndata) ...
    spdiags(ones(ndata,1)* (cosd(I) * sind(D)),0,ndata,ndata) ...
    spdiags(ones(ndata,1)* sind(I),0,ndata,ndata)];

% avg_sens = mean(P*G*spdiags(1./wr,0,mcell,mcell),1)';
% save([work_dir '\avg_sens.dat'],'-ascii','avg_sens');

%% Load model

% Load synthetic model
m = load([work_dir '\' model_sus]);
m(nullcell==0) = 0;

%% Compute data for 3C-TMI 

% data_3C    = G * m ;

bx = Tx*(M*m);%data_3C(1:ndata); 
by = Ty*(M*m);%data_3C(ndata+1:2*ndata); 
bz = Tz*(M*m);%data_3C(2*ndata+1:3*ndata); 


%% Write fwr file with noise
pct_noise   = 0.02;

% Write TMI data
data_TMI = P * [bx;by;bz];
floor       = (pct_noise.*max(abs(data_TMI)));
noise       = floor .*randn(ndata,1);
d_TMI = data_TMI + noise;
wd_TMI = ones(ndata,1); %./ floor;

% Write 3-components file
floor       = (pct_noise.*max(abs(bx)));
noise       = floor .*randn(ndata,1);
bx = bx + noise;
wdx = ones(ndata,1) ./ floor;

floor       = (pct_noise.*max(abs(by)));
noise       = floor .*randn(ndata,1);
by = by + noise;
wdy = ones(ndata,1) ./ floor;

floor       = (pct_noise.*max(abs(bz)));
noise       = floor .*randn(ndata,1);
bz = bz + noise;
wdz = ones(ndata,1) ./ floor;

d_3C = [bx;by;bz];
wd_3C = [wdx;wdy;wdz];


% write_MAG3D_3C(work_dir,'Synthetic_IND_3C_2pc_noise.obs',H,I,Dazm,...
%     obsx,obsy,obsz,bx,by,bz,wdx,wdy,wdz)

% plot_mag3C(obsx,obsy,d_3C,I,D,'Observed 3C-Data data')

write_MAG3D_TMI(work_dir,'Synthetic_IND_TMI_5m.obs',H,I,Dazm,...
    obsx,obsy,obsz,data_TMI,wd_TMI);

plot_TMI(obsx,obsy,data_TMI,d_TMI,wd_TMI,'Observed vs Predicted Magnitude');

