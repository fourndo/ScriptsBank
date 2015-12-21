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
model_sus = 'model.sus';
% model_azm = 'Dual_azm.sus';
% model_dip = 'Dual_dip.sus';
topofile = 'Gaussian.topo';

% Load mesh file and convert to vectors (UBC format)
[xn,yn,zn]=read_UBC_mesh([work_dir '\' meshfile]);
dx = xn(2:end) - xn(1:end-1);

mcell = (length(xn)-1) * (length(yn)-1) * (length(zn)-1);

%% Load observation loc  
[H, I, Dazm, D, obsx, obsy, obsz, ~, ~] = read_MAG3D_obs([work_dir '\' obsfile]);
ndata = length(obsx);
obsx = obsx + randn(ndata,1)*1e-3;
obsy = obsy + randn(ndata,1)*1e-3;
obsz = obsz + randn(ndata,1)*1e-3;

%% Create model magnetization vectors
m_azm = ones(mcell,1)*Dazm;
m_dip = ones(mcell,1)*I;
mv = azmdip_2_xyz(m_azm,m_dip,mcell);

% m_azm = load([work_dir '\' model_azm]);
% m_dip = load([work_dir '\' model_dip]);
% mv = azmdip_2_xyz(m_azm_REM,m_dip_REM,mcell);

M = [spdiags(H * mv(:,1),0,mcell,mcell);spdiags(H * mv(:,2),0,mcell,mcell);spdiags(H * mv(:,3),0,mcell,mcell)];

%% 3D vertical nodal location
[Zn,Xn,Yn] = ndgrid(zn,xn,yn);

% Load topo
topo = read_UBC_topo([work_dir '\' topofile]);

% Create nullcell
[nullcell,tcellID,ztopo_n] = topocheck(Xn,Yn,Zn,topo);

% Get index of active cells
acellID = find(nullcell==1);
ccode{6}='k';ccode{5}='r';ccode{4}='m';ccode{3}='.-y';ccode{2}='.-g';ccode{1}='.-c';

for oo = 0:4
% Set octree levels
OctLev = oo;
% Get nodal discretization for octree levels
acelln = MAG3C_RTC_OctNodes(Xn,Yn,Zn,acellID,ztopo_n,0);
tcelln = MAG3C_RTC_OctNodes(Xn,Yn,Zn,tcellID,ztopo_n,OctLev);

% MAG3C_RTC_write_octreecells(work_dir,'acell_cntr.dat',acelln);
MAG3C_RTC_write_octreecells(work_dir,['tcell_cntr' num2str(oo) '.dat'],tcelln);

% Compute center of mass for octree cells and active cells
t_cntm = MAG3C_RTC_cntmass(tcelln);
% Remove all topocell without cells below topo
tcelln = tcelln(isnan(t_cntm(:,1))==0,:,:);
t_cntm = t_cntm(isnan(t_cntm(:,1))==0,:);

% Compute center mass of active cells
a_cntm = MAG3C_RTC_cntmass(acelln);

% Link topocells (tcell) to closest activecells (acell)
[nlink,alink] = MAG3C_RTC_linkID(a_cntm,t_cntm);

% Fix the octree level used for each obs
jlink = ones(length(obsx),size(t_cntm,1));
% [jlink] = MAG3C_RTC_jlink(obsx, obsy, obsz, t_cntm, OctLev, min(dx));

%% Get depth weighting
% wr = MAG3C_RTC_get_dwr(obsx,obsy,obsz,xn,yn,zn,acellID,acelln,tcelln,alink);
% wr = (dwz).^(1/2);
% wr = wr./(max(wr));
% save([work_dir '\wr.dat'],'-ascii','wr');

% Pre-allocate to store fields
Tx = zeros(ndata,3*mcell);
Ty = zeros(ndata,3*mcell);
Tz = zeros(ndata,3*mcell);

progress = -1;
tic

 figure(1);  title('B-field - RTC')      
for ii = 1:ndata


%     dobsx = xns - obsx(ii);
%     
%     dobsy = yns - obsy(ii);
%     
%     dobsz = obsz(ii) - zns;
    
    % compute kernel for active cells
    [Tx(ii,:),Ty(ii,:),Tz(ii,:)] = MAG3C_RTC_T(obsx(ii),obsy(ii),obsz(ii),mcell,acelln,acellID,tcelln,alink,jlink(ii,:));

    % refine kernel for topocell
%     [Tx(ii,:),Ty(ii,:),Tz(ii,:)] = MAG3C_T_octree(Tx(ii,:),Ty(ii,:),Tz(ii,:), Xn,Yn,Zn,nullcell,topocell,toponode,sprcell,obsx(ii),obsy(ii),obsz(ii));
    
    d_iter = floor(ii/ndata*20);
    if  d_iter > progress

        fprintf('Computed %i pct of data in %8.5f sec\n',d_iter*5,toc)
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
% G       = [Tx;Ty;Tz] * M;


% Create projection matrix for TMI
P = [spdiags(ones(ndata,1)* (cosd(I) * cosd(D)),0,ndata,ndata) ...
    spdiags(ones(ndata,1)* (cosd(I) * sind(D)),0,ndata,ndata) ...
    spdiags(ones(ndata,1)* sind(I),0,ndata,ndata)];

% avg_sens = mean(P*G*spdiags(1./wr,0,mcell,mcell),1)';
% save([work_dir '\avg_sens.dat'],'-ascii','avg_sens');

%% Load model

% Load synthetic model
% m = load([work_dir '\' model_sus]);
m = nullcell*0.01;

%% Compute data for 3C-TMI 

% data_3C    = G * m ;

bx = Tx*(M*m);%data_3C(1:ndata); 
by = Ty*(M*m);%data_3C(ndata+1:2*ndata); 
bz = Tz*(M*m);%data_3C(2*ndata+1:3*ndata); 


%% Write fwr file with noise
pct_noise   = 0.00;

% Write TMI data
data_TMI = P * [bx;by;bz];
uncert       = (pct_noise.*max(abs(data_TMI)));
noise       = uncert .*randn(ndata,1);
d_TMI = data_TMI + noise;
wd_TMI = ones(ndata,1); %./ floor;

% Write 3-components file
uncert       = (pct_noise.*max(abs(bx)));
noise       = uncert .*randn(ndata,1);
bx = bx + noise;
wdx = ones(ndata,1) ./ uncert;

uncert       = (pct_noise.*max(abs(by)));
noise       = uncert .*randn(ndata,1);
by = by + noise;
wdy = ones(ndata,1) ./ uncert;

uncert       = (pct_noise.*max(abs(bz)));
noise       = uncert .*randn(ndata,1);
bz = bz + noise;
wdz = ones(ndata,1) ./ uncert;

d_3C = [bx;by;bz];
wd_3C = [wdx;wdy;wdz];

% semilogy(obsx,sqrt(bx.^2 + by.^2 +bz.^2),ccode{oo+1}); hold on
write_MAG3D_TMI(work_dir,['Synthetic_Octree' num2str(oo) '_5m.obs'],H,I,Dazm,...
    obsx,obsy,obsz,data_TMI,wd_TMI);
end

% write_MAG3D_3C(work_dir,'Synthetic_IND_3C_2pc_noise.obs',H,I,Dazm,...
%     obsx,obsy,obsz,bx,by,bz,wdx,wdy,wdz)

% plot_mag3C(obsx,obsy,d_3C,I,D,'Observed 3C-Data data')


% 
% plot_TMI(obsx,obsy,data_TMI,d_TMI,wd_TMI,'Observed vs Predicted Magnitude');

