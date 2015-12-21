% Generate model and observations for 3D gravity
% Dominique Fournier 2013/01/23
close all
clear all

% addpath C:\Users\dominiquef\Dropbox\Master\INVMAG3D\

addpath ..\FUNC_LIB\;
% Project folders
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Research\Modelling\Topo_adjust\Shell_Sphere';

obsfile = 'Obs_loc.dat';
meshfile = 'Mesh_5m.msh';
model_sus = 'Spheres.sus';
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
% obsx = obsx + randn(ndata,1)*1e-3;
% obsy = obsy + randn(ndata,1)*1e-3;
% obsz = obsz + randn(ndata,1)*1e-3;

%% Create model magnetization vectors
m_azm = ones(mcell,1)*Dazm;
m_dip = ones(mcell,1)*I;
mv = azmdip_2_xyz(m_azm,m_dip);

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
save([work_dir '\nullcell.dat'],'-ascii','nullcell');

% Get index of active cells
acellID = find(nullcell==1);

% for oo = 0:4
% Set octree levels
OctLev = [1 2 3 4];

% Get nodal discretization for octree levels
acelln = MAG3C_RTC_OctNodes(Xn,Yn,Zn,acellID,ztopo_n,0);
tcelln = MAG3C_RTC_OctNodes(Xn,Yn,Zn,tcellID,ztopo_n,OctLev);

% MAG3C_RTC_write_octreecells(work_dir,'acell_cntr.dat',acelln);
% MAG3C_RTC_write_octreecells(work_dir,['tcell_cntr' num2str(oo) '.dat'],tcelln);

% Compute center of mass for octree cells and active cells
t_cntm = MAG3C_RTC_cntmass(tcelln);

% Remove all topocell without cells below topo
keeper = isnan(t_cntm(:,1))==0;
tcelln = tcelln(keeper,:,:);
t_cntm = t_cntm(keeper,:);
tcellID = tcellID(keeper);

% Compute center mass of active cells
a_cntm = MAG3C_RTC_cntmass(acelln);

% Link topocells (tcell) to closest activecells (acell)
[nlink,alink] = MAG3C_RTC_linkID(a_cntm,t_cntm);

% Create projector matrix from active to topocell 
aPt = speye(mcell);
aPt = aPt(acellID(alink),:);
aPt = kron(speye(3),aPt);

save([work_dir '\aPt'],'aPt');

% Fix the octree level used for each obs
% jlink = ones(length(obsx),size(t_cntm,1));
jlink = MAG3C_RTC_jlink(obsx, obsy, obsz, t_cntm, OctLev, min(dx));

%% Get depth weighting
% wr = MAG3C_RTC_get_dwr(obsx,obsy,obsz,xn,yn,zn,acellID,acelln,tcelln,alink);
% wr = (dwz).^(1/2);
% wr = wr./(max(wr));
% save([work_dir '\wr.dat'],'-ascii','wr');

% Pre-allocate to store fields
Tx = zeros(ndata,3*mcell);
Ty = zeros(ndata,3*mcell);
Tz = zeros(ndata,3*mcell);

tTx = zeros(ndata,3*size(tcelln,1));
tTy = zeros(ndata,3*size(tcelln,1));
tTz = zeros(ndata,3*size(tcelln,1));

progress = -1;
tic

      
for ii = 1:ndata


%     dobsx = xns - obsx(ii);
%     
%     dobsy = yns - obsy(ii);
%     
%     dobsz = obsz(ii) - zns;
    
    % compute kernel for active cells
    [Tx(ii,:),Ty(ii,:),Tz(ii,:)] = MAG3C_T(obsx(ii),obsy(ii),obsz(ii),mcell,acelln,acellID);

    % refine kernel for topocell
    [tTx(ii,:),tTy(ii,:),tTz(ii,:)] = MAG3C_RTC_T(obsx(ii),obsy(ii),obsz(ii),size(tcelln,1),tcelln,jlink(ii,:));
    
    d_iter = floor(ii/ndata*20);
    if  d_iter > progress

        fprintf('Computed %i pct of data in %8.5f sec\n',d_iter*5,toc)
        progress = d_iter;

    end
            
end
Tx = sparse(Tx);
Ty = sparse(Ty);
Tz = sparse(Tz);

save([work_dir '\Tx'],'Tx');
save([work_dir '\Ty'],'Ty');
save([work_dir '\Tz'],'Tz');

save([work_dir '\tTx'],'tTx');
save([work_dir '\tTy'],'tTy');
save([work_dir '\tTz'],'tTz');

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
m = load([work_dir '\' model_sus]);
% m = nullcell*0.01;

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

% RTC correction
bx = bx + tTx*(aPt*(M*m));
by = by + tTy*(aPt*(M*m));
bz = bz + tTz*(aPt*(M*m));

data_RTC_TMI = P * [bx;by;bz];

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
write_MAG3D_TMI(work_dir,'Synthetic_RTC_5m.obs',H,I,Dazm,...
    obsx,obsy,obsz,data_RTC_TMI,wd_TMI);

write_MAG3D_TMI(work_dir,'Synthetic_5m.obs',H,I,Dazm,...
    obsx,obsy,obsz,data_TMI,wd_TMI);

% write_MAG3D_3C(work_dir,'Synthetic_IND_3C_2pc_noise.obs',H,I,Dazm,...
%     obsx,obsy,obsz,bx,by,bz,wdx,wdy,wdz)

% plot_mag3C(obsx,obsy,d_3C,I,D,'Observed 3C-Data data')


% 
plot_TMI(obsx,obsy,data_TMI,data_RTC_TMI,wd_TMI,'Observed vs Predicted Magnitude');

