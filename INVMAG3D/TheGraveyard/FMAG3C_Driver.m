% Generate model and observations for 3D gravity
% Dominique Fournier 2013/01/23
close all
clear all

% addpath C:\Users\dominiquef\Dropbox\Master\INVMAG3D\

addpath ..\FUNC_LIB\;
% Project folders
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Documents\Projects\Research\Modelling\Synthetic\Dual_Block';
obsfile     = 'Obs_loc.dat';
meshfile    = 'Mesh_40m.msh';
model_sus   = 'Dual_susc.sus';
% model_azm = 'Dual_azm.sus';
% model_dip = 'Dual_dip.sus';
topofile = [];

% Load mesh file and convert to vectors (UBC format)
[xn,yn,zn]=read_UBC_mesh([work_dir '\' meshfile]);

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
save([work_dir '\nullcell.dat'],'-ascii','nullcell');

% Get index of active cells
cellID = find(nullcell==1);


% Get nodal discretization for octree levels
celln = MAG3C_RTC_OctNodes(Xn,Yn,Zn,cellID,ztopo_n,0);

% Load synthetic model
m = load([work_dir '\' model_sus]);
mcell = length(m);

% Load observation loc  
[H, I, Dazm, D, obsx, obsy, obsz, ~, ~] = read_MAG3D_obs([work_dir '\' obsfile]);
% plot_mag3C(obsx,obsy,oldB,I,D,'Observed 3C-Data data')
ndata = length(obsx);



%% Create model magnetization vectors
m_azm = ones(mcell,1)*Dazm;
m_dip = ones(mcell,1)*I;
mv = azmdip_2_xyz(m_azm,m_dip,mcell);
M = [spdiags(H * mv(:,1),0,mcell,mcell);spdiags(H * mv(:,2),0,mcell,mcell);spdiags(H * mv(:,3),0,mcell,mcell)];

% Compute sensitivities
% Pre-allocate
Tx = zeros(ndata,3*mcell);
Ty = zeros(ndata,3*mcell);
Tz = zeros(ndata,3*mcell);
progress = -1;
tic   
for ii = 1:ndata
   
    % compute kernel for active cells
    [Tx(ii,:),Ty(ii,:),Tz(ii,:)] = MAG3C_T(obsx(ii),obsy(ii),obsz(ii),mcell,celln,cellID);

    d_iter = floor(ii/ndata*20);
    if  d_iter > progress

        fprintf('Computed %i pct of data in %8.5f sec\n',d_iter*5,toc)
        progress = d_iter;

    end
            
end
save([work_dir '\Tx'],'Tx');
save([work_dir '\Ty'],'Ty');
save([work_dir '\Tz'],'Tz');
% load([work_dir '\Tx']);
% load([work_dir '\Ty']);
% load([work_dir '\Tz']);

% Create projection matrix for TMI
P = [spdiags(ones(ndata,1)* (cosd(I) * cosd(D)),0,ndata,ndata) ...
    spdiags(ones(ndata,1)* (cosd(I) * sind(D)),0,ndata,ndata) ...
    spdiags(ones(ndata,1)* sind(I),0,ndata,ndata)];

% avg_sens = mean(P*G*spdiags(1./wr,0,mcell,mcell),1)';
% save([work_dir '\avg_sens.dat'],'-ascii','avg_sens');

%% Compute data for 3C-TMI 

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
wd_TMI = abs(d_TMI)*pct_noise + pct_noise*std(d_TMI); %./ floor;

% Write 3-components file
floor       = (pct_noise.*max(abs(bx)));
noise       = floor .*randn(ndata,1);
bx = bx + noise;
wdx = abs(bx)*pct_noise + pct_noise*std(bx);

floor       = (pct_noise.*max(abs(by)));
noise       = floor .*randn(ndata,1);
by = by + noise;
wdy = abs(by)*pct_noise + pct_noise*std(by);

floor       = (pct_noise.*max(abs(bz)));
noise       = floor .*randn(ndata,1);
bz = bz + noise;
wdz = abs(bz)*pct_noise + pct_noise*std(bz);

d_3C = [bx;by;bz];
wd_3C = [wdx;wdy;wdz];


write_MAG3D_3C(work_dir,'Synthetic_3C_2pc_noise.obs',H,I,Dazm,...
    obsx,obsy,obsz,bx,by,bz,wdx,wdy,wdz)

plot_mag3C(obsx,obsy,d_3C,I,D,'Observed 3C-Data data')

write_MAG3D_TMI(work_dir,'Synthetic_TMI_2pc_noise.obs',H,I,Dazm,...
    obsx,obsy,obsz,data_TMI,wd_TMI);

plot_TMI(obsx,obsy,data_TMI,d_TMI,wd_TMI,'Observed vs Predicted Magnitude');

