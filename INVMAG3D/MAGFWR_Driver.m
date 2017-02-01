% Generate model and observations for 3D gravity
% Dominique Fournier 2013/01/23
close all
clear all

% addpath C:\Users\dominiquef\Dropbox\Master\INVMAG3D\

addpath ..\FUNC_LIB\;
% Project folders
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Shared\KrisDom\Osborne';

inpfile = 'MAG3Cfwr.inp';

[meshfile,obsfile,suscfile,magfile,topofile] = MAG3Cfwd_read_inp([work_dir '\' inpfile]);

% Load mesh file and convert to vectors (UBC format)
[xn,yn,zn]=read_UBC_mesh([work_dir '\' meshfile]);

nx = length(xn)-1;
ny = length(yn)-1;
nz = length(zn)-1;

% Load susceptibility model
m_sus = load([work_dir '\' suscfile]); %m_sus = m_sus;
mcell = length(m_sus);

flr_noise   = 0.0;
floor = 1.0;

%% 3D vertical nodal location
[Zn,Xn,Yn] = ndgrid(zn,xn,yn);

if isempty(topofile)==1
    
%     nxn = length(xn);
%     nyn = length(yn);
%     topo = [reshape(Xn(1,:,:),nxn*nyn,1) reshape(Yn(1,:,:),nxn*nyn,1) reshape(Zn(1,:,:),nxn*nyn,1)+1e-8];
    nullcell = m_sus~=-100;
else
    % Load topo
    topo = read_UBC_topo([work_dir '\' topofile]);
    [nullcell,tcell,~] = topocheck(xn,yn,zn,topo+1e-5);
    
end

% Create nullcell

%save([work_dir '\nullcell.dat'],'-ascii','nullcell');

%## TEMPORARY MOD##
% nullcell = m_sus~=-100;
% % Core out
% m_sus = reshape(m_sus,nz,nx,ny);
% m_sus(1:end-2,3:end-2,3:end-2) = 0;
% m_sus = m_sus(:);

% Get index of active cells
cellID = find(nullcell==1);


% Get nodal discretization for octree levels
%celln = MAG3C_RTC_OctNodes(Xn,Yn,Zn,cellID,ztopo_n,0);
% Create node pair for each cell
znzn = kron(kron(ones(ny,1),ones(nx,1)),[zn(1:end-1)' zn(2:end)']);
xnxn = kron(kron(ones(ny,1),[xn(1:end-1)' xn(2:end)']),ones(nz,1));
ynyn = kron(kron([yn(1:end-1)' yn(2:end)'],ones(nx,1)),ones(nz,1));

celln = [znzn(:,1) xnxn(:,1) ynyn(:,1) znzn(:,2) xnxn(:,2) ynyn(:,2)];
celln = celln(nullcell==1,:);



% Create selector matrix for active cells
X = spdiags(nullcell,0,mcell,mcell);
X = X(nullcell==1,:);

m_sus = X*m_sus;
% Load observation loc  
[H, HI, HD, MI, MD, dtype, obsx, obsy, obsz, d, ~] = read_MAG3D_obs([work_dir '\' obsfile]);

D = mod(450-HD,360);
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
    
    mag_model = [ones(mcell,1)*MD ones(mcell,1)*MI];
    
else
%     mag_model = [ones(mcell,1)*Dazm ones(mcell,1)*I];
    mag_model = load([work_dir '\' magfile]);
    
end

mag_model = X*mag_model;
mcell = sum(nullcell);

%% Create model magnetization vectors
% Azimuth and dip of magnitization
if size(mag_model,2)==2
    
    mag_xyz = azmdip_2_xyz( mag_model(:,1) , mag_model(:,2) );
    
    
    M = [spdiags(H * mag_xyz(:,1),0,mcell,mcell);spdiags(H * mag_xyz(:,2),0,mcell,mcell);spdiags(H * mag_xyz(:,3),0,mcell,mcell)];

else
    
    % Normalize vector to make sure unity
    mag_model = mag_model./kron([1 1 1],sqrt(sum(mag_model.^2,2))+1e-8);
    
    M = [spdiags(H * mag_model(:,1),0,mcell,mcell);spdiags(H * mag_model(:,2),0,mcell,mcell);spdiags(H * mag_model(:,3),0,mcell,mcell)];

   
end

% [p,s,t] = azmdip_2_pst(BD,BI,mcell);


p = [(cosd(HI) * cosd(D)) (cosd(HI) * sind(D)) sind(HI)];

% Compute sensitivities
% Pre-allocate
bx = zeros(ndata,1);
by = zeros(ndata,1);
bz = zeros(ndata,1);
btmi = zeros(ndata,1);

parfor_progress(ndata)
pooljob = parpool(3);
parfor ii = 1:ndata
   
    % compute kernel for active cells
    [Tx,Ty,Tz] = MAG3C_T(obsx(ii),obsy(ii),obsz(ii),celln);

    Gx = Tx * M;
    Gy = Ty * M;
    Gz = Tz * M;
    
    G = p * [Gx;Gy;Gz];
    
    bx(ii) = Gx * m_sus;%data_3C(1:ndata); 
    by(ii) = Gy * m_sus;%data_3C(ndata+1:2*ndata); 
    bz(ii) = Gz * m_sus;%data_3C(2*ndata+1:3*ndata); 
    
    btmi(ii) = G * m_sus;
    
%     parfor_progress;
            
end
delete(pooljob);
parfor_progress(0);
% save([work_dir '\Tx'],'Tx');
% save([work_dir '\Ty'],'Ty');
% save([work_dir '\Tz'],'Tz');
% load([work_dir '\Tx']);
% load([work_dir '\Ty']);
% load([work_dir '\Tz']);

% Create projection matrix for TMI
P = [spdiags(ones(ndata,1)* (cosd(HI) * cosd(D)),0,ndata,ndata) ...
    spdiags(ones(ndata,1)* (cosd(HI) * sind(D)),0,ndata,ndata) ...
    spdiags(ones(ndata,1)* sind(HI),0,ndata,ndata)];

% avg_sens = mean(P*G*spdiags(1./wr,0,mcell,mcell),1)';
% save([work_dir '\avg_sens.dat'],'-ascii','avg_sens');



%% Write fwr file with noise

% rand_noise = pct_noise*randn(ndata,1);


% Write TMI data
% data_TMI = P * [bx;by;bz];

% Amplitude data
lBl = sqrt( bx.^2 + by.^2 + bz.^2 );

% Write 3-components file
% floor       = (pct_noise.*max(abs(bx)));

bx = bx + flr_noise * randn(ndata,1);
wdx = floor*ones(ndata,1);%abs(bx)*pct_noise + pct_noise*std(bx);

% floor       = (pct_noise.*max(abs(by)));
% noise       = floor .*randn(ndata,1);
by = by + flr_noise * randn(ndata,1);
wdy = floor*ones(ndata,1);%abs(by)*pct_noise + pct_noise*std(by);

% floor       = (pct_noise.*max(abs(bz)));
% noise       = floor .*randn(ndata,1);
bz = bz + flr_noise * randn(ndata,1);
wdz = floor*ones(ndata,1);%abs(bz)*pct_noise + pct_noise*std(bz);

d_3C = [bx;by;bz];
wd_3C = [wdx;wdy;wdz];

% Write TMI data
d_TMI = P * [bx;by;bz];
% floor       = (pct_noise.*max(abs(data_TMI)));
% noise       = floor .*rand(ndata,1);
% d_TMI = data_TMI + noise;
wd_TMI = floor*ones(ndata,1);%abs(d_TMI)*pct_noise + pct_noise*std(d_TMI); %./ floor;


write_MAG3D_3C([work_dir '\' obsfile(1:end-4) '_3C.obs'],H,HI,HD,...
    obsx,obsy,obsz,bx,by,bz,wdx,wdy,wdz);

plot_mag3C(obsx,obsy,d_3C,HI,D,'Observed 3C-Data data');

write_MAG3D_TMI([work_dir '\' obsfile(1:end-4) '_TMI.obs'],H,HI,HD,MI,MD,...
    obsx,obsy,obsz,d_TMI,wd_TMI);

plot_TMI(obsx,obsy,btmi,d_TMI,wd_TMI,'Observed vs Predicted Magnitude',false);

% noise       = floor .*randn(ndata,1);
lBl = lBl + flr_noise * randn(ndata,1);
wd = floor*ones(ndata,1);
% lBl = sqrt( bx.^2 + by.^2 + bz.^2 );
% wd = sqrt( wdx.^2 + wdy.^2 + wdz.^2 );

write_MAG3D_TMI([work_dir '\' obsfile(1:end-4) '_lBl.obs'],H,HI,HD,MI,MD,obsx,obsy,obsz,lBl,wd)

