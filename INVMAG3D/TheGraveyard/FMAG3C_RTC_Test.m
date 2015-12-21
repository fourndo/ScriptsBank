% Generate model and observations for 3D gravity
% Dominique Fournier 2013/01/23
close all
clear all

% addpath C:\Users\dominiquef\Dropbox\Master\INVMAG3D\

addpath ..\FUNC_LIB\;
% Project folders
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Documents\Projects\Research\Modelling\Topo_adjust\Sphere_test';

obsfile = 'Obs_loc.dat';
meshfile = 'Unit_3by3.msh';
model_sus = 'Unit_3by3.sus';
% model_azm = 'Dual_azm.sus';
% model_dip = 'Dual_dip.sus';
% topofile = [work_dir '\Gaussian.topo'];

%% Load observation loc  
[H, I, Dazm, D, obsx, obsy, obsz, ~, ~] = read_MAG3D_obs([work_dir '\' obsfile]);
ndata = length(obsx);
obsx = [0:0.1:10]'+randn(ndata,1)*1e-4;
obsy = obsy+randn(ndata,1)*1e-4;
obsz = obsz +randn(ndata,1)*1e-4;

%% Load mesh file and convert to vectors (UBC format)
[xn,yn,zn]=read_UBC_mesh([work_dir '\' meshfile]);
dx = xn(2:end) - xn(1:end-1);

mcell = (length(xn)-1) * (length(yn)-1) * (length(zn)-1);

% 3D vertical nodal location
[Zn,Xn,Yn] = ndgrid(zn,xn,yn);

% Create nullcell
[nullcell,tcellID] = topocheck_unitsphere(Xn,Yn,Zn);

% Get index of active cells
acellID = find(nullcell==1);

%% Create model magnetization vectors
m_azm = ones(mcell,1)*Dazm;
m_dip = ones(mcell,1)*I;
mv = azmdip_2_xyz(m_azm,m_dip,mcell);

M = [spdiags(H * mv(:,1),0,mcell,mcell);spdiags(H * mv(:,2),0,mcell,mcell);spdiags(H * mv(:,3),0,mcell,mcell)];

%% Compute theoritical field from dipole of radius R at location cntr
% Center of sphere [x y z]
cntr_grid = [0 0 0];

% Radius of sphere
R = 1.5;

% Compute 3D and 2D distance from center
r = sqrt( obsx.^2 + obsy.^2 + obsz.^2 );
rxy = sqrt( obsx.^2 + obsy.^2);

phi = asin( rxy ./ r );

theta = 2*pi - acos( -obsx./ rxy );
theta(obsy(:)<0) = 2*pi - theta(obsy(:)<0);

rx = sin(phi) .* cos(theta); 
ry = sin(phi) .* sin(theta);
rz = cos(phi);

% Dipole moment
m = -1/(3) * R^3 *0.1* H  ./ (r.^3);

Mx = mv(1,1); My = mv(1,2) ; Mz = mv(1,3);

Bx = m .* (3* rx .*( Mx * rx + My * ry + Mz * rz ) - Mx);
By = m .* (3* ry .*( Mx * rx + My * ry + Mz * rz ) - My);
Bz = m .* (3* rz .*( Mx * rx + My * ry + Mz * rz ) - Mz);

B = sqrt( Bx.^2 + By.^2 + Bz.^2 );

% Define vector projection for TMI
% Compute TMI
TMI_true = (Bx)*Mx + (By)*My + (Bz)*Mz;

figure(1);semilogy(obsx/1.5,B);hold on; title('\bf|B| theoritical vs calculated')

%% Set octree levels
count = 1;
ccode{6}='k';ccode{5}='r';ccode{4}='m';ccode{3}='.-y';ccode{2}='.-g';ccode{1}='.-c';
% for oo = 0:5
OctLev = [1 2 3 4];

% Get nodal discretization for octree levels
acelln = MAG3C_RTC_OctNodes_unitsphere(Xn,Yn,Zn,acellID,0);
tcelln = MAG3C_RTC_OctNodes_unitsphere(Xn,Yn,Zn,tcellID,OctLev);

% MAG3C_RTC_write_octreecells(work_dir,'acell_cntr.dat',acelln);
% MAG3C_RTC_write_octreecells(work_dir,['tcell_cntr_' num2str(oo) '.dat'],tcelln);

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

% Find the octree level used for each obs unit 
% jlink = ones(length(obsx),size(t_cntm,1));
jlink = MAG3C_RTC_jlink(obsx, obsy, obsz, t_cntm, OctLev, 3*min(dx));

% %% Get depth weighting
% wr = MAG3C_RTC_get_dwr(obsx,obsy,obsz,xn,yn,zn,acellID,acelln,tcelln,alink);
% % wr = (dwz).^(1/2);
% % wr = wr./(max(wr));
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
    
    % compute kernel for active cells
%     [Tx(ii,:),Ty(ii,:),Tz(ii,:)] = MAG3C_RTC_T(obsx(ii),obsy(ii),obsz(ii),mcell,acelln,acellID,tcelln,alink,jlink(ii,:));
    [Tx(ii,:),Ty(ii,:),Tz(ii,:)] = MAG3C_T(obsx(ii),obsy(ii),obsz(ii),mcell,acelln,acellID);

    % refine kernel for topocell
    [tTx(ii,:),tTy(ii,:),tTz(ii,:)] = MAG3C_RTC_T_v2(obsx(ii),obsy(ii),obsz(ii),size(tcelln,1),tcelln,jlink(ii,:));
    
%     d_iter = floor(ii/ndata*20);
%     if  d_iter > progress
% 
%         fprintf('Computed %i pct of data in %8.5f sec\n',d_iter*5,toc)
%         progress = d_iter;
% 
%     end
            
end

% Create projection matrix for TMI
P = [spdiags(ones(ndata,1)* (cosd(I) * cosd(D)),0,ndata,ndata) ...
    spdiags(ones(ndata,1)* (cosd(I) * sind(D)),0,ndata,ndata) ...
    spdiags(ones(ndata,1)* sind(I),0,ndata,ndata)];

%% Load model

% Load synthetic model
m = load([work_dir '\' model_sus]);
m(nullcell==0) = 0;

%% Compute data for 3C-TMI 

% data_3C    = G * m ;

bx = Tx*(M*m);% + tTx * aPt * (M*m);%data_3C(1:ndata); 
by = Ty*(M*m);% + tTy * aPt * (M*m);%data_3C(ndata+1:2*ndata); 
bz = Tz*(M*m);% + tTz * aPt * (M*m);%data_3C(2*ndata+1:3*ndata); 

% TMI_RTC = P*[bx;by;bz];

b = sqrt(bx.^2 + by.^2 + bz.^2);

bx = bx + tTx * aPt * (M*m);%data_3C(1:ndata); 
by = by + tTy * aPt * (M*m);%data_3C(ndata+1:2*ndata); 
bz = bz + tTz * aPt * (M*m);%data_3C(2*ndata+1:3*ndata);

b_RTC = sqrt(bx.^2 + by.^2 + bz.^2);

figure(1);semilogy(obsx/1.5,b,'k');hold on
semilogy(obsx/1.5,b_RTC,'r');hold on

% Measure residual between true and modelled dipole
ratio = (B-b_RTC)./(B);
figure(2);plot(obsx/1.5,ratio,ccode{count});hold on

    
dratio = ratio./max(ratio);
figure(3);plot(obsx/1.5,dratio,ccode{count});hold on

count = count +1;

fprintf('%i cells \t %12.8f sec.\n',sum(sum(nnz(tcelln(:,end,:))))/6,toc);
% end
figure(1);
axis([0 6 0 100])
legend('Dipole','Octree 0','Octree 1','Octree 2','Octree 3','Octree 4','Octree 5')
xlabel('Length scale')
ylabel('\bf |B|')

figure(2);title('\bfNormalized residual r = ( |B|-|b| ) / |B|');hold on
legend('Octree 0','Octree 1','Octree 2','Octree 3','Octree 4','Octree 5')
xlabel('Length scale')
ylabel('\bf r')
xlim([0 6])

figure(3);title('\bfNormalized residual by its maximum r* = r / max(r)');hold on
legend('Octree 0','Octree 1','Octree 2','Octree 3','Octree 4','Octree 5')
xlabel('Length scale')
ylabel('\bf r / max(r)')
xlim([0 6])