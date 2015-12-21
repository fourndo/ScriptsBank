% Generate model and observations for 3D gravity
% Dominique Fournier 2013/01/23
close all
clear all

% addpath C:\Users\dominiquef\Dropbox\Master\INVMAG3D\

addpath ..\FUNC_LIB\;

% Project folders
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Documents\Projects\Research\Modelling\Inversion\Synthetic\Sphere';


%% Compute theoritical field from dipole of radius R at location cntr
R = 5;
Dazm = -45;      % Declinaison of dipole from East
D = mod(450-Dazm,360);
I = 45;     % Inclinaison of dipole
H = 25000;  % Strength of inducing field (nT)
    
% Remenant mag
MDazm = 45;
MI = 45;

% Center of sphere [x y z]
cntr_grid = [0 0 0];
cntr(1,1:3) = [1.5*R 0 0];%
cntr(2,1:3) = [-1.5*R 0 0];%

% Observation points
obsx = [-30:5:-25 -20:2.5:20 25:5:30] + 1e-8;
obsy = [-20:5:-15 -10:2.5:10 15:5:20] + 1e-8;

ndata = length(obsx) * length(obsy);

[obsx,obsy] = ndgrid(obsx,obsy);

obsx = reshape(obsx,ndata,1);
obsy = reshape(obsy,ndata,1);
obsz = ones(ndata,1) * 10 + 1e-8;

% Compute offset between each data point and dipole
d_x{1} = obsx - cntr(1,1) ;
d_y{1} = obsy - cntr(1,2) ;
d_z{1} = obsz - cntr(1,3) ;

r{1} = sqrt( d_x{1}.^2 + d_y{1}.^2 + d_z{1}.^2 );
rxy{1} = sqrt( d_x{1}.^2 + d_y{1}.^2);

phi{1}  = asin( rxy{1} ./ r{1} );

theta{1} = 2*pi - acos( -d_x{1}./ rxy{1} );
theta{1}(d_y{1}(:)<0) = 2*pi - theta{1}(d_y{1}(:)<0);

rx{1} = sin(phi{1}) .* cos(theta{1}); 
ry{1} = sin(phi{1}) .* sin(theta{1});
rz{1} = cos(phi{1});

% Compute offset between each data point and dipole
d_x{2} = obsx - cntr(2,1) ;
d_y{2} = obsy - cntr(2,2) ;
d_z{2} = obsz - cntr(2,3) ;

r{2} = sqrt( d_x{2}.^2 + d_y{2}.^2 + d_z{2}.^2 );
rxy{2} = sqrt( d_x{2}.^2 + d_y{2}.^2);

phi{2}  = asin( rxy{2} ./ r{2} ) ;

theta{2} = 2*pi - acos( -d_x{2}./ rxy{2} ) ;
theta{2}(d_y{2}<0) = 2*pi - theta{2}(d_y{2}<0) ;

rx{2} = sin(phi{2}) .* cos(theta{2}); 
ry{2} = sin(phi{2}) .* sin(theta{2});
rz{2} = cos(phi{2});


B = ones(ndata,1);
wd = ones(ndata,1);


for ii = 1:length(MDazm)
for jj = 1:length(MI)
file_name = ['Forward_D_' num2str(abs(MDazm(ii))) '_I_' num2str(MI(jj)) '.dat'];
% write_MAG3D_3C(work_dir,file_name,H,I(jj),Dazm(ii),obsx,obsy,obsz,B,wd)

% [H, I, Dazm, D, obsx, obsy, obsz, Bx, By, Bz, wdx, wdy, wdz, TMI, magB] = read_MAG3D_3C([work_dir '\magfor3d.mag']);
% Adjust D:azimuth from North to Cartesian
MD(ii) = mod(450-MDazm(ii),360);

m{2} = 1/(3) * R^3 * H  ./ (r{2}.^3);

INDx = cosd(I) * cosd(D) ;
INDy = cosd(I) * sind(D) ;
INDz = sind(I) ;

bx{2} = m{2} .* (3* rx{2} .*( INDx * rx{2} + INDy * ry{2} + INDz * rz{2} ) - INDx);
by{2} = m{2} .* (3* ry{2} .*( INDx * rx{2} + INDy * ry{2} + INDz * rz{2} ) - INDy);
bz{2} = m{2} .* (3* rz{2} .*( INDx * rx{2} + INDy * ry{2} + INDz * rz{2} ) - INDz);

m{1} = 1/(3) * R^3 * H  ./ (r{1}.^3);

REMx = cosd(MI(jj)) * cosd(MD(ii)) ;
REMy = cosd(MI(jj)) * sind(MD(ii)) ;
REMz = sind(MI(jj)) ;

bx{1} = m{1} .* (3* rx{1} .*( REMx * rx{1} + REMy * ry{1} + REMz * rz{1} ) - REMx);
by{1} = m{1} .* (3* ry{1} .*( REMx * rx{1} + REMy * ry{1} + REMz * rz{1} ) - REMy);
bz{1} = m{1} .* (3* rz{1} .*( REMx * rx{1} + REMy * ry{1} + REMz * rz{1} ) - REMz);


% Define vector projection for TMI
% Compute TMI
TMI_true = (bx{1}+bx{2})*INDx + (by{1}+by{2})*INDy + (bz{1}+bz{2})*INDz;

% Compute magnitude of field
magB_true = sqrt( (bx{1}+bx{2}).^2 + (by{1}+by{2}).^2 + (bz{1}+bz{2}).^2 ) ;

% Plot data in 2D
plot_mag3C(obsx,obsy,[bx{1}+bx{2} ; (by{1}+by{2}) ; (bz{1}+bz{2})],I,D,'Theoritical field')

%% Discretize sphere with decreasing cell size and compare 
% predicted vs dipole field
dcell = 0.25;%5*0.5.^[1:6];




% Number of cells in each direction
nx = ceil(2*(R*3 + dcell) / dcell);
ny = ceil(2*(R + dcell) / dcell);

% Create an odd number of cells to get the dipole in the middle
if mod(nx,2)==0
    
    nx = nx+1;
    
end

if mod(ny,2)==0
    
    ny = ny+1;
    
end

% Create mesh
% nx = n; %size(X,1);    %number of cell in X
% ny = n; %size(X,2);    %number of cell in Y
nz = ny; %size(X,3);    %number of cell in Z

dx = ones(nx,1) * dcell;
dy = ones(ny,1) * dcell;
dz = ones(nz,1) * dcell;

mcell = nx*ny*nz;

x0 = cntr_grid(1) - sum(dx(1:floor(nx/2))) - dcell/2;
y0 = cntr_grid(2) - sum(dy(1:floor(ny/2))) - dcell/2;
z0 = cntr_grid(3) + sum(dx(1:floor(nz/2))) + dcell/2;

% Compute cell center location
xx = x0 + cumsum(dx) - dx/2;
yy = y0 + cumsum(dy) - dy/2;
zz = z0 - cumsum(dz) + dz/2;

% Create mesh
[ZZ,XX,YY] = ndgrid(zz,xx,yy);
XX = reshape(XX,mcell,1);
YY = reshape(YY,mcell,1);
ZZ = reshape(ZZ,mcell,1);

% Create model magnetization vectors
m_azm = ones(mcell,1) * Dazm(ii);
m_dip = ones(mcell,1) * I;

% Compute radial distance from dipole #1
rr = sqrt( ( XX-cntr(1,1) ).^2 + ( YY-cntr(1,2) ).^2 + ( ZZ-cntr(1,3) ).^2 );

% Create model
m = zeros(mcell,1);
m(rr<R) = 1;

% Create model magnetization vectors
m_azm(rr<R) = MDazm(ii);
m_dip(rr<R) = MI(jj);

% Create nullcell (all included for this test)
nullcell = ones(mcell,1);

% Compute radial distance from dipole #1
rr = sqrt( ( XX-cntr(2,1) ).^2 + ( YY-cntr(2,2) ).^2 + ( ZZ-cntr(2,3) ).^2 );

% Create model
m(rr<R) = 1;

M = azmdip_2_xyz(m_azm,m_dip,mcell);

save([work_dir '\M_test.dat'],'-ascii','M');

% Compute forward model
[ Bx, By, Bz, TMI, mB, obsx, obsy, obsz ] = FMAG3C( m, M, H, D, I, nullcell, obsx, obsy, obsz, x0, y0, z0, dx, dy, dz);

% wd = ones(ndata,1);
% file_name = ['Forward_D_' num2str(abs(Dazm(ii))) '_I_' num2str(I(jj)) '.dat'];
% write_MAG3D_3C(work_dir,file_name,H,I(jj),MD(ii),obsx,obsy,obsz,Bx,By,Bz,wd,wd,wd)

% Plot data in 2D
plot_mag3C(obsx,obsy,[Bx;By;Bz],I,D,['FMAG3C - D:' num2str(Dazm) ' I:' num2str(I(jj))])
hold on
subplot(2,2,4);plot([-12.5 -12.5 -2.5 -2.5 -12.5],[-5 5 5 -5 -5],'r','LineWidth',2);
hold on
subplot(2,2,4);plot([2.5 2.5 12.5 12.5 2.5],[-5 5 5 -5 -5],'r','LineWidth',2)

figure(1)
hold on
subplot(2,2,4);plot([-12.5 -12.5 -2.5 -2.5 -12.5],[-5 5 5 -5 -5],'r','LineWidth',2);
hold on
subplot(2,2,4);plot([2.5 2.5 12.5 12.5 2.5],[-5 5 5 -5 -5],'r','LineWidth',2)
% Plot residual between theoritical and computed
% plot(n(ii).^3,abs((TMI(ii)-TMI_true)/TMI_true),'*');hold on



end
end
%%
% fprintf('**Analytical Test**\nTensor mesh | | | | |\n')
% fprintf('ncells\t\tCell size\t\t|residual|\n')
% for ii  = 1:length(n)
%     
% fprintf('%i\t\t%3.2e\t\t%3.2e   \n',n(ii)^3,dcell(ii),abs((TMI(ii)-TMI_true)))
%     
% end   
wd = ones(length(TMI_true),1);
write_MAG3D_TMI(work_dir,'Spheres_TMI.obs',H,I,Dazm,obsx,obsy,obsz,TMI_true,wd);
write_MAG3D_3C(work_dir,'Spheres_3C.obs',H,I,Dazm,obsx,obsy,obsz,Bx,By,Bz,wd,wd,wd)
write_UBC_mesh(work_dir,x0,y0,z0,dx,dy,dz);

[H, I, Dazm, D, obsx, obsy, obsz, B, Wd] = read_MAG3D_3C([work_dir '\EquiSource\magfor3d.mag']);
plot_mag3C(obsx,obsy,B,I,D,['Equivalent source - D:' num2str(Dazm) ' I:' num2str(I(jj))]);
hold on
subplot(2,2,4);plot([-12.5 -12.5 -2.5 -2.5 -12.5],[-5 5 5 -5 -5],'r','LineWidth',2);
hold on
subplot(2,2,4);plot([2.5 2.5 12.5 12.5 2.5],[-5 5 5 -5 -5],'r','LineWidth',2)
% save([work_dir '\model.dat'],'-ascii','m');