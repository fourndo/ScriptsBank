% Attempt at modeling the DC expreiment in 2D

clear all
close all


addpath functions
addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\arrow3

work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\MtIsa\Modeling\DCIP2D';

dsep = '\';

meshfile = 'mesh2d_fine.txt';
confile = 'model2d_fine.con';

dvec = 5;
%% User Input - Build mesh and set model parameters

[xn,zn] = read_UBC_mesh_2D([work_dir dsep meshfile]);

zn = zn(15:end);

dx = xn(2:end) - xn(1:end-1); dx = dx';
dz = zn(2:end) - zn(1:end-1); dz = dz';

xc = ( xn(2:end) + xn(1:end-1) ) / 2;
zc = ( zn(2:end) + zn(1:end-1) ) / 2;

nx=length(dx) ; nz=length(dz);  % Number of cells

model = importdata([work_dir dsep confile],' ',1);

m = model.data;

m2d = reshape(m,nx,nz+14);
m2d = m2d(:,15:end);
m = m2d(:);

% m2d = m2d.^0*1e-3;
% m2d(48:58,10:20) = 1e-1;
% m = m2d(:);
[Xn,Zn] = ndgrid(xn,-zn); 
[XX,ZZ] = ndgrid(xc,-zc);
xx = XX(1:dvec:end,dvec:dvec:end); xx = xx(:);
zz = ZZ(1:dvec:end,dvec:dvec:end); zz = zz(:);
% figure;surface(XX,ZZ,ZZ*0,m2d,'EdgeColor','none','FaceAlpha',0.5); hold on
% view([-180 90]);
% axis equal

% set(gca,'XDir','normal')
% set(gca,'ZDir','normal')

uo = 4*pi*10^-7;            % Magnetic permeability of free-space
s_air = 1e-15;              % Air conductivity  
s_background = 1e-0;        % Background conductivity            
s_anomaly = 1e+3;           % Anomalous conductivity
w = 10^3.0;                 % Frequency [hz]

mcell = nx * nz;
nfx = (nx+1) * nz ;
nfz = (nz+1) * nx ;

skindepth = sqrt( 2/ (s_background * uo * w));



%% Averaging operator cell face to center

av = @(n) spdiags (ones(n,1)*[0.5 0.5],[-1,0],n+1,n);

avx = av(nx+1); avx = avx(2:end-1,:);
avz = av(nz+1); avz = avz(2:end-1,:);

Avfcx = kron ( speye (nz) ,  avx);
Avfcz = kron ( avz , speye (nx));

% AVGfc = [Avx spalloc(mcell,nzf,0);spalloc(mcell,nxf,0) Avz];
%% Averaging operator cell center to face
avx = av(nx); avx(1) = 1; avx(end) = 1;
avz = av(nz); avz(1) = 1; avz(end) = 1;

Avcfx = kron ( speye (nz) ,  avx);
Avcfz = kron ( avz , speye (nx));

M = [Avcfx;Avcfz];
%% Create dimension matrix
dxm = dx(1:end-1)/2 + dx(2:end)/2; dxm=[dxm(1);dxm;dxm(end)];
dzm = dz(1:end-1)/2 + dz(2:end)/2; dzm=[dzm(1);dzm;dzm(end)];

% % Create hmid dimensions matrix
% dXm = spdiags(1./dxm,0,nx+2,nx+2);
% dZm = spdiags(1./dzm,0,nz+2,nz+2);

% Create cell dimension matrix
dX = spdiags(1./dx,0,nx,nx); 
dZ = spdiags(1./dz,0,nz,nz); 

% Create over-sized cell dimension matrix with first and last dimension
% repeated for BCs.
dXxl=spdiags([1./dx(1);1./dx;1./dx(end)],0,nx+1,nx+1);
dZxl=spdiags([1./dz(1);1./dz;1./dz(end)],0,nz+1,nz+1);



%% Create primary gradient (center to face)

ddx = @(n) spdiags (ones (n+1,1)*[-1,1],[0,1],n+1,n+2);

d_dx = dXxl * ddx(nx); d_dx=d_dx(:,2:end-1); 
d_dx([1,1])=0; d_dx([end,end])=0;
Dx = kron(speye(nz),d_dx);

d_dz = dZxl * ddx(nz); d_dz=d_dz(:,2:end-1);
d_dz([1,1])=0;
Dz = kron(d_dz,speye(nx));

GRAD = [Dx; Dz];

%% Create primary divergence (face to center)
% First and last column of Dz and Dy are set to 0, since A = 0 on BC
d_dx=dX * ddx(nx-1); %d_dx(end,end)=0;
% d_dx([1,1])=0; d_dx([end,end])=0;
Dx = kron(speye(nz),d_dx);

d_dz=dZ * ddx(nz-1); %d_dz([end,end])=0;
Dz = kron(d_dz,speye(nx));

DIV = [Dx Dz];

%% Laplacian
LAP = DIV*GRAD;

%%
% aviobj = VideoWriter([work_dir '\Efield'],'MPEG-4');
% open(aviobj)

set(figure(1), 'Position', [25 50 900 700])

xloc = [0:5:275];
for jj = 1 : length(xloc)
%% Create boundary condition
s = zeros(nx,nz);
s(116+xloc(jj) ,1) = 10;
% s(18,1) = -10;

%% Solve
phi = (DIV*spdiags((M*(m.^-1)).^-1,0,nfx+nfz,nfx+nfz)*GRAD)\s(:);



% Compute currents in cell center
J = spdiags((M*(m.^-1)).^-1,0,nfx+nfz,nfx+nfz)*GRAD*phi;

Jx = J(1:nfx); Jx = Avfcx*Jx; 
Jz = J(nfx+1:end); Jz = Avfcz*Jz; 

lJl = sqrt(Jx.^2 + Jz.^2);
ljl = spdiags(1./lJl,0,mcell,mcell);

lJl = log10(lJl);

jx = ljl*Jx;
jz = ljl*Jz;

jx = reshape(jx,nx,nz); jx = jx(1:dvec:end,dvec:dvec:end); jx = jx(:);
jz = reshape(jz,nx,nz); jz = jz(1:dvec:end,dvec:dvec:end); jz = jz(:);
lJl = reshape(lJl,nx,nz); lJl = lJl(1:dvec:end,dvec:dvec:end); lJl = lJl(:);

set(gcf,'color','w');
ax1 = axes('Position',[0.075 0.325 .80 .80]);

caxis([min(lJl) max(lJl)]);
colormap(jet);
map=get(gcf,'colormap');
domain=0:1/(size(map,1)-1):1;

domain = (domain*(max(lJl) - min(lJl))*1.1)+min(lJl)*1.05;

cc=interp1(domain,map,lJl);
set(gca,'ColorOrder',cc)

    
set(gca,'XDir','normal')
set(gca,'ZDir','normal')
arrow3([xx zz zz*0],...
    [xx + jx zz - jz jz*0],'0o',500); hold on

axis equal
axis([400 1800 -600 -10]);
% Plot potential
ss = surface(XX-5,ZZ+5,ZZ*0,log10(m2d),'EdgeColor','none');
hold on
colormap(ax1,gray)
caxis([min(log10(m)) 1]);
colorbar
% quiver(XX,ZZ,Jx,Jz)
scatter3(xn(116+xloc(jj))+5,-30,5,300,'vk','filled')
set(gca,'XTickLabel',[])
ylabel('Z (m)')
title('Current density','interpreter','latex')
text(1850,0,'$log10(\sigma) \; (S/m)$','interpreter','latex','HorizontalAlignment','center')
%%
ax2 = axes('Position',[0.075 -.175 .80 .80]);

rho = (LAP*phi);

surface(XX,ZZ,ZZ*0,reshape(rho,nx,nz),'EdgeColor','none'); hold on
set(gca,'XDir','normal')
set(gca,'ZDir','normal')

colormap(ax2,jet)
caxis([-2 2])
colorbar
axis equal
axis([400 1800 -600 -10]);
scatter3(xn(116+xloc(jj))+5,-30,0,300,'vk','filled')
xlabel('East')
ylabel('Z (m)')
title('Charge density ($\nabla^2\;\phi$)','interpreter','latex')
text(1875,0,'$C/(V \cdot m)$','interpreter','latex','HorizontalAlignment','center')
%%
frame = getframe(figure(1));
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
if jj == 1;
  imwrite(imind,cm,[work_dir '\Efield.gif'],'gif', 'Loopcount',inf,'DelayTime',0.2);
else
  imwrite(imind,cm,[work_dir '\Efield.gif'],'gif','WriteMode','append','DelayTime',0.2);
end

% writeVideo(aviobj,frame);

clf(figure(1))
    
end

% close(aviobj);