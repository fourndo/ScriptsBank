% Attempt at modeling the DC expreiment in 2D

clear all
close all


addpath functions
addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\arrow3

work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\MtIsa\Modeling\DCIP2D';

dsep = '\';

meshfile = 'mesh2d_fine.txt';
confile = 'model2d_fine.con';

dvec = 10;
eps0 = 8.85418782e-12;
%% User Input - Build mesh and set model parameters

[xn,zn] = read_UBC_mesh_2D([work_dir dsep meshfile]);

dx = xn(2:end) - xn(1:end-1); dx = dx';
dz = zn(2:end) - zn(1:end-1); dz = dz';

xc = ( xn(2:end) + xn(1:end-1) ) / 2;
zc = ( zn(2:end) + zn(1:end-1) ) / 2;

nx=length(dx) ; nz=length(dz);  % Number of cells

model = importdata([work_dir dsep confile],' ',1);

m = model.data;

m2d = reshape(m,nx,nz);
% m2d = m2d.^0*1e-3;
% m2d(48:58,10:20) = 1e-1;
% m = m2d(:);

[Xn,Zn] = ndgrid(xn,-zn); 
[XX,ZZ] = ndgrid(xc,-zc);
xx = XX(1:dvec:end,1:dvec:end); xx = xx(:);
zz = ZZ(1:dvec:end,1:dvec:end); zz = zz(:);
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
nfz = (nx+1) * nz ;
nfx = (nz+1) * nx ;

skindepth = sqrt( 2/ (s_background * uo * w));

%% Averaging operator cell center to center
av = @(n) spdiags (ones(n,1)*[1 1 1]/3,[-1,0,1],n+1,n);
avx = av(nx); 
avx = avx(1:end-1,:);
avx(1,1:2) = 0.5; 
avx(end,end-1:end) = 0.5;

avz = av(nz); 
avz = avz(1:end-1,:);
avz(1,1:2) = 0.5; 
avz(end,end-1:end) = 0.5;

Avccx = kron ( speye (nz) ,  avx);
Avccz = kron ( avz , speye (nx));

% AVFC = [Avfcx spalloc(mcell,mcell,0);spalloc(mcell,mcell,0) Avfcz];


%% Averaging operator cell edge to center

av = @(n) spdiags (ones(n,1)*[0.5 0.5],[-1,0],n+1,n);

avx = av(nx+1); avx = avx(2:end-1,:);
avz = av(nz+1); avz = avz(2:end-1,:);

Avfcz = kron ( speye (nz) ,  avx);
Avfcx = kron ( avz , speye (nx));


%% Averaging operator cell center to edge
avx = av(nx); avx(1) = 1; avx(end) = 1;
avz = av(nz); avz(1) = 1; avz(end) = 1;

Avcfz = kron ( speye (nz) ,  avx);
Avcfx = kron ( avz , speye (nx));

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



%% Create primary gradient (node to edge)

ddx = @(n) spdiags (ones (n+1,1)*[-1,1],[0,1],n+1,n+2);

d_dx = dX * ddx(nx-1);% d_dx=d_dx(:,2:end-1); 
% d_dx([1,1])=0; d_dx([end,end])=0;
Dx = kron(speye(nz+1),d_dx);

d_dz = dZ * ddx(nz-1); %d_dz=d_dz(:,2:end-1);
% d_dz([1,1])=0; d_dz([end,end])=0;
Dz = kron(d_dz,speye(nx+1));

GRAD = [Dx; Dz];

%% Create primary divergence (edge to node)
% First and last column of Dz and Dy are set to 0, since A = 0 on BC
d_dx=dXxl * ddx(nx); d_dx=d_dx(:,2:end-1);
Dx = kron(speye(nz+1),d_dx);

d_dz=dZxl * ddx(nz); d_dz=d_dz(:,2:end-1);
Dz = kron(d_dz,speye(nx+1));

DIV = [Dx Dz];

%% Laplacian
LAP = DIV*GRAD;

%% Iteration over sources
set(figure(1), 'Position', [25 50 900 700])

xloc = [0:10:260];
ddx = min(dx);
ddz = min(dz);

for jj = 1 : length(xloc)
%% Create boundary condition
s = zeros(nx+1,nz+1);
s(102+xloc(jj) ,1) = 10;
s(410,1) = -10;

%% Solve
phi = (DIV*(spdiags((M*(m.^-1)).^-1,0,nfx+nfz,nfx+nfz)*GRAD))\s(:);



% Compute currents in cell center
E = GRAD*phi;

Jx = spdiags(m,0,mcell,mcell)*Avfcx*E(1:nfx); %Jx = Jx; 
Jz = spdiags(m,0,mcell,mcell)*Avfcz*E(nfx+1:end); %Jz = Jz; 

lJl = sqrt(Jx.^2 + Jz.^2);
ljl = spdiags(1./lJl,0,mcell,mcell);
jx = ljl*Jx;
jz = ljl*Jz;

lJl = sqrt(lJl);

jx = reshape(jx,nx,nz); jx = jx(1:dvec:end,1:dvec:end); jx = jx(:);
jz = reshape(jz,nx,nz); jz = jz(1:dvec:end,1:dvec:end); jz = jz(:);
lJl = reshape(lJl,nx,nz); lJl = lJl(1:dvec:end,1:dvec:end); lJl = lJl(:);

set(gcf,'color','w');
ax1 = axes('Position',[0.075 0.30 .85 .85]);

caxis([min(lJl) max(lJl)/100]);
colormap(jet);
map=get(gcf,'colormap');
domain=0:1/(size(map,1)-1):1;

domain = (domain*(max(lJl) + abs(min(lJl)))*1.1)-abs(min(lJl))*1.05;

cc=interp1(domain,map,lJl);
set(gca,'ColorOrder',cc)
   
set(gca,'XDir','normal')
set(gca,'ZDir','normal')
arrow3([xx-jx zz+jz zz*0],...
    [xx+jx  zz-jz  jz*0],'0o',600); hold on

axis equal
axis([400 1800 -600 50]);
% Plot potential
ss = surface(XX-5,ZZ+5,ZZ*0,log10(m2d),'EdgeColor','none','FaceAlpha',0.5);
hold on
colormap(ax1,gray)
caxis([min(log10(m)) max(log10(m))]);
colorbar
% quiver(XX,ZZ,Jx,Jz)

set(gca,'XTickLabel',[])
ylabel('Depth')
title('Current density','interpreter','latex')
text(1850,10,'$log10(\sigma) \; (S/m)$','interpreter','latex','HorizontalAlignment','center')
%%
ax2 = axes('Position',[0.075 -.175 .86 .86]);

rho = (LAP*phi);

surface(Xn,Zn,Zn*0,reshape(rho,nx+1,nz+1),'EdgeColor','none'); hold on
set(gca,'XDir','normal')
set(gca,'ZDir','normal')

colormap(ax2,jet)
caxis([-1 1])
colorbar
axis equal
axis([400 1800 -600 0]);

xlabel('East')
title('Charge density ($\nabla^2\;\phi$)','interpreter','latex')
text(1875,10,'$C/(V \cdot m)$','interpreter','latex','HorizontalAlignment','center')
%%
frame = getframe(figure(1));
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
if jj == 1;
  imwrite(imind,cm,[work_dir '\Efield.gif'],'gif', 'Loopcount',inf,'DelayTime',0.2);
else
  imwrite(imind,cm,[work_dir '\Efield.gif'],'gif','WriteMode','append','DelayTime',0.2);
end

clf(figure(1))
    
end