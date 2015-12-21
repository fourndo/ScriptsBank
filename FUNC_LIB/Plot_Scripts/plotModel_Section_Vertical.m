% function plotModel_Section(mmat,IDX,lbel)
% Inputs:
%   model - model Matrix
%   padcells - cluster result
clear all
% close all

% padE=57;
% padW=42;
% padN=5;
% padS=5;
% padT=18;
% padB=8;
% 
% Zsection = 52;
% 
% DO27section1 = 12;
% DO27section2 = 12;
% DO18section = 20;
% 
% inp_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\DIGHEM\1D\DIGHEM_1D_inv_3Dmod';
% mod_file = 'Inv_iter9.con';
% meshfile = 'UBC_mesh_small.msh';


padE=12;
padW=12;
padN=30;
padS=4;
padT=0;
padB=4;
Zsection = 10;

DO27section1 = 20;
DO27section2 = 24;
DO18section = 52;

% inp_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\EOSC 556b (2013) EM\TKC_Project\Codes\CLUSTERING\Models';
inp_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\EOSC 556b (2013) EM\TKC_Project\Codes\CLUSTERING\Cluster_models';

mod_file{1} = 'Kmeans_DO27_Grav_Mag_DIGHEM_4Clusters.dat';
mod_file{2} = 'Kmeans_DO27_Grav_Mag_VTEM_4Clusters.dat';
mod_file{3} = 'Kmeans_DO27_Grav_Mag_DIGHEM_VTEM_4Clusters.dat';

meshfile = 'TKC_common_mesh.msh';

% Load point file for mapping
secfile = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\EOSC 556b (2013) EM\TKC_Project\Modelling\Geology_Objects\Section_trace.dat';
% secfile = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\EOSC 556b (2013) EM\TKC_Project\Modelling\Geology_Objects\Density_Estimate_Section.dat';
outline = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\EOSC 556b (2013) EM\TKC_Project\Modelling\Geology_Objects\DO27_Outline.dat';

ndv = -100;

% cpal = 'jet';%[255 255 255;0 71 189;99 209 62;255 204 0;255 0 0;255 153 200]/255;
% cvec = [-5 -4 -3 -2 -1.5 -1];
% crange = -5:1e-3:max(cvec);
% barlabel = '$log10(\sigma) (S/m)$';

cpal = [255 255 255;52 170 220;76 217 100;255 204 0;255 59 48]/255;
% cpal = 1:-0.20:0.2;
% cpal = kron([1 1 1],cpal');
cvec = [0 1 2 3 4];
crange = cvec;
barlabel = 'Rock Unit';

%% Load in section
DO27 = load([secfile]);
Out_grav = load([outline]);

%% Read mesh
[xn,yn,zn] = read_UBC_mesh([inp_dir filesep meshfile]);
dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);

xc = ( xn(1:end-1) + xn(2:end) ) / 2;
yc = ( yn(1:end-1) + yn(2:end) ) / 2;
zc = ( zn(1:end-1) + zn(2:end) ) / 2;

set(figure, 'Position', [50 25 850 550]);
ttle = ['a' 'b' 'c'];
for ii = 1 : 3
%% Load model
m = load([inp_dir '\' mod_file{ii}]);

% temp = m;
% temp(m==2) = 1;
% temp(m==1) = 4;
% temp(m==4) = 3;
% temp(m==3) = 2;
% m =temp;

nullcell = m~=ndv;
% Convert to log
% m = log10(m);

% Reshape to a 2D model
m = reshape(m,nz,nx,ny);

mmax = -100;

%% Plot horizontal section
ax1 = axes('Position',[0.01+(ii-1)*0.3 .5 0.4 0.4]);

m2D = squeeze(reshape(m(Zsection,padW+1:end-padE,padS+1:end-padN),nx-padW-padE,ny-padN-padS));

[XX,YY] = ndgrid(xc(padW+1:end-padE),yc(padS+1:end-padN));
ZZ = ones(size(XX)) * zc(Zsection);

h = surf(XX,YY,ZZ,m2D,'EdgeColor','none'); hold on
view([0 90])
set(gca,'YDir','normal')

if ii == 1
ylabel('$Northing~(m)$', 'interpreter', 'latex','FontSize',14)
else
    set(gca,'YTickLabel',[])
end
% xlabel('$Easting~(m)$', 'interpreter', 'latex','FontSize',14)
set(gca,'XTickLabel',[])
axis equal tight
caxis([min(cvec) max(cvec)])
% colorbar('SouthOutside')
title(['(' ttle(ii) ')']);

if ii == 2
text(mean(xc),yc(padS),'$Depth\;\approx\;100\;m$','interpreter', 'latex','VerticalAlignment','top','HorizontalAlignment','center')
end
% plot outline and section location
plot3(Out_grav(:,1),Out_grav(:,2),Out_grav(:,3),'k','LineWidth',2) 
plot3([xc(padW+1) xc(end-padE)],[yc(DO27section1) yc(DO27section1)],[zc(1) zc(1)],'k')
text(min(XX(:)),yc(DO27section1),zc(Zsection),'A','VerticalAlignment','top','HorizontalAlignment','left')
text(max(XX(:)),yc(DO27section1),zc(Zsection),['A' char(39)],'VerticalAlignment','top','HorizontalAlignment','right')

% plot3([xc(padW+1) xc(end-padE)],[yc(DO27section1) yc(DO27section1)],[zc(1) zc(1)],'k')
% text(min(XX(:)),yc(DO27section1),zc(Zsection),'B','VerticalAlignment','top','HorizontalAlignment','left')
% text(max(XX(:)),yc(DO27section1),zc(Zsection),['B' char(39)],'VerticalAlignment','top','HorizontalAlignment','right')

plot3([xc(padW+1) xc(end-padE)],[yc(DO18section) yc(DO18section)],[zc(1) zc(1)],'k')
% text(min(XX(:)),yc(DO18section),zc(Zsection),'A','VerticalAlignment','top','HorizontalAlignment','left')
% text(max(XX(:)),yc(DO18section),zc(Zsection),['A' char(39)],'VerticalAlignment','top','HorizontalAlignment','right')
ylim([yc(padS+1) 7134125])


% if mmax < max(m2D(:))
%     mmax = max(m2D(:));
% end
%% Plot vertical section DO27
ax2 = axes('Position',[0.09+(ii-1)*0.3 .225 .24 .24]);

m2D = squeeze(reshape(m(padT+1:end-padB,padW+1:end-padE,DO27section1),nz-padT-padB,nx-padE-padW));

[ZZ,XX] = ndgrid(zc(padT+1:end-padB),xc(padW+1:end-padE));
YY = ones(size(XX)) * yc(DO27section1);

h = surf(XX,YY,ZZ,m2D,'EdgeColor','none'); hold on
set(gca,'YDir','normal')

if ii == 1
zlabel('$Depth~(m)$', 'interpreter', 'latex','FontSize',14)
else
    set(gca,'ZTickLabel',[])
end


xlabel('$Easting~(m)$', 'interpreter', 'latex','FontSize',12)
axis equal tight
view([0 0])
caxis([min(cvec) max(cvec)])
% colorbar('SouthOutside')

text(min(XX(:)),YY(1)-10,max(ZZ(:)),'A','VerticalAlignment','bottom','HorizontalAlignment','left')
text(max(XX(:)),YY(1)-10,max(ZZ(:)),['A' char(39)],'VerticalAlignment','bottom','HorizontalAlignment','right')
scatter3(DO27(:,1),DO27(:,2),DO27(:,3),1,'k')



%% Define color map for plotting
% caxis([-8 mmax]);
% cvec = mmax*[0 0.2 0.4 0.6 0.8 1.05];

bb = cpal;
% bb = ([1 1 1;jet]);

colormap(ax1,bb);
colormap(ax2,bb);

% colormap(ax4,bb);

end
%% Add colorbar
axbar = axes('Position',[0.3 -.25 .42 .42]);

% bb = interp1([cvec'],[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,0:1e-4:mmax);

cbar = colorbar(axbar,'NorthOutside');
colormap(axbar,bb);

set(cbar,'Ticks',[0.1 0.3 0.5 0.7 0.9])
set(cbar,'TickLabels',cvec)
set(gca,'Visible','off');
text(0.50,1,barlabel, 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')