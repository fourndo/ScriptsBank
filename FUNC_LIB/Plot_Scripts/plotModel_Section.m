% function plotModel_Section(mmat,IDX,lbel)
% Inputs:
%   model - model Matrix
%   padcells - cluster result
clear all
% close all

addpath ..
% close all

padE=57;
padW=42;
padN=5;
padS=5;
padT=18;
padB=10;

zsection = 10;

DO27section1 = 47;
DO27section2 = 12;
DO18section = 79;

inp_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\DIGHEM\1D';
mod1_file = 'Inv_MOD_iter26.con';
% mod1_file = 'Ref_1em3\Inv_MOD_iter5.con';
mod2_file=[];
meshfile = 'UBC_mesh_small_v2.msh';


% padE=12;
% padW=12;
% padN=4;
% padS=4;
% padT=0;
% padB=4;
% Zsection = 9;
% 
% DO27section1 = 22;
% % DO27section2 = 24;
% DO18section = 54;
% 
% % inp_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\EOSC 556b (2013) EM\TKC_Project\Codes\CLUSTERING\Models';
% inp_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\DIGHEM\1D\Ref_5em4';
% % mod_file = 'dighemcondtensor.con';
% % mod_file='dighem_commonmesh.con';
% % mod_file = 'Kmeans_FULL_Grav_Mag_DIGHEM_5Clusters.dat';
% mod_file = 'Inv_MOD_iter7.con';
% % mod_file = 'Dens_Susc_4Clusters.dat';
% % mod_file = 'jtDenSusCon4Clusters.dat';
% % mod_file = 'Con2Clusters.dat';
% meshfile = 'UBC_mesh_small_v2.msh';

% Load point file for mapping
secfile = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\EOSC 556b (2013) EM\TKC_Project\Modelling\Geology_Objects\Section_trace.dat';
% secfile = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\EOSC 556b (2013) EM\TKC_Project\Modelling\Geology_Objects\Density_Estimate_Section.dat';
outline = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\EOSC 556b (2013) EM\TKC_Project\Modelling\Geology_Objects\DO27_Outline.dat';
hydro = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\EOSC 556b (2013) EM\TKC_Project\Modelling\Geology_Objects\Lake_trace.dat';

ndv = 1e-8;


cvec = [-5 -4 -3 -2];
cpal = colormap(jet);
caxis([min(cvec) max(cvec)])
%[255 255 255;0 71 189;99 209 62;255 204 0;255 0 0;255 153 200]/255;
% bb = interp1([cvec'],[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,0:1e-4:mmax);
crange = size(cpal,1);
barlabel = '$log10(\sigma) (S/m)$';

% cpal = [255 255 255;0 71 189;99 209 62;255 204 0;255 0 0;255 153 200]/255;
% % cpal = 1:-0.20:0.2;
% % cpal = kron([1 1 1],cpal');
% cvec = [0 1 2 3 4 5];
% 
% barlabel = 'Rock Unit';

%% Load in section
DO27 = load([secfile]);
Out_grav = load([outline]);
Out_lake = load([hydro]);

%% Read mesh
[xn,yn,zn] = read_UBC_mesh([inp_dir filesep meshfile]);
dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);

xc = ( xn(1:end-1) + xn(2:end) ) / 2;
yc = ( yn(1:end-1) + yn(2:end) ) / 2;
zc = ( zn(1:end-1) + zn(2:end) ) / 2;

set(figure, 'Position', [50 25 850 550]);

%% Load model
m = load([inp_dir '\' mod1_file]);
m(m==-100) = ndv;
temp = m;
% temp(m==1) = 2;
% temp(m==2) = 5;
% temp(m==3) = 5;
% temp(m==4) = 3;
% temp(m==5) = 3;
m =temp;

% crange = length(unique(m(m~=ndv)));
nullcell = m~=ndv;
% Convert to log

if ~isempty(mod2_file)
    
    m2 = load([inp_dir '\' mod2_file]);

    DOI = abs(log10(m) - log10(m2)) / abs(log10(1e-3) - log10(5e-4));
    DOI(nullcell==0) = 0;
    DOI = DOI/max(DOI);
    
else
    
    DOI = ones(size(nullcell));
    
end

m = log10(m);
m(nullcell==0) = nan;


cpal = colormap(jet);
crange = size(cpal,1);
barlabel = '$log10(\sigma) (S/m)$';

% Reshape to a 2D model
m = reshape(m,nz,nx,ny);
DOI = reshape(DOI,nz,nx,ny);

mmax = -100;

%% Plot horizontal section
ax1 = axes('Position',[-0.15 .1 0.85 0.85]);

m_pad = m(:,padW+1:end-padE,padS+1:end-padN);

[m2D] = get_model_top(m_pad,size(m_pad,2),size(m_pad,3),size(m_pad,1),nan,zsection);
% m2D = squeeze(reshape(m(Zsection,padW+1:end-padE,padS+1:end-padN),nx-padW-padE,ny-padN-padS));

[XX,YY] = ndgrid(xc(padW+1:end-padE),yc(padS+1:end-padN));
ZZ = ones(size(XX)) * zc(zsection);

h = surf(XX,YY,ZZ,m2D,'EdgeColor','none'); hold on
set(h,'alphadata',~isnan(m2D))
view([0 90])
set(gca,'YDir','normal')
ylabel('$Northing~(m)$', 'interpreter', 'latex','FontSize',14)
xlabel('$Easting~(m)$', 'interpreter', 'latex','FontSize',14)

cc = colormap(cpal(1:crange,:));
caxis([min(cvec) max(cvec)])
% colorbar('SouthOutside')
title('$Depth\;\approx\;25\;m$','interpreter', 'latex')

% plot outline and section location
plot3(Out_grav(:,1),Out_grav(:,2),Out_grav(:,3).^0*zc(1),'k--','LineWidth',3) 
plot3(Out_lake(:,1),Out_lake(:,2),Out_lake(:,3).^0*zc(1),'r.','MarkerSize',2) 
plot3([xc(padW+1) xc(end-padE)],[yc(DO27section1) yc(DO27section1)],[zc(1) zc(1)],'w--','LineWidth',2)
text(min(XX(:)),yc(DO27section1),zc(zsection),'B','VerticalAlignment','middle','HorizontalAlignment','right','FontSize',12)
text(max(XX(:)),yc(DO27section1),zc(zsection),['B' char(39)],'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',12)

% plot3([xc(padW+1) xc(end-padE)],[yc(DO27section1) yc(DO27section1)],[zc(1) zc(1)],'k')
% text(min(XX(:)),yc(DO27section1),zc(Zsection),'B','VerticalAlignment','top','HorizontalAlignment','left')
% text(max(XX(:)),yc(DO27section1),zc(Zsection),['B' char(39)],'VerticalAlignment','top','HorizontalAlignment','right')

plot3([xc(padW+1) xc(end-padE)],[yc(DO18section) yc(DO18section)],[zc(1) zc(1)],'w--','LineWidth',2)
text(min(XX(:)),yc(DO18section),zc(zsection),'A','VerticalAlignment','middle','HorizontalAlignment','right','FontSize',12)
text(max(XX(:)),yc(DO18section),zc(zsection),['A' char(39)],'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',12)

axis equal
axis([min(XX(:)) max(XX(:)) min(YY(:)) max(YY(:))])

% if mmax < max(m2D(:))
%     mmax = max(m2D(:));
% end
%% Plot vertical section DO27
ax2 = axes('Position',[0.52 .175 .45 .45]);

m2D = squeeze(reshape(m(padT+1:end-padB,padW+1:end-padE,DO27section1),nz-padT-padB,nx-padE-padW));

doi2D = squeeze(reshape(DOI(padT+1:end-padB,padW+1:end-padE,DO27section1),nz-padT-padB,nx-padE-padW));
doi2D = doi2D/max(doi2D(:));

[ZZ,XX] = ndgrid(zc(padT+1:end-padB),xc(padW+1:end-padE));
YY = ones(size(XX)) * yc(DO27section1);

h = surf(XX,YY,ZZ,m2D,'EdgeColor','none','FaceAlpha','flat'); hold on
if ~isempty(mod2_file)
    
    set(h,'alphadata',~isnan(m2D).*(1-doi2D))

else
    
    set(h,'alphadata',~isnan(m2D))
    
end
set(gca,'YDir','normal')
zlabel('$Depth~(m)$', 'interpreter', 'latex','FontSize',14)
xlabel('$Easting~(m)$', 'interpreter', 'latex','FontSize',12)
axis equal tight
view([0 0])
cc = colormap(cpal(1:crange,:));
caxis([min(cvec) max(cvec)])
% colorbar('SouthOutside')

text(min(XX(:)),YY(1)-10,max(ZZ(:)),'B','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',12)
text(max(XX(:)),YY(1)-10,max(ZZ(:)),['B' char(39)],'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',12)
scatter3(DO27(:,1),DO27(:,2),DO27(:,3),1,'k')
% if mmax < max(m2D(:))
%     mmax = max(m2D(:));
% end

%% Plot vertical section DO27
% ax3 = axes('Position',[0.55 .4 .35 .35]);
% 
% m2D = squeeze(reshape(m(padT+1:end-padB,padW+1:end-padE,DO27section2),nz-padT-padB,nx-padE-padW));
% 
% [ZZ,XX] = ndgrid(zc(padT+1:end-padB),xc(padW+1:end-padE));
% YY = ones(size(XX)) * yc(DO27section2);
% 
% h = surf(XX,YY,ZZ,m2D,'EdgeColor','none'); hold on
% set(gca,'YDir','normal')
% zlabel('$Depth~(m)$', 'interpreter', 'latex','FontSize',12)
% set(gca,'XTickLabel',[])
% axis equal tight
% view([0 0])
% caxis([min(cvec) max(cvec)])
% % colorbar('SouthOutside')
% text(min(XX(:)),YY(1),max(ZZ(:)),'B','VerticalAlignment','top','HorizontalAlignment','left')
% text(max(XX(:)),YY(1),max(ZZ(:)),['B' char(39)],'VerticalAlignment','top','HorizontalAlignment','right')


%% Plot vertical section DO18
ax3 = axes('Position',[0.52 .55 .45 .45]);

m2D = squeeze(reshape(m(padT+1:end-padB,padW+1:end-padE,DO18section),nz-padT-padB,nx-padE-padW));

doi2D = squeeze(reshape(DOI(padT+1:end-padB,padW+1:end-padE,DO27section1),nz-padT-padB,nx-padE-padW));
doi2D = doi2D/max(doi2D(:));

[ZZ,XX] = ndgrid(zc(padT+1:end-padB),xc(padW+1:end-padE));
YY = ones(size(XX)) * yc(DO18section);

h = surf(XX,YY,ZZ,m2D,'EdgeColor','none','FaceAlpha','flat'); hold on
if ~isempty(mod2_file)
    
    set(h,'alphadata',~isnan(m2D).*(1-doi2D))

else
    
    set(h,'alphadata',~isnan(m2D))
    
end
set(gca,'YDir','normal')
zlabel('$Depth~(m)$', 'interpreter', 'latex','FontSize',12)
set(gca,'XTickLabel',[])
% zlabel('$Depth~(m)$', 'interpreter', 'latex','FontSize',14)
% xlabel('$Easting~(m)$', 'interpreter', 'latex','FontSize',14)
cc = colormap(cpal(1:crange,:));
axis equal tight
view([0 0])
caxis([min(cvec) max(cvec)])
% colorbar('SouthOutside')
text(min(XX(:)),YY(1)-1,max(ZZ(:)),'A','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',12)
text(max(XX(:)),YY(1)-1,max(ZZ(:)),['A' char(39)],'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',12)



% if mmax < max(m2D(:))
%     mmax = max(m2D(:));
% end

%% Define color map for plotting
% caxis([-8 mmax]);
% cvec = mmax*[0 0.2 0.4 0.6 0.8 1.05];

% bb = cpal;
% bb = ([1 1 1;jet]);

% colormap(ax1,bb);
% colormap(ax2,bb);
% colormap(ax3,bb);

% colormap(ax4,bb);
%% Add colorbar
axbar = axes('Position',[0.52 -.25 .42 .42]);

% bb = interp1([cvec'],[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,0:1e-4:mmax);

cbar = colorbar(axbar,'NorthOutside');
colormap(axbar,cc);

set(cbar,'Ticks',[1/length(cvec)/2:1/length(cvec):1])
set(cbar,'TickLabels',cvec)
set(gca,'Visible','off');
text(0.50,1,barlabel, 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')