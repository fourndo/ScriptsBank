% Tile_Builder
% 
% Function takes gridded data set and designed a tiling scheme
% with constraints on the overlap and data coverage.
% Used for the MAG3D_Tile code.
%
% INPUT
% meshfile
% Obsfile
% tresh: Cutoff value used for the trend
% 
% OUTPUT
% Grid merged

clear all
close all

addpath '..\.'

%% INPUT VARIABLES
work_dir = 'C:\LC\Private\dominiquef\Projects\4414_Minsim\Modeling\GRAV\Data_Processing';

meshfile = 'ROT_Grid50m.msh';

modfile2 = 'Trout_infill.dat';
modfile1 = 'Grav08_Merged.dat';

dsep = '\';

tresh = 300;

clim = [-60 -40];
%% Load in files
[xn,yn,zn] = read_UBC_mesh([work_dir dsep meshfile]);
dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);

xc = ( xn(2:end) + xn(1:end-1) ) / 2;
yc = ( yn(2:end) + yn(1:end-1) ) / 2;


dx = min(dx);
dy = min(dy);

[Xc,Yc] = ndgrid(xc,yc);


%% Plot mesh and obs
m1 = load([work_dir dsep modfile1]); m1 = reshape(m1,nx,ny); m1(m1==-99999) = nan;
m2 = load([work_dir dsep modfile2]); m2 = reshape(m2,nx,ny); m2(m2==-99999) = nan;

m2(m2==-99999) = nan;
m1(m1==-99999) = nan;

%% Find collocated points on grid
indx = ~isnan(m1) & (m1 < tresh);
A = [Xc(indx) Yc(indx)];
trend1 = (A'*A)\A'*m1(indx);

indx = ~isnan(m2) & (m2 < tresh);
A = [Xc(indx) Yc(indx)];
trend2 = (A'*A)\A'*m2(indx);


% Compute polynomial correction 
out1 = (Yc(~isnan(m2))*trend1(2) + Xc(~isnan(m2))*trend1(1));
out2 = (Yc(~isnan(m2))*trend2(2) + Xc(~isnan(m2))*trend2(1));

corr = out2-out1;

m2_new = m2;
m2_new(~isnan(m2))  = m2(~isnan(m2)) - corr/2;


out_model = m2_new(:);
out_model(isnan(out_model)) = -99999;
save([work_dir dsep 'Model2_corr.dat'],'-ascii','out_model');
% set(figure, 'Position', [50 25 900 600]);
% h = mesh(xc,yc,trend'); hold on
%data = fwr-obs;
%%
set(figure, 'Position', [50 25 900 500]);
axes('Position',[0.45 .2 .75 .75])
h = imagesc(xc,yc,m2_new');hold on
contour(xc,yc,m2_new',-200:100:1800,'k')
set(h,'alphadata',~isnan(m2'))
colormap(jet);
set(gca,'Ydir','normal')
set(gca,'YTickLabel',[])
axis equal tight
% ylabel('North (m)');
xlabel('East (m)');
title('$Leveled\;model$', 'interpreter', 'latex','FontSize',12)
caxis(clim);
% view([60 60])

% set(figure, 'Position', [50 25 900 600]);
axes('Position',[-0.15 .2 .75 .75])
h = imagesc(xc,yc,m1'); hold on
contour(xc,yc,m1',-200:100:1800,'k')
set(h,'alphadata',~isnan(m1'))
colormap(jet);
% mesh(xc,yc,(Yc*trend1(2) + Xc*trend1(1))','FaceColor','k','FaceAlpha',0.25); hold on
set(gca,'Ydir','normal')
axis equal tight
ylabel('North (m)');
xlabel('East (m)');
title('$Model\;A$', 'interpreter', 'latex','FontSize',12)
hold on
% view([60 60])
caxis(clim);

% set(figure, 'Position', [50 25 900 600]);
axes('Position',[0.15 .2 .75 .75])
h = imagesc(xc,yc,m2');hold on
% contour(xc,yc,m2',-200:100:1800,'k')
set(h,'alphadata',~isnan(m2'))
aa = colormap(jet);
% mesh(xc,yc,(Yc*trend2(2) + Xc*trend2(1))','FaceColor','k','FaceAlpha',0.25); hold on
set(gca,'Ydir','normal')
axis equal tight
set(gca,'YTickLabel',[])
% ylabel('North (m)');
xlabel('East (m)');
title('$Model\;B$', 'interpreter', 'latex','FontSize',12)
hold on
% view([60 60])
caxis(clim);
%% Add colorbar
ax = axes('Position',[0.4 0.0 .25 .1]);

cbar = colorbar(ax,'NorthOutside');
colormap(ax,aa);

set(cbar,'Ticks',[0 0.5 1])
set(cbar,'TickLabels',[min(clim) mean(clim) max(clim)])
% caxis([min(m(:)) max(m(:))]);
% colormap(ax,[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0]/255);
% cbar = colorbar(ax,'NorthOutside')
% 
% set(cbar,'Ticks',[0 0.33 0.66 1])
% set(cbar,'TickLabels',[0 0.02 0.05 0.075])

set(gca,'Visible','off');
text(0.50,1.1,'$nT$', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
% text(0.5,1.65,'$(b)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center')

%% Interpolate on transition

msh{1} = get_UBC_mesh([work_dir dsep meshfile]); msh{2} = msh{1};
model{1} = m1(:); model{2} = m2_new(:);

out_model = Tile_Model_Merge_Func(msh,model,msh{1}, nan, 100,'all', 'allpad');
out_model = reshape(out_model,nx,ny);


%%
set(figure, 'Position', [50 25 900 500]);
axes('Position',[0.45 .2 .75 .75])
h = imagesc(xc,yc,out_model');hold on
% contour(xc,yc,out_model',-200:100:1800,'k')
set(h,'alphadata',~isnan(out_model'))
colormap(jet);
set(gca,'Ydir','normal')
set(gca,'YTickLabel',[])
axis equal tight
% ylabel('North (m)');
xlabel('East (m)');
title('$Leveled\;model$', 'interpreter', 'latex','FontSize',12)
caxis(clim);
% view([60 60])

% set(figure, 'Position', [50 25 900 600]);
axes('Position',[-0.15 .2 .75 .75])
h = imagesc(xc,yc,m1'); hold on
% contour(xc,yc,m1',-200:100:1800,'k')
set(h,'alphadata',~isnan(m1'))
colormap(jet);
% mesh(xc,yc,(Yc*trend1(2) + Xc*trend1(1))','FaceColor','k','FaceAlpha',0.25); hold on
set(gca,'Ydir','normal')
axis equal tight
ylabel('North (m)');
xlabel('East (m)');
title('$Model\;A$', 'interpreter', 'latex','FontSize',12)
hold on
% view([60 60])
caxis(clim);

% set(figure, 'Position', [50 25 900 600]);
axes('Position',[0.15 .2 .75 .75])
h = imagesc(xc,yc,m2');hold on
% contour(xc,yc,m2',-200:100:1800,'k')
set(h,'alphadata',~isnan(m2'))
aa = colormap(jet);
% mesh(xc,yc,(Yc*trend2(2) + Xc*trend2(1))','FaceColor','k','FaceAlpha',0.25); hold on
set(gca,'Ydir','normal')
axis equal tight
set(gca,'YTickLabel',[])
% ylabel('North (m)');
xlabel('East (m)');
title('$Model\;B$', 'interpreter', 'latex','FontSize',12)
hold on
% view([60 60])
caxis(clim);
%% Add colorbar
ax = axes('Position',[0.4 0.0 .25 .1]);

cbar = colorbar(ax,'NorthOutside');
colormap(ax,aa);

set(cbar,'Ticks',[0 0.5 1])
set(cbar,'TickLabels',[min(clim) mean(clim) max(clim)])
% caxis([min(m(:)) max(m(:))]);
% colormap(ax,[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0]/255);
% cbar = colorbar(ax,'NorthOutside')
% 
% set(cbar,'Ticks',[0 0.33 0.66 1])
% set(cbar,'TickLabels',[0 0.02 0.05 0.075])

set(gca,'Visible','off');
text(0.50,1.1,'$nT$', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
% text(0.5,1.65,'$(b)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center')

%% Export model
out_model = out_model(:);
out_model(isnan(out_model)) = -99999;
save([work_dir dsep 'Merged.dat'],'-ascii','out_model');


