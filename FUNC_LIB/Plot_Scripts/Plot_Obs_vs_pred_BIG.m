% Tile_Builder
% 
% Function takes in a mesh and data points and designed a tiling scheme
% with constraints on the overlap and data coverage.
% Used for the MAG3D_Tile code.
%
% INPUT
% meshfile
% Obsfile
% Olap
% 
% OUTPUT
% Tile file: SW-NE corner of each tile [x1 y1 x2 y2; ...]

clear all
% close all

addpath '..\.'

%% INPUT VARIABLES
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Paul_Lake\Modeling\Inversion\Misery\l0l2';

% obsfile = '..\..\Obs_Paul_Lake_SUB_5pc_5nT_DETREND.dat';
obsfile = 'Tile_data.dat';
prefile = 'Tile1_RegRem.tmi';
dsep = '\';

lx = 1000;

clim = [-50 75];

%% Plot mesh and obs
set(figure, 'Position', [50 25 900 600]);

[H, I, Dazm, D, Obsx, Obsy, Obsz, data, wd_full] = read_MAG3D_obs([work_dir dsep obsfile]);

xx = min(Obsx):50:max(Obsx);
yy = min(Obsy):50:max(Obsy);
[YY,XX] = ndgrid(yy,xx);

F = scatteredInterpolant(Obsy, Obsx, data ,'natural');

grid_d = F(YY,XX);

F = scatteredInterpolant(Obsy, Obsx, wd_full ,'natural');

grid_wd = F(YY,XX);
%% Flag cell grid too far from obs
flag = zeros(size(YY));

for ii = 1 : length(Obsx)
    
    indx = XX > Obsx(ii) - 250 & XX < Obsx(ii) + 250 & YY > Obsy(ii) - 250 & YY < Obsy(ii) + 250;
    
    flag(indx) = 1;
    
end
% save([work_dir dsep 'Obs_flag'],'flag');
% load([work_dir dsep 'Obs_flag']);

%%
% msh = mesh(Xn,Yn,ones(size(Xn)),'FaceColor','none'); hold on
ax1 = axes('Position',[0.1 .55 .4 .4]);
view([0 90]) 
h = imagesc(xx,yy,grid_d);hold on
set(h,'alphadata',~isnan(grid_d) .* flag)
caxis(ax1,clim)
bb = colormap(jet);
xlim([min(Obsx) max(Obsx)]); rangex = max(Obsx) - min(Obsx); 
ylim([min(Obsy) max(Obsy)]); rangey = max(Obsy) - min(Obsy);
set(gca,'Ydir','normal')
axis equal tight
% set(gca,'XTickLabel',[])
ylabel('North (m)');
xlabel('East (m)');
% colorbar(ax1,'EastOutside');
title('$Observed\;data$', 'interpreter', 'latex','FontSize',12)
hold on
text(5.1e+5,7.16e+6,'$(a)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center')
%% Add color bar
ax = axes('Position',[0.375 .725 .1 .2]);

% bb = interp1([cvec'],[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,0:1e-4:mmax);

cbar = colorbar(ax,'EastOutside');
colormap(ax,jet);
caxis(ax,clim)
% set(cbar,'Ticks',[0 0.2 0.4 0.6 0.8 1])
% set(cbar,'TickLabels',round(cvec*10000)/10000)

set(gca,'Visible','off');
text(1.25,-.1,'$(nT)$', 'interpreter', 'latex','FontSize',10,'HorizontalAlignment','center')


%% Add color bar
ax = axes('Position',[0.85 .725 .1 .2]);

% bb = interp1([cvec'],[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,0:1e-4:mmax);

cbar = colorbar(ax,'EastOutside');
colormap(ax,jet);
caxis(ax,clim)
% set(cbar,'Ticks',[0 0.2 0.4 0.6 0.8 1])
% set(cbar,'TickLabels',round(cvec*10000)/10000)

set(gca,'Visible','off');
text(1.25,-.1,'$(nT)$', 'interpreter', 'latex','FontSize',10,'HorizontalAlignment','center')

%% Load predicted and grid
[~, ~, ~, ~, Obsx, Obsy, Obsz, pred, ~] = read_MAG3D_obs([work_dir dsep prefile]);

F = scatteredInterpolant(Obsy, Obsx, pred ,'natural');

grid_pred = F(YY,XX);

%%
% msh = mesh(Xn,Yn,ones(size(Xn)),'FaceColor','none'); hold on
ax2 = axes('Position',[0.58 .55 .4 .4]);
view([0 90]) 
h = imagesc(xx,yy,grid_pred);hold on
set(h,'alphadata',~isnan(grid_d) .* flag)
caxis(ax2,clim)
bb = colormap(jet);
xlim([min(Obsx) max(Obsx)]); rangex = max(Obsx) - min(Obsx); 
ylim([min(Obsy) max(Obsy)]); rangey = max(Obsy) - min(Obsy);
set(gca,'Ydir','normal')
axis equal tight
set(gca,'YTickLabel',[])
% ylabel('North (m)');
xlabel('East (m)');
% colorbar(ax1,'EastOutside');
title('$Levelled\;data$', 'interpreter', 'latex','FontSize',12)
hold on
text(5.1e+5,7.16e+6,'$(b)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center')
%% Plot residual
ax3 = axes('Position',[0.33 .075 .4 .4]);
view([0 90]) 
h = imagesc(xx,yy,(grid_d - grid_pred));hold on
set(h,'alphadata',~isnan(grid_d) .* flag)
caxis(ax3,[-20 20])
bb = colormap(jet);
xlim([min(Obsx) max(Obsx)]); rangex = max(Obsx) - min(Obsx); 
ylim([min(Obsy) max(Obsy)]); rangey = max(Obsy) - min(Obsy);
set(gca,'Ydir','normal')
axis equal tight
% set(gca,'XTickLabel',[])
ylabel('North (m)');
xlabel('East (m)');
% colorbar(ax1,'EastOutside');
title('$Regional\;Removal$', 'interpreter', 'latex','FontSize',12)
hold on
text(5.1e+5,7.16e+6,'$(c)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center')

%% Add color bar
ax = axes('Position',[0.6 .25 .1 .2]);

% bb = interp1([cvec'],[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,0:1e-4:mmax);

cbar = colorbar(ax,'EastOutside');
colormap(ax,jet);
caxis(ax,[-20 20])
% set(cbar,'Ticks',[0 0.2 0.4 0.6 0.8 1])
% set(cbar,'TickLabels',round(cvec*10000)/10000)

set(gca,'Visible','off');
text(1.25,-.1,'$(nT)$', 'interpreter', 'latex','FontSize',10,'HorizontalAlignment','center')



