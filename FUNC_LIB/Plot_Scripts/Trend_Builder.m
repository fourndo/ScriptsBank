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
close all

addpath '..\.'

%% INPUT VARIABLES
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Paul_Lake\Modeling\Inversion';

meshfile = 'Mesh_50m.msh';

obsfile = 'Obs_Paul_Lake_SUB_1pc_10nT.dat';

dsep = '\';

lx = 1000;

%% Load in files
[xn,yn,zn] = read_UBC_mesh([work_dir dsep meshfile]);
dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);

dx = min(dx);
dy = min(dy);

[Xn,Yn] = ndgrid(xn,yn);


%% Plot mesh and obs
set(figure, 'Position', [50 25 900 600]);

[H, I, Dazm, D, Obsx, Obsy, Obsz, data, wd_full] = read_MAG3D_obs([work_dir dsep obsfile]);

xx = min(Obsx):50:max(Obsx);
yy = min(Obsy):50:max(Obsy);
[YY,XX] = ndgrid(yy,xx);

F = scatteredInterpolant(Obsy, Obsx, data ,'natural');

grid_d = F(YY,XX);

%% Flag cell grid too far from obs
% flag = zeros(size(YY));
% 
% for ii = 1 : length(Obsx)
%     
%     indx = XX > Obsx(ii) - 250 & XX < Obsx(ii) + 250 & YY > Obsy(ii) - 250 & YY < Obsy(ii) + 250;
%     
%     flag(indx) = 1;
%     
% end
% save([work_dir dsep 'Obs_flag'],'flag');
load([work_dir dsep 'Obs_flag']);

%%
% msh = mesh(Xn,Yn,ones(size(Xn)),'FaceColor','none'); hold on
ax1 = axes('Position',[0.1 .55 .4 .4]);
view([0 90]) 
h = imagesc(xx,yy,grid_d);hold on
set(h,'alphadata',~isnan(grid_d) .* flag)
caxis(ax1,[-500 100])
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
text(5.1e+5,7.155e+6,'$(a)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center')
%% Add color bar
ax = axes('Position',[0.475 .6 .1 .3]);

% bb = interp1([cvec'],[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,0:1e-4:mmax);

cbar = colorbar(ax,'EastOutside');
colormap(ax,jet);
caxis(ax,[-500 100])
% set(cbar,'Ticks',[0 0.2 0.4 0.6 0.8 1])
% set(cbar,'TickLabels',round(cvec*10000)/10000)

set(gca,'Visible','off');
text(2,-.1,'$(nT)$', 'interpreter', 'latex','FontSize',10,'HorizontalAlignment','center')


%% Load extra data
load([work_dir '\Pipes43101']);


%% Pick regions of background and compute median value\
xc = 1;
count = 0;
while ~isempty(xc)
    
    count = count + 1;
    [xc,yc] = ginput(1);

    Xnew = [(xc - lx) (xc + lx)];
    Ynew = [(yc - lx) (yc + lx)];

    if isempty(xc)

        
        continue

    end
    
    
    
    aa = patch([Xnew(1);Xnew(1);Xnew(2);Xnew(2)],[Ynew(1);Ynew(2);Ynew(2);Ynew(1)],[1 1 1],'LineWidth',1,'EdgeColor','b');

    alpha(aa,0.5);
    
    % Select observation points within the area
    indx = Obsx > Xnew(1) & Obsx < Xnew(2) & Obsy > Ynew(1) & Obsy < Ynew(2);
    temp = data(indx);
    
    % Remove NaN selections
    indx = isnan(temp)==0;
    
    
    % Compute median value
    temp = median(temp(indx));
    
    target(count,:) = [xc yc temp];
    
    % Plot those points to screen
    scatter(xc,yc,50,temp,'filled');
    
end

caxis([-500 250])




% Compute 1th order polynomial
A = [target(:,1) target(:,2)];
m = (A'*A)\A'*target(:,3);

trend = (Obsy*m(2) + Obsx*m(1));
d_detrend = data - trend;

F = scatteredInterpolant(Obsy, Obsx, trend,'natural');

trend_d = F(YY,XX);


% msh = mesh(Xn,Yn,ones(size(Xn)),'FaceColor','none'); hold on
ax2 = axes('Position',[0.575 .55 .4 .4]);
view([0 90]) 
h = imagesc(xx,yy,trend_d);hold on
set(h,'alphadata',~isnan(trend_d) .* flag)
caxis(ax2,[-500 100])
colormap(jet)
% colorbar
xlim([min(Obsx) max(Obsx)]); rangex = max(Obsx) - min(Obsx); 
ylim([min(Obsy) max(Obsy)]); rangey = max(Obsy) - min(Obsy);
set(gca,'Ydir','normal')
axis equal tight
title('$1^{th}\;order\;polynomial$', 'interpreter', 'latex','FontSize',12)
% colorbar(ax2,'EastOutside');
set(gca,'YTickLabel',[])
xlabel('East (m)');
text(5.1e+5,7.155e+6,'$(b)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center')


%% Plot de-trended data
ax3 = axes('Position',[0.325 0.05 .4 .4]);
view([0 90]) 
h = imagesc(xx,yy,grid_d - trend_d);hold on
set(h,'alphadata',~isnan(grid_d) .* flag)
caxis(ax3,[-250 250])
bb = colormap(jet);
% colorbar
xlim([min(Obsx) max(Obsx)]); rangex = max(Obsx) - min(Obsx); 
ylim([min(Obsy) max(Obsy)]); rangey = max(Obsy) - min(Obsy);
set(gca,'Ydir','normal')
axis equal tight
xlabel('East (m)')
ylabel('North (m)');

% colorbar(ax3,'EastOutside');
title('$De-trended\;data$', 'interpreter', 'latex','FontSize',12)
hold on

text(5.1e+5,7.155e+6,'$(c)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center')
%% Add color bar
ax = axes('Position',[0.7 .1 .1 .3]);

% bb = interp1([cvec'],[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,0:1e-4:mmax);

cbar = colorbar(ax,'EastOutside');
colormap(ax,jet);
caxis(ax,[-250 250])
% set(cbar,'Ticks',[0 0.2 0.4 0.6 0.8 1])
% set(cbar,'TickLabels',round(cvec*10000)/10000)

set(gca,'Visible','off');
text(2,-.1,'$(nT)$', 'interpreter', 'latex','FontSize',10,'HorizontalAlignment','center')


%% Save polynomial and de-trended data
wd_detrend = abs(d_detrend) * 0.02 + 5;
write_MAG3D_TMI([work_dir dsep 'Obs_Paul_Lake_SUB_1pc_10nT_DETREND.dat'],H,I,Dazm,Obsx,Obsy,Obsz,d_detrend,wd_detrend);

write_MAG3D_TMI([work_dir dsep 'Trend.dat'],H,I,Dazm,Obsx,Obsy,Obsz,trend,wd_detrend);


