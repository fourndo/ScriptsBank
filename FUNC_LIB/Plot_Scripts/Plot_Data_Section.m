%% Plot Data

clear all
close all


addpath ..\..\FUNC_LIB;

%% Input Files
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Nut_Cracker\Induced_MAG3C';
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Nut_Cracker\Induced_MAG3C';
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Paul_Lake\Modeling\Inversion\Tile_50m_VOI\ALL_Tiles_l0l2';

% meshfile = '..\Mesh_20m.msh';
meshfile = 'Tile17.msh';

% model = '\REM_l2l2\MAG3D_TMI_iter_4.sus';
% model = '\IND_l0l2_NEW\MAG3D_TMI_iter_12.sus';
% model = 'Tile6_MAI_esus.sus';
% 
% model_true = 'Tile6_MAI_esus.sus';

% norm_vec = 'Lp_vec.txt';
% wmodel = 'Wvec.txt';

obsfile = '..\..\Obs_Paul_Lake_SUB_5pc_5nT_DETREND.dat';

predfile = 'Tile17RegRem.obs';

% Load predicted lp
[~, ~, ~, ~, obsx, obsy, obsz, dpre, wd] = read_MAG3D_obs([work_dir '\' predfile]);

% Load data
[~, ~, ~, ~, Obsx, Obsy, Obsz, data, ~] = read_MAG3D_obs([work_dir '\' obsfile]);

%Indx
indx = Obsx >=min(obsx) & Obsx <=max(obsx) & Obsy >=min(obsy) & Obsy <=max(obsy);

% data = data(indx);

nx = 100;
ny = 100;

xmin = ( min(obsx) );
xmax = ( max(obsx) );
ymin = ( min(obsy) );
ymax = ( max(obsy) );

dx = (( xmax - xmin) / nx);
dy = (( ymax - ymin) / ny);

x = xmin + cumsum(ones(1,nx)*dx);
y = ymin + cumsum(ones(1,ny)*dy);

[X,Y] = meshgrid(x,y);

% res = (data - dpre);

data_interp     = griddata(Obsx(indx), Obsy(indx), data(indx),X,Y,'linear'); 
% data_interp(isnan(data_interp)) = min(data_interp(:))*2;

dpre_interp     = griddata(obsx, obsy, dpre,X,Y,'linear');

% res_interp        = griddata(obsx, obsy, res ,X,Y,'linear'); 
res_interp = data_interp - dpre_interp;

set(figure, 'Position', [50 0 850 600]); 

a = axes('Position',[0.025 .55 .375 .375]);
h =imagesc(x,y,data_interp);hold on
set(h,'alphadata',~isnan(data_interp))
caxis([min(dpre_interp(:)) max(dpre_interp(:))]);
colormap(a,jet);
scatter(obsx,obsy,2,'k.')
set(gca,'YDir','normal')
% xlabel('\bfEasting (m)')
ylabel('$y$', 'interpreter', 'latex','FontSize',14)
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
set(get(gca,'YLabel'),'Rotation',360);
axis([min(x) max(x) min(y) max(y)])
grid on
axis equal tight
title('$(a)$', 'interpreter', 'latex','FontSize',14)
% text(min(xx)-dx*20, mean(yy),'$(a)$', 'interpreter', 'latex','FontSize',14)
plot(x([1 end]),y([end 1]),'k')
text(x(1),y(end),'A' ,'interpreter', 'latex','FontSize',14,'HorizontalAlignment','right','VerticalAlignment','top')
text(x(end),y(1),'$A"$', 'interpreter', 'latex','FontSize',14)

a = axes('Position',[0.625 .55 .375 .375]);
h =imagesc(x,y,dpre_interp);hold on
set(h,'alphadata',~isnan(data_interp))
caxis([min(dpre_interp(:)) max(dpre_interp(:))]);
colormap(a,jet);
scatter(obsx,obsy,2,'k.')
set(gca,'YDir','normal')
% xlabel('\bfEasting (m)')
% ylabel('$y$', 'interpreter', 'latex','FontSize',14)
% xlabel('$x$', 'interpreter', 'latex','FontSize',14)
set(get(gca,'YLabel'),'Rotation',360);
axis([min(x) max(x) min(y) max(y)])
set(gca,'YTickLabel',[])
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
grid on
axis equal tight
title('$(c)$', 'interpreter', 'latex','FontSize',14)
% text(400, mean(yy),'$(a)$', 'interpreter', 'latex','FontSize',14)
plot(x([1 end]),y([end 1]),'k')
text(x(1),y(end),'B' ,'interpreter', 'latex','FontSize',14,'HorizontalAlignment','right','VerticalAlignment','top')
text(x(end),y(1),'$B"$', 'interpreter', 'latex','FontSize',14)


a = axes('Position',[0.325 .55 .375 .375]);
h =imagesc(x,y,res_interp);hold on
set(h,'alphadata',~isnan(res_interp))
caxis([min(res_interp(:)) max(res_interp(:))])
colormap(a,jet);
scatter(obsx,obsy,2,'k.')
set(gca,'YDir','normal')
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
% set(gca,'XTickLabel',[])
% ylabel('$y$', 'interpreter', 'latex','FontSize',14)
% xlabel('$x$', 'interpreter', 'latex','FontSize',14)
set(get(gca,'YLabel'),'Rotation',360);
axis([min(x) max(x) min(y) max(y)])
grid on
axis equal tight
title('$(b)$', 'interpreter', 'latex','FontSize',14)
% text(400, mean(yy),'$(b)$', 'interpreter', 'latex','FontSize',14)
set(gca,'YTickLabel',[])


%% Plot section through the data
data_sec = diag(rot90(fliplr(data_interp)),0);

a = axes('Position',[0.075 .1 .28 .28]);
plot(x,data_sec);
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
ylabel('$\vec B$', 'interpreter', 'latex','FontSize',14)
set(gca,'YTickLabelRotation',90)
text(x(1),75,'A' ,'interpreter', 'latex','FontSize',14,'HorizontalAlignment','right')
text(x(end),75,'$A"$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','right')
grid on

dpre_sec = diag(rot90(fliplr(dpre_interp)),0);

a = axes('Position',[0.675 .1 .28 .28]);
plot(x,dpre_sec);
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
ylabel('$\vec B$', 'interpreter', 'latex','FontSize',14)
set(gca,'YTickLabelRotation',90)
text(x(1),30,'B' ,'interpreter', 'latex','FontSize',14,'HorizontalAlignment','right')
text(x(end),30,'$B"$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','right')
grid on


%% Add colorbars
ax = axes('Position',[0.375 .21 .28 .28]);
cbar = colorbar('NorthOutside');
colormap(ax,jet);
set(cbar,'Ticks',[0.1 0.9])
set(cbar,'TickLabels',round([min(res_interp(:))*1.1 max(res_interp(:))*.9]))
set(gca,'Visible','off');
text(0.5,1.6,'Regional Field (nT)', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
% text(0.5,1.75,'$(c)$', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
%% Add colorbars
ax = axes('Position',[0.075 .21 .28 .28]);
cbar = colorbar('NorthOutside');
colormap(ax,jet);
set(cbar,'Ticks',[0.1 0.9])
set(cbar,'TickLabels',round([min(data_interp(:))*1.1 max(data_interp(:))*.9]))
set(gca,'Visible','off');
text(0.5,1.6,'Observed (nT)', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
% text(0.5,1.75,'$(c)$', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')

%% Add colorbars
ax = axes('Position',[0.675 .21 .28 .28]);
cbar = colorbar('NorthOutside');
colormap(ax,jet);
set(cbar,'Ticks',[0.1 0.9])
set(cbar,'TickLabels',round([min(dpre_interp(:))*1.1 max(dpre_interp(:))*.9]))
set(gca,'Visible','off');
text(0.5,1.6,'Residual (nT)', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
% text(0.5,1.75,'$(c)$', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
