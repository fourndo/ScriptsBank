%% Plot Data

clear all
close all


addpath ..\..\FUNC_LIB;

%% Input Files
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Data\GRAV';

outline = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\EOSC 556b (2013) EM\TKC_Project\Modelling\Geology_Objects\DO27_Outline.dat';
Out_grav = load([outline]);

obsfile = 'newG.grv';

predfile = 'Tile17RegRem.obs';

barlabel = 'mGal';

% Load data
[obsx, obsy, obsz, d, wd] = read_GRAV3D_obs([work_dir '\' obsfile]);

% data = data(indx);

nx = 300;
ny = 300;

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

data_interp     = griddata(obsx, obsy, d,X,Y,'natural');


set(figure, 'Position', [50 0 700 800]); 

a = axes('Position',[0.1 .15 .85 .85]);
h =imagesc(x,y,data_interp);hold on
set(h,'alphadata',~isnan(data_interp))
caxis([min(data_interp(:)) max(data_interp(:))]);
colormap(a,jet);
scatter(obsx,obsy,25,'k.')
set(gca,'YDir','normal')
% xlabel('\bfEasting (m)')
ylabel('$y$', 'interpreter', 'latex','FontSize',14)
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
set(get(gca,'YLabel'),'Rotation',360);
axis([min(x) max(x) min(y) max(y)])
grid on
axis equal tight
% title('$(a)$', 'interpreter', 'latex','FontSize',14)
% text(min(xx)-dx*20, mean(yy),'$(a)$', 'interpreter', 'latex','FontSize',14)

contour(X,Y,data_interp,[-.85 -.85],'k','LineWidth',4)

% plot3(Out_grav(:,1),Out_grav(:,2),Out_grav(:,3),'k','LineWidth',2) 

%% ADD colorbar
axbar = axes('Position',[0.325 0.05 .42 .42]);

% bb = interp1([cvec'],[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,0:1e-4:mmax);

cbar = colorbar(axbar,'SouthOutside');
caxis([min(data_interp(:)) max(data_interp(:))]);
colormap(axbar,jet);

% set(cbar,'Ticks',[1/length(cvec)/2:1/length(cvec):1])
% set(cbar,'TickLabels',cvec)
set(gca,'Visible','off');
text(1.1,-.2,barlabel, 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
