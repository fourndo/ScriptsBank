function plot_TMI(obsx,obsy,data,pred,wd,head)
% Scatter plot of 3-component field mag.
% Also output TMI, and magnitude of field
% Input:
% obsx, obsy: X and Y location of observation points
% d3C: 3-components magnetic fields in nT
% I: Inclinaison of inducing field
% D: Declinaison of inducing field
% head: Header used for plot title

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

r = data - pred;
norm_r = r .* wd;

% Grid data
data_interp     = griddata(obsx, obsy, data,X,Y,'linear'); 
% data_interp     = flipud(data_interp);

pred_interp     = griddata(obsx, obsy, pred,X,Y,'linear');  
% pred_interp     = flipud(pred_interp);

r_interp        = griddata(obsx, obsy, r ,X,Y,'linear');  
% r_interp        = flipud(r_interp);

norm_r_interp   = griddata(obsx, obsy, norm_r ,X,Y,'linear');  
% norm_r_interp   = flipud(norm_r_interp);

set(figure, 'Position', [25 100 1800 900])
subplot(2,2,1);

imagesc(x,y,data_interp);hold on
scatter(obsx,obsy,2,'k.')
set(gca,'YDir','normal')
% xlabel('\bfEasting (m)')
ylabel('\bfNorthing (m)')
caxis([min(data) max(data)]);
colorbar
axis([min(x) max(x) min(y) max(y)])
title(['\bf' head ' - Observed data'])
axis equal

subplot(2,2,2);
% scatter(obsx,obsy,30,By,'filled')
imagesc(x,y,pred_interp);hold on
scatter(obsx,obsy,2,'k.')
set(gca,'YDir','normal')
xlabel('\bfEasting (m)')
ylabel('\bfNorthing (m)')
caxis([min(data) max(data)]);
colorbar
axis([min(x) max(x) min(y) max(y)])
title(['\bf' head ' - Predicted data'])
axis equal

subplot(2,2,3);
% scatter(obsx,obsy,30,Bz,'filled')
imagesc(x,y,r_interp);hold on
scatter(obsx,obsy,2,'k.')
set(gca,'YDir','normal')
xlabel('\bfEasting (m)')
% ylabel('\bfNorthing (m)')
% caxis([-1000 1500]);
colorbar
axis([min(x) max(x) min(y) max(y)])
title(['\bf' head ' - Residual'])
axis equal

subplot(2,2,4);
% scatter(obsx,obsy,30,TMI,'filled')
imagesc(x,y,norm_r_interp);hold on
scatter(obsx,obsy,2,'k.')
set(gca,'YDir','normal')
xlabel('\bfEasting (m)')
ylabel('\bfNorthing (m)')
% caxis([-1100 800]);
colorbar
axis([min(x) max(x) min(y) max(y)])
title(['\bf' head ' - Normalized Residual'])
axis equal