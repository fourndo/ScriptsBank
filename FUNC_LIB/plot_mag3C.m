function plot_mag3C(obsx,obsy,d3C,I,D,head)
% Scatter plot of 3-component field mag.
% Also output TMI, and magnitude of field
% Input:
% obsx, obsy: X and Y location of observation points
% d3C: 3-components magnetic fields in nT
% I: Inclinaison of inducing field
% D: Declinaison of inducing field
% head: Header used for plot title

ndata = length(obsx);

nx = 200;
ny = 200;

xmin = ( min(obsx) );
xmax = ( max(obsx) );
ymin = ( min(obsy) );
ymax = ( max(obsy) );

dx = (( xmax - xmin) / nx);
dy = (( ymax - ymin) / ny);

x = xmin + cumsum(ones(1,nx)*dx);
y = ymin + cumsum(ones(1,ny)*dy);

[X,Y] = meshgrid(x,y);

% Compute TMI and |B|
P = [spdiags(ones(ndata,1)* (cosd(I) * cosd(D)),0,ndata,ndata) ...
    spdiags(ones(ndata,1)* (cosd(I) * sind(D)),0,ndata,ndata) ...
    spdiags(ones(ndata,1)* sind(I),0,ndata,ndata)];

TMI = P*d3C;

Bx = d3C(1:ndata);
By = d3C((ndata+1):2*ndata);
Bz = d3C(2*ndata+1:3*ndata);

magB = sqrt( Bx.^2 + By.^2 + Bz.^2 );

% Grid data
Bx_interp = griddata(obsx, obsy, Bx,X,Y,'linear'); 

By_interp = griddata(obsx, obsy, By,X,Y,'linear');  

Bz_interp = griddata(obsx, obsy, Bz,X,Y,'linear');  

TMI_interp = griddata(obsx, obsy, TMI,X,Y,'linear');  

magB_interp = griddata(obsx, obsy, magB,X,Y,'linear');  


set(figure, 'Position', [25 100 1800 900])
subplot(3,2,1);

imagesc(x,y,Bx_interp);hold on
colormap(hsv);
h = gca;
caxis([min(Bx_interp(:))*1.1 max(Bx_interp(:))])
cmap = jet;
cmap = [ 1 1 1 ; cmap ];
colormap(h,cmap);
scatter(obsx,obsy,2,'k.');axis([xmin xmax ymin ymax]);
% xlabel('\bfEasting (m)')
% ylabel('\bfNorthing (m)')
% caxis([-1000 1500]);
set(gca,'YDir','normal')
cbar = colorbar;
title(['\bf' head ' - Bx'])
axis([min(x) max(x) min(y) max(y)])
axis equal

subplot(3,2,3);
% scatter(obsx,obsy,30,By,'filled')
imagesc(x,y,By_interp);hold on
h = gca;
caxis([min(By_interp(:))*1.1 max(By_interp(:))])
cmap = jet;
cmap = [ 1 1 1 ; cmap ];
colormap(h,cmap);
scatter(obsx,obsy,2,'k.');axis([xmin xmax ymin ymax]);
set(gca,'YDir','normal')
% xlabel('\bfEasting (m)')
ylabel('\bfNorthing (m)')
% caxis([-1000 1500]);
colorbar
title(['\bf' head ' - By'])
axis([min(x) max(x) min(y) max(y)])
axis equal

subplot(3,2,5);
% scatter(obsx,obsy,30,Bz,'filled')
imagesc(x,y,Bz_interp);hold on
h = gca;
caxis([min(Bz_interp(:))*1.1 max(Bz_interp(:))])
cmap = jet;
cmap = [ 1 1 1 ; cmap ];
colormap(h,cmap);
scatter(obsx,obsy,2,'k.');axis([xmin xmax ymin ymax]);
set(gca,'YDir','normal')
xlabel('\bfEasting (m)')
% ylabel('\bfNorthing (m)')
% caxis([-1000 1500]);
colorbar
title(['\bf' head ' - Bz'])
axis([min(x) max(x) min(y) max(y)])
axis equal

subplot(2,2,2);
% scatter(obsx,obsy,30,TMI,'filled')
imagesc(x,y,TMI_interp);hold on
h = gca;
caxis([min(TMI_interp(:))*1.1 max(TMI_interp(:))])
cmap = jet;
cmap = [ 1 1 1 ; cmap ];
colormap(h,cmap);
scatter(obsx,obsy,2,'k.');axis([xmin xmax ymin ymax]);
set(gca,'YDir','normal')
xlabel('\bfEasting (m)')
ylabel('\bfNorthing (m)')
% caxis([-1100 800]);
colorbar
title(['\bf' head ' - TMI'])
axis([min(x) max(x) min(y) max(y)])
axis equal

subplot(2,2,4);
% scatter(obsx,obsy,30,magB,'filled')
imagesc(x,y,magB_interp);hold on
h = gca;
caxis([min(TMI_interp(:)) max(magB_interp(:))])
cmap = jet;
cmap = [ 1 1 1 ; cmap ];
colormap(h,cmap);
scatter(obsx,obsy,2,'k.');axis([xmin xmax ymin ymax]);
set(gca,'YDir','normal')
xlabel('\bfEasting (m)')
ylabel('\bfNorthing (m)')
% caxis([-800 2000]);
colorbar
title(['\bf' head ' - Field Magnitude'])
axis([min(x) max(x) min(y) max(y)])
axis equal