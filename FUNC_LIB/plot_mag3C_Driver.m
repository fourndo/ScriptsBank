% Scatter plot of 3-component field mag.
% Also output TMI, and magnitude of field
% Input:
% obsx, obsy: X and Y location of observation points
% d3C: 3-components magnetic fields in nT
% I: Inclinaison of inducing field
% D: Declinaison of inducing field
% head: Header used for plot title
clear all
close all

%% INPUT FILES
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Equi_source\SingleBlock_Li\Corner_off';
obsfile = 'Obs_loc_3C.obs';
prefile = 'OMES_3C.pre';

%% SCRIPT STARTS HERE
[H, I, Dazm, D, obsx, obsy, obsz, data, wd] = read_MAG3D_obs([work_dir '\' obsfile]);
[~, ~, ~, ~, ~, ~, ~, pred, wd] = read_MAG3D_obs([work_dir '\' prefile]);

d{1} = data;
d{2} = pred;
d{3} = (data - pred) ./wd;

head{1} = 'Obs';
head{2} = 'Pre';
head{3} = 'Norm Res';

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

set(figure, 'Position', [10 30 775 975]);

for ii = 1 : 3
    
TMI_data = P*d{ii};

Bx = d{ii}(1:ndata);
By = d{ii}((ndata+1):2*ndata);
Bz = d{ii}(2*ndata+1:3*ndata);

magB = sqrt( Bx.^2 + By.^2 + Bz.^2 );

% Grid data
Bx_interp = griddata(obsx, obsy, Bx,X,Y,'linear'); 
% Bx_interp = flipud(Bx_interp);
By_interp = griddata(obsx, obsy, By,X,Y,'linear');  
% By_interp = flipud(By_interp);
Bz_interp = griddata(obsx, obsy, Bz,X,Y,'linear');  
% Bz_interp = flipud(Bz_interp);
TMI_interp = griddata(obsx, obsy, TMI_data,X,Y,'linear');  
% TMI_interp = flipud(TMI_interp);
magB_interp = griddata(obsx, obsy, magB,X,Y,'linear');  
% magB_interp = flipud(magB_interp);

if ii == 1
    axes('Position',[0.075 0.8 .2 .2])
    imagesc(x,y,Bx_interp);hold on
    lim_Bx = [min(Bx_interp(:))*1.1 max(Bx_interp(:))];
    caxis(lim_Bx);
    ylabel('\bfBx')
    set(get(gca,'YLabel'),'Rotation',360);
elseif ii==2
    axes('Position',[0.05+(ii-1)/4 0.8 .285 .2])
    imagesc(x,y,Bx_interp);hold on
    caxis(lim_Bx);
    set(gca,'YTickLabel',[]);
    h = colorbar('EastOutside');
    set(get(h,'title'),'string','\bf(nT)');
    
else
    axes('Position',[0.65 0.8 .285 .2])
    imagesc(x,y,Bx_interp);hold on
    set(gca,'YTickLabel',[]);
    h = colorbar('EastOutside');
    set(get(h,'title'),'string','\bf(nT)');
    
end

h = gca;
% scatter(obsx,obsy,2,'k.');axis([xmin xmax ymin ymax]);
% xlabel('\bfEasting (m)')
% ylabel('\bfNorthing (m)')
% caxis([-1000 1500]);
% caxis([min(Bx_interp(:))*1.1 max(Bx_interp(:))])
cmap = jet;
cmap = [ 1 1 1 ; cmap ];
colormap(h,cmap);
axis equal
set(gca,'XTickLabel',[]);
set(gca,'Ydir','normal')
axis([min(x) max(x) min(y) max(y)])
title(['\bf' head{ii}])


if ii == 1
    axes('Position',[0.075 0.6 .2 .2])
    imagesc(x,y,By_interp);hold on
        lim_By = [min(By_interp(:))*1.1 max(By_interp(:))];
    caxis(lim_By);
        ylabel('\bfBy')
    set(get(gca,'YLabel'),'Rotation',360);
elseif ii==2
    axes('Position',[0.05+(ii-1)/4 0.6 .285 .2])
    imagesc(x,y,By_interp);hold on
    caxis(lim_By);
    set(gca,'YTickLabel',[]);
     h = colorbar('EastOutside');
    set(get(h,'title'),'string','\bf(nT)');
    
else
    axes('Position',[0.65 0.6 .285 .2])
    imagesc(x,y,By_interp);hold on
    set(gca,'YTickLabel',[]);
    h = colorbar('EastOutside');
    set(get(h,'title'),'string','\bf (nT)');
    
end
% scatter(obsx,obsy,30,By,'filled')
% scatter(obsx,obsy,2,'k.');axis([xmin xmax ymin ymax]);
% xlabel('\bfEasting (m)')
% ylabel('\bfNorthing (m)')
set(gca,'XTickLabel',[]);
% caxis([-1000 1500]);
h = gca;
% caxis([min(By_interp(:))*1.1 max(By_interp(:))])
cmap = jet;
cmap = [ 1 1 1 ; cmap ];
colormap(h,cmap);
axis equal
set(gca,'Ydir','normal')
axis([min(x) max(x) min(y) max(y)])
hold on


if ii == 1
    axes('Position',[0.075 0.4 .2 .2])
    imagesc(x,y,Bz_interp);hold on
    ylabel('\bfBz')
    lim_Bz = [min(Bz_interp(:))*1.1 max(Bz_interp(:))];
    caxis(lim_Bz);
    set(get(gca,'YLabel'),'Rotation',360);
elseif ii == 2
    axes('Position',[0.05+(ii-1)/4 0.4 .285 .2])
    imagesc(x,y,Bz_interp);hold on
    set(gca,'YTickLabel',[]);
    caxis(lim_Bz);
     h = colorbar('EastOutside');
    set(get(h,'title'),'string','\bf(nT)');
    
else
    axes('Position',[0.65 0.4 .285 .2])
    imagesc(x,y,Bz_interp);hold on
    set(gca,'YTickLabel',[]);
    h = colorbar('EastOutside');
    set(get(h,'title'),'string','\bf(nT)');
    
end
% scatter(obsx,obsy,30,Bz,'filled')
% scatter(obsx,obsy,2,'k.');axis([xmin xmax ymin ymax]);
% xlabel('\bfEasting (m)')
% ylabel('\bfNorthing (m)')
% caxis([-1000 1500]);
h = gca;
% caxis([min(Bz_interp(:))*1.1 max(Bz_interp(:))])
cmap = jet;
cmap = [ 1 1 1 ; cmap ];
colormap(h,cmap);
axis equal
set(gca,'Ydir','normal')
axis([min(x) max(x) min(y) max(y)])
set(gca,'XTickLabel',[]);
hold on


if ii == 1
    axes('Position',[0.075 0.2 .2 .2])
    imagesc(x,y,TMI_interp);hold on
    ylabel('\bfTMI')
    lim_TMI = [min(TMI_interp(:))*1.1 max(TMI_interp(:))];
    caxis(lim_TMI);
    set(get(gca,'YLabel'),'Rotation',360);
elseif ii == 2
    axes('Position',[0.05+(ii-1)/4 0.2 .285 .2])
    imagesc(x,y,TMI_interp);hold on
    set(gca,'YTickLabel',[]);
    caxis(lim_TMI);
     h = colorbar('EastOutside');
    set(get(h,'title'),'string','\bf(nT)');

else
    axes('Position',[0.65 0.2 .285 .2])
    imagesc(x,y,TMI_interp);hold on
    set(gca,'YTickLabel',[]);
    h = colorbar('EastOutside');
    set(get(h,'title'),'string','\bf(nT)');

end
% scatter(obsx,obsy,30,TMI,'filled')
% scatter(obsx,obsy,2,'k.');axis([xmin xmax ymin ymax]);
% xlabel('\bfEasting (m)')
% ylabel('\bfNorthing (m)')
% caxis([-1100 800]);
h = gca;
% caxis([min(TMI_interp(:))*1.1 max(TMI_interp(:))])
cmap = jet;
cmap = [ 1 1 1 ; cmap ];
colormap(h,cmap);
axis equal
set(gca,'Ydir','normal')
axis([min(x) max(x) min(y) max(y)])
set(gca,'XTickLabel',[]);
hold on


if ii == 1
    axes('Position',[0.075 0.0 .2 .2])
    imagesc(x,y,magB_interp);hold on
    lim_magB = [min(TMI_interp(:)) max(magB_interp(:))];
    caxis(lim_magB);
    ylabel('\bf|B|')
    set(get(gca,'YLabel'),'Rotation',360);
elseif ii == 2
    axes('Position',[0.05+(ii-1)/4 0.0 .285 .2])
    imagesc(x,y,magB_interp);hold on
    set(gca,'YTickLabel',[]);
    caxis(lim_magB);
    h = colorbar('EastOutside');
    set(get(h,'title'),'string','\bf(nT)');
    
else
    axes('Position',[0.65 0.0 .285 .2])
    imagesc(x,y,magB_interp);hold on
    caxis([-1 max(magB_interp(:))])
    set(gca,'YTickLabel',[]);
    h = colorbar('EastOutside');
    set(get(h,'title'),'string','\bf(nT)');
    
end
% scatter(obsx,obsy,30,magB,'filled')
% scatter(obsx,obsy,2,'k.');axis([xmin xmax ymin ymax]);
% xlabel('\bfEasting (m)')
% ylabel('\bfNorthing (m)')
% caxis([-800 2000]);
h = gca;
% caxis([min(TMI_interp(:)) max(magB_interp(:))])
cmap = jet;
cmap = [ 1 1 1 ; cmap ];
colormap(h,cmap);
axis equal
set(gca,'Ydir','normal')
axis([min(x) max(x) min(y) max(y)])
set(gca,'XTickLabel',[]);
hold on

end