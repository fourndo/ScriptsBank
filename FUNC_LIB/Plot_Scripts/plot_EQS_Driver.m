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
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Nut_Cracker\Tiled_AMI\Tile1';
obsfile = '..\..\Obs_RAW_GRID_3C.obs';
prefile = 'Tile1_EMS_3C.tmi';
meshfile = 'Tile1.msh';
modfile = 'Tile1_EMS.sus';

%% SCRIPT STARTS HERE
[H, I, Dazm, D, obsx, obsy, obsz, data, wd] = read_MAG3D_obs([work_dir '\' obsfile]);
[~, ~, ~, ~, ~, ~, ~, pred, wd] = read_MAG3D_obs([work_dir '\' prefile]);

d{1} = data;
d{2} = pred;
d{3} = (data - pred) ./wd;

head{1} = 'Obs';
head{2} = 'Pre';
head{3} = 'Norm Res';
ii = 3;
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



%% Load mesh and model

[xn,yn,zn] = read_UBC_mesh([work_dir '\' meshfile]);
xc = (xn(2:end) + xn(1:end-1)) /2;
yc = (yn(2:end) + yn(1:end-1)) /2;
zc = (zn(2:end) + zn(1:end-1)) /2;

dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);

model = load([work_dir '\' modfile]);
model = reshape(model,nz,nx,ny);

% Extract first cell under topo
m = zeros(ny,nx);
for kk = 1 : ny
    
    for jj = 1 : nx
        
        index = model(model(:,jj,kk)~=-100,jj,kk);
        
        m(kk,jj) = index(1);
        
    end
    
end
        
        

%% Plot 
for ii = 1 :3

set(figure, 'Position', [10 30 775 975]);    
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

axes('Position',[0.070 0.45 .525 .525]);
h = gca;
% set(figure(1),'CurrentAxes',h)
imagesc(xc,yc,m);hold on
scatter(obsx,obsy,2,'r.'); hold on
ylabel('\bfY')
xlabel('\bfX')
set(get(gca,'YLabel'),'Rotation',360);

set(gca,'YDir','normal')
set(gca, 'XAxisLocation', 'top')
axis equal
axis([min(xc) max(xc) min(yc) max(yc)])
caxis([min(m(:)) max(m(:))]);
    
%% Add colorbar
axes('Position',[0.070 0.45 .525 .525]);
% colorbar('SouthOutside');
c = colorbar('SouthOutside');
caxis([min(m(:)) max(m(:))]);
set(gca,'Visible','off');
% set(get(c,'title'),'string','$\mathbf{\kappa}\;(SI)$', 'interpreter', 'latex','FontSize',12);
set(get(c,'title'),'string','$\mathbf{\kappa}_{eff}$', 'interpreter', 'latex','FontSize',12);
text(0.5,-0.25,'$(a)$', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center');


%%
axes('Position',[0.675 0.70 .275 .355]);
h = gca;
% set(figure(1),'CurrentAxes',h)
p = imagesc(x,y,Bx_interp);hold on
set(p,'alphadata',~isnan(Bx_interp))
% caxis([min(TMI_interp(:))*1.1 max(TMI_interp(:))])
colormap(h,jet);
% c = colorbar(h,'SouthOutside');
% c.Position(4) = c.Position(4) * 0.5;
% c.Position(2) = c.Position(2) * 0.99;
% set(get(c,'title'),'string','$\mathbf{B}_x\;(nT)$', 'interpreter', 'latex','FontSize',12);
axis equal
% set(h,'XAxisLocation', 'top');
set(h,'YAxisLocation', 'right');
axis([min(x) max(x) min(y) max(y)])
set(h,'YDir','normal')

caxis([min(Bx_interp(:)) max(Bx_interp(:))]);
%% Add colorbar
axes('Position',[0.675 0.68 .275 .225]);
h = gca;
c = colorbar('SouthOutside');
caxis(h,[min(Bx_interp(:)) max(Bx_interp(:))]);
colormap(h,jet)
set(gca,'Visible','off');
% set(get(c,'title'),'string','$\mathbf{\kappa}\;(SI)$', 'interpreter', 'latex','FontSize',12);
set(get(c,'title'),'string','$\mathbf{B}_x\;(nT)$', 'interpreter', 'latex','FontSize',12);
text(0.5,-0.6,'$(b)$', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center');
% caxis([min(Bx_interp(:)) max(Bx_interp(:))]);
%%    
axes('Position',[0.675 0.36 .275 .355]);
h = gca;
% set(figure(1),'CurrentAxes',h)
p = imagesc(x,y,By_interp);hold on
set(p,'alphadata',~isnan(By_interp))
% caxis([min(TMI_interp(:))*1.1 max(TMI_interp(:))])
colormap(h,jet);
% c = colorbar(h,'SouthOutside');
% c.Position(4) = c.Position(4) * 0.5;
% c.Position(2) = c.Position(2) * 0.99;
% set(get(c,'title'),'string','$\mathbf{B}_y\;(nT)$', 'interpreter', 'latex','FontSize',12);
set(h,'YAxisLocation', 'right');
axis equal
axis([min(x) max(x) min(y) max(y)])
set(h,'YDir','normal')
caxis([min(By_interp(:)) max(By_interp(:))]);
%% Add colorbar
axes('Position',[0.675 0.34 .275 .225]);
h = gca;
c = colorbar('SouthOutside');
caxis(h,[min(By_interp(:)) max(By_interp(:))]);
colormap(h,jet)
set(gca,'Visible','off');
% set(get(c,'title'),'string','$\mathbf{\kappa}\;(SI)$', 'interpreter', 'latex','FontSize',12);
set(get(c,'title'),'string','$\mathbf{B}_y\;(nT)$', 'interpreter', 'latex','FontSize',12);
text(0.5,-0.6,'$(c)$', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center');
% caxis([min(Bx_interp(:)) max(Bx_interp(:))]);

%%
axes('Position',[0.675 0.03 .275 .355]);
h = gca;
% set(figure(1),'CurrentAxes',h)
p = imagesc(x,y,Bz_interp);hold on
set(p,'alphadata',~isnan(Bz_interp))
% caxis([min(TMI_interp(:))*1.1 max(TMI_interp(:))])
colormap(h,jet);
set(h,'YAxisLocation', 'right');
% c = colorbar(h,'SouthOutside');
% c.Position(4) = c.Position(4) * 0.5;
% c.Position(2) = c.Position(2) * 0.95;
% set(get(c,'title'),'string','$\mathbf{B}_z\;(nT)$', 'interpreter', 'latex','FontSize',12);
axis equal
axis([min(x) max(x) min(y) max(y)])
% set(gca,'XTickLabel',[]);
set(gca,'YDir','normal')
caxis([min(Bz_interp(:)) max(Bz_interp(:))]);
%% Add colorbar
axes('Position',[0.675 0.025 .275 .225]);
h = gca;
c = colorbar('SouthOutside');
caxis(h,[min(Bz_interp(:)) max(Bz_interp(:))]);
colormap(h,jet)
set(gca,'Visible','off');
% set(get(c,'title'),'string','$\mathbf{\kappa}\;(SI)$', 'interpreter', 'latex','FontSize',12);
set(get(c,'title'),'string','$\mathbf{B}_z\;(nT)$', 'interpreter', 'latex','FontSize',12);
text(0.5,-0.6,'$(d)$', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center');


%%
axes('Position',[0.075 0.03 .275 .355]);
h = gca;
% set(figure(1),'CurrentAxes',h)
p = imagesc(x,y,TMI_interp);hold on
set(p,'alphadata',~isnan(TMI_interp))
% caxis([min(TMI_interp(:))*1.1 max(TMI_interp(:))])
colormap(h,jet);
% set(h,'YTickLabel',[]);
% c = colorbar(h,'SouthOutside');
% c.Position(4) = c.Position(4) * 0.5;
% c.Position(2) = c.Position(2) * 0.95;
% set(get(c,'title'),'string','$\mathbf{TMI}\;(nT)$', 'interpreter', 'latex','FontSize',12);
axis equal
axis([min(x) max(x) min(y) max(y)])
% set(h,'XTickLabel',[]);
set(h,'YDir','normal')
text(-2,-55,'$(e)$', 'interpreter', 'latex','FontSize',12);
hold on

caxis([min(TMI_interp(:)) max(TMI_interp(:))]);
%% Add colorbar
axes('Position',[0.075 0.025 .275 .225]);
h = gca;
c = colorbar('SouthOutside');
caxis(h,[min(TMI_interp(:)) max(TMI_interp(:))]);
colormap(h,jet)
set(gca,'Visible','off');
% set(get(c,'title'),'string','$\mathbf{\kappa}\;(SI)$', 'interpreter', 'latex','FontSize',12);
set(get(c,'title'),'string','$\mathbf{TMI}\;(nT)$', 'interpreter', 'latex','FontSize',12);
text(0.5,-0.6,'$(f)$', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center');

%%
h = axes('Position',[0.375 0.03 .275 .355]);
p = imagesc(x,y,magB_interp);hold on
set(p,'alphadata',~isnan(magB_interp))
% caxis([min(TMI_interp(:))*1.1 max(TMI_interp(:))])
colormap(h,jet);
set(gca,'YTickLabel',[]);
% colorbar(h,'SouthOutside');
% c = colorbar(h,'SouthOutside');
% c.Position(4) = c.Position(4) * 0.5;
% c.Position(2) = c.Position(2) * 0.95;
% set(get(c,'title'),'string','$\mathbf{|B|}\;(nT)$', 'interpreter', 'latex','FontSize',12);
axis equal
axis([min(x) max(x) min(y) max(y)])
% set(gca,'XTickLabel',[]);
set(gca,'YDir','normal')
caxis([min(magB_interp(:)) max(magB_interp(:))]);
%% Add colorbar
axes('Position',[0.375 0.025 .275 .225]);
h = gca;
c = colorbar('SouthOutside');
caxis(h,[min(magB_interp(:)) max(magB_interp(:))]);
colormap(h,jet)
set(gca,'Visible','off');
set(get(c,'title'),'string','$\mathbf{|B|}\;(nT)$', 'interpreter', 'latex','FontSize',12);
text(0.5,-0.6,'$(e)$', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center');

end