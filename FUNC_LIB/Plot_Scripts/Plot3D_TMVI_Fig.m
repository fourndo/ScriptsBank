% Create figure through 3D model
clear all



addpath ..\..\FUNC_LIB;

%% Input Files
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Nut_Cracker\Tiled_AMI\Tile1';
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Nut_Cracker\Tiled_CMI\Tile1';

meshfile = 'Tile1.msh';

obsfile = 'Tile_data.dat';
model_true = '..\..\Effec_sus_20mGrid.sus';


% mag_vecfile = '\l2l2\Mvec_TMVI_iter_.fld';
% predfile = '\l2l2\TMVI_iter_.pre';

mag_vecfile = '\Tile1_MVI.fld';
predfile = '\Tile1_MVI.pre';

% mag_vecfile = '..\..\magvec.fld';
% predfile = '..\Obs_REM_GRID_TMI.obs';

% mag_vecfile = '..\magvec.fld';
obsfile = '..\..\Obs_RAW_REM_GRID_TMI.obs';


zpanel = 5;
ypanel = 8;

padE = 4;
padW = 4;

padN = 4;
padS = 4;

padT = 0;
padB = 4;

iso_cut_surf = 0.01;
iso_cut_vec = 0.01;

mmax = 0.075;

cam_ang = [45 45];

vscale = 1;

% Define cutting planes
nvec(1,1:3) = [0 -0.2 1]; xo = 1000; yo = 850  ; zo = 1420;
nvec(2,1:3) = [0 1 0]; xo(2) = 1000; yo(2) = 730  ; zo(2) = 1440;
nvec(3,1:3) = [1 0 0]; xo(3) = 1000; yo(3) = 850  ; zo(3) = 1440;

nvec = spdiags( 1./ sqrt(sum(nvec.^2,2)) , 0 ,3 ,3) * nvec;

cut_name = ['A','B','C'];

% Color scheme

%% Load in model and plot
set(figure, 'Position', [1000 50 850 650]); 

% Load data
[H, I, Dazm, D, obsx, obsy, obsz, d, wd] = read_MAG3D_obs([work_dir '\' obsfile]);
% H=1000;

[xn,yn,zn] = read_UBC_mesh([work_dir '\' meshfile]);

% Move to local coordinates
obsx = obsx - xn(1);
obsy = obsy - yn(1);
% obsz = obsz - zn(1);

xn = xn - xn(1);
yn = yn - yn(1);
% zn = zn - zn(1);

dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);
xx = (xn(2:end) + xn(1:end-1))/2;   xx = xx(padW+1:end-padE);
yy = (yn(2:end) + yn(1:end-1))/2;   yy = yy(padS+1:end-padN);
zz = (zn(2:end) + zn(1:end-1))/2;   zz = zz(padT+1:end-padB);

[XX,ZZ,YY] = meshgrid(xx,zz,yy); 

% Load mag true
m_true = load([work_dir '\' model_true]);
m_true = reshape(m_true,nz,nx,ny);


% Load magnetization vector
mag_model = load([work_dir '\' mag_vecfile]);

% Load inverted model
% m = load([work_dir '\..\..\Effec_sus_20mGrid.sus']);
m = sqrt(sum(mag_model.^2,2));
m = reshape(m,nz,nx,ny);


mcell = size(mag_model,1);
%% Create model magnetization vectors
% Azimuth and dip of magnitization
if size(mag_model,2)==2
    
    mag_xyz = azmdip_2_xyz( mag_model(:,1)+180 , mag_model(:,2) );
    

else
    
    mag_xyz = mag_model;
%     mag_xyz(mag_xyz(:,3)~=0,1:2)=0;
%     mag_xyz(mag_xyz(:,3)~=0,3)=m(mag_xyz(:,3)~=0);
%     M = [spdiags(m_true(:).*mag_model(:,1),0,mcell,mcell);spdiags(m_true(:).*mag_model(:,2),0,mcell,mcell);spdiags(m_true(:).*mag_model(:,3),0,mcell,mcell)];

   
end

mx = reshape(mag_xyz(:,1),nz,nx,ny).*H;
my = reshape(mag_xyz(:,2),nz,nx,ny).*H;
mz = reshape(-mag_xyz(:,3),nz,nx,ny).*H;

mx = mx( padT+1:end-padB,padW+1:end-padE,padS+1:end-padN );
my = my( padT+1:end-padB,padW+1:end-padE,padS+1:end-padN );
mz = mz( padT+1:end-padB,padW+1:end-padE,padS+1:end-padN );
    
m_true = m_true( padT+1:end-padB,padW+1:end-padE,padS+1:end-padN );
m = m( padT+1:end-padB,padW+1:end-padE,padS+1:end-padN );

for ii = 1 : 2
% Define cutting plane

if ii == 1
    
[XX2D,YY2D] = meshgrid(min(xx):10:max(xx),min(yy):10:max(yy));
ZZ2D = -( nvec(ii,1)*(XX2D - xo(ii)) + nvec(ii,2)*(YY2D  - yo(ii)) - zo(ii) ) / nvec(ii,3) ;


elseif ii == 2
[XX2D,ZZ2D] = meshgrid(min(xx):10:max(xx),min(zz):10:max(zz));
YY2D = ones(size(XX2D)) * yo(ii);

else
[YY2D,ZZ2D] = meshgrid(min(yy):10:max(yy),min(zz):10:max(zz));
XX2D = ones(size(YY2D)) * xo(ii);

end
%%




% Interpolate on cutting plane
% m2D = get_model_top(m1,nx,ny,nz,-100);

% cvec = [min(m(:)) prctile(m(m>0),[10 30 60]) max(m(:))];
% bb = interp1([cvec'],[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0]/255,sort(m(:)));

F = TriScatteredInterp(XX(:),YY(:),ZZ(:),m(:),'nearest'); 
m2D = F(XX2D,YY2D,ZZ2D);

F = TriScatteredInterp(XX(:),YY(:),ZZ(:),mx(:),'nearest'); 
mx2D = F(XX2D,YY2D,ZZ2D);

F = TriScatteredInterp(XX(:),YY(:),ZZ(:),my(:),'nearest'); 
my2D = F(XX2D,YY2D,ZZ2D);

F = TriScatteredInterp(XX(:),YY(:),ZZ(:),mz(:),'nearest'); 
mz2D = F(XX2D,YY2D,ZZ2D);



% Interpolate on cutting plane
% m2D = get_model_top(m1,nx,ny,nz,-100);

F = TriScatteredInterp(XX(:),YY(:),ZZ(:),m_true(:),'nearest'); 
m2D_true = F(XX2D,YY2D,ZZ2D);

% temp = (temp');

if ii == 1
    
    ax1 = axes('Position',[0.59 .45 .35 .35]);
%     h = imagesc(xx,yy,m2D); hold on
    
    scatter(min(XX2D(:)), max(YY2D(:)),'k.'); hold on
    text(min(XX2D(:)), max(YY2D(:)),['\textbf{' cut_name(ii) '}'],'interpreter', 'latex','FontSize',14,'VerticalAlignment','top');
    scatter(max(XX2D(:)), min(YY2D(:)),'k.')
    text(max(XX2D(:)), min(YY2D(:)),['\textbf{' cut_name(ii) '"}'],'interpreter', 'latex','FontSize',14,'HorizontalAlignment','right');
%     title('$Model$','interpreter', 'latex','FontSize',14);

%     qq = quiver(xx,yy,mx2D,my2D,'LineWidth',1,'Color','k','MaxHeadSize',1);
    h = surf(XX2D,YY2D,zeros(size(XX2D)),m2D,'EdgeColor','none'); hold on
    alpha(0.5)
    view([0 90])
    set(gca,'YDir','normal')
    plot_vec(ax1,XX2D,YY2D,[],mx2D,my2D,[],m2D,iso_cut_vec,0,0.5,1.25,3)

    axis([min(XX2D(:)) max(XX2D(:)) min(YY2D(:)) max(YY2D(:))])
    
    set(gca,'XTickLabel',[],'Box','on')
    ylabel('$y$', 'interpreter', 'latex','FontSize',14)
    
    axis equal
    grid on

    set(gca,'YDir','normal')
    hold on
     contour(ax1,XX2D,YY2D,m2D_true,[0.005 0.01],'k');  
     
elseif ii == 2
    
    ax2 = axes('Position',[0.59 .15 .35 .35]);
    
%     h = imagesc(xx,zz,m2D); hold on

    h = surf(XX2D,ZZ2D,zeros(size(XX2D)),m2D,'EdgeColor','none'); hold on
    alpha(0.5)
    view([0 90])
    set(gca,'YDir','normal')
    plot_vec(ax2,XX2D,ZZ2D,[],mx2D,mz2D,[],m2D,iso_cut_vec,0,0.5,1.25,3)
        
    scatter(min(XX2D(:)), max(ZZ2D(:)),'k.'); hold on
    text(min(XX2D(:)), max(ZZ2D(:)),['\textbf{' cut_name(ii) '}'],'interpreter', 'latex','FontSize',14,'VerticalAlignment','top');
    scatter(max(XX2D(:)), min(ZZ2D(:)),'k.')
    text(max(XX2D(:)), min(ZZ2D(:)),['\textbf{' cut_name(ii) '"}'],'interpreter', 'latex','FontSize',14,'HorizontalAlignment','right');
    
    caxis([0 mmax])
    axis([min(XX2D(:)) max(XX2D(:)) min(ZZ2D(:)) max(ZZ2D(:))])
    
    xlabel('$x$', 'interpreter', 'latex','FontSize',14)
    ylabel('$z$', 'interpreter', 'latex','FontSize',14)
    
    axis equal
    grid on

    set(gca,'YDir','normal','Box','on')
    hold on

    
    
%     plot([XX2D(:)-mx2D(:)*100, ZZ2D(:)-mz2D(:)*100],[XX2D(:)+mx2D(:)*100, ZZ2D(:)+mz2D(:)*100]);
    contour(ax2,XX2D,ZZ2D,m2D_true,[0.005 0.01],'k');
%     qq = quiver(xx,zz,mx2D,mz2D,'LineWidth',1,'Color','k','MaxHeadSize',1);
    
end


    colormap(jet);
    caxis([0 mmax])

% set(h,'CData',reshape(bb,size(XX2D,1),size(XX2D,2),3));%colormap(ax,bb);
% colormap(ax,[255 255 255;77 190 238;255 127 0;255 0 0]/255);

% if mmax < max(m2D(:))
%     mmax = max(m2D(:));
%     caxis([0 mmax]);
% end


set(get(gca,'YLabel'),'Rotation',360);

% if ii == 1
%     
%     scatter(obsx,obsy,1,'k')
% 
% end
set(gca, 'YAxisLocation', 'right')

% text(min(xx)-dx(1)*2, mean(yn),'$(b)$','interpreter', 'latex','FontSize',14);



% Add DO27 marker
% scatter3(2000,3035,450,50,'*','r')
% text(1600,3135,'$DO27$','interpreter', 'latex','FontSize',12);
% 
% scatter3(2350,3600,450,50,'*','r')
% text(1950,3700,'$DO18$','interpreter', 'latex','FontSize',12);

end

cvec = mmax*[0 0.2 0.4 0.6 0.8 1.05];
bb = interp1([cvec'],[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,0:1e-4:mmax,'linear');
colormap(ax1,bb);
colormap(ax2,bb);
caxis(ax1,[0 mmax])
caxis(ax2,[0 mmax])
%% Add color bar
%%
ax = axes('Position',[0.59 -0.25 .35 .35]);

% bb = interp1([cvec'],[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0]/255,min(cvec):1e-4:max(cvec));

cbar = colorbar(ax,'NorthOutside');
% caxis([0 mmax]);

colormap(ax,bb);
caxis([0 mmax])
% set(cbar,'Ticks',[0 0.25 0.5 0.75 1])
% set(cbar,'TickLabels',round(cvec*10000)/10000)

set(gca,'Visible','off');
text(0.5,1.5,'$\kappa_{e}$', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
text(0.5,1.65,'$(b)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center')

%%
% ax3 = axes('Position',[0.06 .25 .525 .525]);
% plot_vec(ax3,XX(:),YY(:),ZZ(:),mx(:),my(:),mz(:),m(:),1e-3,40,1.25,4)
% colormap(ax3,bb);
% caxis(ax3,[0 mmax]);
ax = axes('Position',[0.06 .25 .525 .525]);

aa = isosurface(XX,YY,ZZ,m,iso_cut_surf);

% Create color map for depth
bb = interp1([max(aa.vertices(:,3));median(aa.vertices(:,3));min(aa.vertices(:,3))],[215 48 39;255 255 191;69 117 180]/255,aa.vertices(:,3));

view(cam_ang);
p = patch('Faces',aa.faces,'Vertices',aa.vertices,'FaceColor','b','LineStyle','none','FaceColor', 'flat','FaceVertexCData',bb);
ylabel('$y$', 'interpreter', 'latex','FontSize',14)
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
zlabel('$z$', 'interpreter', 'latex','FontSize',14)
hold on
axis equal
set(gca,'DataAspectRatio',[1 1 vscale])
camlight HEADLIGHT 
lighting phong
set(gca,'YDir','normal')
grid on
xlim([min(xx) max(xx)])
ylim([min(yy) max(yy)])
scatter3(obsx,obsy,obsz,1,'k')

% legend(['Iso-surface: ' num2str(iso_cut) ' SI'],'Location','WestOutside')

for ii = 1 : 2
% Define cutting plane

    if ii == 1

    [XX2D,YY2D] = meshgrid(xx,yy);
    ZZ2D = -( nvec(ii,1)*(XX2D - xo(ii)) + nvec(ii,2)*(YY2D  - yo(ii)) - zo(ii) ) / nvec(ii,3) ;

    elseif ii == 2
    [XX2D,ZZ2D] = meshgrid(xx,zz);
    YY2D = ones(size(XX2D)) * yo(ii);

    else
    [YY2D,ZZ2D] = meshgrid(yy,zz);
    XX2D = ones(size(YY2D)) * xo(ii);

    end
    
% Add cutting plane
ss = surf(XX2D,YY2D,ZZ2D); hold on
    if ii == 1
        set(ss,'LineStyle','none','FaceColor',[.25 .25 .25])
    elseif ii == 2
        set(ss,'LineStyle','none','FaceColor',[.25 .25 .25])
    else
        set(ss,'LineStyle','none','FaceColor',[1 1 0])
    end

alpha(ss,.4)

scatter3(min(XX2D(:)), max(YY2D(:)), max(ZZ2D(:)),'k.')
text(min(XX2D(:)), max(YY2D(:)), max(ZZ2D(:)),['\textbf{' cut_name(ii) '}'],'interpreter', 'latex','FontSize',14,'VerticalAlignment','top');
scatter3(max(XX2D(:)), min(YY2D(:)), min(ZZ2D(:)),'k.')
text(max(XX2D(:)), min(YY2D(:)), min(ZZ2D(:)),['\textbf{' cut_name(ii) '"}'],'interpreter', 'latex','FontSize',14,'HorizontalAlignment','right');
end


%% Add colorbars
ax = axes('Position',[0.15 -.25 .35 .35]);

bb = interp1([max(aa.vertices(:,3));median(aa.vertices(:,3));min(aa.vertices(:,3))],[215 48 39;255 255 191;69 117 180]/255,sort(aa.vertices(:,3)));
colormap(ax,bb);
cbar = colorbar('NorthOutside');
set(cbar,'Ticks',[0 0.5 1])
set(cbar,'TickLabels',round([min(aa.vertices(:,3)) median(aa.vertices(:,3)) max(aa.vertices(:,3))]))
set(gca,'Visible','off');
text(0.5,1.5,'$Elevation (m)$', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
text(0.5,1.65,'$(a)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center')




%% Plot Data
% Load data
[~, ~, ~, ~, ~, ~, ~, data, ~] = read_MAG3D_obs([work_dir '\' obsfile]);

% Load predicted lp
[~, ~, ~, ~, ~, ~, ~, dpre, ~] = read_MAG3D_obs([work_dir '\' predfile]);


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

res = (data - dpre)./wd;

data_interp     = griddata(obsx, obsy, data,X,Y,'linear'); 
% data_interp(isnan(data_interp)) = min(data_interp(:))*2;

dpre_interp     = griddata(obsx, obsy, dpre,X,Y,'linear');

res_interp        = griddata(obsx, obsy, res ,X,Y,'linear'); 

set(figure, 'Position', [50 0 850 650]); 

a = axes('Position',[0.075 .65 .28 .28]);
h =imagesc(x,y,data_interp);hold on
set(h,'alphadata',~isnan(data_interp))
caxis([min(data_interp(:)) max(data_interp(:))]);
colormap(a,jet);
scatter(obsx,obsy,2,'k.')
set(gca,'YDir','normal')
% xlabel('\bfEasting (m)')
ylabel('$y$', 'interpreter', 'latex','FontSize',14)
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
set(get(gca,'YLabel'),'Rotation',360);
axis([min(x) max(x) min(y) max(y)])
% set(gca,'XTickLabel',[])
grid on
axis equal
title('$d^{Fwr}$', 'interpreter', 'latex','FontSize',14)
% text(min(xx)-dx*20, mean(yy),'$(a)$', 'interpreter', 'latex','FontSize',14)

a = axes('Position',[0.375 .65 .28 .28]);
h =imagesc(x,y,dpre_interp);hold on
set(h,'alphadata',~isnan(data_interp))
caxis([min(data_interp(:)) max(data_interp(:))]);
colormap(a,jet);
scatter(obsx,obsy,2,'k.')
set(gca,'YDir','normal')
% xlabel('\bfEasting (m)')
% ylabel('$y$', 'interpreter', 'latex','FontSize',14)
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
set(get(gca,'YLabel'),'Rotation',360);
axis([min(x) max(x) min(y) max(y)])
set(gca,'YTickLabel',[])
grid on
axis equal
title('$d^{Obs}$', 'interpreter', 'latex','FontSize',14)
% text(400, mean(yy),'$(a)$', 'interpreter', 'latex','FontSize',14)



a = axes('Position',[0.675 .65 .28 .28]);
h =imagesc(x,y,res_interp);hold on
set(h,'alphadata',~isnan(res_interp))
caxis([-5 5])
colormap(a,jet);
scatter(obsx,obsy,2,'k.')
set(gca,'YDir','normal')
set(gca,'YTickLabel',[])
% ylabel('$y$', 'interpreter', 'latex','FontSize',14)
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
set(get(gca,'YLabel'),'Rotation',360);
axis([min(x) max(x) min(y) max(y)])
grid on
axis equal
title('$Normalized\;Residual$', 'interpreter', 'latex','FontSize',14)
% text(400, mean(yy),'$(b)$', 'interpreter', 'latex','FontSize',14)



%% Add colorbars
ax = axes('Position',[0.675 .25 .28 .28]);
cbar = colorbar('NorthOutside');
colormap(ax,jet);
set(cbar,'Ticks',[0 1])
set(cbar,'TickLabels',round([-5 5]))
set(gca,'Visible','off');
text(0.5,1.5,'${\mathbf{W}_d(d^{fwd} - d^{obs})} \;(nT)$', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
text(0.5,1.75,'$(c)$', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
%% Add colorbars
ax = axes('Position',[0.225 0.25 .28 .28]);
cbar = colorbar('NorthOutside');
colormap(ax,jet);
set(cbar,'Ticks',[0 1])
set(cbar,'TickLabels',round([min(data_interp(:)) max(data_interp(:))]))
set(gca,'Visible','off');
text(0.5,1.5,'$TMI\;(nT)$', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
text(-.11,1.75,'$(a)$', 'interpreter', 'latex','FontSize',12)
text(1.0,1.75,'$(b)$', 'interpreter', 'latex','FontSize',12)