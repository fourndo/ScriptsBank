% MagVector_Calc
% Function inputs a magnetization vector model and extract average
% declination and inclination of anomaly

clear all 


addpath '..\.'
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Paul_Lake\Modeling\Inversion';
topofile = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Paul_Lake\Data\CDED\CDED_076_c_d_250k.dat';
location = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Paul_Lake\Data\Thematic\Canada\Canada_pt.dat';

meshfile = 'Mesh_50m_v2.msh';
amp_file = 'CMI_Tiles\ALL_Tiles_rerun\Merged_Amp_model.amp';
ind_file = 'CMI_Tiles\ALL_Tiles_rerun\Merged_Ind_model.ind';
rem_file = 'CMI_Tiles\ALL_Tiles_rerun\Merged_Rem_model.rem';
fld_file = 'CMI_Tiles\ALL_Tiles_rerun\Merged_M_model.fld';

pipe_file = 'Pipes43101';

tile_file = 'Tiles_75_ALL.dat';

obsfile = 'Obs_Paul_Lake_SUB_5pc_5nT_DETREND.dat';

prefile = 'CMI_Tiles\ALL_Tiles_rerun\Merged_MVI_TMI.pre';
fwrfile = 'CMI_Tiles\ALL_Tiles_rerun\PRED_TMI_ALL.pre';
dsep = '\';


%% Load tiles and deposit name
load([work_dir dsep pipe_file]); pipes = Pipes43101;

%% Load mesh
zpanel = 7;
ypanel = 10;

padE = 4;
padW = 4;

padN = 10;
padS = 10;

padT = 0;
padB = 4;

[xn,yn,zn] = read_UBC_mesh([work_dir dsep meshfile]);
xc = ( xn(2:end) + xn(1:end-1) )/2; nx = length(xc);
yc = ( yn(2:end) + yn(1:end-1) )/2; ny = length(yc);
zc = ( zn(1:end-1) + zn(2:end) )/2; nz = length(zc);

mmax = 0.05;

%% Load model and remove paddings
xc = xc(padW+1:end-padE);
yc = yc(padS+1:end-padN);

%% Plot section

set(figure, 'Position', [50 25 850 850]); 

%% Plot data

[H, I, Dazm, D, Obsx, Obsy, Obsz, data, wd_full] = read_MAG3D_obs([work_dir dsep obsfile]);

xx = min(Obsx):50:max(Obsx);
yy = min(Obsy):50:max(Obsy);
[YY,XX] = ndgrid(yy,xx);

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
F = scatteredInterpolant(Obsy, Obsx, data ,'natural','none');

grid_d = F(YY,XX);

F = scatteredInterpolant(Obsy, Obsx, wd_full ,'natural','none');

grid_wd = F(YY,XX);

% msh = mesh(Xn,Yn,ones(size(Xn)),'FaceColor','none'); hold on
ax1 = axes('Position',[0.35 .65 .4 .4]);
view([0 90]) 
h = imagesc(xx,yy,grid_d);hold on
set(h,'alphadata',~isnan(grid_d) .* flag)
caxis(ax1,[-250 250])
colormap(jet)
xlim([min(Obsx) max(Obsx)]); rangex = max(Obsx) - min(Obsx); 
ylim([min(Obsy) max(Obsy)]); rangey = max(Obsy) - min(Obsy);
set(gca,'Ydir','normal')
axis equal tight
% set(gca,'XTickLabel',[])
ylabel('North (m)');
xlabel('East (m)');
% colorbar(ax1,'EastOutside');
title('$d^{(obs)}$', 'interpreter', 'latex','FontSize',12)
hold on
text(5.45e+5,7.1775e+6,'$(a)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
% scatter(Obsx, Obsy, 1,'k.'); hold on
% view([0 90]) 
% scatter(Obsx, Obsy, 2,data); hold on
% caxis([-250 250])
% colormap(jet)
% 
% axis equal
% axis([510000 548000 7158000 7182000])
% 
% 
% set(gca,'XTick',[515000 529000 543000])
% set(gca,'XTickLabel',[ ])
% % xlabel('East (m)')
% 
% set(gca,'YTick',[7160000 7170000 7180000])
% set(gca,'YTickLabel',['7,160,000';'7,170,000';'7,180,000'])
% set(gca,'YTickLabelRotation',90)
% ylabel('North (m)')
% 
% set(ax4,'Ydir','normal')
% text(511000,7181000,'(a)','interpreter','latex','FontSize',12)



% tiles = load([work_dir dsep tile_file]);
% Xtile = tiles(:,[1 3]);
% Ytile = tiles(:,[2 4]);
% 
% % Plot tiles
% h = patch([Xtile(:,1)';Xtile(:,1)';Xtile(:,2)';Xtile(:,2)'],[Ytile(:,1)';Ytile(:,2)';Ytile(:,2)';Ytile(:,1)'],[0.5 0 0],'LineWidth',1);
% alpha(h,0);

%% Add color bar
ax = axes('Position',[0.48 .425 .1 .2]);

% bb = interp1([cvec'],[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,0:1e-4:mmax);

cbar = colorbar(ax,'EastOutside');
caxis([-250 250])
colormap(jet)
% set(cbar,'Ticks',[0 0.1 0.3 0.5 0.7 1])
% set(cbar,'TickLabels',round(cvec*10000)/10000)

set(gca,'Visible','off');
text(2,-.2,'(nT)', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
% text(0.5,1.65,'$(b)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center')

%% Add color bar
ax = axes('Position',[0.75 .75 .1 .2]);

% bb = interp1([cvec'],[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,0:1e-4:mmax);

cbar = colorbar(ax,'EastOutside');
caxis([-250 250])
colormap(jet)
% set(cbar,'Ticks',[0 0.1 0.3 0.5 0.7 1])
% set(cbar,'TickLabels',round(cvec*10000)/10000)

set(gca,'Visible','off');
text(2,-.2,'(nT)', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
% text(0.5,1.65,'$(b)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center')

%% Plot predicted
% ax1 = axes('Position',[0.10 -.1 .87 .87]);

[H, I, Dazm, D, Obsx, Obsy, Obsz, pred, ~] = read_MAG3D_obs([work_dir dsep prefile]);

F = scatteredInterpolant(Obsy, Obsx, pred,'natural','none');

grid_pre = F(YY,XX);


% msh = mesh(Xn,Yn,ones(size(Xn)),'FaceColor','none'); hold on
ax2 = axes('Position',[0.1 .325 .4 .4]);
view([0 90]) 
h = imagesc(xx,yy,grid_pre);hold on
set(h,'alphadata',~isnan(grid_pre).* flag)
caxis(ax2,[-250 250])
% scatter(Obsx, Obsy, 1,'k.'); hold on
colormap(jet)
% colorbar
xlim([min(Obsx) max(Obsx)]); rangex = max(Obsx) - min(Obsx); 
ylim([min(Obsy) max(Obsy)]); rangey = max(Obsy) - min(Obsy);
set(gca,'Ydir','normal')
axis equal tight
title('$d^{(pre)}$', 'interpreter', 'latex','FontSize',12)
% colorbar(ax2,'EastOutside');
set(gca,'XTickLabel',[])
% xlabel('East (m)');
ylabel('North (m)');
text(5.45e+5,7.1775e+6,'$(b)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')

%% Plot predicted
% ax1 = axes('Position',[0.10 -.1 .87 .87]);

[H, I, Dazm, D, Obsx, Obsy, Obsz, fwrd, ~] = read_MAG3D_obs([work_dir dsep fwrfile]);

F = scatteredInterpolant(Obsy, Obsx, fwrd,'natural','none');

grid_fwr = F(YY,XX);


% msh = mesh(Xn,Yn,ones(size(Xn)),'FaceColor','none'); hold on
ax2 = axes('Position',[0.575 .325 .4 .4]);
view([0 90]) 
h = imagesc(xx,yy,grid_fwr);hold on
set(h,'alphadata',~isnan(grid_fwr).* flag)
caxis(ax2,[-250 250])
% scatter(Obsx, Obsy, 1,'k.'); hold on
colormap(jet)
% colorbar
xlim([min(Obsx) max(Obsx)]); rangex = max(Obsx) - min(Obsx); 
ylim([min(Obsy) max(Obsy)]); rangey = max(Obsy) - min(Obsy);
set(gca,'Ydir','normal')
axis equal tight
title('$d^{(fwr)}$', 'interpreter', 'latex','FontSize',12)
% colorbar(ax2,'EastOutside');
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
% xlabel('East (m)');
text(5.45e+5,7.1775e+6,'$(d)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
%% Plot Residual
ax3 = axes('Position',[0.1 0.00 .4 .4]);
view([0 90]) 

% F = scatteredInterpolant(Obsy, Obsx, (data - pred)./wd_full ,'natural','none');

grid_res = (grid_d - grid_pre)./grid_wd;

h = imagesc(xx,yy,grid_res);hold on
set(h,'alphadata',~isnan(grid_res).* flag)
caxis(ax3,[-5 5])
colormap(jet);
% scatter(Obsx, Obsy, 1,'k.'); hold on
% colorbar
xlim([min(Obsx) max(Obsx)]); rangex = max(Obsx) - min(Obsx); 
ylim([min(Obsy) max(Obsy)]); rangey = max(Obsy) - min(Obsy);
set(gca,'Ydir','normal')
axis equal tight
xlabel('East (m)')
ylabel('North (m)');

% colorbar(ax3,'EastOutside');
title('$\|\mathbf{W_d(d^{(pre)} - d^{(obs)})}\|$', 'interpreter', 'latex','FontSize',12)
hold on
text(5.45e+5,7.1775e+6,'$(c)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')

%% Plot Residual
ax4 = axes('Position',[0.575 0.00 .4 .4]);
view([0 90]) 

% F = scatteredInterpolant(Obsy, Obsx, (data - pred)./wd_full ,'natural','none');

grid_res = (grid_d - grid_fwr)./grid_wd;

h = imagesc(xx,yy,grid_res);hold on
set(h,'alphadata',~isnan(grid_res).* flag)
caxis(ax4,[-5 5])
colormap(jet);
% scatter(Obsx, Obsy, 1,'k.'); hold on
% colorbar
xlim([min(Obsx) max(Obsx)]); rangex = max(Obsx) - min(Obsx); 
ylim([min(Obsy) max(Obsy)]); rangey = max(Obsy) - min(Obsy);
set(gca,'Ydir','normal')
axis equal tight
xlabel('East (m)')
% ylabel('North (m)');
set(gca,'YTickLabel',[])
% colorbar(ax3,'EastOutside');
title('$\|\mathbf{W_d(d^{(fwr)} - d^{(obs)})}\|$', 'interpreter', 'latex','FontSize',12)
hold on
text(5.45e+5,7.1775e+6,'$(e)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
%% Add color bar
ax = axes('Position',[0.48 .1 .1 .2]);

% bb = interp1([cvec'],[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,0:1e-4:mmax);

cbar = colorbar(ax,'EastOutside');
colormap(ax,jet);
caxis(ax,[-5 5])
% set(cbar,'Ticks',[0 0.2 0.4 0.6 0.8 1])
% set(cbar,'TickLabels',round(cvec*10000)/10000)

set(gca,'Visible','off');
text(2,-.1,'$\sigma(nT)$', 'interpreter', 'latex','FontSize',10,'HorizontalAlignment','center')


%% Plot model

% Load magnetization model
set(figure, 'Position', [50 25 1500 400]); 

ax2 = axes('Position',[0.05 0.1 .3 .9]);

m = load([work_dir '\' fld_file]);

Ptmi = [(cosd(-I) * cosd(D)) (cosd(-I) * sind(D)) sind(-I)];

m_ind = (Ptmi * m')';


m_pneg = m_ind; m_pneg(m_pneg<0) = 0;

m_rem = [m(:,1) - Ptmi(1)*m_ind m(:,2) - Ptmi(2)*m_ind m(:,3) - Ptmi(3)*m_ind];

m_ind(m_ind>0) = 0;

m_rem = sqrt(sum(m_rem.^2,2));

%%
m = reshape(abs(m_ind),nz,nx,ny);
m = m( padT+1:end-padB,padW+1:end-padE,padS+1:end-padN );
m2D = squeeze(m(zpanel,:,:));
mesh(xc,yc,m2D','FaceColor','interp',...
   'EdgeColor','none');hold on
view([0 90])

axis equal
axis([510000 548000 7158000 7182000])

set(gca,'XTick',[515000 529000 543000])
set(gca,'XTickLabel',[])


set(gca,'YTick',[7160000 7170000 7180000])
set(gca,'YTickLabel',['7,160,000';'7,170,000';'7,180,000'])
set(gca,'YTickLabelRotation',90)
ylabel('North (m)')

text(511000,7181000,'(a)','interpreter','latex','FontSize',12)

cvec = mmax*[0 0.2 0.4 0.6 0.8 1];
aa = [255 255 255;37 120 238;255 255 0;255 127 0;255 0 0;255 153 200]/255;
bb = interp1([cvec'],aa,0:1e-4:mmax,'linear');
colormap(ax2,bb)
caxis([0 mmax])
grid on

%% Add color bar
ax = axes('Position',[0.25 .6 .1 .3]);

bb = interp1([cvec'],[255 255 255;37 120 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,0:1e-4:mmax);

cbar = colorbar(ax,'EastOutside');
colormap(ax,bb);

set(cbar,'Ticks',[0 0.2 0.4 0.6 0.8 1])
set(cbar,'TickLabels',round(cvec*10000)/10000)

set(gca,'Visible','off');
text(6,-.2,'$\kappa_e (SI)$', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
% text(0.5,1.65,'$(b)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center')

%% Plot ratio
ax2 = axes('Position',[0.365 0.1 .3 .9]);

m_rem = m_rem(:);

% m_rem = load([work_dir '\' rem_file]);
m_rem = reshape(m_rem,nz,nx,ny);
m_rem = m_rem( padT+1:end-padB,padW+1:end-padE,padS+1:end-padN );
% 
m2D_rem = squeeze(m_rem(zpanel,:,:));

% m_ind = load([work_dir '\' ind_file]);
% m_ind = reshape(m_ind,nz,nx,ny);
% m_ind = m_ind( padT+1:end-padB,padW+1:end-padE,padS+1:end-padN );

% m2D_ind = squeeze(m_ind(zpanel,:,:));

m2D_Q = m2D_rem;
% Compute correlation between NRM and kHo
% m2D_Q = ( m2D_rem / max(m2D_rem(:)) ) .* ( m2D_ind / max(m2D_ind(:)) );
% m2D_Q = ( m2D_rem  ) .* ( 1-m2D_ind  );
% m2D_Q(m2D_ind==0 | m2D_ind==-100) = 0;

mesh(xc,yc,m2D_Q','FaceColor','interp',...
   'EdgeColor','none');hold on
view([0 90])

axis equal
axis([510000 548000 7158000 7182000])

set(gca,'XTick',[515000 529000 543000])
set(gca,'XTickLabel',['515,000';'529,000';'543,000'])
xlabel('East (m)')

xlabel('East (m)')

set(gca,'YTickLabel',[])
% set(gca,'YTick',[7160000 7170000 7180000])
% set(gca,'YTickLabel',['7,160,000';'7,170,000';'7,180,000'])
% set(gca,'YTickLabelRotation',90)
% ylabel('North (m)')

text(511000,7181000,'(b)','interpreter','latex','FontSize',12)

mmax = 0.05;
cvec = mmax*[0 0.2 0.4 0.6 0.8 1];
aa = [255 255 255;37 120 238;255 255 0;255 127 0;255 0 0;255 153 200]/255;
bb = interp1([cvec'],aa,0:1e-3:mmax,'linear');
colormap(ax2,bb)
caxis([0 mmax])
grid on



%% Add color bar
ax = axes('Position',[0.565 .58 .1 .3]);

bb = interp1([cvec'],[255 255 255;37 120 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,0:1e-4:mmax);

cbar = colorbar(ax,'EastOutside');
colormap(ax,bb);

set(cbar,'Ticks',[0 0.2 0.4 0.6 0.8 1])
set(cbar,'TickLabels',round(cvec*10000)/10000)

set(gca,'Visible','off');
text(6,-.2,'$\kappa_e (SI)$', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
% text(0.5,1.65,'$(b)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center')

%% Plot ratio
ax3 = axes('Position',[0.68 0.1 .3 .9]);


% m_rem = load([work_dir '\' rem_file]);
m_pneg = reshape(abs(m_pneg),nz,nx,ny);
m_pneg = m_pneg( padT+1:end-padB,padW+1:end-padE,padS+1:end-padN );
% 
m2D_rem = squeeze(m_pneg(zpanel,:,:));

% m_ind = load([work_dir '\' ind_file]);
% m_ind = reshape(m_ind,nz,nx,ny);
% m_ind = m_ind( padT+1:end-padB,padW+1:end-padE,padS+1:end-padN );

% m2D_ind = squeeze(m_ind(zpanel,:,:));

m2D_Q = m2D_rem;
% Compute correlation between NRM and kHo
% m2D_Q = ( m2D_rem / max(m2D_rem(:)) ) .* ( m2D_ind / max(m2D_ind(:)) );
% m2D_Q = ( m2D_rem  ) .* ( 1-m2D_ind  );
% m2D_Q(m2D_ind==0 | m2D_ind==-100) = 0;

mesh(xc,yc,m2D_Q','FaceColor','interp',...
   'EdgeColor','none');hold on
view([0 90])

axis equal
axis([510000 548000 7158000 7182000])

set(gca,'XTick',[515000 529000 543000])
set(gca,'XTickLabel',['515,000';'529,000';'543,000'])
xlabel('East (m)')

% xlabel('East (m)')

set(gca,'YTickLabel',[])
% set(gca,'YTick',[7160000 7170000 7180000])
% set(gca,'YTickLabel',['7,160,000';'7,170,000';'7,180,000'])
% set(gca,'YTickLabelRotation',90)
% ylabel('North (m)')

text(511000,7181000,'(c)','interpreter','latex','FontSize',12)

mmax = 0.05;
cvec = mmax*[0 0.2 0.4 0.6 0.8 1];
aa = [255 255 255;37 120 238;255 255 0;255 127 0;255 0 0;255 153 200]/255;
bb = interp1([cvec'],aa,0:1e-3:mmax,'linear');
colormap(ax3,bb)
caxis([0 mmax])
grid on



%% Add color bar
ax = axes('Position',[0.88 .58 .1 .3]);

bb = interp1([cvec'],[255 255 255;37 120 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,0:1e-4:mmax);

cbar = colorbar(ax,'EastOutside');
colormap(ax,bb);

set(cbar,'Ticks',[0 0.2 0.4 0.6 0.8 1])
set(cbar,'TickLabels',round(cvec*10000)/10000)

set(gca,'Visible','off');
text(6,-.2,'$\kappa_e (SI)$', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
% text(0.5,1.65,'$(b)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center')

%% Plot amplitude on a wide plot
set(figure, 'Position', [50 25 800 600]); 
ax1 = axes('Position',[0.15 0.1 .8 .8]);

m = load([work_dir '\' amp_file]);
m = reshape(m,nz,nx,ny);
m = m( padT+1:end-padB,padW+1:end-padE,padS+1:end-padN );

m2D = squeeze(m(zpanel,:,:));

mesh(xc,yc,m2D','FaceColor','interp',...
   'EdgeColor','none');hold on
view([0 90])

axis equal
axis([510000 548000 7158000 7182000])


set(gca,'XTick',[515000 529000 543000])
set(gca,'XTickLabel',['515,000';'529,000';'543,000'])
xlabel('East (m)')

set(gca,'YTick',[7160000 7170000 7180000])
set(gca,'YTickLabel',['7,160,000';'7,170,000';'7,180,000'])
set(gca,'YTickLabelRotation',90)
ylabel('North (m)')

% text(504000,7170000,'(b)','interpreter','latex','FontSize',12)
mmax = 0.05;
cvec = mmax*[0 0.1 0.3 0.5 0.7 1];
aa = [255 255 255;37 120 238;255 255 0;255 127 0;255 0 0;255 153 200]/255;
bb = interp1([cvec'],aa,0:1e-4:mmax,'linear');
colormap(ax1,bb)
caxis([0 mmax])
grid on

%%
for ii = 1 : size(pipes,1)
    
    scatter3(pipes{ii,5},pipes{ii,6},5,50,[0 0 0])
    
end