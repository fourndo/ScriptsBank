% MagVector_Calc
% Function inputs a magnetization vector model and extract average
% declination and inclination of anomaly

clear all 
close all

addpath '..\.'
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Paul_Lake\Modeling\Inversion';
topofile = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Paul_Lake\Data\CDED\CDED_076_c_d_250k.dat';
location = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Paul_Lake\Data\Thematic\Canada\Canada_pt.dat';

meshfile = 'Mesh_75m.msh';
amp_file = 'Tile_75m_ALL\ALL_Tiles_75m_l0l2\Merged_Amp_model.amp';
rem_file = 'Tile_75m_ALL\ALL_Tiles_75m_l0l2\Merged_Ind_model.ind';
ind_file = 'Tile_75m_ALL\ALL_Tiles_75m_l0l2\Merged_Rem_model.rem';

pipe_file = 'Pipes43101';

obsfile = 'Obs_Paul_Lake_SUB_1pc_5nT.dat';

dsep = '\';


%% Load tiles and deposit name
load([work_dir dsep pipe_file]); pipes = Pipes43101;

%% Load mesh
zpanel = 5;
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

set(figure, 'Position', [50 25 800 1000]); 

%% Plot amplitude
ax1 = axes('Position',[0.11 -.1 .6 .6]);

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

text(504500,7170000,'(c)','interpreter','latex','FontSize',12)

cvec = mmax*[0 0.1 0.3 0.5 0.7 1];
aa = [255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255;
bb = interp1([cvec'],aa,0:1e-4:mmax,'linear');
colormap(ax1,bb)
caxis([0 mmax])
grid on

%% Plot pipe location


for ii = 1 : size(pipes,1)
    
    scatter3(pipes{ii,5},pipes{ii,6},5,50,[0 0 0])
    
end

%% Add color bar
ax = axes('Position',[0.6 .175 .1 .15]);

bb = interp1([cvec'],[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,0:1e-4:mmax);

cbar = colorbar(ax,'EastOutside');
colormap(ax,bb);

set(cbar,'Ticks',[0 0.1 0.3 0.5 0.7 1])
set(cbar,'TickLabels',round(cvec*10000)/10000)

set(gca,'Visible','off');
text(3,-.1,'$\kappa_e (SI)$', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
% text(0.5,1.65,'$(b)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center')

%% Plot remanent
ax2 = axes('Position',[0.11 0.22 .6 .6]);

m = load([work_dir '\' ind_file]);
m = reshape(m,nz,nx,ny);
m = m( padT+1:end-padB,padW+1:end-padE,padS+1:end-padN );

m2D = squeeze(m(zpanel,:,:));

mesh(xc,yc,m2D','FaceColor','interp',...
   'EdgeColor','none');hold on
view([0 90])

axis equal
axis([510000 548000 7158000 7182000])


set(gca,'XTick',[515000 529000 543000])
set(gca,'XTickLabel',[])
% xlabel('East (m)')

set(gca,'YTick',[7160000 7170000 7180000])
set(gca,'YTickLabel',['7,160,000';'7,170,000';'7,180,000'])
set(gca,'YTickLabelRotation',90)
ylabel('North (m)')

text(504500,7170000,'(b)','interpreter','latex','FontSize',12)

cvec = mmax*[0 0.2 0.4 0.6 0.8 1];
aa = [255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255;
bb = interp1([cvec'],aa,0:1e-4:mmax,'linear');
colormap(ax2,bb)
caxis([0 mmax])
grid on

%% Add color bar
ax = axes('Position',[0.6 .495 .1 .15]);

bb = interp1([cvec'],[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,0:1e-4:mmax);

cbar = colorbar(ax,'EastOutside');
colormap(ax,bb);

set(cbar,'Ticks',[0 0.2 0.4 0.6 0.8 1])
set(cbar,'TickLabels',round(cvec*10000)/10000)

set(gca,'Visible','off');
text(3,-.1,'$\kappa_e (SI)$', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
% text(0.5,1.65,'$(b)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center')

%% Plot induced
ax3 = axes('Position',[0.11 0.54 .6 .6]);

m = load([work_dir '\' rem_file]);
m = reshape(m,nz,nx,ny);
m = m( padT+1:end-padB,padW+1:end-padE,padS+1:end-padN );

m2D = squeeze(m(zpanel,:,:));

mesh(xc,yc,m2D','FaceColor','interp',...
   'EdgeColor','none');hold on
view([0 90])

axis equal
axis([510000 548000 7158000 7182000])


set(gca,'XTick',[515000 529000 543000])
set(gca,'XTickLabel',[])
% xlabel('East (m)')

set(gca,'YTick',[7160000 7170000 7180000])
set(gca,'YTickLabel',['7,160,000';'7,170,000';'7,180,000'])
set(gca,'YTickLabelRotation',90)
ylabel('North (m)')

text(504500,7170000,'(a)','interpreter','latex','FontSize',12)

cvec = mmax*[0 0.2 0.4 0.6 0.8 1];
aa = [255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255;
bb = interp1([cvec'],aa,0:1e-4:mmax,'linear');
colormap(ax3,bb)
caxis([0 mmax])
grid on

%% Add color bar
ax = axes('Position',[0.6 .815 .1 .15]);

bb = interp1([cvec'],[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,0:1e-4:mmax);

cbar = colorbar(ax,'EastOutside');
colormap(ax,bb);

set(cbar,'Ticks',[0 0.2 0.4 0.6 0.8 1])
set(cbar,'TickLabels',round(cvec*10000)/10000)

set(gca,'Visible','off');
text(3,-.1,'$\kappa_e (SI)$', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
% text(0.5,1.65,'$(b)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center')

%% Plot Data
[H, avg_I, Dazm, avg_D, Obsx, Obsy, Obsz, data, wd_full] = read_MAG3D_obs([work_dir dsep obsfile]);

% scatter3(Obsx, Obsy , Obsz,1,'k');



% text(518140,7176600,'Ekati','interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top','BackgroundColor','w')
% plot(518140,7176600,'ro','MarkerFaceColor','r','MarkerSize',5)

