% MagVector_Calc
% Function inputs a magnetization vector model and extract average
% declination and inclination of anomaly

clear all 


addpath '..\.'
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Paul_Lake\Modeling\Inversion';
topofile = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Paul_Lake\Data\CDED\CDED_076_c_d_250k.dat';
location = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Paul_Lake\Data\Thematic\Canada\Canada_pt.dat';

meshfile = 'Mesh_50m.msh';
amp_file = 'Tile_AMI_50m\Merged_Amp_model.amp';
rem_file = 'Tile_AMI_50m\Merged_Ind_model.ind';
ind_file = 'Tile_AMI_50m\Merged_Rem_model.rem';
fld_file = 'Tile_AMI_50m\Merged_M_model.fld';

pipe_file = 'Pipes43101';

tile_file = 'Tiles_75_ALL.dat';

obsfile = 'Obs_Paul_Lake_SUB_1pc_5nT.dat';

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

set(figure, 'Position', [50 25 600 800]); 

%% Plot pipe location

ax4 = axes('Position',[0.10 0.35 .87 .87]);

[H, I, Dazm, D, Obsx, Obsy, Obsz, data, wd_full] = read_MAG3D_obs([work_dir dsep obsfile]);

% xx = min(Obsx):50:max(Obsx);
% yy = min(Obsy):50:max(Obsy);
% [YY,XX] = ndgrid(yy,xx);

% F = scatteredInterpolant(Obsy, Obsx, data ,'linear','none');

% grid_d = F(YY,XX);

% msh = mesh(Xn,Yn,ones(size(Xn)),'FaceColor','none'); hold on
view([0 90]) 
scatter(Obsx, Obsy, 2,data); hold on
caxis([-500 250])
colormap(jet)

axis equal
axis([510000 548000 7158000 7182000])


set(gca,'XTick',[515000 529000 543000])
set(gca,'XTickLabel',[ ])
% xlabel('East (m)')

set(gca,'YTick',[7160000 7170000 7180000])
set(gca,'YTickLabel',['7,160,000';'7,170,000';'7,180,000'])
set(gca,'YTickLabelRotation',90)
ylabel('North (m)')

set(ax4,'Ydir','normal')
text(511000,7181000,'(a)','interpreter','latex','FontSize',12)



tiles = load([work_dir dsep tile_file]);
Xtile = tiles(:,[1 3]);
Ytile = tiles(:,[2 4]);

% Plot tiles
h = patch([Xtile(:,1)';Xtile(:,1)';Xtile(:,2)';Xtile(:,2)'],[Ytile(:,1)';Ytile(:,2)';Ytile(:,2)';Ytile(:,1)'],[0.5 0 0],'LineWidth',1);
alpha(h,0);

%% Add color bar
ax = axes('Position',[0.83 .8 .1 .15]);

% bb = interp1([cvec'],[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,0:1e-4:mmax);

cbar = colorbar(ax,'EastOutside');
colormap(ax,jet);
caxis([-500 250])
% set(cbar,'Ticks',[0 0.1 0.3 0.5 0.7 1])
% set(cbar,'TickLabels',round(cvec*10000)/10000)

set(gca,'Visible','off');
text(6,-.2,'(nT)', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
% text(0.5,1.65,'$(b)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center')


%% Plot amplitude
ax1 = axes('Position',[0.10 -.1 .87 .87]);

m = load([work_dir '\' fld_file]);
mx = reshape(m(:,1),nz,nx,ny);
my = reshape(m(:,2),nz,nx,ny);
mz = reshape(m(:,3),nz,nx,ny);

mx = mx( padT+1:end-padB,padW+1:end-padE,padS+1:end-padN );
my = my( padT+1:end-padB,padW+1:end-padE,padS+1:end-padN );
mz = mz( padT+1:end-padB,padW+1:end-padE,padS+1:end-padN );

mx2D = squeeze(mx(zpanel,:,:));
my2D = squeeze(my(zpanel,:,:));
mz2D = squeeze(mz(zpanel,:,:));

M = (mx2D.^2 + my2D.^2 + mz2D.^2) .^ 0.5;

% Compute cluster on field
m = [mx2D(:)./M(:) my2D(:)./M(:) mz2D(:)./M(:)];

cls_m = kmeans(m(M~=0),5);

cls_model = zeros(size(mx2D));
cls_model(M~=0) = cls_m;
% cls_model = reshape(cls_m , size(mx2D));

mesh(xc,yc,cls_model','FaceColor','interp',...
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

text(511000,7181000,'(b)','interpreter','latex','FontSize',12)

cvec = mmax*[0 .2 0.4 0.6 0.8 1];
aa = [255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255;
bb = interp1([cvec'],aa,0:1e-4:mmax,'linear');
colormap(ax1,aa)
caxis([0 mmax])
grid on


% for ii = 1 : size(pipes,1)
%     
%     scatter3(pipes{ii,5},pipes{ii,6},5,50,[0 0 0])
%     
% end

%% Add color bar
ax = axes('Position',[0.83 .35 .1 .15]);

bb = interp1([cvec'],[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,0:1e-4:mmax);

cbar = colorbar(ax,'EastOutside');
colormap(ax,bb);

set(cbar,'Ticks',[0 .2 0.4 0.6 0.8 1])
set(cbar,'TickLabels',0:5)

set(gca,'Visible','off');
text(6,-.2,'$\kappa_e (SI)$', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
% text(0.5,1.65,'$(b)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center')
