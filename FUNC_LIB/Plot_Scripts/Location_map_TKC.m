% MagVector_Calc
% Function inputs a magnetization vector model and extract average
% declination and inclination of anomaly

clear all 
close all

addpath '..\.'
work_dir = 'C:\Users\DominiqueFournier\ownCloud\Research\Paul_Lake\Data';
topofile = 'C:\Users\DominiqueFournier\ownCloud\Research\Paul_Lake\Data\CDED\CDED_076_c_d_250k.dat';
location = 'C:\Users\DominiqueFournier\ownCloud\Research\Paul_Lake\Data\Thematic\Canada\Canada_pt.dat';

pipe_file = 'Ekati\PipesUTM';

obsfile = '..\Modeling\Obs_Paul_Lake_SUB_1pc_10nT_DETREND.dat';

dsep = '\';


%% Load tiles and deposit name
load([work_dir dsep pipe_file]);

%% Load topography and generate contours
% topo = read_UBC_topo([work_dir dsep topofile]);
% topo = load([topofile]);
% 
% % Create topo surface
% x = round(min(topo(:,1)) : 100 : max(topo(:,1)));
% y = round(min(topo(:,2)) : 100 : max(topo(:,2)));
% 
% [X,Y] = ndgrid(x,y);
% 
% T = scatteredInterpolant(topo(:,1),topo(:,2),topo(:,3));
% tgrid = T(X,Y);

load('C:\Users\DominiqueFournier\ownCloud\Research\Paul_Lake\Data\CDED\CDED_076_c_d_250k')
set(figure, 'Position', [50 25 850 650]); 
ax1 = axes('Position',[0.1 .1 .8 .8]);

% load('C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Paul_Lake\Data\CDED\CDED_076_c_d_250k') 
imagesc(X(:,1),Y(1,:),tgrid');hold on
colormap([107 220 255;255 255 255]/255);
caxis([408 415])
contour(X,Y,tgrid,[450:50:550],'Color',[0.7 0.7 0.7]); 
contour(X,Y,tgrid,[416 416],'k'); 
set(gca,'Ydir','normal')

axis([480000 580000 7133000 7205000])

text(518000,7149000,'Lac de Gras','interpreter','latex')
text(486000,7197700,'Exeter Lake','interpreter','latex')

set(gca,'XTick',[490000 520000 560000])
set(gca,'XTickLabel',['490,000';'520,000';'560,000'])
xlabel('East (m)')

set(gca,'YTick',[7150000 7170000 7190000])
set(gca,'YTickLabel',['7,150,000';'7,170,000';'7,190,000'])
set(gca,'YTickLabelRotation',90)
ylabel('North (m)')

grid on
%% Plot pipe location

plot(pipes.X,pipes.Y,'k^','MarkerSize',5,'MarkerFaceColor','k')
% for ii = 1 : size(pipes.NAME,1)
%     
%     text(pipes.X(ii),pipes.Y(ii),pipes.NAME{ii},'interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top')
%     
% end
text(518140,7176600,'Ekati','interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top','BackgroundColor','w')
plot(518140,7176600,'ro','MarkerFaceColor','r','MarkerSize',5)

text(534900,7152500,'Diavik','interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top','BackgroundColor','w')
plot(534900,7152500,'ro','MarkerFaceColor','r','MarkerSize',5)

text(557500,7133500,'DO-27','interpreter','latex','HorizontalAlignment','left','VerticalAlignment','bottom','BackgroundColor','w')
plot(557500,7133500,'ro','MarkerFaceColor','r','MarkerSize',5)

%% Add inset map of Canada
ax2 = axes('Position',[0.69 .69 .20 .20]);

cdn = load([location]);

plot(cdn(:,1),cdn(:,2),'k.','MarkerSize',0.001); hold on

axis([min(cdn(:,1)) max(cdn(:,1)) min(cdn(:,2)) max(cdn(:,2))])
set(ax2,'XTickLabel',[],'YTickLabel',[])

plot(mean(X(:,1)),mean(Y(1,:)),'ro','MarkerFaceColor','r','MarkerSize',5)

%% Create figure with data and known pipes
set(figure, 'Position', [50 25 850 650]); 

ax1 = axes('Position',[0.1 .1 .8 .8]);

% contour(X,Y,tgrid,[450:50:550],'Color',[0.7 0.7 0.7]); 
contour(X,Y,tgrid,[450:50:550],'Color',[0.7 0.7 0.7]); hold on
contour(X,Y,tgrid,[416 416],'k'); 
set(gca,'Ydir','normal')

axis([510000 550000 7155000 7185000])

text(518000,7149000,'Lac de Gras','interpreter','latex')
text(486000,7197700,'Exeter Lake','interpreter','latex')

set(gca,'XTick',[490000 520000 560000])
set(gca,'XTickLabel',['490,000';'520,000';'560,000'])
xlabel('East (m)')

set(gca,'YTick',[7150000 7170000 7190000])
set(gca,'YTickLabel',['7,150,000';'7,170,000';'7,190,000'])
set(gca,'YTickLabelRotation',90)
ylabel('North (m)')
grid on

[H, BI, BD, MI, MD, Obsx, Obsy, Obsz, data, wd_full] = read_MAG3D_obs([work_dir dsep obsfile]);

xx = min(Obsx):50:max(Obsx);
yy = min(Obsy):50:max(Obsy);
[YY,XX] = ndgrid(yy,xx);


% msh = mesh(Xn,Yn,ones(size(Xn)),'FaceColor','none'); hold on
view([0 90]) 
% imagesc(xx,yy,grid_d);hold on
scatter(Obsx, Obsy, 2,data); hold on
caxis([-250 250])
colormap(jet)

%% Plot pipe location
plot(pipes.X,pipes.Y,'k^','MarkerSize',5,'MarkerFaceColor','k')
for ii = 1 : size(pipes.NAME,1)
    
    text(pipes.X(ii),pipes.Y(ii),pipes.NAME{ii},'interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top')
    
end
text(518140,7176600,'Ekati','interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top','BackgroundColor','w')
plot(518140,7176600,'ro','MarkerFaceColor','r','MarkerSize',5)

text(534900,7152500,'Diavik','interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top','BackgroundColor','w')
plot(534900,7152500,'ro','MarkerFaceColor','r','MarkerSize',5)

%%
ax2 = axes('Position',[0.69 .69 .20 .20]);

cdn = load([location]);

plot(cdn(:,1),cdn(:,2),'k.','MarkerSize',0.001); hold on

axis([min(cdn(:,1)) max(cdn(:,1)) min(cdn(:,2)) max(cdn(:,2))])
set(ax2,'XTickLabel',[],'YTickLabel',[])

plot(mean(X(:,1)),mean(Y(1,:)),'ro','MarkerFaceColor','r','MarkerSize',5)