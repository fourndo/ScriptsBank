% MagVector_Calc
% Function inputs a magnetization vector model and extract average
% declination and inclination of anomaly

clear all 
close all


addpath '..\.'
addpath 'C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\m_map';

work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Documents\MATLAB';

pole_file = 'Cande_1995_Normal_Polarity_TABLE.txt';
dsep = '\';

%% Vine and Matthew survey lines
lt1 = (45+[17 19]/60);
ln1 = -[28+27/60 11+29/60];

lt2 = [(5+30/60) 10+10/60];
ln2 = [61+57/60 66+27/60];

lt3 = [49.28,49.28];
ln3 = [-135,-123.1207];  
%% Load polarity file
% pole = load([work_dir dsep pole_file]);
% 
% 
% set(figure, 'Position', [50 25 600 600]); 
% 
% ax1 = axes('Position',[0.14 .15 .8 .1]);
% for ii = 1 : size(pole,1)
%     
%     patch([pole(ii,1) pole(ii,1) pole(ii,2) pole(ii,2)],[0 1 1 0],[0 0 0])
%     hold on
% end
% 
% axis([min(pole(:)) 90 0 1])
% set(gca,'Xdir','reverse','YTickLabel',[])
% text(45,-0.6,'Time (Ma)','HorizontalAlignment','center', 'interpreter', 'latex','FontSize',14)
% text(45,1.2,'Polarity','HorizontalAlignment','center', 'interpreter', 'latex','FontSize',14)
% text(102,0.5,'(b)','HorizontalAlignment','center', 'interpreter', 'latex','FontSize',14)
% 
% ax1 = axes('Position',[0.9 .26 0.025 0.025]);
% patch([1 1 0 0],[0 1 1 0],[0 0 0])
% text(-.2,0.2,'Normal', 'HorizontalAlignment','right','interpreter', 'latex','FontSize',14)
% set(gca,'Visible','off');
% axis equal
% 
% ax1 = axes('Position',[0.15 .26 0.025 0.025]);
% patch([1 1 0 0],[0 1 1 0],[1 1 1])
% text(1.1,0.2,'Reversed', 'interpreter', 'latex','FontSize',14)
% set(gca,'Visible','off');
% axis equal
%% Load world map data
load([work_dir '\Emag.mat']);

%%
lon = unique(Emag(:,1));
lat = unique(Emag(:,2));
% Dazm = squeeze(Emag(3,:,:));
% I = squeeze(Emag(4,:,:));
A = reshape(Emag(:,3),length(lon),length(lat))';
A(A==-99999) = nan;
%%
% limlon = [6000 7000];
% limlat = [2000 3000];

dd = 2;

slat = lat(1:dd:end);
slon = lon(1:dd:end);
sA = A(1:dd:end,1:dd:end);

% lon = round((-90+2/60:2/60:90-2/60)*1000)/1000;
% lat = round((-180+2/60:2/60:180-2/60)*1000)/1000;
% A = zeros(length(lon),length(lat));
% 
% for ii = 1 : size(Emag,1)
% 
%         A(Emag(ii,1)==lon,Emag(ii,2)==lat) = Emag(ii,3);
%         
% end
% A = A  /max(A(:));
% B = abs(A * 50000 + 20000);

% Compute XYZ coordinates from lon-lat
% x = 1.25*cosd(lon).*cosd(lat);
% y = 1.25*sind(lon).*cosd(lat);
% z = 1.25*sind(lat);
% 
% xx = cosd(lon).*cosd(lat);
% yy = sind(lon).*cosd(lat);
% zz = sind(lat);

% Rotate coordinates for mapping
% xyz_temp = [x(:) y(:) z(:)];
% 
% for ii = 1 : size(xyz_temp,1)
%     
% Rx = [1 0 0;
%       0 cosd(rotI) -sind(rotI);
%       0 sind(rotI) cosd(rotI)];
% 
% Rz = [cosd(rotD) -sind(rotD) 0;
%       sind(rotD) cosd(rotD) 0;
%       0 0 1];
%   
% xyz_temp(ii,:) = [Rz * Rx * xyz_temp(ii,:)']'; 
% 
% end
% 
% x(:) = xyz_temp(:,1);
% y(:) = xyz_temp(:,2);
% z(:) = xyz_temp(:,3);


set(figure, 'Position', [50 25 800 800]); 
axes('Position',[0.1 .2 .8 .8]);

cmin = min(sA(sA~=-99999));
cmax = max(sA(sA~=-99999));

ccode = interp1([-205 -100 0 200],[1 1 1;1 0 0;1 1 1;0 0 1],-200:0.1:200);
colormap(jet)
% m_proj('hammer-aitoff','clongitude',-150);
m_proj('lambert','long',[-140 -123],'lat',[46 55]);
% Rather than rearranging the data so its limits match the
% plot I just draw it twice (you can see the join at 180W
% because of the quirks of flat pcolor) (Note that 
% all the global projections have 360 deg ambiguities)
aa=m_pcolor(slon,slat,sA);shading flat;colormap(ccode);
caxis([-500 500])
hold on;
m_pcolor(slon-360,slat,sA);shading flat;colormap(ccode);

m_line(ln1,lt1,'color','k','linewi',3);

m_coast('line','LineWidth',2,'Color',[0.25 0.25 0.25]);
m_grid('xaxis','bottom');

text(-3.5,0,'(a)','HorizontalAlignment','center', 'interpreter', 'latex','FontSize',14)

% set(figure, 'Position', [50 25 800 500]); 
% axes('Position',[0 .05 1 .8]);
% m_proj('lambert','long',[110 160],'lat',[-46 -10]);
% % Rather than rearranging the data so its limits match the
% % plot I just draw it twice (you can see the join at 180W
% % because of the quirks of flat pcolor) (Note that 
% % all the global projections have 360 deg ambiguities)
% m_pcolor(slon,slat,sA);shading flat;colormap(ccode);
% caxis([-500 500])
% hold on;
% m_pcolor(slon-360,slat,sA);shading flat;colormap(ccode);
m_line(-123.1207,49.28,'marker','square','color','k','markersize',16,'linewi',3);
% m_text(-124.1207,50,'Vancouver','FontSize',14,'HorizontalAlignment','right','Color','w')
m_line(ln3,lt3,'color','k','linewi',2,'linestyle','--');
xlabel('Longitude $W^{\circ}$', 'interpreter', 'latex')
ylabel('Latitude $N^{\circ}$', 'interpreter', 'latex')
title('Residual Total Field - West Coast, Canada')
axes('Position',[0.075 .06 .85 .15]);

indx = lon>ln3(1) & lon<ln3(2);

ln3 = [-135,-123.1207];
plot(lon(indx),A(530,indx),'k','LineWidth',2); hold on
plot(-123.1207,49.28,'rs','MarkerSize',16)
text(-123,50,'Vancouver','FontSize',12,'HorizontalAlignment','left','VerticalAlignment','top')
ylim([-500 500])
ylabel('nT')
xlabel('Longitude - $49.28^{\circ}$N', 'interpreter', 'latex')
grid on

% 
% m_coast('line','LineWidth',2,'Color','k');
% m_grid('xaxis','bottom');

%% Plot color bar
ax3 = axes('Position',[0.575 .50 .4 .4]);
cbar = colorbar(ax3,'EastOutside');
set(gca,'Visible','off');
set(cbar,'Ticks',[0.1 0.5 0.9])
set(cbar,'TickLabels',[-500 0 500])
text(1.25,1.,'(nT)', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
text(0.5,1.3,'Anomalous Total Field','HorizontalAlignment','center', 'interpreter', 'latex','FontSize',14)