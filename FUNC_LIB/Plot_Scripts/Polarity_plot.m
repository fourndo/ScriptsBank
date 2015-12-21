% MagVector_Calc
% Function inputs a magnetization vector model and extract average
% declination and inclination of anomaly

clear all 
close all


addpath '..\.'
addpath 'C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\m_map';

work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Paul_Lake\Data';

pole_file = 'Cande_1995_Normal_Polarity_TABLE.txt';
dsep = '\';


%% Load polarity file
pole = load([work_dir dsep pole_file]);


set(figure, 'Position', [50 25 600 600]); 

ax1 = axes('Position',[0.14 .15 .8 .1]);
for ii = 1 : size(pole,1)
    
    patch([pole(ii,1) pole(ii,1) pole(ii,2) pole(ii,2)],[0 1 1 0],[0 0 0])
    hold on
end

axis([min(pole(:)) 90 0 1])
set(gca,'Xdir','reverse','YTickLabel',[])
text(45,-0.6,'Time (Ma)','HorizontalAlignment','center', 'interpreter', 'latex','FontSize',14)
text(45,1.2,'Polarity','HorizontalAlignment','center', 'interpreter', 'latex','FontSize',14)
text(102,0.5,'(b)','HorizontalAlignment','center', 'interpreter', 'latex','FontSize',14)

ax1 = axes('Position',[0.9 .26 0.025 0.025]);
patch([1 1 0 0],[0 1 1 0],[0 0 0])
text(-.2,0.2,'Normal', 'HorizontalAlignment','right','interpreter', 'latex','FontSize',14)
set(gca,'Visible','off');
axis equal

ax1 = axes('Position',[0.15 .26 0.025 0.025]);
patch([1 1 0 0],[0 1 1 0],[1 1 1])
text(1.1,0.2,'Reversed', 'interpreter', 'latex','FontSize',14)
set(gca,'Visible','off');
axis equal
%% Load world map data
load('LonLatDI');
rotD = 0;
rotI = -90;

lon = squeeze(Emag(1,:,:)*180/pi);
lat = squeeze(Emag(2,:,:)*180/pi);
Dazm = squeeze(Emag(3,:,:));
I = squeeze(Emag(4,:,:));
A = squeeze(Emag(5,:,:));

A = A  /max(A(:));
B = abs(A * 50000 + 20000);

% Compute XYZ coordinates from lon-lat
x = 1.25*cosd(lon).*cosd(lat);
y = 1.25*sind(lon).*cosd(lat);
z = 1.25*sind(lat);

xx = cosd(lon).*cosd(lat);
yy = sind(lon).*cosd(lat);
zz = sind(lat);

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



ax2 = axes('Position',[0.125 .30 .85 .85]);

ccode = interp1([min(A(:)) 0 max(A(:))],[1 0 0;1 1 1;0 0 1],min(A(:)):0.1:max(A(:)));
colormap(ccode)
m_proj('hammer-aitoff','clongitude',-150);

% Rather than rearranging the data so its limits match the
% plot I just draw it twice (you can see the join at 180W
% because of the quirks of flat pcolor) (Note that 
% all the global projections have 360 deg ambiguities)
m_pcolor(lon,lat,A);shading flat;colormap(ccode);
caxis([-1.05 0.95])
hold on;
m_pcolor(lon-360,lat,A);shading flat;colormap(ccode);

m_coast('line','LineWidth',0.5,'Color','k');
m_grid('xaxis','middle');

text(-3.5,0,'(a)','HorizontalAlignment','center', 'interpreter', 'latex','FontSize',14)

%% Plot color bar
ax3 = axes('Position',[0.35 .350 .4 .4]);
cbar = colorbar(ax3,'SouthOutside');
set(gca,'Visible','off');
set(cbar,'Ticks',[0.1 0.5 0.9])
set(cbar,'TickLabels',[50000 30000 -50000])
text(1.1,-.3,'(nT)', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')

%% Plot inclination and declination
set(figure, 'Position', [50 25 800 800]); 

ax2 = axes('Position',[0.2 .45 .65 .65]);

ccode = interp1([min(Dazm(:)) 0 max(Dazm(:))],[1 0 0;1 1 1;0 0 1],min(Dazm(:)):0.1:max(Dazm(:)));
colormap(ccode)
m_proj('hammer-aitoff','clongitude',-150);

% Rather than rearranging the data so its limits match the
% plot I just draw it twice (you can see the join at 180W
% because of the quirks of flat pcolor) (Note that 
% all the global projections have 360 deg ambiguities)
m_pcolor(lon,lat,Dazm);shading flat;colormap(ccode);
caxis([-120 120])
hold on;
m_pcolor(lon-360,lat,Dazm);shading flat;colormap(ccode);

m_coast('line','LineWidth',0.5,'Color','k');
m_grid('xaxis','middle');

text(-3.5,0,'(a)','HorizontalAlignment','center', 'interpreter', 'latex','FontSize',14)

%% Plot color bar
ax3 = axes('Position',[0.325 .5 .4 .4]);
cbar = colorbar(ax3,'SouthOutside');
set(gca,'Visible','off');
set(cbar,'Ticks',[0.1 0.5 0.9])
set(cbar,'TickLabels',[-120 0 120])
text(1.1,-.2,'Decl ($^\circ$)', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')

%% Plot inclination and declination
ax1 = axes('Position',[0.2 -.00 .65 .65]);

ccode = interp1([min(I(:)) 0 max(I(:))],[1 0 0;1 1 1;0 0 1],min(I(:)):0.1:max(I(:)));
colormap(ccode)
m_proj('hammer-aitoff','clongitude',-150);

% Rather than rearranging the data so its limits match the
% plot I just draw it twice (you can see the join at 180W
% because of the quirks of flat pcolor) (Note that 
% all the global projections have 360 deg ambiguities)
m_pcolor(lon,lat,I);shading flat;colormap(ccode);
caxis([-90 90])
hold on;
m_pcolor(lon-360,lat,I);shading flat;colormap(ccode);

m_coast('line','LineWidth',0.5,'Color','k');
m_grid('xaxis','middle');

text(-3.5,0,'(a)','HorizontalAlignment','center', 'interpreter', 'latex','FontSize',14)

%% Plot color bar
ax3 = axes('Position',[0.325 0.05 .4 .4]);
cbar = colorbar(ax3,'SouthOutside');
set(gca,'Visible','off');
set(cbar,'Ticks',[0.1 0.5 0.9])
set(cbar,'TickLabels',[-90 0 90])
text(1.1,-.2,'Incl ($^\circ$)', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
%% Test to plot vectors in 3D
% % m_grid('linest','-','xticklabels',[],'yticklabels',[]);
% hold on
% 
% % Compute vector field from D and I

D = mod(450-Dazm,360);
I = I - 90;

vecxy = [B(:).*cosd(D(:)) B(:).*sind(D(:))];
   

vecxyz = [B(:).*cosd(D(:)).*cosd(I(:)) B(:).*sind(D(:)).*cosd(I(:)) B(:).*sind(I(:))];

xyz_temp = [x(:) y(:) z(:)];

sub = zeros(size(x));
sub(1:16:256,1:16:end,1:16:end) =1;
sub = sub==1;

figure;
imagesc(lon(:,1),lat(1,:),I'); hold on
set(gca,'Ydir','Normal')
quiver3(lon(sub),lat(sub),lat(sub).^0*5,vecxyz(sub,1),vecxyz(sub,2),vecxyz(sub,3),'k')

% Rotation matrix
for ii = 1 : length(lon(:))
    
Rx = [1 0 0;
      0 cosd(lat(ii)) -sind(lat(ii));
      0 sind(lat(ii)) cosd(lat(ii))];

Rz = [cosd(lon(ii)) -sind(lon(ii)) 0;
      sind(lon(ii)) cosd(lon(ii)) 0;
      0 0 1];
  
vecxyz(ii,:) = [Rx * Rz * vecxyz(ii,:)']'; 

% Rx = [1 0 0;
%       0 cosd(rotI) -sind(rotI);
%       0 sind(rotI) cosd(rotI)];
% 
% Rz = [cosd(rotD) -sind(rotD) 0;
%       sind(rotD) cosd(rotD) 0;
%       0 0 1];
%   
% xyz_temp(ii,:) = [Rz * Rx * xyz_temp(ii,:)']'; 

end

% x(:) = xyz_temp(:,1);
% y(:) = xyz_temp(:,2);
% z(:) = xyz_temp(:,3);

% cvec = [-180:90:180];
% ccode = interp1(cvec',[77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,D(:));

figure;
h = surf(xx,yy,zz,A,'EdgeColor','none','FaceAlpha',0.75);
hold on
angle([45 45])

q = quiver3(x(sub),y(sub),z(sub),vecxyz(sub,1),vecxyz(sub,2),vecxyz(sub,3),1,'Color','k');
axis equal


%%

% m_proj('ortho','lat',10','long',-110');
% m_coast('patch','k');
% % m_grid('linest','-','xticklabels',[],'yticklabels',[]);
% hold on
% figure
% m_proj('ortho','lat',48','long',-123');
% m_coast('patch','r');
% figure;
% h = surf(x,y,z,D,'EdgeColor','none');
% axis equal
% colormap(jet)
% 
% hold on
% 
% % cvec = [-180:90:180];
% % ccode = interp1(cvec',[77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,D(:));
% sub = zeros(size(x));
% sub(1:6:end,1:6:end,1:6:end) =1;
% sub = sub==1;
% q = quiver3(x(sub)*1.1,y(sub)*1.1,z(sub)*1.1,vecxyz(sub,1),vecxyz(sub,2),vecxyz(sub,3),2);
