% Create dipole field from sphere and plot vectors

clear all 
close all

addpath ..\arrow3

work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\UBC\ResearchDomFournier\General\Figures';

%% Create grid points

xx = -1.25:0.02:1.25; nx = length(xx);
yy = -1.25:0.02:1.25; ny = length(yy);
zz = 0.75;

[XX,YY] = ndgrid(xx,yy);
ZZ = ones(size(XX))*zz;

mcell = nx * ny;

vl = 1e-3;

%% Formula for dipole field in cartesian coordinates
I = 0;
D = 90 ;

dpr = 0.2; % Dipole radius

Ptmi = [(sind(I) * cosd(D)) (sind(I) * sind(D)) cosd(I)];

% bx = @(p,tt,pp,r) 3*p./r.^3.* sin(pp+dphi).*cos(pp+dphi).*cos(tt+dtheta);
% by = @(p,tt,pp,r) 3*p./r.^3.* sin(pp+dphi).*cos(pp+dphi).*sin(tt+dtheta);
% bz = @(p,tt,pp,r) 3*p./r.^3.* (cos(pp+dphi).^2 - 1/3);

% Coordinate-free
b = @(p,uu,rr,r,N) p*spdiags(1./r(:).^3,0,N,N) * ( 3 * spdiags(sum(uu.*rr,2),0,N,N) * rr - uu );

% Compute distances and angles for dipole at the origin
% r2D = sqrt( XX.^2 + YY.^2 );
lrl = sqrt( XX.^2 + YY.^2 + ZZ.^2 );

r = spdiags(1./lrl(:) , 0 , mcell ,mcell ) * [XX(:) YY(:) ZZ(:)];

Rx = [1 0 0;
      0 cosd(I) -sind(I);
      0 sind(I) cosd(I)];

Rz = [cosd(D) -sind(D) 0;
      sind(D) cosd(D) 0;
      0 0 1];
  
u = Rz * Rx * [0;0;1];
u = kron(ones(mcell,1),u');

B = b(-15,u,r,lrl,mcell);
% Compute phi angle
% phi = asin(r2D./r3D+1e-8);

% Compute theta angle
% theta = atan(YY./(XX + 1e-8));

% theta(XX < 0 ) = theta(XX < 0 ) + pi;
% theta(XX >= 0 & YY <= 0 ) = theta(XX >= 0 & YY <= 0 ) + 2*pi;
% Project field onto inducing

Btmi = -Ptmi * B';

% Compute field
Bx = reshape( B(:,1) , nx ,ny );
By = reshape( B(:,2) , nx ,ny );
Bz = reshape( B(:,3) , nx ,ny );
Btmi = reshape( Btmi , nx ,ny );


% Plot frame with dipole sphere
dphi = -pi/2:0.05:pi/2;
dtta = -pi:0.05:pi;

[dP,dT] = ndgrid(dphi,dtta);

% Compute XYZ coordinates from lon-lat
x = dpr*cos(dT).*cos(dP); nx =length(dphi);
y = dpr*sin(dT).*cos(dP); ny =length(dtta);
z = dpr*sin(dP);

mcell =  nx*ny;

lrl = sqrt(x.^2 + y.^2 + z.^2);
r = spdiags(1./lrl(:) , 0 , mcell ,mcell ) * [x(:) y(:) z(:)];

u = Rz * Rx * [0;0;1];
u = kron(ones(mcell,1),u');

b_sph = b(-3,u,r,lrl,mcell);
b_sph = reshape(sum(b_sph.^2,2).^0.5,nx,ny);
% reshape

set(figure, 'Position', [50 25 800 600]);
axes('Position',[-0.1 -.1 1.25 1.25]);
surf(x,y,z,b_sph*-100,'EdgeColor','none','FaceColor','k','FaceAlpha',0.25);
% alpha(0.25)
view([-130 10])
shading INTERP
axis equal
set(gca,'Visible','off');
hold on

% surf(x-0.5,y-0.5,z-0.25,b_sph/10,'EdgeColor','none','FaceColor','k','FaceAlpha',0.1);

% Create surface with b-field on surface
h = surface(XX,YY,ZZ,Btmi,'EdgeColor','none','FaceAlpha',0.75);

caxis([min(Btmi(:))*2 max(Btmi(:))]) 
colormap([0 0 0;jet])

axis([-1.5 1.5 -1.5 1.5 -1.5 1.5])
% Plot axis
% plot3([0 0],[0 0],[-0.5 1.5],'k')
% plot3([-1 1],[0 0],[0 0],'k')
% plot3([0 0],[-1 1],[0 0],'k')

% text(0,0,-0.5,'Z', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','top')

% text(-1,0,0,'X', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','left','VerticalAlignment','middle')

% text(0,-1,0,'Y', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','left','VerticalAlignment','middle')
% Add observation location
% quiver3(-0.85,-0.35,zz,0.85,0.35,-zz,'k','LineWidth',2,'MaxHeadSize',0.25,'AutoScale','off')
% scatter3(-0.85,-0.35,zz+0.01,100,'r^','MarkerFaceColor','w','LineWidth',2)
% text(-0.85,-0.35,zz+0.1,'$\vec b$', 'interpreter', 'latex','Color','w','FontSize',16,'HorizontalAlignment','left','VerticalAlignment','middle')
% text(-0.5,-0.25,zz/2,'$\vec r$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','left','VerticalAlignment','middle')
% text(0,1,0,'Obs', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','right','VerticalAlignment','middle')
%% Seed the sphere with monopoles and trace the field lines
% Plot frame with dipole sphere
% dphi = ([0 20 180]+I)*pi/180 ;
% dtta = [-45 0 45]*pi/180;

dphi = ([-105 -90 -75 0:15:180])*pi/180 ;
dtta = [45]*pi/180;

[dP,dT] = ndgrid(dphi,dtta);

% Compute XYZ coordinates from lon-lat
x = (dpr+0.01)*cos(dT(:)).*cos(dP(:)); nx =length(x);
y = (dpr+0.01)*sin(dT(:)).*cos(dP(:)); ny =length(y);
z = (dpr+0.01)*sin(dP(:));

mcell =  nx;

u = Rz * Rx * [0;0;1]; u = u';
dl = 0.025;

for ii = 1 : mcell
    
    x_in = x(ii);
    y_in = y(ii);
    z_in = z(ii);
    
    lrl = sqrt(x_in.^2 + y_in.^2 + z_in.^2);
    jj = 0;
    
%     plot3(x_in,y_in,z_in, 'ro','MarkerSize',5)
    
    ll = 0;
    flag = 1;
    
    while abs(lrl) > 0.2 && ll < 2 
        jj = jj + 1;
        
        lrl = sqrt(x_in(jj).^2 + y_in(jj).^2 + z_in(jj).^2);
        r = [x_in(jj) y_in(jj) z_in(jj)] / lrl;

        b_temp = b(-1,u,r,lrl,1);
        b_temp = b_temp/norm(b_temp);
        
        x_in = [x_in x_in(end) - flag * dl*b_temp(1)];
        y_in = [y_in y_in(end) - flag * dl*b_temp(2)];
        z_in = [z_in z_in(end) - flag * dl*b_temp(3)];
        
        lrl = sqrt(x_in(end).^2 + y_in(end).^2 + z_in(end).^2);
        
        if jj==1 &&  lrl < 0.2
            
            flag = -1;
            x_in(end) = x_in(jj) - flag * dl*b_temp(1);
            y_in(end) = y_in(jj) - flag * dl*b_temp(2);
            z_in(end) = z_in(jj) - flag * dl*b_temp(3);
            lrl = sqrt(x_in(end).^2 + y_in(end).^2 + z_in(end).^2);
        end
            
        ll = ll + dl;
               
%         if (z_in(end) > zz && z_in(end-1) < zz) || jj == 1
%             
%             quiver3(x_in(end-1),y_in(end-1),z_in(end-1),0.1*b_temp(1),0.1*b_temp(2),0.1*b_temp(3),'k','LineWidth',2,'MaxHeadSize',5);
%         end
        
    end
    
    
    
    if length(x_in)>1 %&& ll > 0.1
    temp = round(length(x_in) / 2);
    dl_temp = 2*[x_in(temp) - x_in(temp+1) y_in(temp) - y_in(temp+1) z_in(temp) - z_in(temp+1)];
%     arrow3([x_in(end-1) y_in(end-1) z_in(end-1)],[x_in(end-1)+vl*b_temp(1) y_in(end-1)+vl*b_temp(2) z_in(end-1)+vl*b_temp(3)],'k',1,2,'cone');

    plot3(x_in,y_in,z_in, 'k--','LineWidth',0.25)
        if flag==-1
            
            arrow3([x_in(temp) y_in(temp) z_in(temp)],[x_in(temp+1)-x_in(temp) y_in(temp+1)-y_in(temp) z_in(temp+1)-z_in(temp)],'k',0.75,1.5,'cone');

        else
            
            arrow3([x_in(temp) y_in(temp) z_in(temp)],[x_in(temp-1)-x_in(temp) y_in(temp-1)-y_in(temp) z_in(temp-1)-z_in(temp)],'k',0.75,1.5,'cone');

        end
    end
    
end

%% Plot amplitude of field for a range of inclination
% Formula for dipole field in cartesian coordinates
% set(figure(2), 'Position', [50 25 1000 800]);
% set(figure(3), 'Position', [50 25 1000 800]);
set(figure(2), 'Position', [50 25 1000 500]);
xx = -1.25:0.05:1.25; nx = length(xx);
yy = -1.25:0.05:1.25; ny = length(yy);
zz = 0.75;

[XX,YY] = ndgrid(xx,yy);
ZZ = ones(size(XX))*zz;

mcell = nx * ny;

peak = zeros(6,1);

b_mat = [];
% figure(2)

for ii = 1 : 24
I = 0 + 15*(ii-1);
D = 90 ;
% bx = @(p,tt,pp,r) 3*p./r.^3.* sin(pp+dphi).*cos(pp+dphi).*cos(tt+dtheta);
% by = @(p,tt,pp,r) 3*p./r.^3.* sin(pp+dphi).*cos(pp+dphi).*sin(tt+dtheta);
% bz = @(p,tt,pp,r) 3*p./r.^3.* (cos(pp+dphi).^2 - 1/3);

% Coordinate-free

% Compute distances and angles for dipole at the origin
% r2D = sqrt( XX.^2 + YY.^2 );

lrl = sqrt( XX.^2 + YY.^2 + ZZ.^2 );

r = spdiags(1./lrl(:) , 0 , mcell ,mcell ) * [XX(:) YY(:) ZZ(:)];

Rx = [1 0 0;
      0 cosd(I) -sind(I);
      0 sind(I) cosd(I)];

Rz = [cosd(D) -sind(D) 0;
      sind(D) cosd(D) 0;
      0 0 1];
  
u = Rz * Rx * [0;0;1];
u = kron(ones(mcell,1),u');

B = b(-15,u,r,lrl,mcell);


lBl = ( sum(B.^2,2) ).^0.5;

peak(ii) = max(lBl);



lBl = reshape(lBl,length(xx),length(yy));

lbl = lBl(:,26);
b_mat = [b_mat lbl];

Bx = reshape(B(:,3),nx,ny);
bx = Bx(:,26);



% figure(3)
% % plot(xx,lBl(26,:));hold on
% maxB = max(lbl);
% plot(xx,lbl,'b');hold on
% text(xx(lbl==maxB),maxB,['I=' num2str(I)],'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',16,'BackgroundColor','w');
% 
% figure(4)
% maxBx = max(bx);
% plot(xx,bx,'b');hold on
% temp = xx(bx==maxBx);
% text(temp(1),maxBx,['I=' num2str(I)],'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',16,'BackgroundColor','w');

% figure(2)
axes('Position',[0.5 0.5 0.35 0.35]);
aa = surf(xx,yy,lBl);
axis tight
zlim([0 80]);
view([45 60]);
drawnow
title('|B|')
grid off

% figure
axes('Position',[0.1 .5 0.35 0.35]);
surf(xx+max(xx),yy,reshape(B(:,1),nx,ny));
axis tight
view([45 60]);
zlim([-75 75]);
title('Bx')
grid off

axes('Position',[0.5 .05 0.35 0.35]);
surf(xx+max(xx),yy,reshape(B(:,2),nx,ny));
axis tight
view([45 60]);
zlim([-75 75]);
title('By')
grid off

axes('Position',[0.1 .05 0.35 0.35]);
surf(xx+max(xx),yy,reshape(B(:,3),nx,ny));
axis tight
view([45 60]);
zlim([-75 75]);
title('Bz')
grid off

% axes('Position',[0.5 .0 0.5 0.5]);
% dpl = [cosd(I) -sind(I);sind(I) cosd(I)]*[0;1];

% arrow3([-dpl(1) -dpl(2) 0],[dpl(1) dpl(2) 0],'k',0.1,.5)
% axis([-5 5 -5 5 -5 5])
% axis equal square
% grid off
% view([45 45]);
% set(gca,'XTickLabel',[])
% subplot(4,7,ii+14)
% imagesc(xx,yy,reshape(B(:,2),nx,ny));
% axis equal tight
% 
% subplot(4,7,ii+21)
% imagesc(xx,yy,reshape(B(:,3),nx,ny));
% axis equal tight


frame = getframe(figure(2));
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
if ii == 1;
  imwrite(imind,cm,[work_dir '\field.gif'],'gif', 'Loopcount',inf,'DelayTime',0.1);
else
  imwrite(imind,cm,[work_dir '\field.gif'],'gif','WriteMode','append','DelayTime',0.1);
end
      
% close(figure(2))
end

figure(3)
plot(xx,max(b_mat'),'r');hold on
axis tight
grid on

figure(4)
axis tight
grid on