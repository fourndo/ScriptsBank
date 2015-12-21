% Create dipole field from sphere and plot vectors

clear all 
close all

addpath ..\arrow3
addpath ..\gridfitdir

work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\UBC\ResearchDomFournier\General\Figures';

% Spherical coordinates e-field
e = @(q,rr,r,N) q*spdiags(1./r(:).^2,0,N,N) * ( rr );

%% Create q source
dphi = linspace(-pi/2,pi/2,10);
dtta = linspace(-pi,pi,10);

[dP,dT] = ndgrid(dphi,dtta);

xx_q = 0.2* cos(dT).*cos(dP); 
yy_q = 0.2* sin(dT).*cos(dP);
zz_q = 0.2* sin(dP);

lrl = sqrt( xx_q.^2 + yy_q.^2 + zz_q.^2 );

nx =size(xx_q,1);
ny =size(xx_q,2);
mcell = nx*ny;

r = spdiags(1./lrl(:) , 0 , mcell ,mcell ) * [xx_q(:) yy_q(:) zz_q(:)];

E_q = e(0.1,r,lrl,mcell);



Ex = reshape( E_q(:,1) , nx ,ny );
Ey = reshape( E_q(:,2) , nx ,ny );
Ez = reshape( E_q(:,3) , nx ,ny );
lEl_q = sqrt( Ex.^2 + Ey.^2 + Ez.^2 );

%% Create random surface and plot e-field from q
    dphi = linspace(-pi/2,pi/2,4);
    dtta = linspace(-pi,pi,4);

    [dp,dt] = ndgrid(dphi,dtta);
    randr = randn(size(dp));
    randr = randr/max(abs(randr(:)))*2.5;
    
for jj = 1 : 75


    if jj<=10
        
        dr = 3.25 - randr;
        
    elseif jj < 20 && jj > 10
        
        dr = 3.25 - randr*(1-(jj-11)/11);

    else
        
        dr = (3.25 - randr*(0))*(1-(jj-20)/55);
        
    end
    dr(end,:) = dr(end); dr(1,:) = dr(1);
    dr(:,end) = dr(:,1);

    % Compute XYZ coordinates from lon-lat

    [xx,yy] = ndgrid(1:size(dp,1),1:size(dp,2));

    [XX,YY] = ndgrid(linspace(1,size(dp,1),100),linspace(1,size(dp,2),100));

    dR = griddata(xx,yy,dr,XX,YY,'natural');

    % F = scatteredInterpolant(xx(:),yy(:),zz(:),'linear','linear');

    dphi = linspace(-pi/2,pi/2,100);
    dtta = linspace(-pi,pi,100);

    [dP,dT] = ndgrid(dphi,dtta);

    % dr = (5 - randn(size(dP))/2);
    % Compute XYZ coordinates from lon-lat
    xx = dR.*cos(dT).*cos(dP); 
    yy = dR.*sin(dT).*cos(dP);
    zz = dR .* sin(dP);



    nx =size(xx,1);
    ny =size(xx,2);


    mcell = nx*ny;


    %% Formula for dipole field in cartesian coordinates
    % bx = @(p,tt,pp,r) 3*p./r.^3.* sin(pp+dphi).*cos(pp+dphi).*cos(tt+dtheta);
    % by = @(p,tt,pp,r) 3*p./r.^3.* sin(pp+dphi).*cos(pp+dphi).*sin(tt+dtheta);
    % bz = @(p,tt,pp,r) 3*p./r.^3.* (cos(pp+dphi).^2 - 1/3);



    % Compute distances and angles for dipole at the origin
    % r2D = sqrt( XX.^2 + YY.^2 );
    lrl = sqrt( xx.^2 + yy.^2 + zz.^2 );

    r = spdiags(1./lrl(:) , 0 , mcell ,mcell ) * [xx(:) yy(:) zz(:)];

    E = e(15,r,lrl,mcell);
    % Compute phi angle
    % phi = asin(r2D./r3D+1e-8);

    % Compute theta angle
    % theta = atan(YY./(XX + 1e-8));

    % theta(XX < 0 ) = theta(XX < 0 ) + pi;
    % theta(XX >= 0 & YY <= 0 ) = theta(XX >= 0 & YY <= 0 ) + 2*pi;
    % Project field onto inducing


    % Compute field
    Ex = reshape( E(:,1) , nx ,ny );
    Ey = reshape( E(:,2) , nx ,ny );
    Ez = reshape( E(:,3) , nx ,ny );
    lEl = sqrt( Ex.^2 + Ey.^2 + Ez.^2 );

    figure(1)
    axes('Position',[-.6 -0.6 2.25 2.25]);
    axis([-7 7 -7 7 -7 7]);
     view([30 45])
     axis off
     set(gcf,'color','w');
    pbaspect([1 1 1])
%     text(0,0,9,'$\oint_{S} \mathbf{e} \cdot \hat{\mathbf{n}} \; \mathrm{d}a = \frac{Q}{ \epsilon_{0} }$','interpreter','latex','FontSize',20,'HorizontalAlignment','center','VerticalAlignment','middle'); hold on

    if jj == 1
        
        h = surface(xx/10,yy/10,zz/10,lEl/10,'FaceAlpha',0);hold on
    
        h = surface(xx_q,yy_q,zz_q,lEl_q);
%         text(1,1,'$Q$','interpreter','latex','FontSize',15)
%         arrow3([0 0 0],[0-1.5 0 0],'r');
        arrow3([0 0 0],[0+1.5 0 0],'r');
%         arrow3([0 0 0],[0 0 0-1.5],'r');
        arrow3([0 0 0],[0 0 0+1.5],'r');
%         arrow3([0 0 0],[0 0-1.5 0],'r');
        arrow3([0 0 0],[0 0+1.5 0],'r');
        
        arrow3([0 0 0],[-1 0 0]/2,'r:',0.25,0.25);
        arrow3([0 0 0],[0 0 -1]/2,'r:',0.25,0.25);
        arrow3([0 0 0],[0 -1 0]/2,'r:',0.25,0.25);
        arrow3([0 0 0],[0 1 1]/2,'r:',0.25,0.25);
        arrow3([0 0 0],[1 0 1]/2,'r:',0.25,0.25);
        arrow3([0 0 0],[1 1 0]/2,'r:',0.25,0.25);
        arrow3([0 0 0],[0 -1 1]/2,'r:',0.25,0.25);
        arrow3([0 0 0],[0 1 -1]/2,'r:',0.25,0.25);
        arrow3([0 0 0],[-1 0 1]/2,'r:',0.25,0.25);
        arrow3([0 0 0],[-1 0 -1]/2,'r:',0.25,0.25);
        arrow3([0 0 0],[-1 1 0]/2,'r:',0.25,0.25);
        arrow3([0 0 0],[-1 -1 0]/2,'r:',0.25,0.25);
        arrow3([0 0 0],[1 0 -1]/2,'r:',0.25,0.25);
        arrow3([0 0 0],[1 -1 0]/2,'r:',0.25,0.25);
        arrow3([0 0 0],[0 -1 -1]/2,'r:',0.25,0.25);
        
        text(0,1.5,0,'$\vec e_y$','interpreter','latex','FontSize',15)
        text(1.5,0,0,'$\vec e_x$','interpreter','latex','FontSize',15)
        text(0,0,1.5,'$\vec e_z$','interpreter','latex','FontSize',15)
        
    elseif jj <= 10
        
        h = surface(xx/10*jj,yy/10*jj,zz/10*jj,lEl/10*jj,'FaceAlpha',0.1*jj);hold on
        h = surface(xx_q,yy_q,zz_q,lEl_q);
        arrow3([0 0 0],[1+jj/2 0 0],'k');
        arrow3([0 0 0],[0 0 1+jj/2],'k');
        arrow3([0 0 0],[0 1+jj/2 0],'k');
        
        arrow3([0 0 0],[-1-jj/2 0 0],'k:',0.25,0.25);
        arrow3([0 0 0],[0 0 -1-jj/2],'k:',0.25,0.25);
        arrow3([0 0 0],[0 -1-jj/2 0],'k:',0.25,0.25);
        arrow3([0 0 0],[0 1+jj/2 1+jj/2]/2,'k:',0.25,0.25);
        arrow3([0 0 0],[1+jj/2 0 1+jj/2]/2,'k:',0.25,0.25);
        arrow3([0 0 0],[1+jj/2 1+jj/2 0]/2,'k:',0.25,0.25);
        arrow3([0 0 0],[0 -1-jj/2 1+jj/2]/2,'k:',0.25,0.25);
        arrow3([0 0 0],[0 1+jj/2 -1-jj/2]/2,'k:',0.25,0.25);
        arrow3([0 0 0],[-1-jj/2 0 1+jj/2]/2,'k:',0.25,0.25);
        arrow3([0 0 0],[-1-jj/2 0 -1-jj/2]/2,'k:',0.25,0.25);
        arrow3([0 0 0],[-1-jj/2 1+jj/2 0]/2,'k:',0.25,0.25);
        arrow3([0 0 0],[-1-jj/2 -1-jj/2 0]/2,'k:',0.25,0.25);
        arrow3([0 0 0],[1+jj/2 0 -1-jj/2]/2,'k:',0.25,0.25);
        arrow3([0 0 0],[1+jj/2 -1-jj/2 0]/2,'k:',0.25,0.25);
        arrow3([0 0 0],[0 -1-jj/2 -1-jj/2]/2,'k:',0.25,0.25);
        
        text(0,1+jj/2,0,'$\vec e_y$','interpreter','latex','FontSize',15)
        text(1+jj/2,0,0,'$\vec e_x$','interpreter','latex','FontSize',15)
        text(0,0,1+jj/2,'$\vec e_z$','interpreter','latex','FontSize',15)
        text(lrl(14,50)*r(1450,1)*1.25,lrl(14,50)*r(1450,2)*1.25,lrl(14,50)*r(1450,3)*1.25,'$S$','interpreter','latex','FontSize',15)
       
        
    elseif jj < 73 && jj > 10
        
       
        h = surface(xx,yy,zz,lEl); hold on
            text(lrl(14,50)*r(1450,1)*1.25,lrl(14,50)*r(1450,2)*1.25,lrl(14,50)*r(1450,3)*1.25,'$S$','interpreter','latex','FontSize',15)
        arrow3([lrl(74,63)*r(6374,1) lrl(74,63)*r(6374,2) lrl(74,63)*r(6374,3)],[lrl(74,63)*r(6374,1)+.5 lrl(74,63)*r(6374,2)+.5 lrl(74,63)*r(6374,3)+.5],'k',[0.5 0.5 0.5]);
        text(lrl(74,63)*r(6374,1)+.5,lrl(74,63)*r(6374,2)+.5,lrl(74,63)*r(6374,3)+.5,'$\hat n$','interpreter','latex','FontSize',15);
         arrow3([0 0 0],[6 0 0],'k'); 
%         arrow3([0 0 0],[-6 0 0],'k');
        arrow3([0 0 0],[0 0 6],'k');
%         arrow3([0 0 0],[0 0 -6],'k');
        arrow3([0 0 0],[0 6 0],'k');
%         arrow3([0 0 0],[0 -6 0],'k');
        
        arrow3([0 0 0],[-6 0 0]/2,'k:',0.25,0.25);
        arrow3([0 0 0],[0 0 -6]/2,'k:',0.25,0.25);
        arrow3([0 0 0],[0 -6 0]/2,'k:',0.25,0.25);
        arrow3([0 0 0],[0 6 6]/2,'k:',0.25,0.25);
        arrow3([0 0 0],[6 0 6]/2,'k:',0.25,0.25);
        arrow3([0 0 0],[6 6 0]/2,'k:',0.25,0.25);
        arrow3([0 0 0],[0 -6 6]/2,'k:',0.25,0.25);
        arrow3([0 0 0],[0 6 -6]/2,'k:',0.25,0.25);
        arrow3([0 0 0],[-6 0 6]/2,'k:',0.25,0.25);
        arrow3([0 0 0],[-6 0 -6]/2,'k:',0.25,0.25);
        arrow3([0 0 0],[-6 6 0]/2,'k:',0.25,0.25);
        arrow3([0 0 0],[-6 -6 0]/2,'k:',0.25,0.25);
        arrow3([0 0 0],[6 0 -6]/2,'k:',0.25,0.25);
        arrow3([0 0 0],[6 -6 0]/2,'k:',0.25,0.25);
        arrow3([0 0 0],[0 -6 -6]/2,'k:',0.25,0.25);
        
        text(0,6,0,'$\vec e_y$','interpreter','latex','FontSize',15)
        text(6,0,0,'$\vec e_x$','interpreter','latex','FontSize',15)
        text(0,0,6,'$\vec e_z$','interpreter','latex','FontSize',15)
    else
        
        h = surface(xx,yy,zz,lEl,'FaceAlpha',1-(jj-20)/55); hold on
        h = surface(xx_q,yy_q,zz_q,lEl_q);
        text(1,1,'$Q$','interpreter','latex','FontSize',15)
%         arrow3([0 0 0],[0-1.5 0 0],'r');
        arrow3([0 0 0],[0+1.5 0 0],'r');
%         arrow3([0 0 0],[0 0 0-1.5],'r');
        arrow3([0 0 0],[0 0 0+1.5],'r');
%         arrow3([0 0 0],[0 0-1.5 0],'r');
        arrow3([0 0 0],[0 0+1.5 0],'r');
        
        arrow3([0 0 0],[-1 0 0]/2,'r:',0.25,0.25);
        arrow3([0 0 0],[0 0 -1]/2,'r:',0.25,0.25);
        arrow3([0 0 0],[0 -1 0]/2,'r:',0.25,0.25);
        arrow3([0 0 0],[0 1 1]/2,'r:',0.25,0.25);
        arrow3([0 0 0],[1 0 1]/2,'r:',0.25,0.25);
        arrow3([0 0 0],[1 1 0]/2,'r:',0.25,0.25);
        arrow3([0 0 0],[0 -1 1]/2,'r:',0.25,0.25);
        arrow3([0 0 0],[0 1 -1]/2,'r:',0.25,0.25);
        arrow3([0 0 0],[-1 0 1]/2,'r:',0.25,0.25);
        arrow3([0 0 0],[-1 0 -1]/2,'r:',0.25,0.25);
        arrow3([0 0 0],[-1 1 0]/2,'r:',0.25,0.25);
        arrow3([0 0 0],[-1 -1 0]/2,'r:',0.25,0.25);
        arrow3([0 0 0],[1 0 -1]/2,'r:',0.25,0.25);
        arrow3([0 0 0],[1 -1 0]/2,'r:',0.25,0.25);
        arrow3([0 0 0],[0 -1 -1]/2,'r:',0.25,0.25);
        
        text(0,1.5,0,'$\vec e_y$','interpreter','latex','FontSize',15)
        text(1.5,0,0,'$\vec e_x$','interpreter','latex','FontSize',15)
        text(0,0,1.5,'$\vec e_z$','interpreter','latex','FontSize',15)
    end
   
    colormap(jet)
    caxis([0 3])
    
    
%     axis equal

    frame = getframe(figure(1));
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if jj == 1;
      imwrite(imind,cm,[work_dir '\Efield.gif'],'gif', 'Loopcount',inf,'DelayTime',0.5);
    elseif jj < 10
          imwrite(imind,cm,[work_dir '\Efield.gif'],'gif','WriteMode','append','DelayTime',0.1);
    elseif jj == 10
              imwrite(imind,cm,[work_dir '\Efield.gif'],'gif','WriteMode','append','DelayTime',0.5);
    elseif jj == 11
        imwrite(imind,cm,[work_dir '\Efield.gif'],'gif','WriteMode','append','DelayTime',0.5);
    elseif jj == 75
        imwrite(imind,cm,[work_dir '\Efield.gif'],'gif','WriteMode','append','DelayTime',1);
    else
      imwrite(imind,cm,[work_dir '\Efield.gif'],'gif','WriteMode','append','DelayTime',0.05);
    end
    
    clf(figure(1))
end

