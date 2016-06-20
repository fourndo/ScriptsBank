% Minimum curvature test
clear all
close all

% Number of samples
ns = 10;


% Define base grid
nx = 20;
ny = 20;
dx = 1;
dy = 1;


xx = cumsum(ones(nx,1) * dx);
yy = cumsum( ones(ny,1) * dy );

[X,Y] = ndgrid(xx,yy);
% Actv cells
actv = randi([0,1],nx*ny,1);

% Create random array of points
% d = randi([-10,10],nx,ny);
d = exp( - ( (X - mean(xx)).^2/50 + (Y - mean(yy)).^2/50 ) );

% Create Random sample points
xrand = xx(randi([3,length(xx)-2],ns,1)) + 0.5;%randn(ns,1)/10;
yrand = yy(randi([3,length(yy)-2],ns,1)) + 0.5;%randn(ns,1)/10;

% xrand = 10.5;
% yrand = 8.5;
% figure;
% surf(X,Y,d,'facealpha',0.5); hold on
% title('True signal')
% colormap(jet)
% colorbar
% 
% scatter3(xrand,yrand,zeros(size(xrand)),50, 'filled','MarkerEdgeColor','k')

%% Create gradient operator
ddx = @(n) spdiags (ones (n+1,1)*[-1,1],[0,1],n-1,n);


d_dx = ddx(nx); %d_dx(end,end-1) = 1;

dX = kron(speye(ny),spdiags(ones(nx,1)/dx,0,nx,nx));

Gx = kron(speye(ny),d_dx);

% Y gradient
d_dy = ddx(ny); %d_dy(end,end-1) = 1;

dY = kron( spdiags(ones(ny,1)/dy,0,ny,ny), speye(nx) );

Gy = kron(d_dy , speye(nx));

% Second order operator
Gxx = Gx'*Gx;
Gyy = Gy'*Gy;

G = (Gxx + Gyy);
% Compute curvature of base mesh
curv = G*d(:);

figure;
surf(X,Y,reshape(curv,size(d)),'facealpha',0.5); hold on
title('Curvature')
colormap(jet)
colorbar

m0 = d(:);
[x,~,~] = CGiter(m0,G'*G,G'*zeros(nx*ny,1));

% Check to see what the final function would look like
figure;
surf(X,Y,d,'facealpha',0.5); hold on
title('True signal')
colormap(jet)
colorbar
scatter3(X(:),Y(:), x,50, 'filled','MarkerEdgeColor','k')

%% Find the nearest nodes on the base grid

gxx = sparse(ns,nx*ny+ns, 3*ns);
gyy = sparse(ns,nx*ny+ns, 3*ns);

x = zeros(ns,1);
for ii = 1 : length(xrand)
    
    delx = xrand(ii) - xx(1);
    dely = yrand(ii) - yy(1);
    
    indx(ii) = ceil(delx/dx);    
    indy(ii) = ceil(dely/dy); 
    
    dx1(ii) = abs(xrand(ii) - xx(indx(ii)));
    dx2(ii) = abs(xrand(ii) - xx(indx(ii)+1));
    
    dy1(ii) = abs(yrand(ii) - yy(indy(ii)));
    dy2(ii) = abs(yrand(ii) - yy(indy(ii)+1));
   
    x(ii) = ( 8 * ( d(indx(ii)+1,indy(ii)) + d(indx(ii)-1,indy(ii)+1) + d(indx(ii),indy(ii)-1) + d(indx(ii),indy(ii)+1)) -...
        2*( d(indx(ii)+1,indy(ii)+1) + d(indx(ii)-1,indy(ii)+1) + d(indx(ii)+1,indy(ii)-1) + d(indx(ii)-1,indy(ii)-1) ) -...
        (d(indx(ii)+2,indy(ii)) + d(indx(ii)-2,indy(ii)) + d(indx(ii),indy(ii)-2) + d(indx(ii),indy(ii)+2) ) ) /20;
    
    gxx(ii,sub2ind([nx,ny],indx(ii),indy(ii))) = 1;
%     gxx(nx*ny+ii,sub2ind([nx,ny],indx(ii)+1,indy(ii))) = 1/dx2(ii)^2;
    gxx(ii,nx*ny+ii) = -1;%/dx1(ii)*dx2(ii);
    
%     RHS(sub2ind([nx,ny],indx(ii),indy(ii))) = -1/dx1(ii)^2;
%     RHS(sub2ind([nx,ny],indx(ii)+1,indy(ii))) = -1/dx2(ii)^2;
    
%     RHS(ii)= -d(sub2ind([nx,ny],indx(ii),indy(ii))) -...
%         d(sub2ind([nx,ny],indx(ii)+1,indy(ii))) -...
%         d(sub2ind([nx,ny],indx(ii),indy(ii))) -...
%         d(sub2ind([nx,ny],indx(ii),indy(ii)+1));
    
    
%     gxx(ii,ii) = -2;%/(dx2(ii)+dx1(ii));
    
%     RHS(sub2ind([nx,ny],indx(ii),indy(ii))) = -1/dy1(ii)^2;
%     RHS(sub2ind([nx,ny],indx(ii),indy(ii))+1) = -1/dy2(ii)^2;
    
    gyy(ii,sub2ind([nx,ny],indx(ii),indy(ii))) = 1;%/dy1(ii)^2;
%     gyy(nx*ny+ii,sub2ind([nx,ny],indx(ii),indy(ii)+1)) = 1/dy2(ii)^2;
    gyy(ii,nx*ny+ii) = -1;%/dy1(ii)*dy2(ii);
    
end

% Gxx = gxx'*gxx;
% Gyy = gyy'*gyy;
% scatter3(xx(indx),yy(indy),zeros(size(xrand)),25,'r', 'filled','MarkerEdgeColor','r')
% scatter3(xx(indx+1),yy(indy+1),zeros(size(xrand)),25,'b', 'filled','MarkerEdgeColor','b')

%% Bring together everything

% Make sure that none of the interpolated points are used
Gx = [ [Gx sparse(size(Gx,1),ns)];gxx];
Gy = [ [Gy sparse(size(Gy,1),ns)];gyy];

G = Gx'*Gx + Gy'*Gy;

% G = gxx+gyy;
%% Compute solution
F = scatteredInterpolant(X(:),Y(:),d(:));
m0 = F(xrand,yrand );
% x = [d(:);m0(:)];
% 
% P = spdiags([zeros(nx*ny,1);ones(ns,1)],0,length(x),length(x));
% 
% %scatter3([X(:);xrand],[Y(:);yrand], m0,50, 'filled','MarkerEdgeColor','k')
% 
RHS = [curv;zeros(ns,1)];
% 

% for ii = 1 : 5
% [x,~,~] = CGiter(x,(G'*G),(G)'*RHS);
% 
% 
% interp_x = x(end-ns+1:end);

%% Iterate 
% x = zeros(ns,1);
% for ii = 1 : length(xrand)
%    
%     x(ii) = (4*d(sub2ind([nx,ny],indx(ii)-1,indy(ii))) +...
%         4*d(sub2ind([nx,ny],indx(ii),indy(ii)+1)) +...
%         4*d(sub2ind([nx,ny],indx(ii)+1,indy(ii))) -...
%         (d(sub2ind([nx,ny],indx(ii)-2,indy(ii))) +...
%         d(sub2ind([nx,ny],indx(ii)+2,indy(ii))) +...
%         d(sub2ind([nx,ny],indx(ii),indy(ii)+2)) +...
%         d(sub2ind([nx,ny],indx(ii)-1,indy(ii)+1)) +...
%         d(sub2ind([nx,ny],indx(ii)+1,indy(ii)+1))))/7 ;
%     
%     %x(ii) = x(ii)/sum(1/dx1(ii)+1/dx2(ii)+1/dy1(ii)+1/dy2(ii));
%     
% end

%% Plot result
figure;
surf(X,Y,d,'facealpha',0.5); hold on
title('True signal')
colormap(jet)
colorbar


%scatter3(xrand,yrand, interp_x,50, 'filled','MarkerEdgeColor','k')
scatter3([xrand],[yrand],x,50, 'filled','MarkerEdgeColor','k')
% scatter3([X(:);xrand],[Y(:);yrand],x,10, 'filled','MarkerEdgeColor','k','MarkerFaceColor','r')
% scatter3(X(indx==0),Y(indx==0),d(indx==0),20,'kv', 'filled')

% figure;
% surf(X, Y,reshape(x(1:end-ns),nx,ny)); hold on
% title('Recovered signal')
% colormap(jet)
% colorbar

% end