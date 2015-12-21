clear all
close all

% Dipole epxeriment

x1 = -0.05;
y1 = 0;
z1 = 0.0;

x2 = 0.05;
y2 = 0;
z2 = 0.0;

dx = .1;
dz = .1;

ObsX = [-2:dx:2];
ObsZ = [2:-dz:-2];

[X,Z] = meshgrid(ObsX,ObsZ);

ndata = length(ObsX);

nx = size(X,1);
nz = size(X,2);

nobs = nx * nz;


nedges = nx*(nz-1) + (nx-1) * nz;

X = reshape(X,nx*nz,1);
Z = reshape(Z,nx*nz,1);

m1 = 1;
m2 = -0;


% for ii = 1: nobs
    
    r1 = ( (X - x1).^2 + ( Z - z1 ).^2 ) .^0.5 + 1e-6;
    r2 = ( (X - x2).^2 + ( Z - z2 ).^2 ) .^0.5 + 1e-6;
    
    delx = (X - x1);
    delz = (Z - z1);
    theta1 = atand( delx ./ delz);
    
    theta1(delz<0) = theta1(delz<0) + 180;
    theta1(delz(delx<0)>0) = theta1(delz(delx<0)>0) + 360;

    delx = (X - x2);
    delz = (Z - z2);
    theta2 = atand( delx ./ delz);
    
    theta2(delz<0) = theta2(delz<0) + 180;
    theta2(delz(delx<0)>0) = theta2(delz(delx<0)>0) + 360;
%     theta2 = atand( (X - x2) ./ (Z - z2));
    
    phi = m1./r1.^0.5 .* cosd(theta1) + m2./r2.^0.5 .* cosd(theta2+45);
%     phi = reshape(phi,nz,nx);
%     r = m1*2*cos(theta1(ii))/r1(ii)^3 + m2*2*cos(theta2(ii))/r2(ii)^3;
%     theta = m1*sin(theta1(ii))/r1(ii)^3 + m2*sin(theta2(ii))/r2(ii)^3;
    
% end
% xx = r1.*sind(theta1); xx = reshape(xx,nx,nz);
% yy = r1.*cosd(theta1); yy = reshape(yy,nx,nz);

ddx = @(n) spdiags (ones (n-1,1)*[1,-1],[0,1],n-1,n);

d_dx = 1/dx*ddx(nx); %d_dx=d_dx(:,2:end-1); d_dx([1,end])=0;
Dz = kron(speye(nz),d_dx);

d_dz = 1/dz*ddx(nz); %d_dz=d_dz(:,2:end-1); d_dz([1,end])=0;
Dx = kron( d_dz ,speye(nx));

GRAD = [Dz ; Dx];

B = GRAD* phi;

%% Averaging for plot
av = @(n) spdiags (ones (n+1,1)*[0.5,0.5],[-1,0],n,n-1);

avcfx = av(nx); avcfx([1,end]) = 1;
avcfz = av(nz); avcfz([1,end]) = 1;

Avcfz = kron ( speye (nz) , avcfx );
Avcfx = kron ( avcfz , speye (nx) );

Ox = sparse ( nx*nz , (nx-1)*nz );
Oz = sparse ( nx*nz , nx*(nz-1)  );

AVGz = [Avcfz Ox];
AVGx = [Oz Avcfx];

Bz = AVGz*B;
Bx = AVGx*B;
Bx2d = reshape(Bx,nz,nx);
Bz2d = reshape(Bz,nz,nx);

figure; quiver(Bx2d,Bz2d);

%%
% Bz = B(1:(nx-1)*nz);
% Bx = B(nx*(nz-1)+1:end);
% Bx2d = reshape(Bx,nz,nx-1);
% Bz2d = reshape(Bz,nz-1,nx);



% Bx = r.*cos(theta);
% By = r.*sin(theta);

% B = [Bx;By];




figure; imagesc(reshape(phi,nx,nz));title('Scalar Potential');caxis([-m1 m1]);colorbar;
figure; imagesc(Bx2d);caxis([-10 10]);colorbar;title('Bx');

figure; imagesc(Bz2d);caxis([-10 10]);colorbar; title('Bz');
%% Curl
ddx = @(n) spdiags (ones (n+1,1)*[1,-1],[0,1],n-1,n);

d_dx = 1/dz*ddx(nx); %d_dy=d_dy(:,2:end-1);
Dzx = kron(speye(nz-1), d_dx );

d_dz = 1/dx*ddx(nz); %d_dx=d_dx(:,2:end-1); 
Dxz = kron(d_dz , speye(nx-1) );

CURL = [-Dxz Dzx];

figure; imagesc(reshape(CURL*B,nx-1,nz-1));colorbar%caxis([-1 1]);
% for jj = 5:10
%     ObsZ = ones(ndata,1)*jj;
% 
%     m1 = [zeros(1,ndata);ones(1,ndata)];
%     m2 = [ones(1,ndata);zeros(1,ndata)];
%     
%     rvec = [(ObsX' - x0) ObsZ-z0];
% 
%     r = (rvec(:,1).^2 + rvec(:,2).^2).^0.5;
% 
%     H1= zeros(ndata,2);
%     H2 = zeros(ndata,2);
%     
%     for ii = 1:ndata
% 
%         H1(ii,:) = rvec(ii,:).*(rvec(ii,:)*m1(:,ii))./(r(ii).^5) - m1(:,ii)'/r(ii).^3;
%         H2(ii,:) = rvec(ii,:).*(rvec(ii,:)*m2(:,ii))./(r(ii).^5) - m2(:,ii)'/r(ii).^3;
%         
%     end
% 
%     figure(1); plot(ObsX,H1,ObsX,(H1(:,1).^2 + H1(:,2).^2).^0.5,'c:'); hold on
%     figure(2); plot(ObsX,H2,ObsX,(H2(:,1).^2 + H2(:,2).^2).^0.5,'r:'); hold on
% end