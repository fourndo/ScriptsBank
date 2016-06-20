% Minimum curvature test
clear all
close all

%% 1D
% Define base grid
nx = 20;
dx = 1;
xx = cumsum(ones(nx,1)*dx);

x = randn(nx,1);

% Build system
A = zeros(nx,nx);

for ii = 1 : nx
    
    A(ii,:) = abs(xx(ii) - xx).^3;
    
end

% Find coefficients
w = A\x;
m0 = ones(size(A,1),1)*1e-4;
%[w,~,~] = CGiter(m0,A,x);
[w,~] = BiCG_STAB(A,x,m0);

% Compute new solution for random locations
nrand = 100;
xrand = min(xx) + (max(xx)-min(xx)).*rand(nrand,1);
xrand = sort(xrand);

for ii = 1 : nrand
    
    d(ii) = sum(w.*abs(xx-xrand(ii)).^3);
    
end

% Plot result
figure;
plot(xx,x); hold on
plot(xrand,d,'r')

%% 2D
% Define base grid
nx = 20;
ny = 20;

dx = 1;
dy = 1;

xx = cumsum(ones(nx,1)*dx);
yy = cumsum(ones(ny,1)*dy);

[X,Y] = ndgrid(xx,yy);

x = exp( - ( (X - mean(xx)).^2/50 + (Y - mean(yy)).^2/50 ) );

% Build system
A = zeros(nx*ny,nx*ny);

count=1;
for ii = 1 : ny
    
    for jj = 1 : nx
    
        r = (xx(ii) - X(:)).^2 + (yy(jj) - Y(:)).^2 ;
        A(count,:) = r.*(log(sqrt(r+ 1e-8)) - 1 );

        count = count + 1;
        
    end
    
end

% Find coefficients
w = A\x(:);
m0 = ones(size(A,1),1)*1e-4;
%[w,~,~] = CGiter(m0,A,x);
% [w,~] = BiCG_STAB(A,x(:),m0);

% Compute new solution for random locations
nrand = 100;
xrand = min(xx) + (max(xx)-min(xx)).*rand(nrand,1);
yrand = min(xx) + (max(xx)-min(xx)).*rand(nrand,1);

for ii = 1 : nrand
    
    r = (xrand(ii) - X(:)).^2 + (yrand(ii) - Y(:)).^2;
    d(ii) = sum(w.*r.*(log(sqrt(r)) - 1 ));
    
end

% Plot result
figure;
surf(X,Y,x); hold on
scatter3(xrand,yrand,d,'r')