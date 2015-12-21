function [wr] = MAG3C_RTC_get_dwr(obsx,obsy,obsz,xn,yn,zn,acellID,acelln,tcelln,linkat)
% Function [dwz] = compdepthw(x0, y0, z0, dx, dy, dz, nullcell)
% Computes the depth depth weighting for MAG3C
% INPUT
% z0        : Z Coordinates of grid origin
% dx, dy, dz: Cell size dimension
% nullcell  : Topomodel in vector binary format
% 
% OUTPUT
% dwz: depth weighting vector

dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);

% Pre-allocate memory
wr = zeros(1,nz*nx*ny);
awr = zeros(length(acellID),1);

% Number of observation location
ndata = length(obsx);
      
R0 = min([min(dx) min(dy) min(dz)])/4;
p = 1/sqrt(3);

z1 = acelln(:,1,1);
z2 = acelln(:,1,4);
x1 = acelln(:,1,2);
x2 = acelln(:,1,5);
y1 = acelln(:,1,3);
y2 = acelln(:,1,6);

Zm = (z2+z1)/2;
Xm = (x2+x1)/2;
Ym = (y2+y1)/2;

dZ = (z1-z2);
dX = (x2-x1);
dY = (y2-y1);

V = dZ.*dX.*dY;
dX = dX/2;
dY = dY/2;
dZ = dZ/2;
V2 = dZ.*dX.*dY;

fprintf('Start computing RTC distance weighting\n')
progress = -1;
tic
for dd = 1 : ndata
    
        
    nx1 = (Xm - dX * p - obsx(dd)).^2;
    nx2 = (Xm + dX * p - obsx(dd)).^2;

    ny1 = (Ym - dY * p - obsy(dd)).^2;
    ny2 = (Ym + dY * p - obsy(dd)).^2;

    nz1 = (Zm - dZ * p - obsz(dd)).^2;
    nz2 = (Zm + dZ * p - obsz(dd)).^2;

    temp = (sqrt(nx1 + ny1 + nz1) + R0).^-3 + ...
        (sqrt(nx1 + ny1 + nz2) + R0).^-3 +...
        (sqrt(nx2 + ny1 + nz1) + R0).^-3 +...
        (sqrt(nx2 + ny1 + nz2) + R0).^-3 +...
        (sqrt(nx1 + ny2 + nz1) + R0).^-3 +...
        (sqrt(nx1 + ny2 + nz2) + R0).^-3 +...
        (sqrt(nx2 + ny2 + nz1) + R0).^-3 +...
        (sqrt(nx2 + ny2 + nz2) + R0).^-3;
    
    temp = temp.*V2;
    % Loop topocells
    for jj = 1 : size(tcelln,1)
        
        tz1 = tcelln(jj,end,1:6:end);
        tz2 = tcelln(jj,end,4:6:end);
        tx1 = tcelln(jj,end,2:6:end);
        tx2 = tcelln(jj,end,5:6:end);
        ty1 = tcelln(jj,end,3:6:end);
        ty2 = tcelln(jj,end,6:6:end);
        
        tZm = (tz2+tz1)/2;
        tXm = (tx2+tx1)/2;
        tYm = (ty2+ty1)/2;

        tdZ = (tz1-tz2);
        tdX = (tx2-tx1);
        tdY = (ty2-ty1);

        tV = tdZ.*tdX.*tdY;

        if dd == 1
            
            V(linkat(jj)) = V(linkat(jj))+sum(tV);
            
        end

        tdX = tdX/2;
        tdY = tdY/2;
        tdZ = tdZ/2;
        tV2 = tdZ.*tdX.*tdY;

        nx1 = (tXm - tdX * p - obsx(dd)).^2;
        nx2 = (tXm + tdX * p - obsx(dd)).^2;

        ny1 = (tYm - tdY * p - obsy(dd)).^2;
        ny2 = (tYm + tdY * p - obsy(dd)).^2;

        nz1 = (tZm - tdZ * p - obsz(dd)).^2;
        nz2 = (tZm + tdZ * p - obsz(dd)).^2;

        ttemp = (sqrt(nx1 + ny1 + nz1) + R0).^-3 + ...
            (sqrt(nx1 + ny1 + nz2) + R0).^-3 +...
            (sqrt(nx2 + ny1 + nz1) + R0).^-3 +...
            (sqrt(nx2 + ny1 + nz2) + R0).^-3 +...
            (sqrt(nx1 + ny2 + nz1) + R0).^-3 +...
            (sqrt(nx1 + ny2 + nz2) + R0).^-3 +...
            (sqrt(nx2 + ny2 + nz1) + R0).^-3 +...
            (sqrt(nx2 + ny2 + nz2) + R0).^-3;
        
        ttemp = ttemp.*tV2;
        
        temp(linkat(jj)) = temp(linkat(jj))+sum(ttemp);
    
    end
     
    
    awr = awr + temp.^2;

    d_iter = floor(dd/ndata*100);
    if  d_iter > progress

        fprintf('Computed %i pct of data in %8.5f sec\n',d_iter,toc)
        progress = d_iter;

    end
    
end

wr(acellID) = (awr.^(1/2))./V;
wr = reshape(wr,nx*ny*nz,1);
wr = (wr./(max(wr))).^(1/2);
    
    
    




