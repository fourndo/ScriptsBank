function [wr] = get_wr(obsx, obsy, obsz, D, I, xn, yn, zn, nullcell, FLAG)
% Function [dwz] = compdepthw(x0, y0, z0, dx, dy, dz, nullcell)
% Computes the depth depth weighting for MAG3C
% INPUT
% z0        : Z Coordinates of grid origin
% dx, dy, dz: Cell size dimension
% nullcell  : Topomodel in vector binary format
% 
% OUTPUT
% dwz: depth weighting vector

dx = xn(2:end) - xn(1:end-1);
dy = yn(2:end) - yn(1:end-1);
dz = zn(1:end-1) - zn(2:end);

nx = length(dx);
ny = length(dy);
nz = length(dz);

z0 = 4*min(dz);
wr = ones(nz,nx,ny)*1e-8;

ndata = length(obsx);

% 3D cell center location
[Zn,Xn,Yn] = ndgrid(zn,xn,yn);

% 3D cell size
[dZ,dX,dY] = ndgrid(dz,dx,dy);

% Volume of all cells
V = dZ.*dX.*dY;
dX = dX/2;
dY = dY/2;
dZ = dZ/2;
V2 = dZ.*dX.*dY;
switch FLAG
    case 'DEPTH'

M = [cosd(I) * cosd(D) ;
     cosd(I) * sind(D) ;
     sind(I)];


nullcell = reshape(nullcell, nz,nx,ny);

% Get average observation height
avgz = mean(obsz);

% Compute sensitivity for a column of cells
dodx(1) = -dx(1)/2 ; dodx(2) = dx(1)/2 ;
dody(1) = -dy(1)/2 ; dody(2) = dy(1)/2 ;
dodz = avgz - zn;
sens = zeros(nz,1);

for kk = 1 : nz
    
     
    [txx,txz,tyx,tyy,tyz] = MAG3C_T_cell( dodx, dody, [dodz(kk) dodz(kk+1)]);
    sens(kk) = M'*[txx tyx txz;tyx tyy tyz;txz tyz -(txx+tyy)]*M;

end
sens = sens / max(sens);
figure;plot(sens);hold on
% Find the best fitting z0
dz0new = 999;
dz0old = 1000;
z2 = zn(1) - zn( 2 :end );
z1 = zn(1) - zn( 1 : end-1 );
while dz0new < dz0old
    

        
    temp = -( (z2 + z0).^-2 - (z1 + z0).^-2 ) ./ (dz);
%     temp = (temp).^(1/2);
    temp = temp./(max(temp));
    
    dz0old = dz0new;
    dz0new = norm(sens-temp');
    
    z0 = z0 * 0.90;
    plot(temp,'r');hold on;
    
end

for jj = 1 : ny

    for ii = 1: nx

        acell = find(nullcell(:,ii,jj)==1);
        
        zcol = dz(acell);
        
        z2 = zn(1) - zn( (acell(1)+1) :end);
        z1 = zn(1) - zn( acell(1) : end-1 );
        
        wr(nullcell(:,ii,jj)==1,ii,jj) = ...
            -( (z2 + z0).^-2 - (z1 + z0).^-2 ) ./ (zcol);

    end

end
      

    wr = reshape(wr,nx*ny*nz,1);
    wr = (wr).^(1/2);
    wr = wr./(max(wr));

    case 'DISTANCE'
        
    R0 = min([min(dx) min(dy) min(dz)])/4;
    p = 1/sqrt(3);
    
    % Create cell center location
    Zm = (Zn(1:end-1,1:end-1,1:end-1) + Zn(2:end,2:end,2:end))/2;
    Xm = (Xn(1:end-1,1:end-1,1:end-1) + Xn(2:end,2:end,2:end))/2;
    Ym = (Yn(1:end-1,1:end-1,1:end-1) + Yn(2:end,2:end,2:end))/2;
    
    
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
    
    wr = wr + (V2.*temp).^2;
     
end
    
    wr = (wr.^(1/2))./V;
    wr = reshape(wr,nx*ny*nz,1);
    wr = (wr./(max(wr(nullcell==1)))).^(1/2);
    
    
    otherwise
        
        fprintf('Flag should be "depth" or "distance"\n');
        
end




