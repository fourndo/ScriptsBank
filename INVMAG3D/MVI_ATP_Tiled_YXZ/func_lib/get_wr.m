function [wr] = get_wr(obsx, obsy, obsz, D, I, xn, yn, zn, nullcell, FLAG, pwr, R0)
% Function [wr] = get_wr(x0, y0, z0, dx, dy, dz, nullcell, FLAG, pwr, R0)
% Computes the depth or distance weighting for PF Inversion
% INPUT
% obsx, obys, obsz          : Observation locations
% D, I                      : Inducing field Declination and Inclination
% xn, yn, zn                : Mesh node locations
% nullcell                  : Active cells [0:air, 1: ground]
% FLAG                      : Either 'DISTANCE' | 'DEPTH'
% pwr, R0                   : Exponent and treshold value 
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
        
%     R0 = min([min(dx) min(dy) min(dz)])*0;
%     R0 = 1;
    p = 1/sqrt(3);
    
    % Create cell center location
    Zm = (Zn(1:end-1,1:end-1,1:end-1) + Zn(2:end,2:end,2:end))/2;
    Xm = (Xn(1:end-1,1:end-1,1:end-1) + Xn(2:end,2:end,2:end))/2;
    Ym = (Yn(1:end-1,1:end-1,1:end-1) + Yn(2:end,2:end,2:end))/2;
    
    progress= -1;
    tic
    
    for dd = 1 : ndata
        nx1 = (Xm - dX * p - obsx(dd)).^2;
        nx2 = (Xm + dX * p - obsx(dd)).^2;

        ny1 = (Ym - dY * p - obsy(dd)).^2;
        ny2 = (Ym + dY * p - obsy(dd)).^2;

        nz1 = (Zm - dZ * p - obsz(dd)).^2;
        nz2 = (Zm + dZ * p - obsz(dd)).^2;

        R1 = sqrt(nx1 + ny1 + nz1);
        R2 = sqrt(nx1 + ny1 + nz2);
        R3 = sqrt(nx2 + ny1 + nz1);
        R4 = sqrt(nx2 + ny1 + nz2);
        R5 = sqrt(nx1 + ny2 + nz1);
        R6 = sqrt(nx1 + ny2 + nz2);
        R7 = sqrt(nx2 + ny2 + nz1);
        R8 = sqrt(nx2 + ny2 + nz2);
        
        temp = (R1 + R0).^-pwr  + ...
            (R2 + R0).^-pwr  +...
            (R3 + R0).^-pwr  +...
            (R4 + R0).^-pwr  +...
            (R5 + R0).^-pwr  +...
            (R6 + R0).^-pwr  +...
            (R7 + R0).^-pwr  +...
            (R8 + R0).^-pwr  ;

        wr = wr + (V2.*temp).^2;
        
        d_iter = floor(dd/ndata*20);
        if  d_iter > progress

            fprintf('Computed %i pct of data in %8.5f sec\r',d_iter*5,toc)
            progress = d_iter;

        end
    

    end
    
    wr = (wr.^(1/2))./V;
    wr = reshape(wr,nx*ny*nz,1);
    wr = (wr./(max(wr(nullcell==1)))).^(1/2);
    
    
    otherwise
        
        fprintf('Flag should be "Depth" or "DISTANCE"\n');
        
end




