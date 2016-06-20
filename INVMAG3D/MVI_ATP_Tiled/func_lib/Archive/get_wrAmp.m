function [wr] = get_wrAmp(obsx, obsy, obsz, d, D, I, xn, yn, zn, nullcell, FLAG, R, R0)
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

mcell = nx * ny* nz;
z0 = 4*min(dz);
wr = ones(nz*nx*ny,1)*1e-8;

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

d = d / norm(d);
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
    Zm = (Zn(1:end-1,1:end-1,1:end-1) + Zn(2:end,2:end,2:end))/2; Zm = Zm(:);
    Xm = (Xn(1:end-1,1:end-1,1:end-1) + Xn(2:end,2:end,2:end))/2; Xm = Xm(:);
    Ym = (Yn(1:end-1,1:end-1,1:end-1) + Yn(2:end,2:end,2:end))/2; Ym = Ym(:);
    
    dX = dX(:);
    dY = dY(:);
    dZ = dZ(:);
    
    progress= -1;
    tic
    
    % TEST -> ADD LINEAR WEIGTHING IN DIRECTION OF MAGNETIZATION
    
    Ptmi = [(cosd(abs(I)) * cosd(D)) (cosd(abs(I)) * sind(D)) sind(abs(I))];

    wzc = kron(kron(ones(ny,1),ones(nx,1)),( zn(1:end-1)' + zn(2:end)' )/2);
    wxc = kron(kron(ones(ny,1),( xn(1:end-1)' + xn(2:end)' )/2),ones(nz,1));
    wyc = kron(kron(( yn(1:end-1)' + yn(2:end)' )/2,ones(nx,1)),ones(nz,1));

    wm = Ptmi * [wxc';wyc';wzc'];
    
%     ii = find(max(wm));
    %
%     rM = sqrt( (wzc - wzc(ii)).^2 + (wyc - wyc(ii)).^2 + (wxc - wxc(ii)).^2 );
    
    wm = wm - min(wm);
    wm = 1-wm / max(wm) + 0.1;
% 
% 
%     wM = reshape(wM,nz,nx,ny);

%     wr = wr(:)./sqrt(wM);

    for dd = 1 : ndata
        dx1 = ((Xm - dX * p - obsx(dd)))';
        dx2 = ((Xm + dX * p - obsx(dd)))';

        dy1 = ((Ym - dY * p - obsy(dd)))';
        dy2 = ((Ym + dY * p - obsy(dd)))';

        dz1 = ((Zm - dZ * p - obsz(dd)))';
        dz2 = ((Zm + dZ * p - obsz(dd)))';

        R1 = sqrt(dx1.^2 + dy1.^2 + dz1.^2);
        R2 = sqrt(dx1.^2 + dy1.^2 + dz2.^2);
        R3 = sqrt(dx2.^2 + dy1.^2 + dz1.^2);
        R4 = sqrt(dx2.^2 + dy1.^2 + dz2.^2);
        R5 = sqrt(dx1.^2 + dy2.^2 + dz1.^2);
        R6 = sqrt(dx1.^2 + dy2.^2 + dz2.^2);
        R7 = sqrt(dx2.^2 + dy2.^2 + dz1.^2);
        R8 = sqrt(dx2.^2 + dy2.^2 + dz2.^2);
        
%         R= R';
        
%         rr = ([dx1;dy1;dz1]) * spdiags(1./R1',0,mcell,mcell);
%         
%         wm = abs(Ptmi * rr);
%         wm = wm - min(wm);
%         wm = wm/max(wm) + 0.01;
%         r2 = ([dx1;dy1;dz2])* spdiags(1./R2',0,mcell,mcell);
%         r3 = ([dx2;dy1;dz1])* spdiags(1./R3',0,mcell,mcell);
%         r4 = ([dx2;dy1;dz2])* spdiags(1./R4',0,mcell,mcell);
%         r5 = ([dx1;dy2;dz1])* spdiags(1./R5',0,mcell,mcell);
%         r6 = ([dx1;dy2;dz2])* spdiags(1./R6',0,mcell,mcell);
%         r7 = ([dx2;dy2;dz1])* spdiags(1./R7',0,mcell,mcell);
%         r8 = ([dx2;dy2;dz2])* spdiags(1./R8',0,mcell,mcell);
        
        temp = (R1 + R0).^-R  + ...
            (R2 + R0).^-R  +...
            (R3 + R0).^-R  +...
            (R4 + R0).^-R  +...
            (R5 + R0).^-R  +...
            (R6 + R0).^-R  +...
            (R7 + R0).^-R  +...
            (R8 + R0).^-R  ;

        wr = wr + (  V2(:).*temp' ).^2;
        
        d_iter = floor(dd/ndata*20);
        if  d_iter > progress

            fprintf('Computed %i pct of data in %8.5f sec\n',d_iter*5,toc)
            progress = d_iter;

        end
    

    end
    
    wr = (wr.^(1/2))./V(:);
%     wr = reshape(wr,nx*ny*nz,1);
    wr = (wr./(max(wr(nullcell==1))).*wm').^(1/2);
    
    case 'CENTER'
        
        % Create cell center location
        Zm = (Zn(1:end-1,1:end-1,1:end-1) + Zn(2:end,2:end,2:end))/2;
        Xm = (Xn(1:end-1,1:end-1,1:end-1) + Xn(2:end,2:end,2:end))/2;
        Ym = (Yn(1:end-1,1:end-1,1:end-1) + Yn(2:end,2:end,2:end))/2;
        
        cntrx = mean(xn);
        cntry = mean(yn);
        cntrz = mean(zn);
        
%         wrx = (Xm - cntrx).^2; wrx = wrx ./ max([wrx(:);1e-8]);
%         wry = (Ym - cntry).^2; wry = wry ./ max([wry(:);1e-8]);
%         wrz = (Zm - cntrz).^2; wrz = wrz ./ max([wrz(:);1e-8]);
%         wr = (wrx + wry + wrz) .^ (3/2);
%         wr = reshape(wr,nx*ny*nz,1);
%         wr = (wr./(max(wr(nullcell==1)))).^(1/2);
        
    % Compute rotated elliptical weights
    wr = (( (Xm - cntrx) * cosd(D) + (Ym - cntry) * sind(D) + 1e-1 ).^2 ./ (min(dx)).^2 +...
    ( (Xm - cntrx) * sind(D) - (Ym - cntry) * cosd(D) + 1e-1 ).^2 ./ (min(dy)/2).^2) .^ (0.5) ;

    wr = wr ./ max(wr(:));
    wr = (wr ) .^ (1/4);
    
    wr = reshape(wr,nx*ny*nz,1);
    wr = wr./(max(wr(nullcell==1)));
    
    otherwise
        
        fprintf('Flag should be "DEPTH" or "DISTANCE"\n');
        
end




