function [Tx,Ty,Tz] = MAG3C_RTC_T_v2(obsx,obsy,obsz,mcell,tcelln,jlink)
% TMAG3C(X0, Y0, Z0, dX, dY, dZ, X , Y , Z , H, I, D)
%
% NEED TO UPDATE HEADER
% Forward operator for the magnetic response from an orthogonal mesh 
% due to an inducing field H for a given location (single row of G)
% Reference: Sharma (1963) Computation of magnetic anomaly
% Inputs:
% xo, yo, zo : Coordinates of the South-West-Top corner of the mesh
% dx, dy, dz : Vectors of cell sizes for the East, North and Vertical axis
% obsx,  obsy,  obsz  : Coordinates of the observation points 
%      H     : Magnitude of inducing field (nT)
%      I     : Inclinaison of H (degrees from horizontal)
%      D     : Declinaison of H (degrees from North)
% nullcell   : 1D vector for active1 and inactive (0) cells
% Last update: February 17th, 2014

% nx = length(dobsx)-1;
% ny = length(dobsy)-1;
% nz = length(dobsz)-1;

% mcell = nx*ny*nz;
e = 1e-12;
% Pre-allocate to store fields
Tx = zeros(1,3*mcell);
Ty = zeros(1,3*mcell);
Tz = zeros(1,3*mcell);
  
% RTC step
for ii = 1 : size(tcelln,1)

    ind = jlink(ii);
    if ind ==0
        continue
    end
    temp = tcelln(ii,ind, 1 : nnz(tcelln(ii,ind,:)) );

    z1 = temp(1:6:end);
    z2 = temp(4:6:end);
    x2 = temp(2:6:end);
    x1 = temp(5:6:end);
    y2 = temp(3:6:end);
    y1 = temp(6:6:end);

    dy2 = y2 - obsy ;
    dy1 = y1 - obsy ;

    dx2 = x2 - obsx ;
    dx1 = x1 - obsx ;
 
    dz1 = obsz - z1 ;
    dz2 = obsz - z2 ;

    R1 = ( dy2.^2 + dx2.^2 );
    R2 = ( dy2.^2 + dx1.^2 );
    R3 = ( dy1.^2 + dx2.^2 );
    R4 = ( dy1.^2 + dx1.^2 );

    arg1 = sqrt( dz2.^2 + R2 );
    arg2 = sqrt( dz2.^2 + R1 );
    arg3 = sqrt( dz1.^2 + R1 );
    arg4 = sqrt( dz1.^2 + R2 );
    arg5 = sqrt( dz2.^2 + R3 );
    arg6 = sqrt( dz2.^2 + R4 );
    arg7 = sqrt( dz1.^2 + R4 );
    arg8 = sqrt( dz1.^2 + R3 );

    Ty(1,mcell+ii) = sum(atan( dx1 .* dz2 ./ ( dy2 .* arg1 ) ) +...
    - atan( dx2 .* dz2 ./ ( dy2 .* arg2 ) ) +...
      atan( dx2 .* dz1 ./ ( dy2 .* arg3 ) ) +...
    - atan( dx1 .* dz1 ./ ( dy2 .* arg4 ) ) +...
      atan( dx2 .* dz2 ./ ( dy1 .* arg5 ) ) +...
    - atan( dx1 .* dz2 ./ ( dy1 .* arg6 ) ) +...
      atan( dx1 .* dz1 ./ ( dy1 .* arg7 ) ) +...
    - atan( dx2 .* dz1 ./ ( dy1 .* arg8 ) ));

    Ty(1,ii) = sum(log( ( dz2 + arg2 ) ./ (dz1 + arg3 ) ) +...
        -log( ( dz2 + arg1 ) ./ (dz1 + arg4 ) ) +...
         log( ( dz2 + arg6 ) ./ (dz1 + arg7 ) ) +...
        -log( ( dz2 + arg5 ) ./ (dz1 + arg8 ) ));

    Tx(1,ii) = sum(atan( dy1 .* dz2 ./ ( dx2 .* arg5 ) ) +...
        - atan( dy2 .* dz2 ./ ( dx2 .* arg2 ) ) +...
          atan( dy2 .* dz1 ./ ( dx2 .* arg3 ) ) +...
        - atan( dy1 .* dz1 ./ ( dx2 .* arg8 ) ) +...
          atan( dy2 .* dz2 ./ ( dx1 .* arg1 ) ) +...
        - atan( dy1 .* dz2 ./ ( dx1 .* arg6 ) ) +...
          atan( dy1 .* dz1 ./ ( dx1 .* arg7 ) ) +...
        - atan( dy2 .* dz1 ./ ( dx1 .* arg4 ) ));

        R1 = (dy2.^2 + dz1.^2);
        R2 = (dy2.^2 + dz2.^2);
        R3 = (dy1.^2 + dz1.^2);
        R4 = (dy1.^2 + dz2.^2);

        Ty(1,2*mcell+ii) = sum(log( ( dx1 + sqrt( dx1.^2 + R1 ) ) ./ (dx2 + sqrt( dx2.^2 + R1 ) ) ) +...
                -log( ( dx1 + sqrt( dx1.^2 + R2 ) ) ./ (dx2 + sqrt( dx2.^2 + R2 ) ) ) +...
                 log( ( dx1 + sqrt( dx1.^2 + R4 ) ) ./ (dx2 + sqrt( dx2.^2 + R4 ) ) ) +...
                -log( ( dx1 + sqrt( dx1.^2 + R3 ) ) ./ (dx2 + sqrt( dx2.^2 + R3 ) ) ) );

        R1 = (dx2.^2 + dz1.^2);
        R2 = (dx2.^2 + dz2.^2);
        R3 = (dx1.^2 + dz1.^2);
        R4 = (dx1.^2 + dz2.^2);

        Tx(1,2*mcell+ii) = sum( log( ( dy1 + sqrt( dy1.^2 + R1 ) ) ./ (dy2 + sqrt( dy2.^2 + R1 ) ) ) +...
                -log( ( dy1 + sqrt( dy1.^2 + R2 ) ) ./ (dy2 + sqrt( dy2.^2 + R2 ) ) ) +...
                 log( ( dy1 + sqrt( dy1.^2 + R4 ) ) ./ (dy2 + sqrt( dy2.^2 + R4 ) ) ) +...
                -log( ( dy1 + sqrt( dy1.^2 + R3 ) ) ./ (dy2 + sqrt( dy2.^2 + R3 ) ) ) );

        Tz(1,2*mcell+ii) = -( Ty(1,mcell+ii) + Tx(1,ii) );
        Tz(1,mcell+ii) = Ty(1,2*mcell+ii);
        Tx(1,mcell+ii) = Ty(1,ii);
        Tz(1,ii) = Tx(1,2*mcell+ii);
        
end

Tx = Tx/(4*pi);
Ty = Ty/(4*pi);
Tz = Tz/(4*pi);


