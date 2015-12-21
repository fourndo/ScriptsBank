function [Tx,Ty,Tz] = MAG3C_T_row(dobsx,dobsy,dobsz)
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
% nullcell   : 1D vector for active(1) and inactive (0) cells
% Last update: February 17th, 2014

nx = length(dobsx)-1;
ny = length(dobsy)-1;
nz = length(dobsz)-1;

mcell = nx*ny*nz;

% Pre-allocate to store fields
Tx = zeros(1,3*mcell);
Ty = zeros(1,3*mcell);
Tz = zeros(1,3*mcell);
  

nn = 1;

for jj = 1 : ny

%     %First compute de location of the center of the cell
%     Y = y0 + sum(dy(1:jj)) - dy(jj) /2;
% 
%     %Then the distance between the cell corners and the observation
    dy(2) = dobsy(jj+1) ;
    dy(1) = dobsy(jj) ;

   for ii = 1 : nx

%         X = x0 + sum(dx(1:ii)) - dx(ii) /2;
% 
        dx(2) = dobsx(ii+1) ;
        dx(1) = dobsx(ii) ;

       for kk = 1 : nz            

%             Z = z0 - sum(dz(1:kk)) + dz(kk) /2;
% 
            dz(1) = dobsz(kk) ;
            dz(2) = dobsz(kk+1) ;

            R1 = ( dy(2)^2 + dx(2)^2 );
            R2 = ( dy(2)^2 + dx(1)^2 );
            R3 = ( dy(1)^2 + dx(2)^2 );
            R4 = ( dy(1)^2 + dx(1)^2 );

            arg1 = sqrt( dz(2)^2 + R2 );
            arg2 = sqrt( dz(2)^2 + R1 );
            arg3 = sqrt( dz(1)^2 + R1 );
            arg4 = sqrt( dz(1)^2 + R2 );
            arg5 = sqrt( dz(2)^2 + R3 );
            arg6 = sqrt( dz(2)^2 + R4 );
            arg7 = sqrt( dz(1)^2 + R4 );
            arg8 = sqrt( dz(1)^2 + R3 );

            Ty(1,mcell+nn) = atan( dx(1) * dz(2) / ( dy(2) * arg1 ) ) +...
                - atan( dx(2) * dz(2) / ( dy(2) * arg2 ) ) +...
                  atan( dx(2) * dz(1) / ( dy(2) * arg3 ) ) +...
                - atan( dx(1) * dz(1) / ( dy(2) * arg4 ) ) +...
                  atan( dx(2) * dz(2) / ( dy(1) * arg5 ) ) +...
                - atan( dx(1) * dz(2) / ( dy(1) * arg6 ) ) +...
                  atan( dx(1) * dz(1) / ( dy(1) * arg7 ) ) +...
                - atan( dx(2) * dz(1) / ( dy(1) * arg8 ) );

            Ty(1,nn) = log( ( dz(2) + arg2 ) / (dz(1) + arg3 ) ) +...
                    -log( ( dz(2) + arg1 ) / (dz(1) + arg4 ) ) +...
                     log( ( dz(2) + arg6 ) / (dz(1) + arg7 ) ) +...
                    -log( ( dz(2) + arg5 ) / (dz(1) + arg8 ) );


            Tx(1,nn) = atan( dy(1) * dz(2) / ( dx(2) * arg5 ) ) +...
                - atan( dy(2) * dz(2) / ( dx(2) * arg2 ) ) +...
                  atan( dy(2) * dz(1) / ( dx(2) * arg3 ) ) +...
                - atan( dy(1) * dz(1) / ( dx(2) * arg8 ) ) +...
                  atan( dy(2) * dz(2) / ( dx(1) * arg1 ) ) +...
                - atan( dy(1) * dz(2) / ( dx(1) * arg6 ) ) +...
                  atan( dy(1) * dz(1) / ( dx(1) * arg7 ) ) +...
                - atan( dy(2) * dz(1) / ( dx(1) * arg4 ) );

            R1 = (dy(2)^2 + dz(1)^2);
            R2 = (dy(2)^2 + dz(2)^2);
            R3 = (dy(1)^2 + dz(1)^2);
            R4 = (dy(1)^2 + dz(2)^2);

            Ty(1,2*mcell+nn) = log( ( dx(1) + sqrt( dx(1)^2 + R1 ) ) / (dx(2) + sqrt( dx(2)^2 + R1 ) ) ) +...
                    -log( ( dx(1) + sqrt( dx(1)^2 + R2 ) ) / (dx(2) + sqrt( dx(2)^2 + R2 ) ) ) +...
                     log( ( dx(1) + sqrt( dx(1)^2 + R4 ) ) / (dx(2) + sqrt( dx(2)^2 + R4 ) ) ) +...
                    -log( ( dx(1) + sqrt( dx(1)^2 + R3 ) ) / (dx(2) + sqrt( dx(2)^2 + R3 ) ) );

            R1 = (dx(2)^2 + dz(1)^2);
            R2 = (dx(2)^2 + dz(2)^2);
            R3 = (dx(1)^2 + dz(1)^2);
            R4 = (dx(1)^2 + dz(2)^2);

            Tx(1,2*mcell+nn) = log( ( dy(1) + sqrt( dy(1)^2 + R1 ) ) / (dy(2) + sqrt( dy(2)^2 + R1 ) ) ) +...
                    -log( ( dy(1) + sqrt( dy(1)^2 + R2 ) ) / (dy(2) + sqrt( dy(2)^2 + R2 ) ) ) +...
                     log( ( dy(1) + sqrt( dy(1)^2 + R4 ) ) / (dy(2) + sqrt( dy(2)^2 + R4 ) ) ) +...
                    -log( ( dy(1) + sqrt( dy(1)^2 + R3 ) ) / (dy(2) + sqrt( dy(2)^2 + R3 ) ) );

            Tz(1,2*mcell+nn) = -( Ty(1,mcell+nn) + Tx(1,nn) );
            Tz(1,mcell+nn) = Ty(1,2*mcell+nn);
            Tx(1,mcell+nn) = Ty(1,nn);
            Tz(1,nn) = Tx(1,2*mcell+nn);

            nn = nn+1;

       end

   end

end

Tx = Tx/(4*pi);
Ty = Ty/(4*pi);
Tz = Tz/(4*pi);


