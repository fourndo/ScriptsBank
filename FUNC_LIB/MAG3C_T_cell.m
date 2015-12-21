function [txx,txz,tyx,tyy,tyz] = MAG3C_T_cell( dx, dy, dz)
% MAG3C_T_cell(x,y,z)
%
% NEED TO UPDATE HEADER
% Forward operator for the magnetic response from an orthogonal mesh 
% due to an inducing field H for a given location (single row of G)
% Reference: Sharma (1963) Computation of magnetic anomaly
% Inputs:
% x, y, z : Coordinates of the cell corners, 
% (1)=south/west/top, (2)=north,east,bottom
%
% obsx,  obsy,  obsz  : Coordinates of the observation points 
%
% Last update: May 26th, 2014
% Written by: D.Fournier

% Pre-allocate to store kernel

       
%Then the distance between the cell corners and the observation
% dy(2) = y(2) - obsy;
% dy(1) = y(1) - obsy;
% 
% dx(2) = x(2) - obsx;
% dx(1) = x(1) - obsx;
% 
% dz(1) = obsz - z(1) ;
% dz(2) = obsz - z(2) ;

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

tyy = atan( dx(1) * dz(2) / ( dy(2) * arg1 ) ) +...
    - atan( dx(2) * dz(2) / ( dy(2) * arg2 ) ) +...
      atan( dx(2) * dz(1) / ( dy(2) * arg3 ) ) +...
    - atan( dx(1) * dz(1) / ( dy(2) * arg4 ) ) +...
      atan( dx(2) * dz(2) / ( dy(1) * arg5 ) ) +...
    - atan( dx(1) * dz(2) / ( dy(1) * arg6 ) ) +...
      atan( dx(1) * dz(1) / ( dy(1) * arg7 ) ) +...
    - atan( dx(2) * dz(1) / ( dy(1) * arg8 ) );

tyx = log( ( dz(2) + arg2 ) / (dz(1) + arg3 ) ) +...
        -log( ( dz(2) + arg1 ) / (dz(1) + arg4 ) ) +...
         log( ( dz(2) + arg6 ) / (dz(1) + arg7 ) ) +...
        -log( ( dz(2) + arg5 ) / (dz(1) + arg8 ) );


txx = atan( dy(1) * dz(2) / ( dx(2) * arg5 ) ) +...
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

tyz = log( ( dx(1) + sqrt( dx(1)^2 + R1 ) ) / (dx(2) + sqrt( dx(2)^2 + R1 ) ) ) +...
        -log( ( dx(1) + sqrt( dx(1)^2 + R2 ) ) / (dx(2) + sqrt( dx(2)^2 + R2 ) ) ) +...
         log( ( dx(1) + sqrt( dx(1)^2 + R4 ) ) / (dx(2) + sqrt( dx(2)^2 + R4 ) ) ) +...
        -log( ( dx(1) + sqrt( dx(1)^2 + R3 ) ) / (dx(2) + sqrt( dx(2)^2 + R3 ) ) );

R1 = (dx(2)^2 + dz(1)^2);
R2 = (dx(2)^2 + dz(2)^2);
R3 = (dx(1)^2 + dz(1)^2);
R4 = (dx(1)^2 + dz(2)^2);

txz = log( ( dy(1) + sqrt( dy(1)^2 + R1 ) ) / (dy(2) + sqrt( dy(2)^2 + R1 ) ) ) +...
        -log( ( dy(1) + sqrt( dy(1)^2 + R2 ) ) / (dy(2) + sqrt( dy(2)^2 + R2 ) ) ) +...
         log( ( dy(1) + sqrt( dy(1)^2 + R4 ) ) / (dy(2) + sqrt( dy(2)^2 + R4 ) ) ) +...
        -log( ( dy(1) + sqrt( dy(1)^2 + R3 ) ) / (dy(2) + sqrt( dy(2)^2 + R3 ) ) );

% Tz(3) = -( Ty(1,2) + Tx(1,1) );
% Tz(2) = Ty(1,3);
% Tx(2) = Ty(1,1);
% Tz(1) = Tx(1,3);
