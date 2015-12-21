function [Tx,Ty,Tz,wr] = TMAG3C( obsx, obsy, obsz, x0, y0, z0, dx, dy, dz)
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

nx = length(dx);
ny = length(dy);
nz = length(dz);

ndata = size( obsx , 1 );

R0 = min( [min(dx) min(dy) min(dz)] )/4;

mcell = nx*ny*nz;


% Pre-allocate to store fields
Txx = zeros(ndata,mcell);
Txy = zeros(ndata,mcell);
Txz = zeros(ndata,mcell);

Tyx = zeros(ndata,mcell);
Tyy = zeros(ndata,mcell);
Tyz = zeros(ndata,mcell);

Tzx = zeros(ndata,mcell);
Tzy = zeros(ndata,mcell);
Tzz = zeros(ndata,mcell);

wr = zeros(mcell,1);

progress = -1;
tic
for dd = 1 : ndata
    
    count = 1;
    nn = 1;
    
    for jj = 1 : ny
        
        %First compute de location of the center of the cell
        Y = y0 + sum(dy(1:jj)) - dy(jj) /2;
        
        %Then the distance between the cell corners and the observation
        y(1) = ( ( Y + dy(jj) /2 )-obsy(dd)   ) ;
        y(2) = ( ( Y - dy(jj) /2 )-obsy(dd)   ) ;

       for ii = 1 : nx
           
            X = x0 + sum(dx(1:ii)) - dx(ii) /2;

            x(1) = ( ( X + dx(ii) /2 )-obsx(dd)   ) ;
            x(2) = ( ( X - dx(ii) /2 )-obsx(dd)   ) ;

           for kk = 1 : nz            

                Z = z0 - sum(dz(1:kk)) + dz(kk) /2;

                z(1) = ( obsz(dd) - ( Z + dz(kk) /2 ) ) ;
                z(2) = ( obsz(dd) - ( Z - dz(kk) /2 ) ) ;
                
                R1 = ( y(1)^2 + x(1)^2 );
                R2 = ( y(1)^2 + x(2)^2 );
                R3 = ( y(2)^2 + x(1)^2 );
                R4 = ( y(2)^2 + x(2)^2 );
                
                arg1 = sqrt( z(2)^2 + R2 );
                arg2 = sqrt( z(2)^2 + R1 );
                arg3 = sqrt( z(1)^2 + R1 );
                arg4 = sqrt( z(1)^2 + R2 );
                arg5 = sqrt( z(2)^2 + R3 );
                arg6 = sqrt( z(2)^2 + R4 );
                arg7 = sqrt( z(1)^2 + R4 );
                arg8 = sqrt( z(1)^2 + R3 );
                
                Ty(dd,mcell+nn) = atan( x(2) * z(2) / ( y(1) * arg1 ) ) +...
                    - atan( x(1) * z(2) / ( y(1) * arg2 ) ) +...
                      atan( x(1) * z(1) / ( y(1) * arg3 ) ) +...
                    - atan( x(2) * z(1) / ( y(1) * arg4 ) ) +...
                      atan( x(1) * z(2) / ( y(2) * arg5 ) ) +...
                    - atan( x(2) * z(2) / ( y(2) * arg6 ) ) +...
                      atan( x(2) * z(1) / ( y(2) * arg7 ) ) +...
                    - atan( x(1) * z(1) / ( y(2) * arg8 ) );
                
                Ty(dd,nn) = log( ( z(2) + arg2 ) / (z(1) + arg3 ) ) +...
                        -log( ( z(2) + arg1 ) / (z(1) + arg4 ) ) +...
                         log( ( z(2) + arg6 ) / (z(1) + arg7 ) ) +...
                        -log( ( z(2) + arg5 ) / (z(1) + arg8 ) );
                                
                
                Tx(dd,nn) = atan( y(2) * z(2) / ( x(1) * arg5 ) ) +...
                    - atan( y(1) * z(2) / ( x(1) * arg2 ) ) +...
                      atan( y(1) * z(1) / ( x(1) * arg3 ) ) +...
                    - atan( y(2) * z(1) / ( x(1) * arg8 ) ) +...
                      atan( y(1) * z(2) / ( x(2) * arg1 ) ) +...
                    - atan( y(2) * z(2) / ( x(2) * arg6 ) ) +...
                      atan( y(2) * z(1) / ( x(2) * arg7 ) ) +...
                    - atan( y(1) * z(1) / ( x(2) * arg4 ) );
                                
                R1 = (y(1)^2 + z(1)^2);
                R2 = (y(1)^2 + z(2)^2);
                R3 = (y(2)^2 + z(1)^2);
                R4 = (y(2)^2 + z(2)^2);
                
                Ty(dd,2*mcell+nn) = log( ( x(2) + sqrt( x(2)^2 + R1 ) ) / (x(1) + sqrt( x(1)^2 + R1 ) ) ) +...
                        -log( ( x(2) + sqrt( x(2)^2 + R2 ) ) / (x(1) + sqrt( x(1)^2 + R2 ) ) ) +...
                         log( ( x(2) + sqrt( x(2)^2 + R4 ) ) / (x(1) + sqrt( x(1)^2 + R4 ) ) ) +...
                        -log( ( x(2) + sqrt( x(2)^2 + R3 ) ) / (x(1) + sqrt( x(1)^2 + R3 ) ) );
                               
                R1 = (x(1)^2 + z(1)^2);
                R2 = (x(1)^2 + z(2)^2);
                R3 = (x(2)^2 + z(1)^2);
                R4 = (x(2)^2 + z(2)^2);
                
                Tx(dd,2*mcell+nn) = log( ( y(2) + sqrt( y(2)^2 + R1 ) ) / (y(1) + sqrt( y(1)^2 + R1 ) ) ) +...
                        -log( ( y(2) + sqrt( y(2)^2 + R2 ) ) / (y(1) + sqrt( y(1)^2 + R2 ) ) ) +...
                         log( ( y(2) + sqrt( y(2)^2 + R4 ) ) / (y(1) + sqrt( y(1)^2 + R4 ) ) ) +...
                        -log( ( y(2) + sqrt( y(2)^2 + R3 ) ) / (y(1) + sqrt( y(1)^2 + R3 ) ) );
                
                Tz(dd,2*mcell+nn) = -( Ty(dd,mcell+nn) + Tx(dd,nn) );
                Tz(dd,mcell+nn) = Ty(dd,2*mcell+nn);
                Tx(dd,mcell+nn) = Ty(dd,nn);
                Tz(dd,nn) = Tx(dd,2*mcell+nn);
                
                
%                 if dd==1
                % Compute volume of cell
                V = dx(ii) * dy(jj) * dz(kk) ;

                % Distance to center of cell
                R = ((obsx(dd) - X) ^ 2 + (obsy(dd) - Y)^2 + (obsz(dd) - Z)^2)^0.5;

                % Compute the distance weighting
                wr(nn) = wr(nn) + ( V / (R  + R0)^ 3  ) ^ 2 ;
%                 end 
                
                nn = nn+1;
                 
                
                count = count + 1;
                
           end

       end

    end
    
    d_iter = floor(dd/ndata*20);
    if  d_iter > progress
        
        fprintf('Computed %i pct of data in %8.5f sec\n',d_iter*5,toc)
        progress = d_iter;
        
        tic
         
    end
    
end

Tx = Tx/(4*pi);
Ty = Ty/(4*pi);
Tz = Tz/(4*pi);

fprintf('MAG Forward modelling completed\n');

