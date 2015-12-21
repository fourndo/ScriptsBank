function [ Bx, By, Bz, TMI, magB, obsx, obsy, obsz ] = FMAG3C( m, M, H, D, I, nullcell, obsx, obsy, obsz, x0, y0, z0, dx, dy, dz)
% MAG3D_FWR(X0, Y0, Z0, dX, dY, dZ, X , Y , Z , H, I, D)
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
% (NOT FINISHED July 20, 2013)

% [H, ~, ~, obsx, obsy, obsz, ~, ~] = read_MAG3D_obs(datafile);

% Define vector projection for TMI
Fx = cosd(I) * cosd(D);
Fy = cosd(I) * sind(D);
Fz = sind(I);

nx = length(dx);
ny = length(dy);
nz = length(dz);

ndata = size( obsx , 1 );

% R0 = min( [min(dx) min(dy) min(dz)] ) / 4;

mcell = nx * ny * nz;
u0 = 4 * pi * 10^-7;

% Pre-allocate to store fields
Bx = zeros(ndata,1);
By = zeros(ndata,1);
Bz = zeros(ndata,1);
TMI = zeros(ndata,1);
magB = zeros(ndata,1);

Hx = H * M(:,1);
Hy = H * M(:,2);
Hz = H * M(:,3);

% Used for uniform induced magnetization
% a1 = Hx * Fy + Hy * Fx;
% a2 = ( Hx * Fz + Hz * Fx );
% a3 = ( Hy * Fz + Hz * Fy );
% a4 = Hx * Fx;
% a5 = Hy * Fy;
% a6 = Hz * Fz;
progress = -1;
tic
for dd = 1 : ndata
    
    %Pre-allocate memory
    % V=zeros(mcell,1);
    Gx = zeros(1,mcell);
    Gy = zeros(1,mcell);
    Gz = zeros(1,mcell);
    
    count = 1;
    
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

                %If not an air cell
                if nullcell(count) == 1
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
                
                T(2,2) = atan( x(2) * z(2) / ( y(1) * arg1 ) ) +...
                    - atan( x(1) * z(2) / ( y(1) * arg2 ) ) +...
                      atan( x(1) * z(1) / ( y(1) * arg3 ) ) +...
                    - atan( x(2) * z(1) / ( y(1) * arg4 ) ) +...
                      atan( x(1) * z(2) / ( y(2) * arg5 ) ) +...
                    - atan( x(2) * z(2) / ( y(2) * arg6 ) ) +...
                      atan( x(2) * z(1) / ( y(2) * arg7 ) ) +...
                    - atan( x(1) * z(1) / ( y(2) * arg8 ) );
                
                T(2,1) = log( ( z(2) + arg2 ) / (z(1) + arg3 ) ) +...
                        -log( ( z(2) + arg1 ) / (z(1) + arg4 ) ) +...
                         log( ( z(2) + arg6 ) / (z(1) + arg7 ) ) +...
                        -log( ( z(2) + arg5 ) / (z(1) + arg8 ) );
                                
                
                T(1,1) = atan( y(2) * z(2) / ( x(1) * arg5 ) ) +...
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
                
                T(2,3) = log( ( x(2) + sqrt( x(2)^2 + R1 ) ) / (x(1) + sqrt( x(1)^2 + R1 ) ) ) +...
                        -log( ( x(2) + sqrt( x(2)^2 + R2 ) ) / (x(1) + sqrt( x(1)^2 + R2 ) ) ) +...
                         log( ( x(2) + sqrt( x(2)^2 + R4 ) ) / (x(1) + sqrt( x(1)^2 + R4 ) ) ) +...
                        -log( ( x(2) + sqrt( x(2)^2 + R3 ) ) / (x(1) + sqrt( x(1)^2 + R3 ) ) );
                               
                R1 = (x(1)^2 + z(1)^2);
                R2 = (x(1)^2 + z(2)^2);
                R3 = (x(2)^2 + z(1)^2);
                R4 = (x(2)^2 + z(2)^2);
                
                T(1,3) = log( ( y(2) + sqrt( y(2)^2 + R1 ) ) / (y(1) + sqrt( y(1)^2 + R1 ) ) ) +...
                        -log( ( y(2) + sqrt( y(2)^2 + R2 ) ) / (y(1) + sqrt( y(1)^2 + R2 ) ) ) +...
                         log( ( y(2) + sqrt( y(2)^2 + R4 ) ) / (y(1) + sqrt( y(1)^2 + R4 ) ) ) +...
                        -log( ( y(2) + sqrt( y(2)^2 + R3 ) ) / (y(1) + sqrt( y(1)^2 + R3 ) ) );
                
                T(3,3) = -( T(2,2) + T(1,1) );
                T(3,2) = T(2,3);
                T(1,2) = T(2,1);
                T(3,1) = T(1,3);
                
                T = T/(4*pi);
                    % Repeat computation for top and bottom
                   

                    %Compute volume of cell
%                     V = ;

                    Gx(count) = (Hx(count)* T(1,1) + Hy(count)* T(1,2) + Hz(count)* T(1,3) ) ;%*...
%                         u0 * dx(ii) * dy(jj) * dz(kk) ;
                    
                    Gy(count) = (Hx(count)* T(2,1) + Hy(count)* T(2,2) + Hz(count)* T(2,3) ) ;%*...
%                         u0 * dx(ii) * dy(jj) * dz(kk) ;
                    
                    Gz(count) = (Hx(count)* T(3,1) + Hy(count)* T(3,2) + Hz(count)* T(3,3) ) ;%*...
%                         u0 * dx(ii) * dy(jj) * dz(kk) ;
                    % Distance to center of cell
%                     R= ((Obsx - X) ^ 2 + (Obsy - Y)^2 + (Obsz - Z)^2)^0.5;

                    % Compute the distance weighting
%                     Wr(count) = ( V(count) / (R + R0) ^ 3 ) ^ 2 ;

                end
               
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
    
    % Compute forward data
    Bx(dd) = Gx * m;
    By(dd) = Gy * m;
    Bz(dd) = Gz * m;
    
    % Compute TMI
    TMI(dd) = Bx(dd)*Fx + By(dd)*Fy + Bz(dd)*Fz;
    
    % Compute magnitude of field
    magB(dd) = sqrt( Bx(dd)^2 + By(dd)^2 + Bz(dd)^2 ) ;
    
end

fprintf('MAG Forward modelling completed\n');

