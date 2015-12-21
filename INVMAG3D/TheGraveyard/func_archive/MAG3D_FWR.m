function [ d, obsx, obsy, obsz ] = MAG3D_FWR( m, M, datafile, x0, y0, z0, dx, dy, dz)
% MAG3D_FWR(X0, Y0, Z0, dX, dY, dZ, X , Y , Z , H, I, D)
%
% NEED TO UPDATE HEADER
% Forward operator for the magnetic response from an orthogonal mesh 
% due to an inducing field H for a given location (single row of G)
% Reference: Blakely (1996) Potential Theory in Gravity and Magnetism
% Inputs:
% xo, yo, zo : Coordinates of the South-West-Top corner of the mesh
% dx, dy, dz : Vectors of cell sizes for the East, North and Vertical axis
% obsx,  obsy,  obsz  : Coordinates of the observation points 
%      H     : Magnitude of inducing field (nT)
%      I     : Inclinaison of H (degrees from horizontal)
%      D     : Declinaison of H (degrees from North)
% nullcell   : 1D vector for active(1) and inactive (0) cells
% (NOT FINISHED July 20, 2013)

addpath C:\Users\dominiquef\Dropbox\Master\Miscellaneous\
addpath functions

addpath C:\Users\dominiquef\Dropbox\Master\Miscellaneous\

[H, I, D, obsx, obsy, obsz, ~, ~] = read_MAG3D_obs(datafile);

nx = length(dx);
ny = length(dy);
nz = length(dz);

ndata = size( obsx , 1 );

% R0 = min( [min(dx) min(dy) min(dz)] ) / 4;

mcell = nx * ny * nz;
u0 = 4 * pi * 10^-7;


%Pre-allocate memory
% V=zeros(mcell,1);
G = zeros(1,mcell);
d = zeros(ndata,1);
% Define vector projection
Fx = cosd(I) * cosd(D-90);
Fy = cosd(I) * sind(D-90);
Fz = sind(I);

% Define orientation of magnetization.
% Assume that all fields are induced for now...
% Hx = cosd(I) * cosd(D-90);
% Hy = cosd(I) * sind(D-90);
% Hz = sind(I);

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
for dd = 1 : ndata
    
    count = 1;
    
    for jj = 1 : ny
        %First compute de location of the center of the cell
        Y = y0 + sum(dy(1:jj)) - dy(jj) /2;
        %Then the distance between the cell corners and the observation
        y(2) = (obsy(dd) - ( Y + dy(jj) /2 ) ) ;
        y(1) = (obsy(dd) - ( Y - dy(jj) /2 ) )  ;

       for ii = 1 : nx
            X = x0 + sum(dx(1:ii)) - dx(ii) /2;

            x(2) = (obsx(dd) - ( X + dx(ii) /2 ) ) ;
            x(1) = (obsx(dd) -  ( X - dx(ii) /2 ) ) ;

           for kk = 1 : nz            

                Z = z0 - sum(dz(1:kk)) + dz(kk) /2;

                z(1) = (obsz(dd) - ( Z + dz(kk) /2 ) ) ;
                z(2) = (obsz(dd) - ( Z - dz(kk) /2 ) ) ;

%                  if nullcell(kk,ii,jj) ~= 0

                    % Repeat computation for top and bottom
                    for zz= 1:2

                        h2 = z(zz)^2;
                        s1 = (-1)^(zz+1);

                        for xx = 1:2

                            x2 = x(xx)^2;

                            for yy = 1:2

                                s2 = (-1) ^ ( xx + yy );
                                r2 = ( x2 + y(yy) ^ 2 + h2 );
                                r = sqrt( r2 );  
                                rh = r * z(zz);

                                xy = x(xx) * y(yy);

                                arg1 = ( r - x(xx) ) / ( r + x(xx) );
                                arg2 = ( r - y(yy) ) / ( r + y(yy) );
                                arg3 = x2 + rh + h2;
                                arg4 = r2 + rh - x2;

                                tlog = ( Hy(count) * Fz + Hz(count) * Fy ) * log(arg1)/2 +...
                                    ( Hx(count) * Fz + Hz(count) * Fx ) * log(arg2)/2 -...
                                    (Hx(count) * Fy + Hy(count) * Fx) * log( r + z(zz) );

                                tatan = -(Hx(count) * Fx) * atan( xy / arg3 ) + ...
                                    -(Hy(count) * Fy) * atan( xy / arg4 ) +...
                                    (Hz(count) * Fz) * atan( xy / rh );

                                G(count) = G(count) + s1 * s2 *...
                                    ( tlog + tatan );

                            end

                        end

                    end

                    %Compute volume of cell
%                     V = ;

                    G(count) = G(count) * u0 * dx(ii) * dy(jj) * dz(kk) ;
                    % Distance to center of cell
%                     R= ((Obsx - X) ^ 2 + (Obsy - Y)^2 + (Obsz - Z)^2)^0.5;



                    % Compute the distance weighting
%                     Wr(count) = ( V(count) / (R + R0) ^ 3 ) ^ 2 ;

                    count = count + 1;

%                   end                    

           end

       end

    end
    
    % Compute forward data
    d(dd) = G * m;
    
end




