function [G,Wr,V] = Fwr_Mag_Integral_DEV(mcell, x0, y0, z0, dx, dy, dz,...
                           Obsx, Obsy, Obsz, H, I, D, nullcell)
% FwdMag(X0, Y0, Z0, dX, dY, dZ, X , Y , Z , H, I, D)
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

nx =length(dx);
ny =length(dy);
nz =length(dz);

R0 = min( [min(dx) min(dy) min(dz)] ) / 4;
% mcell = nX * nY * nZ;
% uo = 4 * pi * 10^-7;

%Pre-allocate memory
Wr=zeros(mcell,1);
V=zeros(mcell,1);
G=zeros(1,mcell);

% Unit vector
unit = 1 / sqrt(3);

% Define vector projection
Fx = cosd(I) * cosd(D) * unit;
Fy = cosd(I) * sind(D) * unit;
Fz = sind(I) * unit;

% Define orientation of magnetization.
% Assume that all fields are induced for now...
Hx = cosd(I) * cosd(D) * unit;
Hy = cosd(I) * sind(D) * unit;
Hz = sind(I) * unit;

a(1,2) = Hx * Fy + Hy * Fx;
a(1,3) = ( Hx * Fz + Hz * Fx ) / 2;
a(2,3) = ( Hy * Fz + Hz * Fy ) / 2;

count = 1;
for jj = 1 : ny
    %First compute de location of the center of the cell
    Y = y0 + sum(dy(1:jj)) - dy(jj) /2;
    %Then the distance between the cell corners and the observation
    y(1) = ( Obsy - Y + dy(jj) /2 ) ;
    y(2) = ( Obsy - Y - dy(jj) /2 )  ;
    
   for ii = 1 : nx
        X = x0 + sum(dx(1:ii)) - dx(ii) /2;

        x(1) = ( Obsx - X + dx(jj) /2 )  ;
        x(2) = ( Obsx - X - dx(jj) /2 )  ;
        
       for kk = 1 : nz            
                
            Z = z0 - sum(dz(1:kk)) + dz(kk) /2;

            z(1) = ( Obsz - Z - dz(kk) /2 )  ;
            z(2) = ( Obsz - Z + dz(kk) /2 )  ;
                
             if nullcell(kk,ii,jj) ~= 0
                 
                % Compute distance to each corners
                for zz= [2 1]
                    
                    for xx = [2 1]
                        
                        for yy = [2 1]
                            
                            s = (-1)^(xx+yy+zz);
                            r = (x(xx) ^ 2 + y(yy) ^ 2 + z(zz) ^ 2) ^ (0.50);
                                                       
                            G(count) = s * a(2,3) * log( ( r - x(xx) ) / ( r + x(xx) ) ) +...
                                s * a(1,3) * log( ( r - y(yy) ) / ( r + y(yy) ) ) +...
                                -s * a(1,2) * log( r + z(zz)) +...
                                -s * Hx*Fx * atan( x(xx) * y(yy) / ( x(xx)^2 + r * z(zz) + z(zz)^2 ) ) +...
                                -s * Hy*Fy * atan( x(xx) * y(yy) / ( -x(xx)^2 + r * z(zz) + r^2 ) ) +...
                                s * Hz*Fz * atan( x(xx) * y(yy) / ( r * z(zz) ) ) ;
                            
                        end
                        
                    end
                    
                end
                
                G(count) = G(count) * H * 1e-7 * 1e+9;
                % Distance to center of cell
                R= ((Obsx - X) ^ 2 + (Obsy - Y)^2 + (Obsz - Z)^2)^0.5;
            
                %Compute volume of cell
                V(count) = dx(ii)*dy(jj)*dz(kk);

                % Compute the distance weighting
                Wr(count) = ( V(count) / (R + R0) ^ 2 ) ^ 2 ;
                
                count = count + 1;
                
              end                    
        
       end
       
   end
   
end

clear model;
end