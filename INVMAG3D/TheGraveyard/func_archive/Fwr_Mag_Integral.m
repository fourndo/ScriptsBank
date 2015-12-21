function [G,Wr,V] = Fwr_Mag_Integral(mcell, X0, Y0, Z0, dX, dY, dZ,...
                           ObsX, ObsY, ObsZ, H, I, D, nullcell)
% FwdMag(X0, Y0, Z0, dX, dY, dZ, X , Y , Z , H, I, D)
% Forward operator for the magnetic response from an orthogonal mesh 
% due to an inducing field H for a given location (single row of G)
% Reference: Rao & Babu (1993) Computer & Geosciences 19,6 : 781-801
% Inputs:
% Xo, Yo, Zo : Coordinates of the South-West-Top corner of the mesh
% dX, dY, dZ : Vectors of cell sizes for the East, North and Vertical axis
% X,  Y,  Z  : Coordinates of the observation points 
%      H     : Magnitude of inducing field (nT)
%      I     : Inclinaison of H (degrees from horizontal)
%      D     : Declinaison of H (degrees from North)
% nullcell   : 1D vector for active(1) and inactive (0) cells
% (NOT FINISHED July 20, 2013)

nX =length(dX);
nY =length(dY);
nZ =length(dZ);

R0 = min( [min(dX) min(dY) min(dZ)] ) / 4;
% mcell = nX * nY * nZ;
uo = 4 * pi * 10^-7;

%Pre-allocate memory
Wr=zeros(mcell,1);
V=zeros(mcell,1);
G=zeros(1,mcell);

count = 1;
p=cosd(I)*cosd(D);
q=cosd(I)*sind(D);
r=sind(I);

L=cosd(I)*cosd(D);
M=cosd(I)*sind(D);
N=sind(I);

G1=H*(M*r+N*q);
G2=H*(L*r+N*p);
G3=H*(L*q+M*p);
G4=H*(N*r-M*q);
G5=H*(N*r-L*p);

for jj = 1 : nY
    %First compute de location of the center of the cell
    Y = Y0 + sum(dY(1:jj)) - dY(jj) /2;
    %Then the distance between the cell corners and the observation
    dy(1) = ( ObsY - Y + dY(jj) /2 ) ;
    dy(2) = ( ObsY - Y - dY(jj) /2 )  ;
    
   for ii = 1 : nX
        X = X0 + sum(dX(1:ii)) - dX(ii) /2;

        dx(1) = ( ObsX - X + dX(jj) /2 )  ;
        dx(2) = ( ObsX - X - dX(jj) /2 )  ;
        
       for kk = 1 : nZ            
                
            Z = Z0 - sum(dZ(1:kk)) + dZ(kk) /2;

            dz(1) = ( ObsZ - Z - dZ(kk) /2 )  ;
            dz(2) = ( ObsZ - Z + dZ(kk) /2 )  ;
                
             if nullcell(kk,ii,jj) ~= 0
                % Compute distance to each corners
                counter = 1;
                for bb= 1:2
                    for aa= 1:2
                        for cc= 1:2
                            r(counter) = (dx(aa) ^ 2 + dy(bb) ^ 2 + dz(cc) ^ 2) ^ (0.50);
                            counter = counter+1;
                        end
                    end
                end

                F(1) = log( ( (r(2) + dx(1)) * (r(3) + dx(2)) *...
                    (r(5) + dx(1)) * (r(8) + dx(2)) ) /...
                    ( (r(1) + dx(1)) * (r(4) + dx(2)) *...
                    (r(6) + dx(1)) * (r(7) + dx(2)) ) ) ;

                F(2) = log( ( (r(2) + dy(1)) * (r(3) + dy(1)) *...
                    (r(5) + dy(2)) * (r(8) + dy(2)) ) /...
                    ( (r(1) + dy(1)) * (r(4) + dy(1)) *...
                    (r(6) + dy(2)) * (r(7) + dy(2)) ) ) ;

                F(3) = log( ( (r(2) + dz(2)) * (r(3) + dz(1)) *...
                    (r(5) + dz(1)) * (r(8) + dz(2)) ) /...
                    ( (r(1) + dz(1)) * (r(4) + dz(2)) *...
                    (r(6) + dz(2)) * (r(7) + dz(1)) ) ) ;

                F(4) = atan(dx(2) * dz(2) / r(8) / dy(2) ) -...
                    atan(dx(1) * dz(2) / r(6) / dy(2) ) -...
                    atan(dx(2) * dz(2) / r(4) / dy(1) ) +...
                    atan(dx(1) * dz(2) / r(2) / dy(1) ) -...
                    atan(dx(2) * dz(1) / r(7) / dy(2) ) +...
                    atan(dx(1) * dz(1) / r(5) / dy(2) ) +...
                    atan(dx(2) * dz(1) / r(3) / dy(1) ) -...
                    atan(dx(1) * dz(1) / r(1) / dy(1) );

                F(5) = atan(dy(2) * dz(2) / r(8) / dx(2) ) -...
                    atan(dy(2) * dz(2) / r(6) / dx(1) ) -...
                    atan(dy(1) * dz(2) / r(4) / dx(2) ) +...
                    atan(dy(1) * dz(2) / r(2) / dx(1) ) -...
                    atan(dy(2) * dz(1) / r(7) / dx(2) ) +...
                    atan(dy(2) * dz(1) / r(5) / dx(1) ) +...
                    atan(dy(1) * dz(1) / r(3) / dx(2) ) -...
                    atan(dy(1) * dz(1) / r(1) / dx(1) );

                G(count) = G1*F(1) + G2*F(2) + G3*F(3) + G4*F(4) + G5*F(5);
                
 
                
                % Distance to center of cell
                R= ((ObsX - X) ^ 2 + (ObsY - Y)^2 + (ObsZ - Z)^2)^0.5;
            
                %Compute volume of cell
                V(count) = dX(ii)*dY(jj)*dZ(kk);

                % Compute the distance weighting
                Wr(count) = ( V(count) / (R + R0) ^ 2 ) ^ 2 ;
                
                count = count + 1;
                
              end                    
        
       end
       
   end
   
end

clear model;
end