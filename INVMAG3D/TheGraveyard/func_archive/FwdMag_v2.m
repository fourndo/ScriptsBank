function [H,Wr,V] = FwdMag(mcell, X0, Y0, Z0, dX, dY, dZ,...
                           ObsX, ObsY, ObsZ, H0, I, D, nullcell)
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
uo = 1;%4 * pi * 10^-7;

%Pre-allocate memory
Wr=zeros(mcell,1);
V=zeros(mcell,1);
H=zeros(1,mcell);

count = 1;
p=cosd(I)*cosd(D);
q=cosd(I)*sind(D);
r=sind(I);

J = zeros(3,1);
J(1)=uo*H0*(p);
J(2)=uo*H0*(q);
J(3)=uo*H0*(r);

% Pre-defined function
R1 = @(r1,r2) r1(1)^2 + r2(1)^2;
R2 = @(r1,r2) r1(1)^2 + r2(2)^2;
R3 = @(r1,r2) r1(2)^2 + r2(1)^2;
R4 = @(r1,r2) r1(2)^2 + r2(2)^2;
Ln = @(dk,R,r1,r2) log( (dk(2) + (dk(2)^2 + R(r1,r2))^0.5) / (dk(1) + (dk(1)^2 + R(r1,r2))^0.5) );
Atn = @(di,dj,dk,i1,i2,i3,i4,R,r1,r2) atan( (dj(i1) * dk(i1) ) / ( di(i3) * ( dk(i4)^2 + R(r1,r2) )^0.5 ) );

T = zeros(3,3);
for jj = 1 : nY
    %First compute de location of the center of the cell
    Y = Y0 + sum(dY(1:jj)) - dY(jj) /2;
    %Then the distance between the cell corners and the observation
    dy(1) = ( ObsY - (Y + dY(jj) /2) ) ;
    dy(2) = ( ObsY - (Y - dY(jj) /2) )  ;
    
   for ii = 1 : nX
        X = X0 + sum(dX(1:ii)) - dX(ii) /2;

        dx(1) = ( ObsX - (X + dX(jj) /2) )  ;
        dx(2) = ( ObsX - (X - dX(jj) /2) )  ;
        
       for kk = 1 : nZ            
                
            Z = Z0 - sum(dZ(1:kk)) + dZ(kk) /2;

            dz(1) = ( ObsZ - (Z - dZ(kk) /2) )  ;
            dz(2) = ( ObsZ - (Z + dZ(kk) /2) )  ;
                
             if nullcell(kk,ii,jj) ~= 0
                % Compute distance to each corners
                T(2,2) = Atn(dx,dy,dz,2,2,1,2,R2,dx,dy) -...
                    Atn(dx,dy,dz,1,2,1,2,R1,dx,dy) +...
                    Atn(dx,dy,dz,1,1,1,1,R1,dx,dy) -...
                    Atn(dx,dy,dz,2,1,1,1,R2,dx,dy) +...
                    Atn(dx,dy,dz,1,2,2,2,R3,dx,dy) -...
                    Atn(dx,dy,dz,2,2,2,2,R4,dx,dy) +...
                    Atn(dx,dy,dz,2,1,2,1,R4,dx,dy) -...
                    Atn(dx,dy,dz,1,1,2,1,R3,dx,dy);

                T(1,1) = Atn(dy,dx,dz,2,2,1,2,R2,dy,dx) -...
                    Atn(dy,dx,dz,1,2,1,2,R1,dy,dx) +...
                    Atn(dy,dx,dz,1,1,1,1,R1,dy,dx) -...
                    Atn(dy,dx,dz,2,1,1,1,R2,dy,dx) +...
                    Atn(dy,dx,dz,1,2,2,2,R3,dy,dx) -...
                    Atn(dy,dx,dz,2,2,2,2,R4,dy,dx) +...
                    Atn(dy,dx,dz,2,1,2,1,R4,dy,dx) -...
                    Atn(dy,dx,dz,1,1,2,1,R3,dy,dx);
               
                
                T(1,2) = Ln(dz,R1,dx,dy) -...
                    Ln(dz,R2,dx,dy) +...
                    Ln(dz,R4,dx,dy) -...
                    Ln(dz,R3,dx,dy);
                
                T(1,3) = Ln(dy,R1,dx,dz) -...
                    Ln(dy,R2,dx,dz) +...
                    Ln(dy,R4,dx,dz) -...
                    Ln(dy,R3,dx,dz);
                
                T(2,3) = Ln(dx,R1,dy,dz) -...
                    Ln(dx,R2,dy,dz) +...
                    Ln(dx,R4,dy,dz) -...
                    Ln(dx,R3,dy,dz);
                
                T(3,2) = T(2,3);
                
                T(2,1) = T(1,2);
                
                T(3,1) = T(1,3);
                
                T(3,3) = -( T(2,2) + T(1,1) );
                
                
                Hm = (T * J);
                
                H(count) = Hm(1);
 
                
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