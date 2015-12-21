function [G,Wr,V] = forwardGrav(nX, nY, nZ, X0, Y0, Z0, dX, dY, dZ, ObsX,ObsY,ObsZ)
%Create forward operator

NewtG=6.6738e-3;

R0 = min( [min(dX) min(dY) min(dZ)] ) / 4;
mcell = nX * nY * nZ;

%Pre-allocate for weighting matrix
Wr=zeros(mcell,1);
V=zeros(mcell,1);

%Pre-allocate space for one row of forward operator
G=zeros(1,nX * nY * nZ);

count = 1;

for jj = 1 : nY
    %First compute de location of the center of the cell
    Y = Y0 + sum(dY(1:jj)) - dY(jj) /2;
    %Then the distance between the cell and the observation
    dy(1) = ( ObsY - Y + dY(jj) /2 ) ;
    dy(2) = ( ObsY - Y - dY(jj) /2 )  ;
    
   for ii = 1 : nX
        X = X0 + sum(dX(1:ii)) - dX(ii) /2;
% dX = ( X - ObsX ) ^2;

        dx(1) = ( ObsX - X + dX(jj) /2 )  ;
        dx(2) = ( ObsX - X - dX(jj) /2 )  ;
        
       for kk = 1: nZ
            Z = Z0 - sum(dZ(1:kk)) + dZ(kk) /2;
% dZ = ( Z - ObsZ ) ^2;
            dz(1) = ( ObsZ - Z - dZ(kk) /2 )  ;
            dz(2) = ( ObsZ - Z + dZ(kk) /2 )  ;
            
            % Compute contribution from each corners
            for aa= 1:2
                for bb= 1:2
                    for cc= 1:2
                        r = (dx(aa) ^ 2 + dy(bb) ^ 2 + dz(cc) ^ 2) ^ (0.50);
                        
                       G(count) = G(count) - NewtG * ...
                           (-1) ^ aa * (-1) ^ bb * (-1) ^ cc * ...
                           (dx(aa) * log ( dy(bb) + r ) + ...
                           dy(bb) * log ( dx(aa) + r ) - ...
                           dz(cc) * atan ( dx(aa) * dy(bb) / ( dz(cc) * r )));
                    end
                end
            end
%             G(count)=-G(count);
            R= ((ObsX - X) ^ 2 + (ObsY - Y)^2 + (ObsZ - Z)^2)^0.5;
            
            %Compute volume of cell
            V(count) = dX(ii)*dY(jj)*dZ(kk);
            
            % Compute the distance weighting
            Wr(count) = ( V(count) / (R + R0) ^ 2 ) ^ 2 ;
                        
            % Compute the forward operator
            
        count = count + 1;
       end
   end
end

clear model;
end