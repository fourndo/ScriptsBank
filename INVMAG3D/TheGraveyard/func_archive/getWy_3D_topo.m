function [Wy,Vy]=getWy_3D_topo(mcell,dX,dY,dZ,nullcell)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

nX = length(dX);
nY = length(dY);
nZ = length(dZ);

% Wx = zeros(mcell,n*(m-1));
Wy = sparse(mcell - nX*nZ,mcell);
vy = zeros(mcell - nX*nZ,1);
count=1;
for jj=1:nY-1
    for ii=1:nX
        for kk=1:nZ
        
                if nullcell(kk,ii,jj)== 0 
                    
                %count=count+1;
                
                elseif nullcell(kk,ii,jj+1)== 0

                count=count+1;
                
                else 
                    
                dYhalf = ( dY(jj) + dY(jj+1) ) /2;
                
                Wy(count,count+nZ*nX) = sqrt(1 / dYhalf^2 ); 
                Wy(count,count) = -sqrt(1 / dYhalf^2 );

                vy(count) = sqrt( dX(ii) * dZ(kk) * dYhalf );
                
                count=count+1;
                end
        end
    end
end

Vy = spdiags(vy,0,count-1,count-1);
Wy = Wy(1:count-1,:);
end