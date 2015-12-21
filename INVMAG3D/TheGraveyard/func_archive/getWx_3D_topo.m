function [Wx,Vx]=getWx_3D_topo(mcell,dX,dY,dZ,nullcell)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

nX = length(dX);
nY = length(dY);
nZ = length(dZ);

% Wx = zeros(mcell,n*(m-1));
Wx = sparse(mcell - nY*nZ, mcell);
vx = zeros(mcell - nY*nZ,1);
count=1;
skip =0;


for jj=1:nY
    for ii=1:nX
        for kk=1:nZ
        
            if (ii == nX) 
                
            skip=skip+1;
            else
                
                if nullcell(kk,ii,jj)== 0 
                    
                %count=count+1;
                
                elseif nullcell(kk,ii+1,jj)== 0

                count=count+1;
                
                else
                    
                dXhalf = (dX(ii+1) + dX(ii)) / 2;
                
                Wx(count,count+skip+nZ) = sqrt(1 / dXhalf^2); 
                Wx(count,count+skip) = -sqrt(1 / dXhalf^2);
                
                vx(count) = sqrt( dY(jj) * dZ(kk) * dXhalf );
                
                count=count+1;
                
                end
                
            end
            
        end
    end
end

Vx = spdiags(vx,0,count-1,count-1);
Wx = Wx(1:count-1,:);

end