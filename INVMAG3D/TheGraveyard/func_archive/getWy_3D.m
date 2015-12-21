function [Wy,Vy]=getWy_3D(mcell,dX,dY,dZ,nullcell)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

nX = length(dX);
nY = length(dY);
nZ = length(dZ);

mcelly = mcell-sum(sum(nullcell(:,:,end)~=0));

Wy = sparse(mcelly,mcell);
vy = zeros(mcelly,1);
count=1;
for jj=1:nY-1
    for ii=1:nX
        for kk=1:nZ
        
            if (jj ~= nY) && nullcell(kk,ii,jj)~= 0
                
                    
                dYhalf = (dY(jj+1) + dY(jj)) / 2;
                
                ind1 = sub2ind(size(nullcell),kk,ii,jj);
                ind2 = sub2ind(size(nullcell),kk,ii,jj+1);  
                
                Wy(count,nullcell(ind1)) = sqrt(1 / dYhalf^2); 
                Wy(count,nullcell(ind2)) = -sqrt(1 / dYhalf^2);
                
                vy(count) = sqrt( dX(ii) * dZ(kk) * dYhalf );
                
                count=count+1;
                
                
                
            end
        end
    end
end

Vy = spdiags(vy,0,mcelly,mcelly);

end