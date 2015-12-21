function [Wz,Vz]=getWz_3D(mcell,dX,dY,dZ,nullcell)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

nX = length(dX);
nY = length(dY);
nZ = length(dZ);

mcellz = mcell-sum(sum(nullcell(end,:,:)~=0));

Wz = sparse(mcellz,mcell);
vz = zeros(mcellz,1);
count=1;

for jj=1:nY
    for ii=1:nX
        for kk=1:nZ
            
            if (kk ~= nZ) && nullcell(kk,ii,jj)~= 0
                
                    
                dZhalf = (dZ(kk+1) + dZ(kk)) / 2;
                
                ind1 = sub2ind(size(nullcell),kk,ii,jj);
                ind2 = sub2ind(size(nullcell),kk+1,ii,jj);  
                
                Wz(count,nullcell(ind1)) = sqrt(1 / dZhalf^2); 
                Wz(count,nullcell(ind2)) = -sqrt(1 / dZhalf^2);
                
                vz(count) = sqrt( dX(ii) * dY(jj) * dZhalf );
                
                count=count+1;
                                
            end

        end
    end
end

Vz = spdiags(vz,0,mcellz,mcellz);

% Wz = sparse(wz);

end
