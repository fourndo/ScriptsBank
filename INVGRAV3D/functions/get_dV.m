function [dV] = get_dV(nX, nY, nZ, dX, dY, dZ)

dV = zeros(nX*nY*nZ,1);

count=1;
skip = 0;
%% Compute the derivative terms ( WxtWx + WytWy + WztWz ) * m
for jj = 1 : nY
    
    for ii = 1 : nX
        
        for kk = 1 : nZ
            
            dV(count) = sqrt( dZ(kk) * dY(jj) * dX(ii) );
            count=count+1;
                      
        end
    end
end

dV = spdiags(dV,0,nX*nY*nZ,nX*nY*nZ);