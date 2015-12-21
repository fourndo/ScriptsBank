function [Ws]=getWs_3D(mcell,dX,dY,dZ,mnull)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

nX = length(dX);
nY = length(dY);
nZ = length(dZ);

Ws = sparse(mcell,mcell);

count=1;
for jj=1:nY
    for ii=1:nX
        for kk=1:nZ

                if mnull(count)==0
                    
                    count=count+1;
                    
                else
                    
                    Ws(count,count) = sqrt(dZ(kk)*dX(ii)*dY(jj)); 
                    count = count+1;          
                end
                
        end
    end
end

% Wy = sparse(wy);

end