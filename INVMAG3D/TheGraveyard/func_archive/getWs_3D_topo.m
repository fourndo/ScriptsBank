function [Ws]=getWs_3D_topo(mcell,dX,dY,dZ,nullcell)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

nX = length(dX);
nY = length(dY);
nZ = length(dZ);

ws = zeros(mcell,1);

count=1;
for jj=1:nY
    for ii=1:nX
        for kk=1:nZ

                if nullcell(kk,ii,jj)~=0
                                                          
                    ws(count) = sqrt(dZ(kk)*dX(ii)*dY(jj)); 
                    count=count+1;        
                end
                
        end
    end
end

Ws = spdiags(ws,0,mcell,mcell);

end