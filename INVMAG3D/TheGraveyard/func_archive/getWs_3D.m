function [Ws,V]=getWs_3D(mcell,dX,dY,dZ,nullcell)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

nX = length(dX);
nY = length(dY);
nZ = length(dZ);

% Wx = zeros(mcell,n*(m-1));
ws = zeros(mcell,1);
v = zeros(mcell,1);
count=1;
for jj=1:nY
    for ii=1:nX
        for kk=1:nZ

                if nullcell(count)==1
                                                          
                    ws(count) = sqrt(dZ(kk)*dX(ii)*dY(jj)); 
                    v(count) = dZ(kk)*dX(ii)*dY(jj);
                    
                end
              count=count+1;  
        end
    end
end

Ws = spdiags(ws,0,mcell,mcell);
V = spdiags(v,0,mcell,mcell);
end