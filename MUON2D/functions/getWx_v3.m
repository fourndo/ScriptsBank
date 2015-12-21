function [Wx]=getWx_v3(mcell,nX,nZ,dX,dZ)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

% Wx = zeros(mcell,n*(m-1));
wx = zeros(nZ*(nX-1),mcell);
count=1;
for ii=1:nX
    
    for jj=1:nZ
        
        if (ii < nX)
            wx(count,count+nZ) = sqrt(dZ(jj)/dX(ii+1)); 
            wx(count,count) = -sqrt(dZ(jj)/dX(ii));

        end
        count=count+1;
    end
    
end

Wx=sparse(wx);
end