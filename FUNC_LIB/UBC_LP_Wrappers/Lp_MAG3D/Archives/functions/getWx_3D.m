function [Wx]=getWx_3D(mcell,dX,dY,dZ,mnull)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

nX = length(dX);
nY = length(dY);
nZ = length(dZ);

Wx = zeros(mcell,(nX-1)*nY*nZ);
% Wx = sparse((nX-1)*nY*nZ,mcell);
count=1;
skip =0;


for jj=1:nY
    for ii=1:nX
        for kk=1:nZ
        
            if (ii == nX) %|| mnull(kk,ii,jj)== 0 || mnull(kk,ii+1,jj)== 0
                
            skip=skip+1;
            else
                
%                 if mnull(kk,ii,jj)== 0 || mnull(kk,ii+1,jj)== 0
%                     
%                 count=count+1;
%                 
%                 else  
                Wx(count,count+skip+nZ) = sqrt(dZ(kk)*dY(jj)/dX(ii+1)); 
                Wx(count,count+skip) = -sqrt(dZ(kk)*dY(jj)/dX(ii));
                count=count+1;
%                 end
            end
            
        end
    end
end

% Wx = sparse(wx);

end