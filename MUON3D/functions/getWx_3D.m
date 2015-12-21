function [Wx]=getWx_3D(mcell,dX,dY,dZ,mnull)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

nX = length(dX);
nY = length(dY);
nZ = length(dZ);

% Wx = zeros(mcell,n*(m-1));
Wx = sparse((nX-1)*nY*nZ,mcell);
count=1;
skip =0;


for jj=1:nY
    for ii=1:nX
        for kk=1:nZ
        
            if (ii == nX) %|| mnull(count+skip)== 0 || mnull(count+skip+nz)== 0
                
            skip=skip+1;
            else
                
                if mnull(count+skip)== 0 || mnull(count+skip+nZ)== 0
                    
                count=count+1;
                
                else
                dXhalf = (dX(ii+1) + dX(ii)) / 2;
                Wx(count,count+skip+nZ) = sqrt(1 / dXhalf^2); 
                Wx(count,count+skip) = -sqrt(1 / dXhalf^2);
                count=count+1;
                end
            end
            
        end
    end
end

% Wx = sparse(wx);

end