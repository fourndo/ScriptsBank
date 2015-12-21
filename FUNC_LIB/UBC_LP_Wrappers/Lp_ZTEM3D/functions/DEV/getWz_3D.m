function [Wz]=getWz_3D(mcell,dX,dY,dZ,mnull)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

nX = length(dX);
nY = length(dY);
nZ = length(dZ);

Wz = sparse(nX*(nZ-1)*nY,mcell);
count=1;
skip =0;
for jj=1:nY
    for ii=1:nX
        for kk=1:nZ
            if (kk == nZ) %|| mnull(kk,ii,jj)== 0 || mnull(kk+1,ii,jj)== 0
                
                skip=skip+1;
                
                
            else
                                
                if mnull(kk,ii,jj)== 0 || mnull(kk+1,ii,jj)== 0
                    
                count=count+1;
                
                else   
                Wz(count,count+skip) = -sqrt(dX(ii)*dY(jj)/dZ(kk));
                Wz(count,count+skip+1) = sqrt(dX(ii)*dY(jj)/dZ(kk+1));
                count=count+1;
                               
                end
                
            end
        
        end
    end
end

% Wz = sparse(wz);

end
