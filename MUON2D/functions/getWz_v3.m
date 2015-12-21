function [Wz]=getWz_v3(mcell,nX,nZ,dX,dZ)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1


wz = zeros(nX*(nZ-1),mcell);
count=1;
skip =0;

    for ii=1:nX
        for kk=1:nZ
            if (kk == nZ)
                
                skip=skip+1;
                
            else
                
                
                wz(count,count+skip) = -sqrt(dX(ii)/dZ(kk));
                wz(count,count+skip+1) = sqrt(dX(ii)/dZ(kk+1));
                count=count+1;

                
                
            end
        
        end
    end


Wz = sparse(wz);

end