function [Ws]=getWs(mcell,nX,nZ,dX,dZ)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1


Ws = sparse(mcell,mcell);
count=2;
for ii=1:nX
    
    for jj=1:nZ
        
    Ws(count,count)= sqrt(dX(ii)*dZ(jj));
    count=count+1;
    
    end
    
end