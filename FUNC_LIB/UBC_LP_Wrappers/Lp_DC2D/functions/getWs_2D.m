function [Ws]=getWs_2D(mcell,dX,dZ,mnull)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

nX = length(dX);
nZ = length(dZ);

% Wx = zeros(mcell,n*(m-1));
Ws = zeros(mcell,1);

count=1;
    for ii=1:nX
        for kk=1:nZ

                if mnull(count)==0
                    
                    Ws(count)=0;
                    
                else
                    
                    Ws(count) = sqrt(dZ(kk)*dX(ii)); 
                             
                end
              count=count+1;  
        end
    end


