function [Ws]=getWs_3D(mcell,dX,dY,dZ)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

nX = length(dX);
nY = length(dY);
nZ = length(dZ);

% Wx = zeros(mcell,n*(m-1));
Ws = speye(mcell,mcell);

count=1;
for jj=1:nY
    for ii=1:nX
        for kk=1:nZ

%                 if mnull(count)==0
%                     
%                     Ws(count)=0;
%                     
%                 else
                    
                    Ws(count,count) = sqrt(dZ(kk)*dX(ii)*dY(jj)); 
                             
%                 end
              count=count+1;  
        end
    end
end



end