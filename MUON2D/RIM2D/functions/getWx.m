function [Wx]=getWx(mcell,nX,nZ,dX,dZ)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

% Wx = zeros(mcell,n*(m-1));
% Wx = zeros(nZ*(nX-1)+1,mcell);
% count=2;
% for ii=1:nX
%     
%     for jj=1:nZ
%         
%         if (ii < nX)
%             Wx(count,count+nZ) = sqrt(dZ(jj)/dX(ii+1)); 
%             Wx(count,count) = -sqrt(dZ(jj)/dX(ii));
% 
%         end
%         count=count+1;
%     end
% end

Wx = sparse((nX-1)*nZ+1,mcell);
count=2;
skip =0;



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
            dXhalf = (dX(ii+1) + dX(ii)) / 2;
            Wx(count,count+skip+nZ) = sqrt(1 / dXhalf^2); 
            Wx(count,count+skip) = -sqrt(1 / dXhalf^2);
            count=count+1;
%                 end
        end

    end
end
