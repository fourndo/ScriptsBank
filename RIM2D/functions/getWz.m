function [Wz]=getWz(mcell,nX,nZ,dX,dZ)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1


% Wz = zeros(mcell,mcell);
% count=2;
% for ii=1:nX
%     
%     for jj=1:nZ
%         
%         if (jj < nZ)
%             Wz(count,count+1) = sqrt(dZ(jj+1)/dX(ii)); 
%             Wz(count,count) = -sqrt(dZ(jj)/dX(ii));
%         end
%         
%         count=count+1;
%     end
% end

Wz = sparse(nX*(nZ-1)+1,mcell);
count=2;
skip =0;

for ii=1:nX
    for kk=1:nZ
        if (kk == nZ) %|| mnull(kk,ii,jj)== 0 || mnull(kk+1,ii,jj)== 0

            skip=skip+1;


        else

%                 if mnull(kk,ii,jj)== 0 || mnull(kk+1,ii,jj)== 0
%                     
%                 count=count+1;
%                 
%                 else   
            dZhalf = ( dZ(kk) + dZ(kk+1) ) / 2;
            Wz(count,count+skip) = - sqrt(1 / dZhalf^2 );
            Wz(count,count+skip+1) = sqrt(1 / dZhalf^2 );
            count=count+1;

%                 end

        end

    end
end
