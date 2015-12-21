function [mnull]=getnull(z0,dX,dY,dZ,ZZ)
% Compute cells volume

nX = length(dX);
nY = length(dY);
nZ = length(dZ);


% volume=zeros(nX*nY*nZ,1);
% 
% count=1;
% for ii=1:nY
%     
%     for jj=1:nX
%         
%         for kk=1:nZ
%             
%             volume(count)= dY(ii) * dX(jj) * dZ(kk);
%             
%             count= count+1;
%         end
%         
%     end
%     
% end

%% Assign a null value (0) to cells above topo

mnull=ones(nZ,nX,nY);
z = z0;


for ii=1:nX
    
    for jj= 1:nY
        count=1;
        z=z0;
        
        while z > ZZ(ii,jj)

            mnull(count,ii,jj)= 0;
            z = z - dZ(count);
            count=count+1;
            
        end
    end
    
end



end