function [UBC_model] = make_UBC_model(model,nX,nY,nZ)
% Reshape model from Eldad's (x,y,z) to UBC (z,x,y)

UBC_model=zeros(nZ,nX,nY);

count =1;

for ii = 1:nZ
    
    for jj = 1:nY
         
        for kk = 1:nX
            
            
            UBC_model(ii,kk,jj) = model(count);
            count=count+1;
            
        end
        
    end
    
end

UBC_model=reshape(UBC_model,nX*nY*nZ,1);