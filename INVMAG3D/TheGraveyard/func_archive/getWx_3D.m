function [Wx,Vx]=getWx_3D(mcell,dX,dY,dZ,nullcell)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

nX = length(dX);
nY = length(dY);
nZ = length(dZ);

mcellx = mcell-sum(sum(nullcell(:,end,:)~=0));
% Wx = zeros(mcell,n*(m-1));
Wx = sparse(mcellx,mcell);
vx = zeros(mcellx,1);
count=1;
skip =0;


for jj=1:nY
    for ii=1:nX
        for kk=1:nZ
        
            if (ii ~= nX) && nullcell(kk,ii,jj)~= 0
                
                    
                dXhalf = (dX(ii+1) + dX(ii)) / 2;
                
                ind1 = sub2ind(size(nullcell),kk,ii,jj);
                ind2 = sub2ind(size(nullcell),kk,ii+1,jj);  
                
                Wx(count,nullcell(ind1)) = sqrt(1 / dXhalf^2); 
                Wx(count,nullcell(ind2)) = -sqrt(1 / dXhalf^2);
                
                vx(count) = sqrt( dY(jj) * dZ(kk) * dXhalf );
                
                count=count+1;
                
                
                
            end
            
        end
    end
end

Vx = spdiags(vx,0,mcellx,mcellx);
% Wx = sparse(wx);

end