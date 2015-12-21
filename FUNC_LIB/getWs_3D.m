function [ Ws , v ] = getWs_3D(dx,dy,dz,nullcell)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

nx = length(dx);
ny = length(dy);
nz = length(dz);

nullcell = reshape(nullcell,nz,nx,ny);

mcell = sum(nullcell(:));

v = zeros(mcell,1);
count=1;

for jj=1:ny
    
    for ii=1:nx
        
        for kk=1:nz
        
            if nullcell(kk,ii,jj)~= 0
                

                v(count) = sqrt( dx(ii) * dy(jj) * dz(kk) );

                count=count+1;

                
            end
            
        end
    end
end

Ws = speye(mcell);

end