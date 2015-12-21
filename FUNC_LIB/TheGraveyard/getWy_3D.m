function [ Wy , Vy ] = getWy_3D(dx,dy,dz,nullcell)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

nx = length(dx);
ny = length(dy);
nz = length(dz);

index = cumsum(nullcell);
index = reshape(index,nz,nx,ny);
nullcell = reshape(nullcell,nz,nx,ny);

mcell = sum(nullcell(:));
mcelly = mcell-sum(sum(nullcell(:,:,end)));

Wy = spalloc(mcelly,mcell,2*mcell);
vy = zeros(mcelly,1);
count=1;

for jj=1:ny
    
    for ii=1:nx
        
        for kk=1:nz
        
            if (jj ~= ny) && nullcell(kk,ii,jj)~= 0
                
                if nullcell(kk,ii,jj+1)~= 0


                    dyhalf = (dy(jj+1) + dy(jj)) / 2;

                    Wy(count,index(kk,ii,jj)) = (1 / dy(jj)); 
                    Wy(count,index(kk,ii,jj+1)) = -( 1 / dy(jj+1));

                    vy(count) = sqrt( dx(ii) * dz(kk) * dy(jj) );
                    count=count+1;

                else

                    count=count+1;

                end
                
            end
            
        end
    end
end

Vy = spdiags(vy,0,mcelly,mcelly);
