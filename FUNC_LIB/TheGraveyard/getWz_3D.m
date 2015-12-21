function [ Wz , Vz ] = getWz_3D(dx,dy,dz,nullcell)
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
mcellz = mcell-sum(sum(nullcell(end,:,:)));

Wz = spalloc(mcellz,mcell,2*mcell);
vz = zeros(mcellz,1);
count=1;

for jj=1:ny
    
    for ii=1:nx
        
        for kk=1:nz
        
            if (kk ~= nz) && nullcell(kk,ii,jj)~= 0
                
                if nullcell(kk+1,ii,jj)~= 0


                    dzhalf = (dz(kk+1) + dz(kk)) / 2;

                    Wz(count,index(kk,ii,jj)) = ( 1  / dz(kk)); 
                    Wz(count,index(kk+1,ii,jj)) = -( 1 / dz(kk+1));

                    vz(count) = sqrt( dx(ii) * dy(jj) * dz(kk) );
                    count=count+1;

                else

                    count=count+1;

                end
                
            end
            
        end
    end
end

Vz = spdiags(vz,0,mcellz,mcellz);