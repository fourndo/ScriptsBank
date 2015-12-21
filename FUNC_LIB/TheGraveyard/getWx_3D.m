function [ Wx , Vx ]=getWx_3D(dx,dy,dz,nullcell)
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
mcellx = mcell-sum(sum(nullcell(:,end,:)));

Wx = spalloc(mcellx,mcell,2*mcell);
vx = zeros(mcellx,1);
count=1;

for jj=1:ny
    
    for ii=1:nx
        
        for kk=1:nz
        
            if (ii ~= nx) && nullcell(kk,ii,jj)~= 0
                
                if nullcell(kk,ii+1,jj)~= 0


                    dXhalf = (dx(ii+1) + dx(ii)) / 2;

                    Wx(count,index(kk,ii,jj)) = ( 1 / dx(ii)); 
                    Wx(count,index(kk,ii+1,jj)) = -( 1 / dx(ii+1));

                    vx(count) = sqrt( dy(jj) * dz(kk) * dx(ii) );
                    count=count+1;

                else

                    count=count+1;

                end
                
            end
            
        end
    end
end

Vx = spdiags(vx,0,mcellx,mcellx);

end