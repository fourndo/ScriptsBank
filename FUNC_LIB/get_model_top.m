function [m2D] = get_model_top(m3D,nx,ny,nz,ndv,dz)


m3D = reshape(m3D,nz,nx,ny);
m2D = ones(nx,ny) * ndv;

for ii = 1:nx;
    for jj = 1:ny
        
       indx = find(m3D(:,ii,jj)~=ndv & ~isnan(m3D(:,ii,jj)));
       
       m2D(ii,jj) = m3D(indx(dz),ii,jj);
       
    end
end
