function [ Wx , Wy, Wz, Vx, Vy, Vz ]=get_GRAD_op3D_v4(dx,dy,dz,nullcell)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

nx = length(dx);
ny = length(dy);
nz = length(dz);

% Compute half-cell distance 
dxm = sqrt( dx(1:end-1)/2 + dx(2:end)/2 ) ;
dym = sqrt( dy(1:end-1)/2 + dy(2:end)/2 ) ;
dzm = sqrt( dz(1:end-1)/2 + dz(2:end)/2 ) ;%dzm = [dzm(1) dzm];

% Square root dimensions for dV
dx = sqrt(dx(:));
dy = sqrt(dy(:));
dz = sqrt(dz(:));

mcell = nx*ny*nz;

nullcell = reshape(nullcell,nz,nx,ny);

% Generate Gradient Operators and volume matrices
mcellx = (nx-1)*ny*nz;
mcelly = nx*(ny-1)*nz;
mcellz = nx*ny*(nz-1);

Wx = spalloc(mcellx,mcell,2*mcellx);
Wy = spalloc(mcelly,mcell,2*mcelly);
Wz = spalloc(mcellz,mcell,2*mcellz);

Vx = speye(mcellx);
Vy = speye(mcelly);
Vz = speye(mcellz);

countx = 1;
county = 1;
countz = 1;
count = 1;

for jj = 1 : ny
    
    for ii = 1 : nx
        
        for kk = 1 : nz
            
            if ii < nx
               
                if nullcell(kk,ii,jj)==1 && nullcell(kk,ii+1,jj)==1
                    
                    Wx(countx,count) = -1;
                    Wx(countx,count+nz) = 1;
                    
                elseif nullcell(kk,ii,jj)==1 && nullcell(kk,ii+1,jj)==0 && ii > 1
                    
                    Wx(countx,count) = 1;
                    Wx(countx,count-nz) = -1;
                    
                end
                
                Vx(countx,countx) = dy(jj)*dz(kk)/dxm(ii);
                countx = countx + 1;
                
            end
            
            if jj < ny
               
                if nullcell(kk,ii,jj)==1 && nullcell(kk,ii,jj+1)==1
                    
                    Wy(county,count) = -1;
                    Wy(county,count+nx*nz) = 1;
                    
                elseif nullcell(kk,ii,jj)==1 && nullcell(kk,ii,jj+1)==0 && ii > 1
                    
                    Wy(county,count) = 1;
                    Wy(county,count-nx*nz) = -1;
                    
                end
                
                Vy(county,county) = dx(ii)*dz(kk)/dym(jj);
                county = county + 1;
                
            end
            
            if kk < nz
               
                if nullcell(kk,ii,jj)==1 && nullcell(kk+1,ii,jj)==1
                    
                    Wz(countz,count) = -1;
                    Wz(countz,count+1) = 1;
                    
                elseif nullcell(kk,ii,jj)==1 && nullcell(kk+1,ii,jj)==0 && ii > 1
                    
                    Wz(countz,count) = 1;
                    Wz(countz,count-1) = -1;
                    
                end
                
                Vz(countz,countz) = dx(ii)*dy(jj)/dzm(kk);
                countz = countz + 1;
                
            end
                
            count = count + 1;
                    
         
                    
        end
        
    end
    
end

end