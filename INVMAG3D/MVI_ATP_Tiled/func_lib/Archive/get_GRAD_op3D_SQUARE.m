function [ Ws, Gx , Gy, Gz, V, Vx, Vy, Vz ]=get_GRAD_op3D_SQUARE(dx,dy,dz,nullcell,X)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

nx = length(dx);
ny = length(dy);
nz = length(dz);

% Compute half-cell distance 
dxm =  dx(1:end-1)/2 + dx(2:end)/2 ;
dym =  dy(1:end-1)/2 + dy(2:end)/2 ;
dzm = dz(1:end-1)/2 + dz(2:end)/2 ;%dzm = [dzm(1) dzm];

% Square root dimensions for dV
dx = sqrt(dx(:));
dy = sqrt(dy(:));
dz = sqrt(dz(:));

mcell = nx*ny*nz;

nullcell = reshape(nullcell,nz,nx,ny);

% Generate Gradient Operators and volume matrices
% mcellx = (nx-1)*ny*nz;
% mcelly = nx*(ny-1)*nz;
% mcellz = nx*ny*(nz-1);

Gx = spalloc(mcell,mcell,2*mcell);
Gy = spalloc(mcell,mcell,2*mcell);
Gz = spalloc(mcell,mcell,2*mcell);

Vx = speye(mcell);
Vy = speye(mcell);
Vz = speye(mcell);

countx = 1;
county = 1;
countz = 1;
count = 1;

for jj = 1 : ny
    
    for ii = 1 : nx
        
        for kk = 1 : nz
            
            if ii < nx
               
                if nullcell(kk,ii,jj)==1 && nullcell(kk,ii+1,jj)==1
                    
                    Gx(countx,count) = -1;
                    Gx(countx,count+nz) = 1;
                    
                elseif nullcell(kk,ii,jj)==1 && nullcell(kk,ii+1,jj)==0 && ii > 1
                    
                    if nullcell(kk,ii-1,jj)==1 
                    
                        Gx(countx,count) = -1;
                        Gx(countx,count-nz) = 1;
                    
                    end
                    
                end
                
                Vx(countx,countx) = dx(ii)*dy(jj)*dz(kk)/dxm(ii);
                countx = countx + 1;
                
            else
                
                if nullcell(kk,ii-1,jj)==1 && nullcell(kk,ii,jj)==1
                    
                    Gx(countx,count) = 1;
                    Gx(countx,count-nz) = -1;
                                      
                end
                
                Vx(countx,countx) = dx(ii)*dy(jj)*dz(kk)/dxm(ii-1);
                countx = countx + 1;
                
            end
            
            if jj < ny
               
                if nullcell(kk,ii,jj)==1 && nullcell(kk,ii,jj+1)==1
                    
                    Gy(county,count) = -1;
                    Gy(county,count+nx*nz) = 1;
                    
                elseif nullcell(kk,ii,jj)==1 && nullcell(kk,ii,jj+1)==0 && jj > 1
                    
                    Gy(county,count) = -1;
                    Gy(county,count-nx*nz) = 1;
                    
                end
                
                Vy(county,county) = dx(ii)*dy(jj)*dz(kk)/dym(jj);
                county = county + 1;
            
            else
                
                if nullcell(kk,ii,jj-1)==1 && nullcell(kk,ii,jj)==1
                    
                    Gy(county,count) = 1;
                    Gy(county,count-nx*nz) = -1;
                                      
                end
                
                Vy(county,county) = dx(ii)*dy(jj)*dz(kk)/dym(jj-1);
                county = county + 1;
                    
            end
            
            if kk < nz
               
                if nullcell(kk,ii,jj)==1 && nullcell(kk+1,ii,jj)==1
                    
                    Gz(countz,count) = -1;
                    Gz(countz,count+1) = 1;
                    
                elseif nullcell(kk,ii,jj)==1 && nullcell(kk+1,ii,jj)==0 && zz > 1
                    
                    Gz(countz,count) = -1;
                    Gz(countz,count-1) = 1;
                    
                end
                
                Vz(countz,countz) = dx(ii)*dy(jj)*dz(kk)/dzm(kk);
                countz = countz + 1;
                
            elseif nz~=1
                
                if nullcell(kk-1,ii,jj)==1 && nullcell(kk,ii,jj)==1
                    
                    Gz(countz,count) = 1;
                    Gz(countz,count-1) = -1;
                                      
                end
                
                Vz(countz,countz) = dx(ii)*dy(jj)*dz(kk)/dzm(kk-1);
                countz = countz + 1;
                        
            end
                
            count = count + 1;
                    
         
                    
        end
        
    end
    
end


% Create hmid dimensions matrix
dX = kron( kron( ones(ny,1) , dx(:) ), ones(nz,1) );
dY = kron( kron( dy(:)  , ones(nx,1) ), ones(nz,1));
dZ = kron( kron( ones(ny,1) , ones(nx,1) ), dz(:) );

v = dX .* dY .* dZ;

Ws = X * speye(mcell) * X';

V = X * spdiags( v , 0 , mcell,mcell) * X';

Gx = X * Gx * X';
Gy = X * Gy * X';
Gz = X * Gz * X';

Vx = X * Vx * X';
Vy = X * Vy * X';
Vz = X * Vz * X';


end