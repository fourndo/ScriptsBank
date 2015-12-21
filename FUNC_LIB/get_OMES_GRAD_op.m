function [Ws, Wx, Wy, V, Vx, Vy ]=get_OMES_GRAD_op(dx,dy,dz,OMESid)
%Build the derivative weighting matrix for 2D OMES model


nx = length(dx);
ny = length(dy);
nz = length(dz);

% Convert from 3D to 2D mesh
[~,I,J] = ind2sub([nz nx ny],OMESid);

% Keep only the horizontal address
nullcell = zeros(nx,ny);
nullcell(I,J) = 1;

% Compute half-cell distance 
% dxm = sqrt( dx(1:end-1)/2 + dx(2:end)/2 ) ;
% dym = sqrt( dy(1:end-1)/2 + dy(2:end)/2 ) ;

% Square root dimensions for dV
dx = sqrt(dx(:));
dy = sqrt(dy(:));
% dz = sqrt(dz(:));

mcell = nx*ny*3;

Wx = spalloc(mcell,mcell,2*mcell);
Wy = spalloc(mcell,mcell,2*mcell);

vx = zeros(mcell);
vy = zeros(mcell);
v = zeros(mcell);

% countx = 1;
% county = 1;
count = 1;
flag= 0;
for jj = 1 : ny
    
    for ii = 1 : nx
            
        for kk = 1 : 3
            
            countx = 1;
            dindx  = 1;
            
            if nullcell(ii,jj)==1
               
                if ii==nx
                    
                    countx = -1;
                    dindx = -1;
                    
                end
                
                while nullcell(ii+countx,jj)==0 && ii+countx <= nx
                    
                    % If reach the end of the mesh, increment the other way
                    if ii+countx == nx && nullcell(ii+countx,jj)==0
                        
                        flag = 1;
                        break
                        countx = -1;
                        dindx = -1;
                        
                    end
                    
                    % If only cell on this row, then move to next cell
                    if ii+countx < 1
                        
                        flag = 1;
                        break
                        
                    end
                    
                    countx = countx + dindx;
                    
                end
                
                if flag==1
                   
                    flag = 0;
                    break
                    
                end
                
                if nullcell(ii+countx,jj)==1
                    
                    Wx(count,count) = -dindx;
                    Wx(count,count+countx) = dindx;
                    
                    dxm = sum( dx(ii:dindx:ii+countx)) - (dx(ii) + dx(ii+countx)) / 2 ;
                    vx(count) = dy(jj)/dxm;
                
                end
                
                
            
%             elseif ii == nx && nullcell(ii,jj)==1
%                 
%                 
%                 while nullcell(ii-countx,jj)==0 && ii-countx <= 1
%                     
%                     countx = countx + 1;
%                     
%                 end
%                 
%                 if nullcell(ii-countx,jj)==1
%                     
%                     Wx(count,count) = 1;
%                     Wx(count,count-countx) = -1;
%                     
%                 end
                
%                 dxm = sum( dx(ii:-1:ii-countx)) - (dx(ii) + dx(ii-countx)) / 2 ;
%                 vx(count) = dy(jj)/dxm;
                
            
            
            
            county = 1;
            dindy  = 1;
               
            if jj==ny

                county = -1;
                dindy = -1;

            end

            while nullcell(ii,jj+county)==0 && jj+county <= ny           

                % If reach the end of the mesh, increment the other way
                if jj+county == ny && nullcell(ii,jj+county)==0

                    flag = 1;
                    break
                    county = -1;
                    dindy = -1;

                end

                % If only cell on this row, then move to next cell
                if jj+county < 1

                    flag = 1;
                    break

                end
                
                county = county + dindy;

            end

            if flag==1

                flag = 0;
                break

            end

            if nullcell(ii,jj+county)==1

                Wy(count,count) = -dindy;
                Wy(count,count+3*nx*county) = dindy;
               
                dym = sum( dy(jj:dindy:jj+county)) - (dy(jj) + dy(jj+county)) / 2 ;
                vy(count) = dx(ii)/dym;
            
            end

            
            
            v(count) = dx(ii)*dy(jj);
            
            end
            
            count = count + 1;
        end   
                      
    end

end
    
Ws = spdiags( kron(nullcell(:),ones(3,1)) ,0,mcell,mcell);
V = spdiags(v,0,mcell,mcell);
Vx = spdiags(vx,0,mcell,mcell);
Vy = spdiags(vy,0,mcell,mcell);
