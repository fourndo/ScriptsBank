function [ Wx , Wz, Vx, Vz ]=get_GRAD_op2D(dx,dz,nullcell)
%Build the derivative weighting matrix
%Takes care of the fact that model is a m-by-n matrix converted to a
%m*n-by-1

nx = length(dx);
nz = length(dz);

% Compute half-cell distance 
dxm = sqrt( dx(1:end-1)/2 + dx(2:end)/2 ) ;
dzm = sqrt( dz(1:end-1)/2 + dz(2:end)/2 ) ;%dzm = [dzm(1) dzm];

% Square root dimensions for dV
dx = sqrt(dx(:));
dz = sqrt(dz(:));

mcell = nx*nz;

nullcell = reshape(nullcell,nz,nx);

Wx = spalloc(mcell,mcell,2*mcell);
Wz = spalloc(mcell,mcell,2*mcell);

Vx = speye(mcell);
Vz = speye(mcell);

countx = 1;
countz = 1;
count = 1;

  
for ii = 1 : nx

    for kk = 1 : nz

        if ii < nx

            if nullcell(kk,ii)==1 && nullcell(kk,ii+1)==1

                Wx(countx,count) = -1;
                Wx(countx,count+nz) = 1;

            elseif nullcell(kk,ii)==1 && nullcell(kk,ii+1)==0 && ii > 1

                Wx(countx,count) = -1;
                Wx(countx,count-nz) = 1;

            end

            Vx(countx,countx) = dz(kk)/dxm(ii);
            countx = countx + 1;
            
        else
            
            Wx(countx,count) = -1;
            Wx(countx,count-nz) = 1;

            Vx(countx,countx) = dz(kk)/dxm(end);
            countx = countx + 1;
            
        end

        

        if kk < nz

            if nullcell(kk,ii)==1 && nullcell(kk+1,ii)==1

                Wz(countz,count) = -1;
                Wz(countz,count+1) = 1;

            elseif nullcell(kk,ii)==1 && nullcell(kk+1,ii)==0 && kk > 1

                Wz(countz,count) = -1;
                Wz(countz,count-1) = 1;

            end

            Vz(countz,countz) = dx(ii)/dzm(kk);
            countz = countz + 1;
            
        else
            
            Wz(countz,count) = -1;
            Wz(countz,count-1) = 1;

            Vz(countz,countz) = dx(ii)/dzm(end);
            countz = countz + 1;
            
        end

        count = count + 1;



    end

end
    

end