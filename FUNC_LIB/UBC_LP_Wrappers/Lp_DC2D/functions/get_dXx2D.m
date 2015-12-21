function [dXx] = get_dXx2D(m, nX, nZ, dX, dZ,topo_model)

dXx = zeros((nX-1)*nZ,1);

count=1;
skip = 0;
%% Compute the derivative terms ( WxtWx + WytWy + WztWz ) * m
   
for ii = 1 : nX

    for kk = 1 : nZ
        if (ii == nX) %|| mnull(kk,ii,jj)== 0 || mnull(kk,ii+1,jj)== 0

        skip=skip+1;
        else

            if topo_model(count+skip)==0 || topo_model(count+skip+nZ)==0

                dXx(count)=0;
                count=count+1;

            else
                dXmid = (dX(ii) + dX(ii+1))/2;
                dXx(count) = sqrt( dZ(kk) * dXmid);
                count=count+1;

            end

        end

    end
end


