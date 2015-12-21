function [dmz] = comp_gradzm(m, nX, nY, nZ, dX, dY, dZ ,topo_model)



%% Compute the derivative terms ( WxtWx + WytWy + WztWz ) * m
dmz = zeros(nX* nY* (nZ-1), 1);
count = 1;
skip = 0;
for jj = 1 : nY
    
    for ii = 1 : nX
        
        for kk = 1 : nZ
            if (kk == nZ) %|| mnull(kk,ii,jj)== 0 || mnull(kk+1,ii,jj)== 0

                skip=skip+1;


            else
                if topo_model(count+skip)==0 || topo_model(count+skip+1)==0
                    dmz(count)=0;
                    count=count+1;
                else
                dZmid = (dZ(kk) + dZ(kk+1))/2;   
                dmz(count) = (-m(count+skip) +...
                    m(count + skip + 1) )/dZmid;

                count=count+1;
                end
            end
        end
    end
end

