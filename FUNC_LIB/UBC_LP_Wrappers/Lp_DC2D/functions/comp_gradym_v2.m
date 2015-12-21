function [dmy] = comp_gradym(m, nX, nY, nZ, dX, dY, dZ,topo_model,Wr)



%% Compute the derivative terms ( WxtWx + WytWy + WztWz ) * m
dmy = zeros( nX * (nY -1) * nZ, 1);
count=1;
skip = 0;
for jj = 1 : nY
    
    for ii = 1 : nX
        
        for kk = 1 : nZ
            
            if (jj == nY) %|| mnull(kk,ii,jj)== 0 || mnull(kk,ii+1,jj)== 0

                skip=skip+1;
            else 
                if topo_model(count+skip)==0 || topo_model(count+skip+nZ)==0

                        dmy(count)=0;
                        count=count+1;

                else
                        dYmid = (dY(jj) + dY(jj+1))/2;
                        dmy(count) = (-m(count+skip)*Wr(count+skip) + ...
                            m(count + skip + nZ * nX)*Wr(count + skip + nZ * nX) ) * (1 / dYmid );

                        count=count+1;
               
                end
                
            end
            
        end
    end
end
