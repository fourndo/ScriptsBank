function [dmx] = comp_gradxm(m, nX, nY, nZ, dX, dY, dZ,topo_model)

dmx = zeros((nX-1)*nY*nZ,1);

count=1;
skip = 0;
%% Compute the derivative terms ( WxtWx + WytWy + WztWz ) * m
for jj = 1 : nY
    
    for ii = 1 : nX
        
        for kk = 1 : nZ
            if (ii == nX) %|| mnull(kk,ii,jj)== 0 || mnull(kk,ii+1,jj)== 0
                
            skip=skip+1;
            else
                
                if topo_model(count+skip)==0 || topo_model(count+skip+nZ)==0
                    
                    dmx(count)=0;
                    count=count+1;
                    
                else
                    dXmid = (dX(ii) + dX(ii+1))/2;
                    dmx(count) = (-m(count+skip) * sqrt(dZ(kk)*dY(jj)/dX(ii)) +...
                        m(count+skip+nZ) * sqrt(dZ(kk)*dY(jj)/dX(ii+1)) ) ;
                    count=count+1;
                    
                end
                
            end
            
        end
    end
end

