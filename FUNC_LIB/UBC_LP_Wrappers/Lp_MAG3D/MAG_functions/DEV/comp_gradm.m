function [LHS] = comp_gradm(m, nX, nY, nZ, dX, dY, dZ)

mcell = nX * nY * nZ;
LHS=zeros( mcell, 1);
dmx = zeros(mcell,1);

count=1;

%% Compute the derivative terms ( WxtWx + WytWy + WztWz ) * m
for jj = 1 : nY
    
    for ii = 1 : nX
        
        for kk = 1 : nZ
            if ii == 1
                dmx(count) = m(count) * dX(ii) * dY(jj) * dZ(kk) - ...
                    m(count + nZ) * dX(ii + 1) * dY(jj) * dZ(kk);
                
            elseif ii==nX
                    dmx(count) = -m(count - nZ) * dX(ii - 1) * dY(jj) * dZ(kk) + ...
                    2*m(count) * dX(ii) * dY(jj) * dZ(kk);
            
            else
                dmx(count) = -m(count - nZ) * dX(ii - 1)* dY(jj) * dZ(kk) + ...
                    2 * m(count) * dX(ii) * dY(jj) * dZ(kk) - ...
                    m(count + nZ) * dX(ii + 1) * dY(jj) * dZ(kk);
            end
%             dmx(count)=dmx(count)*Wr(count);
            count=count+1;
        end
    end
end

dmy = zeros( mcell, 1);
count=1;
for jj = 1 : nY
    
    for ii = 1 : nX
        
        for kk = 1 : nZ
            if jj == 1
                dmy(count) = m(count) * dX(ii) * dY(jj) * dZ(kk) - ...
                    m(count + nZ * nX) * dX(ii) * dY(jj + 1) * dZ(kk);
                
            elseif jj==nY
                    dmy(count) = -m(count - nZ*nX) * dX(ii) * dY(jj - 1) * dZ(kk) + ...
                    2*m(count) * dX(ii) * dY(jj) * dZ(kk);
            
            else
                dmy(count) = -m(count - nZ * nX) * dX(ii )* dY(jj - 1) * dZ(kk) + ...
                    2 * m(count) * dX(ii) * dY(jj) * dZ(kk) - ...
                    m(count + nZ* nX) * dX(ii) * dY(jj + 1) * dZ(kk);
            end
%             dmy(count)=dmy(count)*Wr(count);
            count=count+1;
        end
    end
end


dmz = zeros(mcell, 1);
count = 1;

for jj = 1 : nY
    
    for ii = 1 : nX
        
        for kk = 1 : nZ
            if kk == 1
                dmz(count) = m(count) * dX(ii) * dY(jj) * dZ(kk) - ...
                    m(count + 1) * dX(ii) * dY(jj) * dZ(kk + 1);
                
            elseif kk==nZ
                    dmz(count) = -m(count - 1) * dX(ii) * dY(jj) * dZ(kk - 1) + ...
                    2*m(count) * dX(ii) * dY(jj) * dZ(kk);
            
            else
                dmz(count) = -m(count - 1) * dX(ii )* dY(jj) * dZ(kk - 1) + ...
                    2 * m(count) * dX(ii) * dY(jj) * dZ(kk) - ...
                    m(count + 1) * dX(ii) * dY(jj) * dZ(kk + 1);
            end
%             dmz(count)=dmz(count)*Wr(count);
            count=count+1;
        end
    end
end

LHS= dmx + dmx + dmy;