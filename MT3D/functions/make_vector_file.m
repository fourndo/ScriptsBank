function make_vector_file(field,XYZo,dX,dY,dZ,argin)
% Create XYZ pointset with vector attribute <i,j,k>

nX = length(dX);
nY = length(dY);
nZ = length(dZ);

field = reshape(field,nX*nY*nZ,3);

% Compute mid-cell mesh
dxm = dX(1:end-1)/2 + dX(2:end)/2; %dxm=[dxm(1);dxm;dxm(end)];
dym = dY(1:end-1)/2 + dY(2:end)/2; %dym=[dym(1);dym;dym(end)];
dzm = dZ(1:end-1)/2 + dZ(2:end)/2; %dzm=[dzm(1);dzm;dzm(end)];

% Create mid-cell coordinates
xm=[(XYZo(1) + dX(1)/2);(XYZo(1) + dX(1)/2) + cumsum(dxm)];
ym=[(XYZo(2) + dY(1)/2);(XYZo(2) + dY(1)/2) + cumsum(dym)];
zm=[(XYZo(3) - dZ(1)/2);(XYZo(3) - dZ(1)/2) - cumsum(dzm)];


% Number of faces in each directions


vector_file = zeros(nX*nY*nZ,6);

count = 1;

    
    for kk = 1 : nZ
        
        for jj = 1 : nY
            
            for ii = 1 : nX
                
                vector_file(count,1) = xm(ii);
                vector_file(count,2) = ym(jj);
                vector_file(count,3) = zm(kk);
                vector_file(count,4) = field(count,1);
                vector_file(count,5) = field(count,2);
                vector_file(count,6) = field(count,3);
                
                count = count +1;
                
            end
            
        end
        
    end
    
save(argin,'-ascii','vector_file');