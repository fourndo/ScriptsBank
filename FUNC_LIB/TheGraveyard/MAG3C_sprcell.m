function [sprcell] = MAG3C_sprcell(dx,dy,dz,nullcell,topocell)
% Find the closest active cell for each topocell

% Pre-allocate
sprcell = zeros(length(topocell),1);

% Create cell center array
xc = cumsum(dx) - dx/2;
yc = cumsum(dy) - dy/2;
zc = cumsum(dz) - dz/2;

nullcell = reshape(nullcell,length(dz),length(dx),length(dy));

[Zc,Xc,Yc] = ndgrid(zc,xc,yc);

for ii = 1 : length(topocell)
    
    % Compute distance from topocell to all active cells
    r = sqrt( (Xc(topocell(ii)) - Xc ).^2 +...
        ( Yc(topocell(ii)) - Yc ).^2 +...
        ( Zc(topocell(ii)) - Zc + 1e-4 ).^2);
    
    near = min(r(nullcell==1));
    
    temp = find((r==near)&(nullcell==1));
    
    sprcell(ii) = temp(1);
    
end

