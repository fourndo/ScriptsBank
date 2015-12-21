function [L] = checkloop(XYZ,flag)
% Inputs
% XYZ   : [x(:) y(:) z(:)] location file for the loop
% flag  : 'clockw' or 'counterclockw'
%
% Ouput
% L     : XYZ coordinates of the densified and ordered transmitter loop
%
% Mira Geoscience
% Written by: D. Fournier
% Last update: June 26th, 2014
 
switch flag 
    case 'clockwise'
        
        rot = -1;
        
    otherwise
        
        rot = 1;
        
end

% Remove replicates
tester = sum(XYZ,2);
[arg,AA,CC] = unique(tester);

XYZ = XYZ(AA,:);

nn = size(XYZ,1);

% Pre-allocate memory
L = zeros(nn,3);

% Find the xy centroid of the loop and 0 azmuth vector
cx = mean(XYZ(:,1));
cy = mean(XYZ(:,2));
r0 = [1 0];

% Find a starting point (first entry)
L(1,:) = XYZ(1,:);

% Loop over the points and find neighbour in direction
% then densify if dl > dx
figure(1);
plot(L(1,1)-cx,L(1,2)-cy,'*'); hold on
plot(0,0,'ko'); hold on

for ii = 2 : nn
    
    % Order nodes as function of distance from current
    dX = L(ii-1,1) - XYZ(:,1);
    dY = L(ii-1,2) - XYZ(:,2);

    R = sqrt(dX.^2 + dY.^2);
    
    [arg,index] = sort(R);
    
    ri(1) = XYZ(index(1),1) - cx;
    ri(2) = XYZ(index(1),2) - cy;
%     quiver(0,0,ri(1),ri(2)); hold on
    % Look at the first option and keep if in the right direction
    % Otherwise, go down the list of closest nodes until in finds one
    for jj = 2 : nn
        
        rj(1) = XYZ(index(jj),1) - cx;
        rj(2) = XYZ(index(jj),2) - cy;

        % Compute dot product to find angle
        ang(1) = acosd( ( ri * r0' ) / norm(ri) );
        
        if ri(2) < 0
            ang(1) = 360 - ang(1);  % If bottom half circle
            
        end
        
        ang(2) = acosd( ( rj * r0' ) / norm(rj) );
        if rj(2) < 0
            ang(2) = 360 - ang(2);
        end
        
        % Check when crossing the origin
        if sign(rj(2)) ~= sign(ri(2)) && ri(1)>0
            ang(ang>180) = ang(ang>180) - 360; 
        end
        
        % If the closest point 
        if sign(ang(2) - ang(1)) == rot
            
            L(ii,:) = XYZ(index(jj),:);
            quiver(ri(1),ri(2),rj(1)-ri(1),rj(2)-ri(2)); hold on
            break
            
        end
        
    end
    
end

% Last check to make sure the loop is closed
if (L(end,1)~= L(1,1)) || (L(end,2)~= L(1,2))
    
    L = [L;L(1,:)];
    quiver(rj(1),rj(2),L(end,1)-cx-rj(1),L(end,2)-cy-rj(2)); hold on
    
end
plot(L(end,1)-cx,L(end,2)-cy,'ro'); hold on
legend('Start','Centroid','Current');title('\bfTx loop')