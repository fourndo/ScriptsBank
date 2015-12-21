function [jlink] = MAG3C_RTC_jlink(obsx,obsy,obsz,t_cntm,OctLev,Ls)
% Function [link] = MAG3C_RTC_tlink(obsx,obsy,obsz,t_cntm)
% Takes in center of mass information from subfuntion MAG3C_RTC_cntmass 
% and find the octree break down for given obersvation points
%
% INPUT
% obsx,obsy,obsz: Observation location [X(:) Y(:) Z(:)]
%
% t_cntm: active and topocell center of mass information
% *_cntm (mcell-by-4): [Rx Ry Rz V] Coordinate of center of mass
% and volume of everycell. Center of mass of topocell has been computed
% for cells below topo with highest octree level
%
% Lx: Length scale for the number of cells away to be included in the
% highest octree level. If more octree levels are 
%
% OUTPUT
% tlink (nobs-by-tlink): Integer for the octree level for each jth 
% observation and topocell   
%
%

% Pre-allocate
jlink = zeros(size(obsx,1),size(t_cntm,1)); % Number of associated tcell

if isempty(jlink)==1

    jlink = zeros(size(obsx,1),1);
    fprintf('Program as detected no topographic cells\n')
    return

end


for jj = 1 : length(obsx)

    % Compute distance from topocell to all active cells
    dx = t_cntm(:,1) - obsx(jj,1);
    dy = t_cntm(:,2) - obsy(jj,1);
    dz = t_cntm(:,3) - obsz(jj,1);
    
    % Compute distance between current topo cell and all active cells
    R = sqrt( dx.^2 + dy.^2 + dz.^2 );

    for ii = 1:length(OctLev)

        % Assign level to cells
        jlink( jj, R < (sqrt(3) * Ls * (OctLev(end) - ii + 1) ) ) = ii;

    end
   
end

% % Re-order the links for the sequence of active cells
% count = 1;
% for jj = 1 : size(a_cntm,1);
%     
%     index = find(linkup==jj);
%     linkat(count:count+nlink(jj)-1 ) = index;
%     count = count + nlink(jj);
% end
