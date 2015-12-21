function [nlink,linkup] = MAG3C_RTC_linkID(a_cntm,t_cntm)
% Function [link] = MAG3C_RTC_tlink(obsx,obsy,obsz,t_cntm)
% Takes in center of mass information from subfuntion MAG3C_RTC_cntmass 
% and find the octree break down for given obersvation points
%
% INPUT
% 
% t_cntm: active and topocell center of mass information
% *_cntm (mcell-by-4): [Rx Ry Rz V] Coordinate of center of mass
% and volume of everycell. Center of mass of topocell has been computed
% for cells below topo with highest octree level
%
% OUTPUT
% nlink (m_acell-by-1): Integer for the number of tcell associated to each active cell  
% linkID (m_tcell-by-1): Integer for the ID of tcells 
%
% example:
% nlink = [0 2 1 3 0 ...]
% linkID = [34 35 45 56 57 58 ...]
%
% The first active cell has no children tcell, the second actice cell
% has tcell(34:35)

% Pre-allocate
nlink = zeros(size(a_cntm,1),1); % Number of associated tcell
% linkat = zeros(size(a_cntm,1),1) ; % link active to index of t_cells
linkup = zeros(size(t_cntm,1),1) ; % link tcells to active

for ii = 1 : size(t_cntm,1);
    
    % Compute distance from topocell to all active cells
    Rx = t_cntm(ii,1) - a_cntm(:,1);
    Ry = t_cntm(ii,2) - a_cntm(:,2);
    Rz = t_cntm(ii,3) - a_cntm(:,3);
    
    % Compute distance between current topo cell and all active cells
    R = ( Rx.^2 + Ry.^2 + Rz.^2 );
    
    % Find closest
    [~,index] = min(R);
    
    % Store number of tcell and ID
    nlink(index) = nlink(index) + 1;
    linkup(ii) = index ;
    
end

% % Re-order the links for the sequence of active cells
% count = 1;
% for jj = 1 : size(a_cntm,1);
%     
%     index = find(linkup==jj);
%     linkat(count:count+nlink(jj)-1 ) = index;
%     count = count + nlink(jj);
% end
