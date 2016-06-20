function [M] = azmdip_2_xyz(m_azm,m_dip)
% Function converts an input model with azimuth and dip of magnetization
% and converts to x, y, z components for every cell
% M: 3-by-n array

mcell = length(m_azm);

M = zeros(mcell,3);

% Unit vector
% for ii = 1 : mcell

    % Modify azimuth from North to Cartesian
D = mod(450-m_azm,360);
I = m_dip;

M(:,1) = cosd(I) .* cosd(D) ;
M(:,2) = cosd(I) .* sind(D) ;
M(:,3) = sind(I) ;


% end
