function [p,s,t] = azmdip_2_pst(azm,dip,mcell)

Rz = @(x)   [cosd(x) -sind(x) 0;
            sind(x) cosd(x) 0;
            0 0 1];
                    
Rx = @(x)   [ 1 0 0;
            0 cosd(x) -sind(x);
            0 sind(x) cosd(x)];
             
azm = -azm;
% T = TMI * [Tx;Ty;Tz];

% Inducing direction
yvec = [0;1;0];
pvec = Rz(azm) * Rx(dip) * yvec;

% m_azm = ones(mcell,1)*Dazm;
% m_dip = ones(mcell,1)*I;
% mv = H * azmdip_2_xyz(m_azm,m_dip);

p = [pvec(1) * speye(mcell,mcell);pvec(2) * speye(mcell,mcell);pvec(3) * speye(mcell,mcell)];

% Perpendicular S
zvec = [0;0;1];
svec = Rz(azm) * Rx(dip) * zvec;

% m_azm = ones(mcell,1)*(Dazm);
% m_dip = ones(mcell,1)*(I-90);
% mv = H * azmdip_2_xyz(m_azm,m_dip);

s = [svec(1) * speye(mcell,mcell);svec(2) * speye(mcell,mcell);svec(3) * speye(mcell,mcell)];

% Third orthogonal direction T
% Perpendicular S
xvec = [1;0;0];
tvec = Rz(azm) * Rx(dip) * xvec;

t = [tvec(1) * speye(mcell,mcell);tvec(2) * speye(mcell,mcell);tvec(3) * speye(mcell,mcell)];
