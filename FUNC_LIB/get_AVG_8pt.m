function [A] = get_AVG_8pt(dx,dy,dz)
% Create 8-point averaging operator for tensor mesh elevated to the power
% (pwr)
% Larger the power, the wider the averaging
%
% Last update: 2015-01-13
% Written by: D.Fournier

%% FOR DEV ONLY
% % % Create small example for testing
% clear all
% close all
% 
% dx = ones(6,1);
% dy = ones(7,1);
% dz = ones(8,1);
% 
% nx = length(dx);
% ny = length(dy);
% nz = length(dz);
% 
% m = zeros(nz,nx,ny);
% m(3,3,3) = 1;
% m = m(:);

%% SCRIPT START HERE

nx = length(dx);
ny = length(dy);
nz = length(dz);

av = @(n) spdiags (ones (n,1)*[0.25,0.5,0.25],[-1,0,1],n,n);

avcx = av(nx); avcx(1,1:2) = 0.5; avcx(end,end-1:end) = 0.5;
avcy = av(ny); avcy(1,1:2) = 0.5; avcy(end,end-1:end) = 0.5;
avcz = av(nz); avcz(1,1:2) = 0.5; avcz(end,end-1:end) = 0.5;

Avcx = kron ( kron ( speye (ny) , avcx), speye (nz));
Avcy = kron ( kron ( avcy , speye (nx)), speye (nz));
Avcz = kron ( kron ( speye (ny) , speye (nx)), avcz);

% A = Avcx;
A = (Avcx * Avcy * Avcz);
% for ii = 1 : pwr
%     A = A*A;
% end
