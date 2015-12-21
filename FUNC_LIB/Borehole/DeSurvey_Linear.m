% Function de-survey hole
% Load in XYZ,Azm,Dip of borehole
% Ouput XYZ location of samples along hole at specified depths

close all 
clear all

work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\EOSC 556b (2013) EM\TKC_Project\Modelling\Phys_Prop_Desurveyed';

data_file = 'Hole_XYZ_AZD.txt';

dsep = '\';
%% Load data and operators

data = load([work_dir dsep data_file]);

% Vector direction
m = @(a,t,p) [a.*cos(p).*cos(t);...
    a.*cos(p).*sin(t);...
    a.*sin(p)];

%% De-survey
ndata = size(data,1);

d_out = zeros(ndata,3);

for ii = 1 : size(data,1)
    
    % Get azm and dip and convert to Cartesian
    azm = mod(450-data(ii,4),360) * pi / 180;
        
    dip = data(ii,5) * pi / 180;
        
    % Compute vector direction
    v = m(1,azm,dip);
    
    dxyz = v * data(ii,6);
    
    % Compute XYZ of sample
    d_out(ii,:) = data(ii,1:3) + dxyz';
    
end

save([work_dir dsep 'SampleXYZ.dat'],'-ascii','d_out');