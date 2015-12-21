% Function de-survey hole
% Load in XYZ,Azm,Dip of borehole
% Ouput XYZ location of samples along hole at specified depths

close all 
clear all

addpath '..\'

work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\EOSC 556b (2013) EM\TKC_Project\Modelling\Phys_Prop_Desurveyed';

datafile = 'XYZPhys';
meshfile = 'TKC_common_mesh.msh';
nullfile = 'nullcell.dat';

dsep = '\';

rmax = 25;
%% Load data model space
load([work_dir dsep datafile]);
data = table2array(XYZPhys(:,4:end));
xyz = table2array(XYZPhys(:,1:3));

var_id = XYZPhys.Properties.VariableNames;
% Load mesh
[xn,yn,zn]=read_UBC_mesh([work_dir '\' meshfile]);

% Load nullcell
nullcell = load([work_dir dsep nullfile]);

%% Create mesh grid
x = (xn(2:end) + xn(1:end-1)) / 2;
y = (yn(2:end) + yn(1:end-1)) / 2;
z = (zn(2:end) + zn(1:end-1)) / 2;

[Z,X,Y] = ndgrid(z,x,y);



%% Interpolate onto grid
nvar = size(data,2);

for ii = 1 : nvar
    
    % Select all non-nan
    indx = ~isnan(data(:,ii));
    
    %% Create buffer zone around sample points
    buffer = zeros(size(Z));

    for jj = 1 : size(xyz,1)

        if indx(jj) == 1
            
            R = sqrt( (Z - xyz(jj,3) ).^2 + (X - xyz(jj,1) ).^2 + (Y - xyz(jj,2) ).^2);

            buffer(R < rmax) = 1;
            
        end
        
    end

    %% Create interp function
    F = scatteredInterpolant(xyz(indx,1),xyz(indx,2),xyz(indx,3),data(indx,ii),'nearest');
    
    % Interpolate onto cell-center location
    model = F(X,Y,Z);
    
    % Cut topography and write to file
    if ii == 7
        
        model(nullcell==0 | buffer(:) == 0) = 1e-8;
        model(nullcell==1 & buffer(:) == 1) = 1./model(nullcell==1 & buffer(:) == 1);
        
    else
        
        model(nullcell==0 | buffer(:) == 0) = -100;
        
    end
    
    model = model(:);
    
    save([work_dir dsep 'PhysProp_' var_id{ii+3} '_interp.dat'],'-ascii','model');
    
end