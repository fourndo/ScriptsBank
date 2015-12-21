% Compute the sum over each column of a UBC grid
%
% Created by: D. Fournier
% Date: August 23, 2012

clear all
close all

% Set no data value
% Usual NDV (-1(Mag), -100(Grav), -99999(Gocad), 1e-008(DC))
ndv=-1;

DELIMITER = ' ';
HEADERLINES = 0;

% Import the mesh file
meshfile = importdata('C:\PVK_Projects\3815_Capstone_MtMist_grav\Modeling\Analysis\Mesh_25m_300mdepth.msh', DELIMITER, HEADERLINES);

X=meshfile(1,1);    %Number of cells in X-dir
Y=meshfile(1,2);    %Number of cells in Y-dir
Z=meshfile(1,3);    %Number of cells in Z-dir

% Import the model

UBC_model = importdata('C:\PVK_Projects\3815_Capstone_MtMist_grav\Modeling\Analysis\lldenll_llsuscll_model.dat');

model_3D=reshape(UBC_model,Z,Y,X);

% Create new 2D grid for the computed sum of each column
model_2D=ones(Y,X)*ndv;

for ii=1:X
    for jj=1:Y
        column=model_3D(:,jj,ii);
        model_2D(jj,ii)=sum(column(column~=ndv));
    end
end

% Write back to UBC format

UBC_model_2D=reshape(model_2D,Y*X,1);

save('Column_sum.dat','-ascii','UBC_model_2D')