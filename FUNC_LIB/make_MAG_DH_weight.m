% FUNCTION make_MAG_DH_weight
clear all
close all


% work_dir    = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\AeroTEMinv\1D';
work_dir    = 'C:\LC\Private\dominiquef\Projects\Prospective\RealEagleExploration';
meshfile    = 'MESH12_REVISED.msh';
topofile    = 'topo2.XYZ';
modfile = 'Dist_r3_weight.dat';

%% Make nullcell from topo and mesh
[xn,yn,zn]=read_UBC_mesh([work_dir '\' meshfile]);
[Zn,Xn,Yn] = ndgrid(zn,xn,yn);

nx = length(xn) - 1 ;
ny = length(yn) - 1 ;
nz = length(zn) - 1 ;

% Load topography
topo = read_UBC_topo([work_dir '\' topofile]);

% Create discretize topogrphy
[nullcell,tcellID,ztopo_n] = topocheck(Xn,Yn,Zn,topo);
save([work_dir '\nullcell.dat'],'-ascii','nullcell');

load([work_dir '\nullcell.dat']);


%% Load model and normalize
model = load([work_dir '\' modfile]);

N = sum(nullcell);

scale = N / sum(model(nullcell==1));

ws = scale * model;

% Create blanck wx, wy , wz
wx = ones( (nx-1) * ny * nz , 1 );

wy = ones( nx * (ny-1) * nz , 1 );

wz = ones( nx * ny * (nz-1) , 1 );

W = [ws;wx;wy;wz]

save([work_dir '\Wsxyz_scaled.dat'],'-ascii','W');

