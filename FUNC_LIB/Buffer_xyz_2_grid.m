%% Topocheck_DRIVER
clear all
close all

work_dir = 'C:\LC\Private\dominiquef\Projects\3796_AGIC_Research\HyperCube\Processing\MtDore_files_for_Hypercube\Targeting_2015\Regions';

meshfile = 'Mesh_100m_2D.msh';
xyz_file = 'tsites_ALL.dat';

radius = 500;
%% SCRIPT START HERE
% Load mesh file and convert to vectors (UBC format)
[xn,yn,zn]=read_UBC_mesh([work_dir '\' meshfile]);

nx = ( length(xn) - 1 );
ny = ( length(yn) - 1 );
nz = ( length(zn) - 1 );

mcell = nx*ny*nz;

% Get cell center values
xx = ( xn(1:end-1) + xn(2:end) ) / 2;
yy = ( yn(1:end-1) + yn(2:end) ) / 2;
zz = ( zn(1:end-1) + zn(2:end) ) / 2;

[ZZ,XX,YY] = ndgrid(zz,xx,yy);
ZZ = ZZ(:);
XX = XX(:);
YY = YY(:);

% Load data points
data = load([work_dir '\' xyz_file]);
nstn = size(data,1);
stn_id = 1 : nstn;

buff = zeros(nstn,mcell);

% Find buffer zones for each point
for ii = 1 : nstn
    
    R = sqrt(( data(ii,1) - XX ).^2 + ( data(ii,2) - YY ).^2);
    temp = R < radius;
    buff(ii,temp) = R(R < radius);
    
end
    
denom = sum(buff,1);

zones = zeros(mcell,1);
% Discard overlap and choose closest zone
for ii = 1 : mcell
    
    % Take all other zones then current one
%     temp = (sum(buff(stn_id ~= ii,:)*-1) + buff(stn_id == ii,:))./denom;
    if sum(buff(:,ii))~=0
    idx = buff(:,ii)==min(buff(buff(:,ii)~=0,ii));
    idx = find(idx);
    zones(ii) = idx(1);
    end
    
end

save([work_dir '\buffer_zones_ALL.dat'],'-ascii','zones');

%% Create random sampling zones
nsub = round(2*size(buff,1)/3);

for jj = 1 : 10
    
    target = randperm(nstn);
    target = target(1:nsub);
    z_temp = zeros(mcell,1);
    
    for kk = 1 : nsub
        
        z_temp(zones==target(kk)) = target(kk);
        
    end
    
    save([work_dir '\buffer_zones_SUB' num2str(jj) '.dat'],'-ascii','z_temp');
    
end
