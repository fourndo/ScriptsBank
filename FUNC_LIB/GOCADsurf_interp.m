 
% Read in GOCAD surface file and assign to elevation based on cloud of
% points
clear all
close all

work_dir = 'C:\LC\Private\dominiquef\Projects\4239_Kaminak_Coffee_Mag\Modeling';
ts_file = 'CDED_115j10_to_15_NAD83_UTM7N_surf.ts';
topofile = 'CDED_115j10_pt.dat';

% Load topography(XYZ format)
topo = load([work_dir '\' topofile]);

% Set limits of topo surface and specs
dx = 75;
dy = 75;
xmax = 622000;
xmin = 563500;
ymin = 6956900;
ymax = 6983900;

x = xmin:dx:xmax;
y = ymin:dy:ymax;

nx = length(x);
ny = length(y);

[X,Y]=ndgrid(x,y);

nc = size(X,1) * size(X,2);
index = 1:nc;

% Interpolate topo onto grid
% F = scatteredInterpolant(topo(:,1),topo(:,2),topo(:,3));


% topo_grid = F(reshape(X,nc,1),reshape(Y,nc,1));
load([work_dir '\topo_grid']);


% Write to file surface file
fid = fopen([work_dir '\CDED_Interp_75m.ts'],'w');

fprintf(fid,'GOCAD TSurf 1\n'); 
fprintf(fid,'HEADER {\n');
fprintf(fid,'name:CDED_Interp_75m\n');
fprintf(fid,'}\n');
fprintf(fid,'TFACE\n');

% Write vertices
for ii = 1 : nc
    
    fprintf(fid,'VRTX %i %9.2f %9.2f %9.2f\n',ii,X(ii),Y(ii),topo_grid(ii));
    
end

% Write triangles
X = reshape(X,nx,ny);
Y = reshape(Y,nx,ny);
index = reshape(index,nx,ny);
topo_grid = reshape(topo_grid,nx,ny);

for jj = 1 : ny-1
    
    for ii = 1 : nx-1
        
        fprintf(fid,'TRGL %i %i %i\n',index(ii,jj),index(ii+1,jj),index(ii,jj+1)); 
        
        fprintf(fid,'TRGL %i %i %i\n',index(ii+1,jj),index(ii,jj+1),index(ii+1,jj+1));  
            
    end
    
end
fprintf(fid,'END\n');

fclose(fid);