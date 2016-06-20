% Test script tp compute the distance between grid points and a triangle

clear all

addpath ..\.
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Desktop\Test_tria_2_grid';
meshfile = 'Mesh_2p5m.msh';
tsfile = 'Crown.ts';

[vrtx,trgl] = read_GOCAD_surf([work_dir '\' tsfile]);

% Create a grid of points
% nx = 20;
% ny = 20;
% nz = 20;
% 
% x0 = 0;
% y0 = 0;
% z0 = 0;
% 
% dx = ones(nx,1)*10/nx;
% dy = ones(ny,1)*10/nx;
% dz = ones(nz,1)*10/nx;
% 
% xc = x0 + cumsum(dx) - dx/2;
% yc = y0 + cumsum(dy) - dy/2;
% zc = z0 - cumsum(dz) + dz/2;

[xn,yn,zn] = read_UBC_mesh([work_dir '\' meshfile]);

zc = ( zn(2:end) + zn(1:end-1) ) / 2;
yc = ( yn(2:end) + yn(1:end-1) ) / 2;
xc = ( xn(2:end) + xn(1:end-1) ) / 2;

[Zc,Xc,Yc] = ndgrid(zc,xc,yc);

gridcc = [Xc(:) Yc(:) Zc(:)];
nc = size(gridcc,1);

% Empty index for nullcell
nullcell = zeros(nc,1);

% Create list of vertices (x,y,z)
%vrtx = [2 2 -2;8 3 -2;4 5 -5;6 8 -3];

% Create link between vertices
%tri = [1 2 3;3 2 4];


tic
for ii = 1 : (size(trgl,1))
    
    P1 = vrtx(trgl(ii,1),:);
    P2 = vrtx(trgl(ii,2),:);
    P3 = vrtx(trgl(ii,3),:);
    
    % Compute normal
    n = cross((P2 - P1),(P3 - P1));
    lnl = sqrt(sum(n.^2,2));

    % Distance P0P1
    P0P1 = gridcc - repmat(P1,nc,1);
    lP0P1l = sqrt( sum(P0P1.^2,2));

    % Angle between normal and grid points
    cosa = ( P0P1 * n') ./ (lP0P1l *lnl);

    % Length of vector P0P0'
    lvec0l = lP0P1l .*cosa;

    nullcell = nullcell + (lvec0l>=0)*1;

    % Full vector
    %veco = -spdiags(lvec0l,0,nc,nc)*repmat(n/lnl,nc,1);

end

save([work_dir '\vec0.dat'],'-ascii','nullcell')

% nullcell = logical(1 - all(nullcell,2));
% 
% toc 
% 
% % Reshape and only keep cells with all 8 nodes below topo
% ztest = reshape(nullcell,length(zn),length(xn),length(yn));
% 
% ztest = ztest(1:end-1,:,:) & ztest(2:end,:,:);
% xtest = ztest(:,1:end-1,:) & ztest(:,2:end,:);
% ytest = xtest(:,:,1:end-1) & xtest(:,:,2:end);
% 
% nullcell = double(ytest(:));
% 
% %write_UBC_mesh(work_dir,'\Mesh.txt',x0,y0,z0,dx,dy,dz)
% save([work_dir '\vec0.dat'],'-ascii','nullcell')

%% Write triangle file
fid = fopen([work_dir '\triangle.dat'], 'w');

for ii = 1 : (size(trgl,1))
    
    fprintf(fid,'%f \t %f \t %f\n',vrtx(trgl(ii,1),:));
    fprintf(fid,'%f \t %f \t %f\n',vrtx(trgl(ii,2),:));
    fprintf(fid,'%f \t %f \t %f\n',vrtx(trgl(ii,3),:));
    
end

fclose(fid);
