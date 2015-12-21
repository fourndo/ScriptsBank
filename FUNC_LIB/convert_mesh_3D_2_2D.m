% Convert UBC3D mesh to UBC2D mesh
clear all
close all

addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\

work_dir = 'C:\LC\Private\dominiquef\Projects\4253_UBC_GIF_Testing\Codes\DCIP2D';
mesh = 'Mesh_3D.msh';
model = 'Synthetic_3d.chg';
% ndv = -100;
mesh = get_UBC_mesh([work_dir '\' mesh]);
ndv=-99999;

    
% Check orientation
if mesh(1,1)==1
    nx = mesh(1,2);
    dx = mesh(4,1:nx);
    xo = mesh(2,2);
    fprintf('Mesh is oriented NS\n')
    
elseif mesh(1,2)==1
    nx = mesh(1,1);
    dx = mesh(3,1:nx);
    xo = mesh(2,1);
    fprintf('Mesh is oriented EW\n')
    
else
    
    fprintf('Mesh is 3D, the program will stop.\n')
    
end

% Find the minimum number of non-zero in each colum
m = load([work_dir '\' model]);
m = reshape(m,mesh(1,3),mesh(1,1),mesh(1,2));
m(m==-100)=0;
nullcell = m;
nullcell(m~=ndv) = 1;
nullcell(m==ndv) = 0;

nz = min(min(sum(nullcell,1)));
dz = [mesh(5,1:nz)];% (mesh(5,nz)*(1.4 .^[1:10]))] ;
zo = mesh(2,3);


fid = fopen([work_dir '\Mesh_2D.msh'],'w');
% Write horizontal component
fprintf(fid,'%i\n',nx);

for ii = 1:nx
    
    if ii == 1
        
        fprintf(fid,'%12.6f %12.6f 1\n',xo,xo+sum(dx(1:ii)));
        
    else
        
        fprintf(fid,'%12.6f 1\n',xo+sum(dx(1:ii)));
        
    end
    
end

% Write vertical component
fprintf(fid,'\n%i\n',nz);

for ii = 1:nz
    
    if ii == 1
        
        fprintf(fid,'%12.6f %12.6f 1\n',0,sum(dz(1:ii)));
        
    else
        
        fprintf(fid,'%12.6f 1\n',sum(dz(1:ii)));
        
    end
    
end

fclose(fid);

%% Convert model
mod2d = zeros(nz,nx);
% Only grab the top nz ground cells
for jj = 1 : nx
    
    colm = m(nullcell(:,jj)==1,jj);
    mod2d(1:(nz),jj)=colm(1:(nz));
   
end

fid = fopen([work_dir '\Model_2D.dat'],'w');
fprintf(fid,'%i %i \n', nx,nz); 
count = 1;

for ii = 1 : nz
    
    for jj = 1 : nx
        
       fprintf(fid,'%12.6e ', mod2d(ii,jj));       
       
    end
    
        fprintf(fid,'\n');
        
end

fclose(fid);
       