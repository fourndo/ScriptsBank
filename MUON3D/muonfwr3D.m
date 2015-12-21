 %getpath
%Compute N number of randomely generated ray paths.
%First determine if ray intersects a receiver
%If NO, ray path is only plotted.
%If YES, ray path populate a data matrix, then plotted

close all
clear all

addpath data/
addpath functions/

[meshfile]=get_UBCmesh('Mesh_10m_nopad.msh');

dx = meshfile(3,meshfile(3,:)~=0)';
dy = meshfile(4,meshfile(4,:)~=0)';
dz = meshfile(5,meshfile(5,:)~=0)';

nx = meshfile(1,1); %size(X,1);    %number of cell in X
ny = meshfile(1,2); %size(X,2);    %number of cell in Y
nz = meshfile(1,3); %size(X,3);    %number of cell in Z

x0 = meshfile(2,1);
y0 = meshfile(2,2);
z0 = meshfile(2,3);

nx = length(dx);
ny = length(dy);
nz = length(dz);

mcell = nx * ny * nz;

% m = importdata('Model_20m.den');
topo = importdata('Topo_Gaussian.topo');
surf = topo.data;
obs_loc = importdata('Obs_loc.obs');


% Discretize topo on mesh
% Create cell center array
xx = zeros(nx,1);
yy = zeros(ny,1);

for ii=1:nx
    xx(ii)= x0 + sum(dx(1:ii)) - dx(ii)/2;
end

for ii=1:ny
    yy(ii)= y0 + sum(dy(1:ii)) - dy(ii)/2;
end

[YY,XX] = meshgrid(yy,xx);

ZZ = zeros(nx,ny);

%% Using topo file
% for ii=1:size(surf,1);
%     
%     distoX = abs(surf(ii,1) - XX);
%     findx = (min(min(distoX))==distoX);
%             
%     
%     distoY = abs(surf(ii,2) - YY);
%     findy = (min(min(distoY))==distoY);
%     
%     
%     match = (findx.*findy)==1;
%     
%     ZZ(match)= surf(ii,3);
%     
%     
% end
%% Using UBC topocheck file

topocheck = ones(mcell,1);%importdata('topo_model.txt');
topocheck = reshape(topocheck,nz,nx,ny);

for ii = 1:nx
    
    for jj = 1:ny
        
        ZZ(ii,jj) = z0;
        kk = 1;
        
        while topocheck(kk,ii,jj)==0
            
            ZZ(ii,jj) = ZZ(ii,jj)- dz(kk);
            kk = kk+1;
            
        end
        
    end
    
end

fid = fopen('data/Ray_Collar.dat','w');
fprintf(fid,'ID X Y Z Length\n');

fid2 = fopen('data/Ray_Survey.dat','w');
fprintf(fid2,'ID Dist Azmth Dip\n');

% Create angle range
d_angl= 10*pi/180;
phi=0:(d_angl):(20*pi/180);
theta=[0:-(d_angl):-pi d_angl:d_angl:(pi-d_angl)];

%Create ray paths
ndata = length(phi)*length(theta)*size(obs_loc,1);

g=zeros(1,mcell);

G=sparse(ndata,mcell);

count=0;
vcount=0;
for ii=1:size(obs_loc,1)
    
    for jj= 1:length(phi)
        
        for kk = 1:length(theta)
    
        count = count+1;
        vcount = vcount+2;
        
        [g,xout,yout,zout]=getG_3D(obs_loc(ii,:),phi(jj),theta(kk),x0,y0,z0,dx,dy,dz,ZZ);
        
        G(count,:) = sparse(g);

        fprintf(fid,'%i %12.8f %12.8f %12.8f %12.8f\n',count,obs_loc(ii,1),obs_loc(ii,2),obs_loc(ii,3),500);
        fprintf(fid2,'%i %12.8f %12.8f %12.8f\n',count,0,theta(kk)*180/pi,90-phi(jj)*180/pi);
        fprintf(fid2,'%i %12.8f %12.8f %12.8f\n',count,500,theta(kk)*180/pi,90-phi(jj)*180/pi);
        

        end
    end
end

fclose(fid);
fclose(fid2);

% Create cell volume matrix
mnull = getnull(z0,dx,dy,dz,ZZ);
mnull = reshape(mnull,mcell,1);

%Create derivative matrices
Wx=getWx_3D(mcell,dx,dy,dz,mnull);
Wy=getWy_3D(mcell,dx,dy,dz,mnull);
Wz=getWz_3D(mcell,dx,dy,dz,mnull);
Ws=getWs_3D(mcell,dx,dy,dz,mnull);

save ('data/kernel_10m','G');
save ('data/mnull_10m','mnull');%save ('data/Topo_mnull.dat','-ascii','mnull');
save ('data/Wx_10m','Wx');
save ('data/Wy_10m','Wy');
save ('data/Wz_10m','Wz');
save ('data/Ws_10m','Ws');
