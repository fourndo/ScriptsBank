 %getpath
%Compute N number of randomely generated ray paths.
%First determine if ray intersects a receiver
%If NO, ray path is only plotted.
%If YES, ray path populate a data matrix, then plotted

close all
clear all
addpath data
addpath functions

load data


% Build mesh
dX = ones(1,100)*2;
dZ = ones(1,60)*5;


X0 = -20;
Z0 = 300;

% length of land surveyed
nX = length(dX);
nZ = length(dZ);

mcell=nX*nZ;
ndata=size(data,1);
d=zeros(ndata,1);

RxX=160;
TxX=0;
%Generate ray path for every tx-tx pairs
counter=1;

for ii=1:size(data,1)
    
    
        
        rangeX = RxX - TxX;
        rangeZ = data(ii,1) - data(ii,2);
        
        angl_in = atan(rangeZ/rangeX);
        
        R(ii)= (rangeX^2 + rangeZ^2) ^(0.5);
        % Modified data
        d(ii)=-log( (abs(data(ii,4)) * R(ii) ) / cos(angl_in)^2);
        G(ii,:) = compG(X0,Z0,nX,nZ,dX,dZ,RxX,Z0-data(ii,1),TxX,Z0-data(ii,2));
        residual(ii)=sum(G(counter,2:end))-R(ii);
        
        
    
end

%Create scaling matrix
noise_est=abs(1./(0.02*d));
Wd=eye(ndata);
for ii=1:ndata
    Wd(ii,ii)=(max(R)-min(R))/R(ii);
end


%Create mesh group
mesh{1}=[X0 Z0];
mesh{2}=dX;
mesh{3}=dZ;

save ('data/d','d')
save ('data/kernel','G');
% save ('data/model','m');
save ('data/mesh','mesh');
save ('data/Wd','Wd');

%Write UBC mesh file
fid=fopen('mesh.dat','w');
fprintf(fid,'%d %d %d\n',nX,1,nZ);


fid=fopen('mesh.dat','a');
fprintf(fid,'%d %d %d\n',mesh{1}(1),0,mesh{1}(2));

for ii=1:nX
    
fprintf(fid,'%d ',dX(ii));
end

fprintf(fid,'\n');
fprintf(fid,'5');
fprintf(fid,'\n');
for ii=1:nZ
    
fprintf(fid,'%d ',dZ(ii));
end


AA=sum(G);
figure;imagesc(reshape(AA(2:end),nZ,nX))