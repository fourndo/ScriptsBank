 %getpath
%Compute N number of randomely generated ray paths.
%First determine if ray intersects a receiver
%If NO, ray path is only plotted.
%If YES, ray path populate a data matrix, then plotted

close all
clear all
addpath data
addpath functions





% For now all the cells have dimension 1x1x1
dX = ones(1,20);
% dY = [2 1.5 1.25 ones(1,6) 1.25 1.5 2];
dZ = ones(1,30);

X0 = 0;
Y0 = 0;
Z0 = 300;

%length of land surveyed
nX = length(dX);
nZ = length(dZ);


mcell=nX*nZ;

%Densitiy 2D model
background = 0;


model=zeros(nZ,nX);

% model(30:end,:)=-0.05;

anomaly = 0.1;
target = [14 10];
[model] = genModel(nX, nZ, model, anomaly, target);

anomaly = -0.2;
target = [14 20];
[model] = genModel(nX, nZ, model, anomaly, target);

anomaly = 0.3;
target = [8 20];
[model] = genModel(nX, nZ, model, anomaly, target);

anomaly = -.4;
target = [8 10];
[model] = genModel(nX, nZ, model, anomaly, target);



m=zeros(mcell+1,1);
m(2:end)=model(:);

figure(1)
imagesc(model)
caxis([-0.5 0.5])

hold on

%Create transmiter, receiver stations

Rx.Z=1:1:25;%[3:2:58 3:2:58];
Rx.X=ones(1,length(Rx.Z))*19;%[ones(1,length(Rx.Z)/2)*14 ones(1,length(Rx.Z)/2)*1];

nRx=length(Rx.Z);


Tx.Z=2.5:2.5:25;%[2:2:58 2:2:58];
Tx.X=ones(1,length(Tx.Z))*2;%[ones(1,length(Tx.Z)/2)*14 ones(1,length(Tx.Z)/2)*1];

nTx=length(Tx.Z);

ndata= nRx * nTx;
%Set maximum field of view angle (in radian)
maxFOV=2 * pi * 45 / 360;

%Generate ray path for every tx-tx pairs
counter=1;

for ii=1:nTx
    
    
    for jj=1:nRx
        
        rangeX = Rx.X(jj) - Tx.X(ii);
        rangeZ = Rx.Z(jj) - Tx.Z(ii);
        
        angl_in = atan(rangeZ/rangeX);
        
        if angl_in < maxFOV && angl_in > -maxFOV
        R(counter)= (rangeX^2 + rangeZ^2) ^(0.5);
        G(counter,:) = compG(X0,Z0,nX,nZ,dX,dZ,Rx.X(jj),Z0-Rx.Z(jj),Tx.X(ii),Z0-Tx.Z(ii));
        
        TxZ(counter) = Tx.Z(ii);
        RxZ(counter) = Rx.Z(jj);
        counter=counter+1;
        figure(1)
        plot([Tx.X(ii) Rx.X(jj)],[Tx.Z(ii) Rx.Z(jj)],'k:')
        hold on
        end
        
    end
end
figure(1)
pbaspect([0.5 1 1])

%Create data matrix
data = G * m ;

ndata=length(data);
%Corrupt with 5% random noise
% d = awgn(data,-12.5);
% noise = ( (data.*.05) .* randn(length(data),1) );
pct_noise = 0.1;
noise = (pct_noise.*max(abs(data))).*randn(ndata,1);
d = data + noise;
% 
% Wd=eye(ndata);
% for ii=1:ndata
%     Wd(ii,ii)=(max(R)-min(R))/R(ii);
% end

%Create mesh group
mesh{1}=[X0 Y0 Z0];
mesh{2}=dX;
mesh{3}=dZ;

save ('data/d','d');
save ('data/kernel','G');
save ('data/model','m');
save ('data/mesh','mesh');
save ('data/Tx','TxZ');
save ('data/Rx','RxZ');
% save ('data/Wd','Wd');
