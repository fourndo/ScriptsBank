 %getpath
%Compute N number of randomely generated ray paths.
%First determine if ray intersects a receiver
%If NO, ray path is only plotted.
%If YES, ray path populate a data matrix, then plotted

close all
clear all

addpath data/
addpath functions/

mesh = importdata('Mesh_5m_nopad.msh');
m = importdata('Model.den');
topo = importdata('Topo_Gaussian.topo');

%length of land surveyed
% Xmax=60;
% Ymax=60;
% Zmax=20;

Xo=0;
Yo=0;
Zo=0;

dX=ones(1,60);
dY=ones(1,60);
dZ=ones(1,20);

nX=length(dX);
nY=length(dY);
nZ=length(dZ);

mesh{1}=[nX nY nZ];
mesh{2}=[Xo Yo Zo];
mesh{3}=dX;
mesh{4}=dY;
mesh{5}=dZ;

[Y,X,Z]=meshgrid(dY,dX,dZ);

model=zeros(nZ,nX,nY);
%Define target size [X Y lenght(x) length(y)]
% target=[30 10 2 2];

%Generate anomalies here
% Gaussian
% model=0.5*exp(-((X-26).^2+(Y-10).^2)/10) - 0.5*exp(-((X-34).^2+(Y-10).^2)/10);
% 
% Block model
% model(9:12,34:37)= 0.5;
% 
% model(9:12,26:29)= -0.5;
% 20x60 Three anomalies: 2 blocks + 1 gradient block
% load model_A 

% 20x60 Two anomalies: 1 Gaussian + 1 block (same line integral through
% center)
% load model_B   


% model=zeros(nX*nZ,1);
% model(:)=rho(:)-background;

m=zeros( nX * nY * nZ ,1);
m(:)=model(:);

% figure (1)
% imagesc (reshape(m,nZ,nX));
% caxis([-0.4 0.6])
% hold on

%Create receivers
% receiver=[20:2:40];
% plot(receiver,(zeros(1,length(receiver))+Zmax/dz),'g^',...
%     'LineWidth',2,'MarkerEdgeColor','k',...
%     'MarkerFaceColor','y','MarkerSize',10);
% hold on

%Commands to Plot data
DIPrange=30;

%Create ray paths
ndata=300;
G=zeros(ndata,nX*nZ);
count=0;
while count<ndata
    Tx_Z=Zo; %Flat topo for now
    Tx_X = abs(rand(1))*Xmax;
    dip = rand(1)*DIPrange * (-1)^count;
    %length=sqrt(Zmax*tand(dip)^2+Zmax^2);
    hit=round(Tx_X+Zmax*tand(dip));
 
    if sum(receiver==hit)==1
        
        plot([Tx_X Tx_X+Zmax*tand(dip)],[Tx_Z Zmax],'k:');
        hold on
        
        count=count+1;
        TX(count)=Tx_X;
        DIP(count)=dip;
        
        G(count,:)=getG(Xo,Zo,nX,nZ,dX,dZ,Tx_X,Tx_Z,dip*pi/180);
        
        % G2(count,:)=getG2(Xmax,Zmax,hit(1),Zo,dx,dz,dip);

    end
    
end

set(gca,'YDir','reverse')

%Create data matrix 
data=G*(m);

ndata=length(data);
%Corrupt with 5% random noise
% d = awgn(data,-12.5);
noise = (max(abs(data))).*randn(ndata,1);
% d = data + noise;
model=reshape(model,nZ*nX,1);

save ('data/kernel','G');
save ('data/model','model');
save ('data/data','data');
save ('data/mesh','mesh');
save ('data/noise','noise');
% save ('data/DIP','DIP');
% save ('data/TX','TX');

% save ('data/DIPrange','DIPrange');
% save('data/d_DipX','d_DipX');
% save ('kernel2','G2');

%Map data as a function of receiver and angle
% [data_plot]=interp_data(data);
% [d_plot]=interp_data(d);
% Plotting original vs corrupted data
% figure (2)
% ylim([-5 5])
% plot(1:length(data),d,1:length(data),data)
    
% figure (1)
% xlabel('\bfEasting (m)')
% ylabel('\bfDepth (m)')
% colorbar;
% 
% figure;imagesc(receiver,ang,data_plot)
% colorbar;
% xlabel('\bfReceiver position')
% ylabel('\bfDip Angle')
% axis equal
% axis tight
% 
% figure;imagesc(receiver,ang,d_plot)
% colorbar;
% xlabel('\bfReceiver position')
% ylabel('\bfDip Angle')
% axis equal
% axis tight