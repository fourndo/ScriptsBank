 %getpath
%Compute N number of randomely generated ray paths.
%First determine if ray intersects a receiver
%If NO, ray path is only plotted.
%If YES, ray path populate a data matrix, then plotted

close all
clear all

addpath ../data/
addpath functions/

%length of land surveyed

Xo=0;
Yo=0;
Zo=0;
dz=1;
dx=1;

dX=ones(1,60)*dx;
dZ=ones(1,20)*dz;
nX=length(dX);
nZ=length(dZ);

mesh{1}=[nX nZ];
mesh{2}=[Xo Zo];
mesh{3}=dX;
mesh{4}=dZ;

mcell = nX * nZ;

for ii=1:nX
    xvec(ii)=sum(dX(1:ii));
end

for ii=1:nZ
    zvec(ii)=sum(dZ(1:ii));
end

Xmax = Xo + xvec(end);
Zmax = Zo - zvec(end);

[X,Y]=meshgrid(xvec,zvec);

model=zeros(nZ,nX);
%Define target size [X Y lenght(x) length(y)]
% target=[30 10 2 2];

%Generate anomalies here
% Gaussian
model=0.5*exp(-((X-26*dx).^2+(Y-zvec(end)/2).^2)/(10*dx)) - 0.5*exp(-((X-32*dx).^2+(Y-zvec(end)/2).^2)/(10*dx));
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

m=zeros( mcell ,1);
m(:)=model(:);

figure (1)
imagesc (reshape(m,nZ,nX));
caxis([-0.4 0.6])
hold on

%Create receivers
obs_loc(:,1) = [20:2:40]'*dx;
obs_loc(:,2) = Zmax+dz/2;

% plot(receiver,(zeros(1,length(receiver))+Zmax/dz),'g^',...
%     'LineWidth',2,'MarkerEdgeColor','k',...
%     'MarkerFaceColor','y','MarkerSize',10);
% hold on

%Commands to Plot data
DIPrange=20;

d_angl= 10*pi/180;
phi=[(-DIPrange*pi/180):d_angl:d_angl 0:(d_angl):(DIPrange*pi/180)];
% theta=[0:-(d_angl):-pi d_angl:d_angl:(pi-d_angl)];

% Descritized topography (all flat for now)
ZZ= ones(1,nX) * Zo;

%Create ray paths
ndata = length(phi)*size(obs_loc,1);

g=zeros(ndata,mcell);

count=0;
for ii=1:size(obs_loc,1)
    
    for jj= 1:length(phi)
        
    
        count = count+1;
        
        g(count,:)=getG_2D(obs_loc(ii,:),phi(jj),Xo,Zo,dX,dZ,ZZ);
        

    end
end

G = sparse(g);


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