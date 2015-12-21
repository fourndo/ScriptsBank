%Inverting for compact model - Iterative method
%Dominique Fournier
%October 15, 2011

clear all
close all

addpath data/
addpath functions/
addpath output/

load d;
load kernel;
load model;
load mesh;
load Tx;
load Rx;
load Wd;

% d=d-mean(d);

dX=mesh{2};
dZ=mesh{3};
nX=length(dX);
nZ=length(dZ);


% colorbar

mcell=nX*nZ;
ndata=size(G,1);


%Modify data and model to substract background density
% model=zeros(nZ,nX);
% model(:)=m(:);

m0=zeros(mcell,1);

%% Compute compact model using iterative method


%Compact function (I for first iteration)
WctWc=eye(mcell+1,mcell+1);
WctWc(1)=0;

%Initial model
model_comp=ones(mcell+1,1)*0e-1;
      
alphaX= 1e+0;      %Smoothness
alphaZ= 5e+0;      %Smoothness
alphaS= 0.0e+0;      %Smallness
alphaC= 0.0;      %Compact term (0 for first iteration)
gamma=1.0;                 %Exponent for compact term
epsilon=1.0e-11;               %Small term in compact function

beta_in=1e+0;   %Trade-off parameter
count=10;
iter=500;
lambda=4.0e-1;

%Create derivative matrices
Wx=getWx2D(mcell,nX,nZ,dX,dZ);
Wz=getWz2D(mcell,nX,nZ,dX,dZ);
WxtWx=(Wx')*Wx;
WztWz=(Wz')*Wz;


%Weighted on reference model
WtW= alphaX * WxtWx + alphaZ * WztWz;

clear  Wx Wz WxtWx WztWz

%Constant terms
GtG=G'* Wd'*Wd*G;
RHS= (G')* Wd' * Wd * d;

for oo=1:length(beta_in)

    count=1;
    iter=500;
    model_comp=ones(mcell+1,1)*1e-4;


    
    while count<=5 || iter>=10

    beta= beta_in;% / 1.01^count;      %Smoothness
    
    phi_xyz(count)= model_comp' * WtW * model_comp;
    phi_c(count)= model_comp' * WctWc * model_comp;

        if count>1
            
%             beta(count)=beta_in/(1.5*count);%*phi_d(count)/phi_d(1);
%               alphaX=0;  
            for ii=2:mcell+1

                WctWc(ii,ii)=1./(model_comp(ii).^2+epsilon).^(gamma);
                alphaC(count)=lambda * phi_xyz(count)/phi_c(count);%Weight on compactness
    %             alphaC(count)=1/max(WctWc*model_comp);
            end
        end



    modfunc=alphaC(count)*WctWc+alphaX*WtW;

%     beta(oo)=beta(oo)/1.1;



    % CG solver
    A= (GtG + beta(oo) * (modfunc) );
    
    [model_comp,iter]= conjgrad(model_comp,A,RHS);
    
%     if count==1
%      figure(1)
%     imagesc(reshape(model_comp(2:end),nZ,nX))
%     caxis([-0.1 0.1]) 
%     pbaspect([0.5 1 1])
%     title('Smooth model')
%     end
    
    d_pred=G*(model_comp);
    
figure(oo)
    imagesc(reshape(model_comp(3*nZ+2:end-3*nZ),nZ,nX-6))
%     caxis([-0.1 0])
    pbaspect([0.5 1 1])
    colorbar


    count=count+1;
    
    
        if count>100
            break
        end
    end

    
    
    figure(oo+10)
    scatter(RX_bin, TX_bin, 50, d_pred, 'filled');
    caxis([-10 -5])
    
    % model_out=reshape(model_comp(2:end),nZ,nX);
    filename=sprintf('model_beta%4.2f_lambda%d.dat',beta_in(oo),lambda);
    save(filename,'-ascii','model_comp')
    
    % RMS(oo)= (sum (( reshape(model,mcell,1) - model_comp(2:end) ) .^2) / mcell) ^(0.5);
    phi_m(oo) = (model_comp)' * (modfunc) * (model_comp);
    phi_d(oo) = sum( (G * model_comp - d) .^2);
end

% set(gca,'YTick',1.5:nZ/15:nZ)
% set(gca,'YTickLabel',10:20:300)
% set(gca,'XTick',0:nX/10:nX)
% set(gca,'XTickLabel',-20:5:180)
% xlabel('Distance (m)')
% ylabel('Depth (m)')
% pbaspect([0.5 1 1])

% title('\bfRIM2D_v1')
% colormap('jet')
% xlabel('\bfx-axis')
% ylabel('\bfDepth')
% axis([1 nX 1 nZ])
% set(gca,'YDir','reverse')

% colorbar
% figure; imagesc(reshape(m(2:end),nZ,nX))

    figure
    scatter(RX_bin, TX_bin, 50, d, 'filled');
    caxis([-10 -5])
% figure;loglog(phi_m,phi_d);

% figure;plot(lambda,RMS)
figure;
imagesc(reshape(model,nZ,nX))
% caxis([-0.1 1.1])
% colorbar
% figure(2)
% imagesc(reshape(model_comp,nZ,nX))
% caxis([-0.1 1.1])
% title('\bfRIM2D_v1')
% colormap('jet')
% xlabel('\bfx-axis')
% ylabel('\bfDepth')
% axis([1 nX 1 nZ])
% set(gca,'YDir','reverse')
% colorbar

% AA=sum(G);
% figure;imagesc(reshape(AA(2:end),nZ,nX))

