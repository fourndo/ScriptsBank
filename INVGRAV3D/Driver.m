% Generate model and observations for 3D gravity
% Dominique Fournier 2013/01/23
close all
clear all

addpath data
addpath functions


iter_max = 15;

% For now all the cells have dimension 1x1x1
dX = [6 3 1.5 1.25 ones(1,14) 1.25 1.5 3 6]*50;
dY = [6 3 1.5 1.25 ones(1,14) 1.25 1.5 3 6]*50;
dZ = [ones(1,10) 1.25 1.5 3 6]*50;

nX= length(dX);
nY= length(dY);
nZ= length(dZ);

mcell= nX * nY * nZ;

X0 = 0;
Y0 = 0;
Z0 = 0;


%3D density contrast model

figure (1)

[model] = zeros(nZ,nX,nY);

model(3:5,10:12,10:12)=1.0;

topo_model = ones(mcell,1);

% [xx,yy,zz]=meshgrid(0:2:6,0:2:6,0:2:6);
% fx = exp(-((xx-3).^2+(yy-3).^2+(zz-3).^2)/15);
% model(3:6,12:15,12:15)=fx;
load model_intrusive
model = model_intrusive;

% Vectorized model
m=reshape(model,mcell,1);

save('model','model');
save('data/model.den','-ascii','m');

%Create data points
% Center observations for simplicity
cellx = cumsum(dX); celly = cumsum(dY); nx = length(cellx); ny = length(celly);
[ObsX, ObsY, ObsZ] = meshgrid(cellx((floor(nx/2-6)):(ceil(nx/2+6)))+0.5,celly((floor(ny/2-6)):(ceil(ny/2+6)))+0.5,2);

ndata = size(ObsX, 2) * size(ObsY, 1);

count = 1;

% Initiate G matrix
G=zeros(ndata,mcell);

% Compute depth weigthing matrix
% mode 0: distance weigthing , 1: depth weighting
% pow: Power of the exponentiel decay (default 2 for grav, 3 for mag)
Wr=zeros(1,mcell);
mode=0;
pow=2;

%% Compute forward operator (G) and distance weighting
for ii=1:ndata;
       [G(count,:),wr,V] = forwardGrav(nX, nY, nZ, X0, Y0, Z0,...
            dX, dY, dZ, ObsX(ii), ObsY(ii), ObsZ(ii));
        
        Wr = Wr + wr;
count = count + 1;
end

% Square root for the sum of the squares
% Plus another square root of result because inside the objective function,
Wr=Wr.^(1/2);

% Normalize depth weighting with the largest value

Wr = Wr./(V);

Wr = Wr./(max(Wr));

Wr=Wr.^(1/2);

% Wr=Wr.^(1/2);

IWr = spdiags(1./Wr',0,mcell,mcell);
Wr = spdiags(Wr',0,mcell,mcell);


%Create data matrix
data = G * m ;

%% Generate data and noise
%Corrupt with 5% random noise
pct_noise = 0.05;
noise = (pct_noise.*max(abs(data))).*randn(ndata,1);
save ('noise','noise')
Wd = diag(1./abs(noise));
d = data + noise;


% Normalize d and G by standard deviation:

G = Wd * G * IWr;
GtG = G'*G;
% d = data;
% G = G * IWr;
target = sum((Wd*(data - d)).^2);

d = Wd * d;

RHS=(G')*d;

% save ('original','data');
% save('data/data.mat','data');
% save('data/kernel.mat','G');
% save('data/model.mat','m');
% save('data/Wr.mat','Wr');
% save ('kernel2','G2');

d_obs = reshape(Wd\d, size(ObsX,2), size(ObsY,1));


figure (1)
imagesc(d_obs)
xlabel('\bfEasting (m)')
ylabel('\bfNorthing (m)')


%% Inversion
Lx = 4*min(dX);
Ly = 4*min(dY);
Lz = 4*min(dZ);


alphaC = 1;
alphaX = Lx^2 * alphaC;
alphaY = Lx^2 * alphaC;
alphaZ = Lx^2 * alphaC;


Wx = getWx_3D(mcell,dX,dY,dZ);
Wy = getWy_3D(mcell,dX,dY,dZ);
Wz = getWz_3D(mcell,dX,dY,dZ);
Ws = getWs_3D(mcell,dX,dY,dZ);

% dXx = get_dXx(nX, nY, nZ, dX, dY, dZ,topo_model);
% dYy = get_dYy(nX, nY, nZ, dX, dY, dZ,topo_model);
% dZz = get_dZz(nX, nY, nZ, dX, dY, dZ,topo_model);
% dV = get_dV(nX, nY, nZ, dX, dY, dZ);

pvec= [0.0 1.0 2.0]; 
qvec= [0.0 1.0 2.0]; 
lvec= [0.5 1.0 1.5];

nl=length(lvec);
nq=length(qvec);
np=length(pvec);

finalmodels = zeros(length(pvec),length(qvec),length(lvec),mcell);
finalphid = zeros(length(pvec),length(qvec),length(lvec));


for pp=3%1:length(pvec)
    
for qq=1%1:length(qvec)   
    
for ll=3%1:length(lvec)
 
beta = 5e-1;
fprintf('Iteration %i of %i.\n',sub2ind([nl,nq,np],ll,qq,pp),np*nq*nl); 
invmod=ones(mcell,1)*1e-6; 
    for ii=1:iter_max



            if ii==1                 %First iteration

                modfunc= phim_3D(invmod,Wx,Wy,Wz,Ws,Wr,alphaX,alphaY,alphaZ,alphaC,2,2,0,ii);

            else              

                modfunc= phim_3D(invmod,Wx,Wy,Wz,Ws,Wr,alphaX,alphaY,alphaZ,alphaC,pvec(pp),qvec(qq),lvec(ll),ii);

            end
%     modfunc = (alphaX*(Wx*Ws)'*(Wx*Ws) + alphaY*(Wy*Ws)'*(Wy*Ws) + alphaZ*(Wz*Ws)'*(Wz*Ws) + alphaC*(Ws)'*(Ws));

    % CG solver
    A=(GtG + beta(end) * ((modfunc)));

%     [invmod]=CG_Solver(invmod,A,RHS);
    invmod=A\RHS;
       close (figure(3))
     phi_d(ii) = sum((G*(invmod)-d).^2);
%      phi_m(ii) = invmod' * modfunc * invmod
%     figure(2)
%     d_pred = reshape(Wd\G*invmod, size(ObsX,2), size(ObsY,1));
%     imagesc(d_pred);
    
    
        if (phi_d(ii) <= target) && (ii > 2)
           break;
        else
           if phi_d(end) < target*2
              beta = 0.75*beta;
           else
              beta = 0.5*beta;
           end
        end

        fprintf('Iteration %i of %i completed, phi_d: %f\n', ii , iter_max,phi_d(ii) )
         model_out = invmod;
%     if pp==1 || pp==length(pvec) || ll==1 || ll==length(lvec) || qq==1 || qq==length(qvec)
    save(['results/lambda' num2str(lvec(ll)) 'q' num2str(qvec(qq))  'p' num2str(pvec(pp)) '.den'],'-ascii','model_out')
    
%     RegInvMat = (A) \ ( G'*(Wd'*Wd) );

% R = ( sum(RegInvMat.^2,2) ).^(-0.5);
% save('R_l1_pct5_v2.dat','-ascii','R')
% ResDenMat = ( norm(d).^-1 ) * spdiags( R , 0 , mcell , mcell );

    end

%% Computing model resolution matrix

% MRM = (GtG + beta*(Wx'*Wx + Wy'*Wy + Wz'*Wz ))\GtG;
% 


% MRM_top=MRM(:,3223);
% MRM_bottom=MRM(:,3225);
% MRM_low=MRM(:,3530);
% MRM_high=MRM(:,3220);
% 
% save('R_l1.dat','-ascii','R')
% save('MRM_bottom_l1.dat','-ascii','MRM_bottom')
% save('MRM_low_l1.dat','-ascii','MRM_low')
% save('MRM_high_l1.dat','-ascii','MRM_high')

% Rall = sum(abs(MRM),2);
% save(['results/Rall_l' num2str(lvec(ll)) '_lq' num2str(qvec(qq))  '_lp' num2str(pvec(pp)) '.dat'],'-ascii','Rall')

%%
% imagesc(reshape(invmod,nZ,nX));
% caxis([-0.4 0.6])
%     model_out = IWr*invmod;
% %     if pp==1 || pp==length(pvec) || ll==1 || ll==length(lvec) || qq==1 || qq==length(qvec)
%     save(['data/lambda' num2str(lvec(ll)) 'q' num2str(qvec(qq))  'p' num2str(pvec(pp)) '.den'],'-ascii','model_out')
%     end
    
    finalphid(pp,qq,ll) = phi_d(ii);
    finalmodels(pp,qq,ll,:) = (IWr)*model_out;

    
    
end
         
         
end

end

save('3DMAG_lplq.mat','finalmodels','lvec','pvec','qvec','finalphid');


%%
%Create UBC mesh file
fid1=fopen('data/UBC_mesh.msh', 'w');
fprintf(fid1, '%i %i %i\n', nX, nY, nZ);
fprintf(fid1, '%i %i %i\n', X0, Y0, Z0);

for jj=1:nX
    fprintf(fid1, '%4.2f ', dX(jj));    
end

fprintf(fid1,'\n');

for ii=1:nY
           fprintf(fid1,'%4.2f ', dY(ii));
end

fprintf(fid1, '\n');

for kk=1 : nZ
       fprintf(fid1, '%4.2f ', dZ(kk));
end

fclose(fid1);

%Create UBC observation file
count = 1;
fid3 = fopen('data/UBC_obs.obs','w');
fprintf(fid3,'%i\n',ndata);

for ii=1:ndata

    fprintf(fid3,'%4.2f %4.2f %4.2f %e\n',ObsX(ii),ObsY(ii),ObsZ(ii),Wd(ii,ii)\d(ii));
    count = count + 1;

end
fclose(fid3);

figure(2)
d_pred = reshape(Wd\G*invmod, size(ObsX,2), size(ObsY,1));
imagesc(d_pred);