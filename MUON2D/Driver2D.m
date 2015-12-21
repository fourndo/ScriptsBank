close all
clear all

addpath data/
addpath functions/

load model
load kernel
load noise
load mesh
load data

dX=mesh{3};
dZ=mesh{4};
nX=length(dX);
nZ=length(dZ);

% 2D problem
mcell=length(model);
ndata=size(G,1);

%% Create Noisey data:

% Generate raw data
data= G * model;

pct_noise = 0.1;
noise = (pct_noise.*max(abs(data))).*randn(ndata,1);
Wd = diag(1./abs(noise));
d = data + noise;

% save('noise','noise')
%% Sensitivity Matrices:
% Normalize d and G by standard deviation:
target = sum((Wd*(data - d)).^2);
G = Wd*G;
GtG = G'*G;
d = Wd*d;

phid = @(x) sum((G*x - d).^2);


%Create derivative matrices
Wx=getWx_v3(mcell,nX,nZ,dX,dZ);
Wz=getWz_v3(mcell,nX,nZ,dX,dZ);


Lx = 1;%min(dX); % Length Scale
Lz = 1;%min(dZ); % Length Scale

% Compute dimension vectors
[ws,wx,wz]=getdim(mcell,dX,dZ);

alphaX = Lx^2;
alphaZ = Lz^2;


lvec= [0.0 1.0 1.5 2.0 2.5 3.0 3.5 4.0];
pvec= [0 0.25 0.5 0.75 1.0 1.25 1.5 1.75 2]; 
qvec= [0 0.25 0.5 0.75 1.0 1.25 1.5 1.75 2]; 

nl=length(lvec);
nq=length(qvec);
np=length(pvec);

tic
for pp=9%1:length(pvec)
    
for qq=3%1:length(qvec)   
    
for ll=4%1:length(lvec)
 
    beta = 5e+2;
     invmod = ones(size(model))*1e-4;
     fprintf('Iteration %i of %i.\n',sub2ind([nl,nq,np],ll,qq,pp),np*nq*nl); 
    for ii=1:100
   

            if ii==1                 %First iteration

                lambda = 0.0;

                modfunc= phim_v3(invmod,2,Wx,Wz,ws,wx,wz,alphaX,alphaZ,lambda,2,ii);
                
            else
                
                lambda = lvec(ll);

                modfunc= phim_v3(invmod,pvec(pp),Wx,Wz,ws,wx,wz,alphaX,alphaZ,lambda,qvec(qq),ii);
            end


    % CG solver
    A=(GtG + beta(end)*(modfunc));
    RHS=(G')*d;
    % [invmod,iter]=CG_solver(invmod,A,RHS);
    invmod=A\RHS;

    phi_d(ii) = phid(invmod);


        if (phi_d(ii) <= target) && (ii > 3)
           break;
        else
           if phi_d < target*2
              beta = 0.85*beta;
           else
              beta = 0.5*beta;
           end
        end
    
% figure(1);imagesc(reshape(invmod,20,60))
% caxis([0 0.75])
    end



         finalphid(pp,qq,ll) = phi_d(end);
         finalmodels(pp,qq,ll,:) = invmod;

end

end

end

fprintf('Results took %6.2f minutes.\n',toc/60);
save('2Dresults_model_E_scaleALL.mat','finalmodels','lvec','pvec','qvec','finalphid');
figure;imagesc(reshape(invmod,nZ,nX))
