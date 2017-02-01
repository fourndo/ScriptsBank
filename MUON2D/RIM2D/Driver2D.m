close all
clear all

addpath data/
addpath functions/
addpath output/

load d;
load kernel;
load model;
load mesh;
load Tx;
load Rx;
% load Wd;

dX=mesh{2};
dZ=mesh{3};
nX=length(dX);
nZ=length(dZ);

% 2D problem
mcell=size(G,2);
ndata=size(G,1);

% data=G*m;

%% Create Noisey data:
data=G*m;

pct_noise = 0.1;
% noise = (pct_noise.*max(abs(data))).*randn(ndata,1);
% save ('noise','noise')
load noise
Wd = spdiags(1./abs(noise),0,ndata,ndata);
d = data + noise;

[Rx_interp,Tx_interp] = meshgrid (min(RxZ):max(RxZ),min(TxZ):max(TxZ));
dobs_interp = griddata(RxZ,TxZ,d,Rx_interp,Tx_interp);
figure(2);imagesc(((dobs_interp)))


target = sum((Wd*(data - d)).^2);

%% Sensitivity Matrices:
% Normalize d and G by standard deviation:

G = Wd*G;
GtG = G'*G;
d = Wd*d;


%Create derivative matrices
Wx= getWx(mcell,nX,nZ,dX,dZ);
Wz= getWz(mcell,nX,nZ,dX,dZ);
Ws= getWs(mcell,nX,nZ,dX,dZ);

%Initial model
invmod=ones(mcell,1)*1e-6;
                
lvec=[0.25 1.0 1.25 1.5 2];%10.^[-1:4];
pvec=[0 0.5 1.0 1.5 2]; 
qvec=[0 0.5 1.0 1.5 2]; 

nl=length(lvec);
np=length(pvec);
nq=length(qvec);

finalmodels = zeros(np,nq,nl,mcell);
finalphid = zeros(np,nq,nl);



Lx = 4*min(dX);
Lz = 4*min(dZ);

alphaC = 1;
alphaX = Lx^2 * alphaC;
alphaZ = Lz^2 * alphaC;

itermax = 40;

for pp=3%:np
    
    for qq=3%:nq   
    
        for ll=2%:nl
 
beta = 5e+2;
invmod=ones(mcell,1)*1e-6; 
fprintf('Iteration %i of %i.\n',sub2ind([nl,nq,np],ll,qq,pp),np*nq*nl);

            for ii=1:itermax

                % dphi_d=1;    
                count=1;

                % while count<3 



                    if ii==1                 %First iteration

                        modfunc= phim(invmod,Wx,Wz,Ws,alphaX,alphaZ,alphaC,2,2,1,1);

                    else

                        modfunc= phim(invmod,Wx,Wz,Ws,alphaX,alphaZ,alphaC,pvec(pp),qvec(qq),lvec(ll),1);

                    end

                modfunc(1)=0;
                % CG solver
                A=(GtG + beta*(modfunc));
                RHS=(G')*d;
                % [invmod,iter]=CG_solver(invmod,A,RHS);
                invmod=A\RHS;

                phi_d(ii) = sum(((G*(invmod)-d)).^2);


                    d_pred=Wd\G*(invmod);

                    figure(1)
                    imagesc(reshape(invmod(3*nZ+2:end-3*nZ),nZ,nX-6))
                    caxis([-0.5 0.5])
                    pbaspect([0.5 1 1])
                    colorbar
                %     sprintf('Step %i of iteration %i',count,ii)

                            if (phi_d(ii) <= target) && (ii > 2)
                               break;
                            else
                               if phi_d(ii) < target*2
                                  beta = 0.85*beta;
                               elseif phi_d(ii) < target*1
                                  beta = 0.75*beta;
                               else
                                  beta = 0.5*beta;
                               end
                            end
                % count=count+1;
                % end
                %     figure(ll+10)
                %     scatter(RX_bin, TX_bin, 50, d_pred, 'filled');
                %     caxis([-10 -5])

                    % model_out=reshape(model_comp(2:end),nZ,nX);
                %     filename=sprintf('model_beta%4.2f_lambda%d.dat',beta_in(oo),lambda);
                %     save(filename,'-ascii','invmod')

                    % RMS(oo)= (sum (( reshape(model,mcell,1) - model_comp(2:end) ) .^2) / mcell) ^(0.5);
                %     phi_m(ll) = (invmod)' * (modfunc) * (invmod);
                %     phi_d = sum( (G * invmod - d) .^2);

                fprintf('Iteration %i of %i completed, phi_d: %f\n', ii , itermax,phi_d(ii) )
                
            end


         finalphid(pp,qq,ll) = phi_d(end);
         finalmodels(pp,qq,ll,:) = invmod;
        end

    end

end
% ylim([0 0.1])
% for kk=1:size(RMS,1)
% alphaC_min_RMS(kk) = alphaC_in(RMS(kk,:)==min(RMS(kk,:)));
% end
save('2Dresults.mat','finalmodels','lvec','pvec','qvec','finalphid');

% figure;
% scatter(RxZ, TxZ, 50, d_pred, 'filled');
% caxis([-10 -5])

dpred_interp = griddata(RxZ,TxZ,d_pred,Rx_interp,Tx_interp);
figure;imagesc(((dpred_interp)))
      
draw_interp = griddata(RxZ,TxZ,data,Rx_interp,Tx_interp);
figure;imagesc(((draw_interp)));title('Raw data')

residual = (dpred_interp - dobs_interp);
figure;imagesc(((residual)));caxis([-.25 .25])

figure; imagesc(reshape(m(3*nZ+2:end-3*nZ),nZ,nX-6));
colormap(jet)
caxis([-0.5 0.5])
pbaspect([0.5 1 1])
colorbar


% figure;plot(pct_noise,alphaC_min_RMS)
% xlabel('% Noise')
% ylabel('\alphaC for min(RMS)')
% hold on

% Compute trend line
% K=[pct_noise'.^0 pct_noise'.^1 pct_noise'.^2  pct_noise'.^3];
% func=(K'*K)\K'*alphaC_min_RMS';
% 
% plot(pct_noise,func(1) + func(2)*pct_noise + func(3)*pct_noise.^2 + func(4)*pct_noise.^3,'g')
% BB=sort(invmod);
% wctwc=1./(BB.^2 + delta) .^p;
% wcm=alphaC(end)*wctwc.*BB;
% figure(2)
% plot(BB,wcm,'*')
% xlim([0 max(model)])
% ylabel('\bf\alphac*WctWc * model')
% xlabel('\bfmodel')
% axis([0 0.5 0 max(wcm)])
% title('\bfCompact term as a function of model')