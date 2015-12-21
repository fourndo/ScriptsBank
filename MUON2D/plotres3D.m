clear all
close all

% Final results is in np,nq,nl,mcell format
%% Plot in this format:
%           Q
%          ^
%         /
%        /
%       / _ _ _ _ p
%       |
%       |
%       |
%       |
%       lambda (+ up)
%
%% Load results and set up axis vectors
addpath data/ sliceomatic/ 
addpath results/
load model_D
load 2Dresults_model_E_scaleALL
np = length(pvec);
nq = length(qvec);
nl = length(lvec);


%% Calculate model norm
L2 = zeros(np,nq,nl);
RMS = zeros(np,nq,nl);
L1 = zeros(np,nq,nl);
Lfro = zeros(np,nq,nl);
Linf = zeros(np,nq,nl);
RMS = zeros(np,nq,nl);
phid = zeros(np,nq,nl);
for ip = 1:np
   for iq = 1:nq
      for il = 1:nl
         rec = finalmodels(ip,iq,il,:);
         rec = rec(:);
         % Norms:
         L2(ip,iq,il) = norm((rec - model),2);
         RMS(ip,iq,il) = sqrt(sum(rec - model).^2/length(model));
         L1(ip,iq,il) = norm((rec - model),1);
         Lfro(ip,iq,il) = norm((rec - model),'fro');
         Linf(ip,iq,il) = norm((rec - model),Inf);
         % Phid:
         phid(ip,iq,il) = finalphid(ip,iq,il,end);
      end
   end
end

[p,q,l] = meshgrid(pvec,qvec,lvec);

%% Plot results:
figure(1)
subplot(2,3,1)
slice(p,q,l,phid,1,1,1)
title('\Phi_d');
xlabel('Q');
zlabel('\lambda');
ylabel('P');
axis equal
axis tight
colorbar

subplot(2,3,2)
slice(p,q,l,L2,1,1,1)
% caxis([0.8 3])
title('L-2 norm');
xlabel('Q');
zlabel('\lambda');
ylabel('P');
axis equal
axis tight
colorbar

subplot(2,3,3)
slice(p,q,l,L1,1,1,1)
title('L-1 norm');
xlabel('Q');
zlabel('\lambda');
ylabel('P');
axis equal
axis tight
colorbar

subplot(2,3,4)
slice(p,q,l,RMS,1,1,1)
title('RMS model error');
xlabel('Q');
zlabel('\lambda');
ylabel('P');
axis equal
axis tight
colorbar

subplot(2,3,5)
slice(p,q,l,Linf,1,1,1)
title('L-\infty norm');
xlabel('Q');
zlabel('\lambda');
ylabel('P');
axis equal
axis tight
colorbar

subplot(2,3,6)
slice(p,q,l,Lfro,1,1,1)
title('Frobenius norm');
xlabel('Q');
zlabel('\lambda');
ylabel('P');
axis equal
axis tight
colorbar

rec = finalmodels(1,1,1,:);
rec = rec(:);
figure;imagesc(reshape(rec,20,60))
caxis([-.5 0.5])
title(['\bflp=0 lq=0 \lambda=' num2str(lvec(1)) ' l2norm:' num2str(L2(1,1,1))])

rec = finalmodels(1,1,end,:);
rec = rec(:);
figure;imagesc(reshape(rec,20,60))
caxis([-.5 0.5])
title(['\bflp=0 lq=0 \lambda=' num2str(lvec(end)) ' l2norm:' num2str(L2(1,1,end))])

rec = finalmodels(end,1,1,:);
rec = rec(:);
figure;imagesc(reshape(rec,20,60))
caxis([-.5 0.5])
title(['\bflp=2 lq=0 \lambda=' num2str(lvec(1)) ' l2norm:' num2str(L2(end,1,1))])

rec = finalmodels(end,1,end,:);
rec = rec(:);
figure;imagesc(reshape(rec,20,60))
caxis([-.5 0.5])
title(['\bflp=2 lq=0 \lambda=' num2str(lvec(end)) ' l2norm:' num2str(L2(end,1,end))])

rec = finalmodels(end,end,1,:);
rec = rec(:);
figure;imagesc(reshape(rec,20,60))
caxis([-.5 0.5])
title(['\bflp=2 lq=2 \lambda=' num2str(lvec(1)) ' l2norm:' num2str(L2(end,end,1))])

rec = finalmodels(end,end,end,:);
rec = rec(:);
figure;imagesc(reshape(rec,20,60))
caxis([-.5 0.5])
title(['\bflp=2 lq=2 \lambda=' num2str(lvec(end)) ' l2norm:' num2str(L2(end,end,end))])

rec = finalmodels(1,end,1,:);
rec = rec(:);
figure;imagesc(reshape(rec,20,60))
caxis([-.5 0.5])
title(['\bflp=0 lq=2 \lambda=' num2str(lvec(1)) ' l2norm:' num2str(L2(1,end,1))])

rec = finalmodels(1,end,end,:);
rec = rec(:);
figure;imagesc(reshape(rec,20,60))
caxis([-.5 0.5])
title(['\bflp=0 lq=2\lambda=' num2str(lvec(end)) ' l2norm:' num2str(L2(1,end,end))])

minRMS=find(RMS==min(min(min(RMS))));


minL2=find(L2==min(min(min(L2))));
[iip,iiq,iil]=ind2sub(size(L2),minL2);

rec = finalmodels(iip,iiq,iil,:);
rec = rec(:);
figure;imagesc(reshape(rec,20,60))
caxis([-.5 0.5])
title(['\bfMIN L2 lp=' num2str(pvec(iip)) ' lq= ' num2str(qvec(iiq)) '\lambda= ' num2str(lvec(iil)) ' l2norm:' num2str(L2(1,end,1))])

[iip,iiq,iil]=ind2sub(size(L2),minRMS);
rec = finalmodels(iip,iiq,iil,:);
rec = rec(:);
figure;imagesc(reshape(rec,20,60))
caxis([-.5 0.5])
title(['\bfMIN RMSlp=' num2str(pvec(iip)) ' lq= ' num2str(qvec(iiq)) '\lambda= ' num2str(lvec(iil)) ' l2norm:' num2str(L2(1,end,1))])

figure;imagesc(reshape(model,20,60))
% caxis([-.5 0.5])
title('Model_C')
colorbar
%% Can also use slicomatic:
sliceomatic(L2,pvec,qvec,lvec);caxis([0.8 3])