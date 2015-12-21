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
load model
load noise
d = model + noise;

% model=m;
load ModelB_1D_Beta1e2.mat
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
         phid(ip,iq,il) = finalphid(iq,ip,il);
         
      end
   end
end

L2vec=zeros(1,np*nq*nl);
L2vec(:)=L2(:);
AA=sort(L2vec);
AA(81)
AA(162)

[p,q,l] = meshgrid(pvec,pvec,lvec);

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
caxis([0.25 0.6])
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

mcell=length(model);
rec = finalmodels(1,1,1,:);
rec = rec(:);
figure;plot(1:mcell,rec,'LineWidth',2)
hold on
plot(1:mcell,model,'g',1:mcell,d,'r*')
title(['\bfL2 =',num2str(L2(1,1,1)),' \bflp=0 lq=0 \lambda=0.5'])
ylim([-0.8 0.8])

rec = finalmodels(1,1,9,:);
rec = rec(:);
figure;plot(1:mcell,rec,'LineWidth',2)
hold on
plot(1:mcell,model,'g',1:mcell,d,'r*')
title(['\bfL2 =',num2str(L2(1,1,9)),' \bflp=0 lq=0 \lambda=2'])
ylim([-0.8 0.8])

rec = finalmodels(9,1,1,:);
rec = rec(:);
figure;plot(1:mcell,rec,'LineWidth',2)
hold on
plot(1:mcell,model,'g',1:mcell,d,'r*')
title(['\bfL2 =',num2str(L2(9,1,1)),' \bflp=2 lq=0 \lambda=0.5'])
ylim([-0.8 0.8])

rec = finalmodels(9,1,9,:);
rec = rec(:);
figure;plot(1:mcell,rec,'LineWidth',2)
hold on
plot(1:mcell,model,'g',1:mcell,d,'r*')
title(['\bfL2 =',num2str(L2(9,1,9)),' \bflp=2 lq=0 \lambda=2'])
ylim([-0.8 0.8])

rec = finalmodels(9,9,1,:);
rec = rec(:);
figure;plot(1:mcell,rec,'LineWidth',2)
hold on
plot(1:mcell,model,'g',1:mcell,d,'r*')
title(['\bfL2 =',num2str(L2(9,9,1)),' \bflp=2 lq=2 \lambda=0.5'])
ylim([-0.8 0.8])

rec = finalmodels(9,9,9,:);
rec = rec(:);
figure;plot(1:mcell,rec,'LineWidth',2)
hold on
plot(1:mcell,model,'g',1:mcell,d,'r*')
title(['\bfL2 =',num2str(L2(9,9,9)),' \bflp=2 lq=2 \lambda=2'])
ylim([-0.8 0.8])

rec = finalmodels(1,9,1,:);
rec = rec(:);
figure;plot(1:mcell,rec,'LineWidth',2)
hold on
plot(1:mcell,model,'g',1:mcell,d,'r*')
title(['\bfL2 =',num2str(L2(1,9,1)),' \bflp=0 lq=2 \lambda=0.5'])
ylim([-0.8 0.8])

rec = finalmodels(1,9,9,:);
rec = rec(:);
figure;plot(1:mcell,rec,'LineWidth',2)
hold on
plot(1:mcell,model,'g',1:mcell,d,'r*')
title(['\bfL2 =',num2str(L2(1,9,9)),' \bflp=0 lq=2 \lambda=2'])
ylim([-0.8 0.8])

minRMS=find(RMS==min(min(min(RMS))));
[iip,iiq,iil]=ind2sub(size(RMS),minRMS);
rec = finalmodels(iip,iiq,iil,:);
rec = rec(:);
figure;plot(1:mcell,rec,'LineWidth',2)
hold on
plot(1:mcell,model,'g',1:mcell,d,'r*')
% ID=['Smallest RMS p=',num2str(pvec(iip)),', q = ',num2str(qvec(iiq)),', lambda= ',num2str(lvec(iil))'];
title(['Smallest RMS p=',num2str(pvec(iip)),', q = ',num2str(qvec(iiq)),', lambda= ',num2str(lvec(iil))])

ylim([-0.8 0.8])

minL2=find(L2==min(min(min(L2))));
[iip,iiq,iil]=ind2sub(size(L2),minL2);
rec = finalmodels(iip(1),iiq(1),iil(1),:);
rec = rec(:);
figure;plot(1:mcell,rec,'LineWidth',2)
hold on
plot(1:mcell,model,'g',1:mcell,d,'r*')
title(['Smallest L2',num2str(L2(1,9,9)), 'p=',num2str(pvec(iip(2))),', q = ',num2str(qvec(iiq(2))),', lambda= ',num2str(lvec(iil(2)))])

ylim([-0.8 0.8])

fid1=fopen('UBC_L2_modelD.dat', 'w');
for iq = 1:nq
   for ip = 1:np
      for il = 1:nl
         
            fprintf(fid1,'%8.5f\n', L2(ip,iq,il));
      end
   end
end
fclose(fid1);
%% Can also use slicomatic:
sliceomatic(L2,pvec,qvec,lvec);