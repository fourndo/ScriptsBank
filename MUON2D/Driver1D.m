%% House cleaning/set up
close all
clear all

addpath data/
addpath functions/

%% Load correct model/data:
% load model_A
% load data_A
load model_C
load data_C

load kernel

%% Create problem
n=size(G,1);
% dX=ones(1,60);
% dZ=ones(1,20);
dX = 1;
dZ = 1;
nX=60;
nZ=20;

% 2D problem
mcell=length(model);
% Create unit model of ones
unit = ones(size(model));
% ndat = length(d);

%  Create derivative matrices:
% Wz:
Wz = zeros((nX*(nZ-1)),mcell);

p = 0;
for ii=1:(nX*(nZ-1))
   jj = ii + p;
   if (rem(ii,nZ-1)~=0)
      Wz(ii,jj) = -1;
      Wz(ii,jj+1) = 1;
   else
      p=p+1;
      Wz(ii,jj) = -1;
      Wz(ii,jj+1) = 1;
   end
end
% Wx:
Wx = zeros(nZ*(nX-1),mcell);
for ii=1:(nZ*(nX-1))
   Wx(ii,ii) = -1;
   Wx(ii,ii+nZ) = 1;
end

%Ws:
Ws = eye(mcell);

%% Create Noisey data:
data = G*model;
ndat = length(data); 
pct_noise = 0.1;
noise = (pct_noise.*max(abs(data))).*randn(ndat,1);
Wd = diag(1./abs(noise));
d = data + noise;
% example = 'data_B';
% filename = ['data' filesep example '.mat'];
% save(filename,'d','Wd','data');

%% Sensitivity Matrices:
% Normalize d and G by standard deviation:
target = sum((Wd*(data - d)).^2);
G = Wd*G;
GtG = G'*G;
d = Wd*d;
phid = @(x) sum((G*x - d).^2);

%% Length scales to make more realistic (4 cell widths)
% alphaX = 1; % Smoothness term
% alphaZ = 1; % Smoothness term
% alphaS = 1; % Smallness term
Lx = 4*dX; % Length Scale
Lz = 4*dZ; % Length Scale

alphaS = 1.0;
alphaX = Lx^2 * alphaS;
alphaZ = Lz^2 * alphaS;
% Update Ws so that Ws'Ws = alphaS*Sw'*Ws
Ws = sqrt(alphaS).*Ws;
%% For only 1 run:
% p = 1;
% q = 2;
% lambda = 1;

%% Set up for-loop variables
pvec = 1.4;%linspace(0,2,9);
qvec = 0.6;%linspace(0,2,9);
lvec = 1.7;%linspace(0.5,2,10);
np = length(pvec);
nq = length(qvec);
nl = length(lvec);
finalmodels = zeros(np,nq,nl,mcell);
finalphid = zeros(np,nq,nl);
tic;
% --- For Testing purposes: ---
% figure(1);
% figure(2);
% drawnow;
% -----------------------------
for ip = 1:np;
   p = pvec(ip);
   for iq = 1:nq;
      q = qvec(iq);
      for il = 1:nl;
         lambda = lvec(il);
         %% Inversion:
         beta = 20.0;
         invmod = ones(size(model))*1e-4;
         fprintf('Iteration %i of %i.\n',sub2ind([nl,nq,np],il,iq,ip),np*nq*nl);
         for ii=1:40
            % CG solver
            if (ii == 1)
               modfunc = phim(unit,2,Wx,Wz,alphaX,alphaZ);
            else
               % Stop here: check Rs in q by going into function
               modfunc = phim(invmod,p,Wx,Wz,alphaX,alphaZ,lambda,Ws,q);
            end
            A = (GtG + (beta*modfunc));
            RHS = (G')*d;
                        
            % invmod = pcg(A,RHS,1e-5,500);
            invmod = A\RHS;
            % [invmod,iter] = CG_solver(invmod,A,RHS);
            
            phi_d = phid(invmod);
            phi_m = invmod'*modfunc*invmod;
            phi = phi_d + beta*phi_m;
            
            if (phi_d <= target) && (ii > 2)
               break;
            else
               if phi_d < target*2
                  beta = 0.85*beta;
               else
                  beta = 0.5*beta;
               end
            end
            
            %             phi_d(ii) = phid(invmod);
            %             phi_m(ii) = invmod'*modfunc*invmod;
            %             phi(ii) = phi_d(ii) + beta*phi_m(ii);
            %             if (phi_d(ii) <= target) && (ii > 2)
            %                break;
            %             else
            %                beta = beta/2;
            %             end
         end
         finalphid(ip,iq,il) = phi_d;
         finalmodels(ip,iq,il,:) = invmod;
% --- For testing purposes: ---
%          figure(2);
%          subplot(np,nq,sub2ind([nq,np],iq,ip));            
%          plot(phi_d(2:end),'k-');
%          str2 = sprintf('P=%3.1f; Q=%3.1f; Phid=%7.2f',p,q,phi_d(end));
%          title(str2);
%          drawnow;
%          
%          figure(1)
%          subplot(np,nq,sub2ind([nq,np],iq,ip));
%          imagesc(reshape(invmod,nZ,nX));
%          caxis([-0.4 0.6])
%          str2={'\phi d' phi_d(end)};
%          str1={'lq' q};
%          str3={'lp' p};
%          str4={'beta' beta};
%          text(45,13,str1)
%          text(45,18,str2)
%          text(2,13,str3)
%          text(2,18,str4)
%          drawnow;
%          clear phi_d;
% ------------------------------
      end
   end
end

%%
fprintf('Results took %6.2f minutes.\n',toc/60);
save('2Dresults_model_A.mat','finalmodels','lvec','pvec','qvec','finalphid');

%% Display results
% figure;
% title('\Phi_m,\Phi_d,\Phi per iteration');
% semilogy((2:ii),phi_d(2:end),'.-r');
% hold on;
% semilogy((2:ii),phi_m(2:end),'.-k');
% semilogy((2:ii),phi(2:end),'.-g');
% xlabel('Iteration');
% legend('\Phi_d','\Phi_m','\Phi');

% test1 = (Wz'*Wz)*model;
% test2 = (Wx'*Wx)*model;
% figure;
% subplot(311)
% imagesc(reshape(test1,nZ,nX));
% title('Z-derivative');
% caxis([-0.4 0.6])
% axis equal;
% axis tight;
% subplot(312)
% imagesc(reshape(model,nZ,nX));
% caxis([-0.4 0.6])
% title('True');
% axis equal;
% axis tight;
% subplot(313)
% imagesc(reshape(test2,nZ,nX));
% title('X-derivative');
% caxis([-0.4 0.6])
% axis equal;
% axis tight;
