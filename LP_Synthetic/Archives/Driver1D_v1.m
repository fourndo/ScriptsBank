close all
clear all

addpath data\
addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\Lp_norm_codes\Dual_norm_1D\functions
addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\Lp_norm_codes\Dual_norm_1D\sliceomatic\
% Set up 1D problem

%% INPUT VARIABLES

% Model space
nx = 300;

% Kernel
nk = 100; %number of frequencies

decay = -0.20;

basis = 0.5;

% Data noise
amp_pct = 0.02;
floor_pct = 0.02;

%% SCRIPT STARTS HERE
x = (pi/nx : pi/nx : pi) * nx; x = x(:);

mcell=length(x);
dx=ones(1,mcell) * abs((min(x)-max(x))/mcell);

% Create mdoel: Square and tanh functions
model = 0.5* exp(-abs((x*5/nx-pi*3.75).^2)).*tanh((x*10/nx-pi*3.75)); 

model(50:100)=0.25;

model=model(:);


figure;plot(x,model,'LineWidth',2);axis([x(1) x(end) -0.1 1]);hold on
grid on

% Generate kernel functions
G = zeros(nk,mcell);
wr = zeros(mcell,1);
figure;
for ii = 1 : nk 
    
    G(ii,1:mcell) = exp( decay * x  / ( nx * pi )  * (ii-1)/3) .*cos (basis * 2 * x / (nx) * (ii-1)/3);

    wr = wr + exp( decay * x  / ( nx * pi )  * (ii-1));
    
    plot(G(ii,:));hold on

end
grid on

wr = (wr / max(wr)).^0.5;

IWr = spdiags(1./wr,0,mcell,mcell);
Wr = spdiags(wr,0,mcell,mcell);
Wrx = spdiags(wr(1:end-1),0,mcell-1,mcell-1);

figure(1); plot(x,wr,'r');
legend('Model','Sensitivity weighting');
%% Generate raw data and corrupt with noise

% Total number of data
ndata = nk;

% Compute data
data= G * model;

% Corrupt data with Gaussian noise
rand_noise = rand(ndata,1);
save ('noise','rand_noise')
load('noise');

amp_noise = amp_pct * rand_noise .* abs(data);
floor_noise = floor_pct * rand_noise * max(abs(data));

noise = amp_noise + floor_noise;

d = data + noise;

wd =  abs(d)*0.05 + 0.1 * std(data);

figure; 
plot(data);
hold on
errorbar(d,wd,'r*')
legend('Data','Data+noise')
grid on

% Create uncertainty weighting matrix
Wd = spdiags(1./wd,0,ndata,ndata);

% Normalize d and G by standard deviation:
G = Wd * G ;

target = ndata;

% Weight data with uncertainty
d = Wd * d;

%Create weight matrices
[Wx,Vx] = getWx1D(nx,dx);
[Ws] = speye(mcell);
V = spdiags(sqrt(dx'),0,nx,nx);


% Global constant
as = 1.0 / min(dx) ^2;  %Smallnest term
ax = 1.0;               %Smoothness term                     

% Lp-norm parameters
pvec = 0:0.5:2;
qvec = 0:0.5:2;
lvec = 0.1:0.4:1.9;

epsilon = 1e-8;     %Small term in compact function

% Store all the final models
models_out = zeros(length(pvec),length(qvec),length(lvec),mcell);

counter = 1;  

objfunc = @(m,phi,b) sum( ( G * m - d ).^2 ) + ( m' * b * phi * m );

% Iterate over all the combination of norms as specified by the vectors
% pvec, qvec and lvec.

for ll= 1:length(lvec)

    for pp = 1:length(pvec)
        
        for qq = 1:length(qvec)
            
            % Message prompt
            head = ['lp: ' num2str(pvec(pp)) ' lq: ' num2str(qvec(qq)) ' psi: ' num2str(lvec(ll))];
            fprintf('Starting lp inversion %s\n',...
                head)
            
            invmod      = ones(mcell,1)*2e-1;       % Initial model       
                        
            phi_d       = sum((G*invmod - d).^2);   % Initial misfit
            scale_x     = ax; % Scale for gradient norm
            scale_s     = as;
                       % Active cell
            count=0; % Initiate iteration count 
                  
            while phi_d(end) > target*1.05
                
                    count=count+1;

                    if count==1   %First iteration

                        WxtWx = ( Vx * Wx)' * ( Vx * Wx );
                        WstWs = ( V * Ws)' * ( V * Ws);
                        
                        Rs = speye(mcell);
                        Rx = speye(size(Wx,1));
                    
                        phim_start =  ax * WxtWx + as * WstWs;
                        
                        phim = phim_start; 
                        
                        
                        
                        % Initial beta trace(G'G) / trace(phim)
                        beta = sum(sum(G.^2,1)) / sum(diag(phim,0).^2) * 1e+2;                
                        lambda = beta(count);
                        
                        phi(count) = objfunc(invmod,phim,beta(count));
                        
                    else

                        Rs = spdiags( 1./ ( abs(invmod) .^( 2-pvec(pp) ) + epsilon ) .^0.5 ,0,mcell,mcell);
                        Rx = spdiags( 1./ ( abs(Wx * invmod) .^( 2-qvec(qq) ) + epsilon ) .^0.5 ,0,mcell-1,mcell-1);

                        WxtRxWx = ( Vx * Rx * Wx)' * ( Vx * Rx * Wx);
                        WstRsWs = ( V * Rs * Ws)' * ( V * Rs * Ws);

                        scale_s = 0.5 * as *...
                            ( invmod' * phim_start * invmod ) /...
                            ( as * invmod' * WstRsWs * invmod );
                        
                        scale_x = 0.5 * ax *...
                            ( invmod' * phim_start * invmod ) /...
                            ( ax * invmod' * WxtRxWx * invmod );
                        
                        scale_s = lvec(ll) * scale_s;
                        scale_x = (2.0 - lvec(ll)) * scale_x;
                        
                        phim =  scale_s * WstRsWs +...
                                scale_x * WxtRxWx;
                        
                        scale(count) = ( invmod'* phim_start * invmod ) / (invmod'*phim*invmod ) ;
       
                        lambda(count) = beta(count) * scale(count);
                        
                    end

                
 
%                 eigval = eig(G'*G + lambda(count) * (phim) ); 
%                 condnum = max(abs(eigval)) / min(abs(eigval));
%                 
%                 figure(4); plot(eigval,'*');
%                 fprintf('Condition number: %f\n',condnum);
                fprintf('\n# # # # # #\n');
                fprintf('BETA ITER: \t %i  \nbeta: \t %8.5e \n',count,beta(count));

                ncg = 0;
                npcg = 1;
                rddm = 1;
                ddm = [1 1];

            while npcg < 10 && rddm > 1e-3

                A= [ G  ;...
                    sqrt( lambda(count) * scale_s ) * ( Wr * V * Rs * Ws ) ;...
                    sqrt( lambda(count) * scale_x ) * ( Wrx * Vx * Rx * Wx  )] ;

                b = [- (G *invmod - d) ; ...
                    -sqrt( lambda(count) * scale_s ) * ( Wr * V * Rs * Ws * invmod) ;...
                    -sqrt( lambda(count) * scale_x ) * ( Wrx * Vx * Rx * Wx  * invmod )] ;
% 
%                 A= [ G  ;...
%                     sqrt( lambda(count) * scale_s ) * (  V * Rs * Ws ) ;...
%                     sqrt( lambda(count) * scale_x ) * (  Vx * Rx * Wx  )] ;
% 
%                 b = [- (G *invmod - d) ; ...
%                     -sqrt( lambda(count) * scale_s ) * (  V * Rs * Ws * invmod) ;...
%                     -sqrt( lambda(count) * scale_x ) * (  Vx * Rx * Wx  * invmod )] ;

                %% Gauss-Newton step
                dm = zeros(mcell,1);
                [dm,r,iter] = CGLSQ( dm, A , b );

                %% Step length, line search                
                ncg = ncg+iter; % Record the number of CG iterations

                gamma = 2;

                % Initialise phi^k
                phi_temp = 0;   
                while (phi_temp > phi(end) || gamma == 2) && gamma>1e-4

                    gamma = 0.5 * gamma;

                    gdm = gamma * dm;

                    ddm(2) = norm(dm);

                    m_temp = invmod + gdm;

                    phi_temp = objfunc(m_temp,phim,lambda(count));

                end

                
                
                
                if npcg == 1

                    rddm = 1;
                    ddm(1) = ddm(2);

                else

                    rddm = ddm(2)/ddm(1);

                end
% 
%                 phi_d(count) = sum((G*(m_temp)-d).^2);
%                 
%                 if phi_d(count) < target *0.95
%                     
%                     beta(count) = beta(count) * 1.1;
%                     lambda(count) = beta(count) * scale(count);
%                     ncg = 0;
%                     npcg = 1;
%                     rddm = 1;
%                     ddm = [1 1];
%                     continue
%                     
%                 end


                % Update model
                invmod = invmod + dm;

                fprintf('GN iter %i |g| rel:\t\t %8.5e\n',npcg,rddm);
                npcg = npcg + 1;

            end                

%% Save iteration details
            phi_d(count) = sum((G*(invmod)-d).^2);
            phi_m(count) = (invmod)'*(phim)*(invmod);
            phi(count) = objfunc(invmod,phim,lambda(count));          
            

            if phi_d(count) < target * 2
                
                beta(count+1) = 0.5 * beta(count);
                
            else
                
                beta(count+1) = 0.5 * beta(count);
                
            end
           
%             set(figure(3+counter), 'Position', [50 200 1000 500])
%             plot(1:mcell,model); hold on
%             plot(1:mcell, invmod,'r','LineWidth',3); hold off
%             axis([0 mcell -.1 0.6]);
            fprintf('---------->')
            fprintf(' misfit:\t %8.5e ',phi_d(count))
            fprintf('<----------\n')

            
            end
            
%             figure(3); 
%             plot((G*(invmod)).*wd,'go');

            % Load other results and plot
%             load('m_IWR');
%             load('m_phim');
%             figure(3+counter); hold on
%             plot(m_phim,'--','LineWidth',3);
%             xlabel('\bfZ');
%             ylabel('\bfm(Z)');
%             legend('True','\bf Wr || m ||_0 + Wr || \Delta m ||_0','\bf || Wr m ||_0 + || \Delta (Wr m) ||_0','Location','NorthWest')

            % Remove depth weighting
%             invmod = IWr * invmod;
            
%             title(['\bfInversion: ' head])
%             figure; plot(1:mcell,model,1:mcell,invmod,'r','LineWidth',2);axis([0 mcell -.1 1])
%             figure; plot(eigval,'*');
%             figure; plotyy(1:count,phi_d,1:count,phi_m,@semilogy);
            RMS(pp,qq,ll)=(sum((invmod-model).^2)/length(data))^(0.5);
            linf(pp,qq,ll) = norm(model-invmod,'inf');
            l2(pp,qq,ll) = norm(model-invmod,2);
            l1(pp,qq,ll) = norm(model-invmod,1);
            phid(pp,qq,ll) = phi_d(count);
            models_out(pp,qq,ll,:) = invmod(:);
            
%             figure; plotyy(1:count,phi_d,1:count,phi_m);

            fprintf('End of lp inversion. Number of iterations: %i\n',count)
            fprintf('Final model misfit: %8.3e\n\n',norm(model-invmod,2))
            
            counter = counter+1;
            
        end
    end
end

%% Plot result 

% figure;slice(pvec,qvec,lvec,phid(:,:,:),[0.5 1 1.5 2.0],[0.5 1 1.5 2.0],0.5);%caxis([50 200])
% figure;slice(pvec,qvec,lvec,l1(:,:,:,nn),[0.5 1 1.5 2.0],[0.5 1 1.5 2.0],0.5);%caxis([0.05 4])
% figure;slice(pvec,qvec,lvec,l2(:,:,:),[0.5 1 1.5 2.0],[0.5 1 1.5 2.0],0.5);caxis([0.05 4]);title(['l2 residual for model' num2str(nn)])
% L2=l2(:,:,:,nn);L2=L2(:);save(['modelD' num2str(nn)],'-ascii','L2')