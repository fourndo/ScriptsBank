close all
clear all

addpath data\
addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\Lp_norm_codes\Dual_norm_1D\functions
addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\Lp_norm_codes\Dual_norm_1D\sliceomatic\
addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB
% Set up 1D problem

%% INPUT VARIABLES

% Model space
nx = 200;

% Kernel
nk = 40; %number of frequencies

decay = -0.20;

basis = 0.20;

% Data noise
amp_pct = 0.02;
floor_pct = 0.02;

% Lp-norm parameters
pvec = 0%0:1:2;
qvec = 0%0:1:2;
lvec = 1.0%0.1:0.4:1.9;

%% SCRIPT STARTS HERE
%% % Generate kernel functions and depth weighting
x = (0 : 1/(nx) : 1) ;
mcell=nx;
dx=ones(1,mcell) * abs((min(x)-max(x))/mcell);


set(figure, 'Position', [50 200 750 750]);
G = zeros(nk,mcell);
wr = zeros(1,mcell);
for ii = 1 : nk 
    
    b = basis * 2*pi* (ii) ;
    a = decay * (ii) ;
    
    G(ii,1:mcell) = exp( a * x(2:end) ) .* (a/(a^2+b^2)*cos (b * x(2:end)) + b/(a^2+b^2) *sin (b * x(2:end))) -...
        exp( a * x(1:end-1) ) .* (a/(a^2+b^2)*cos (b * x(1:end-1)) + b/(a^2+b^2) *sin (b * x(1:end-1)));

    G(ii,1:mcell) = G(ii,1:mcell) / max(G(ii,1:mcell));
    
    wr = wr + ((exp( a*x(2:end))-exp( a*x(1:end-1)))/a).^2 ;
    plot(x(2:end),G(ii,:));hold on

end


wr = wr.^0.5;
wr = (wr / max(wr)).^0.5;

IWr = spdiags(1./wr',0,mcell,mcell);
Wr = spdiags(wr',0,mcell,mcell);
Wrx = spdiags(wr(1:end-1)',0,mcell-1,mcell-1);

% figure(1); plot(wr','r--','LineWidth',3);
axis square;grid on
% temp = legend('\bfg(z) = e^{-jz} cos(j 2\pi z)','\bfDepth Weighting (Wr)');
ylabel('\bfAmplitude');
xlabel('\bfZ');
% a=get(temp,'children');
% set(a(2),'LineWidth',3,'Color',[1 0 0],'LineStyle','--'); 
%% Generate model

x = (1/nx : 1/(nx) : 1) ;
% Create mdoel: Square and tanh functions
model = 0.4* exp(-abs(((x-0.75) ).^2)*nx); 

model(30:70)=0.20;

model=model(:);


set(figure, 'Position', [50 200 750 750]);
plot(x,model,'LineWidth',2);
axis([0 1 -.025 0.425]);
xlabel('\bf Z');
ylabel('\bf Model m(Z)')
axis square
grid on


%% FORWARD MODELING
% Generate data and corrupt with noise

data= G * model;

% Corrupt data with Gaussian noise
% rand_noise = rand(nk,1);
% save ('noise','rand_noise')
load('noise');

amp_noise = amp_pct * rand_noise .* abs(data);
floor_noise = floor_pct * rand_noise * max(abs(data));

noise = amp_noise + floor_noise;

d = data + noise;

wd =  abs(d)*0.05 + 0.05 * std(data);

set(figure, 'Position', [50 200 750 750])
plot(data,'LineWidth',2);
hold on
errorbar(d,wd,'r.')
legend('Raw data','Noisy data & Uncertainties')
grid on
xlabel('\bf Frequency ( j )');
ylabel('\bf Data');
% Create uncertainty weighting matrix
Wd = spdiags(1./wd,0,nk,nk);

% Normalize d and G by standard deviation:
G = Wd * G ;

target = nk;

% Weight data with uncertainty
d = Wd * d;

%% INVERSION
% Create inversion parameters
[Wx,Vx] = getWx1D(nx,dx);
[Ws] = speye(mcell);
V = spdiags(sqrt(dx'),0,nx,nx);


% Global constant
as = 1.0 / min(dx) ^2;  %Smallnest term
ax = 1.0;               %Smoothness term                     

epsilon = 1e-8;     %Small term in compact function

% Store all the final models
RMS = zeros(length(pvec),length(qvec),length(lvec));
linf = zeros(length(pvec),length(qvec),length(lvec));
l2 = zeros(length(pvec),length(qvec),length(lvec));
l1 = zeros(length(pvec),length(qvec),length(lvec));
misfit = zeros(length(pvec),length(qvec),length(lvec));
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
            
            invmod      = ones(mcell,1)*1e-2;       % Initial model       
                        
            phi_d       = sum((G*invmod - d).^2);   % Initial misfit
            scale_x     = ax; % Scale for gradient norm
            scale_s     = as;
                       % Active cell
            count=0; % Initiate iteration count 
            switcher = 0;
            lp_count= 0;
            cg_iter = 0;
            while switcher ~= 3
                
                    count=count+1;

                    delta(count) = 1;
                    if switcher == 0   %First iteration
                        

                        WxtWx = ( Vx * Wx)' * ( Vx * Wx );
                        WstWs = ( V * Ws)' * ( V * Ws);
                        
                        Rs = speye(mcell);
                        Rx = speye(size(Wx,1));
                    
                        phim_start =  ax * WxtWx + as * WstWs;
                        
                        phim = phim_start; 
                        
                        
                        
                        % Initial beta trace(G'G) / trace(phim)
                        if count==1
                            % Initial beta trace(G'G) / trace(phim) 
                            beta = sum(sum(G.^2,1)) / sum(diag(phim,0).^2) * 1e+5;
                            
                        end
                                        
                        lambda(count) = beta(count);
                        
                        phi(count) = objfunc(invmod,phim,beta(count));
                        
                    else
                        
                        fprintf('\n# # LP-LQ ITER# #');
                        lp_count = lp_count+1;
                        
                        if cg_iter(end) < mcell && switcher == 1 && delta(count-1) > 1e-9
                            
                            delta(count) = 10^-(lp_count/2);
                            
                        else
                            
                            switcher = 2;
                            delta(count) = delta(count-1);
                            
                        end

                        Rs = spdiags( 1./ ( abs(invmod) .^( 2-pvec(pp) ) + delta(count) ) .^0.5 ,0,mcell,mcell);
                        Rx = spdiags( 1./ ( abs(Wx * invmod) .^( 2-qvec(qq) ) + delta(count) ) .^0.5 ,0,mcell-1,mcell-1);

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

                ncg = 0;        % Count total CG steps
                npcg = 1;       % Count number of projected gradient steps
                rddm = 1;       % Mesure the relative change in rel|dm|
                ddm = [1 1];    % Mesure of change in model update |dm|
                m_in = invmod;  % Save initial model for backtracking steps
                rdm = [];
                
            while npcg < 10 && rddm > 1e-3

                A= [ G  ;...
                    sqrt( lambda(count) * scale_s ) * ( Wr * V * Rs * Ws ) ;...
                    sqrt( lambda(count) * scale_x ) * ( Wrx * Vx * Rx * Wx  )] ;

                b = [- (G *invmod - d) ; ...
                    -sqrt( lambda(count) * scale_s ) * ( Wr * V * Rs * Ws * invmod) ;...
                    -sqrt( lambda(count) * scale_x ) * ( Wrx * Vx * Rx * Wx  * invmod )] ;

                %% Gauss-Newton step
                dm = zeros(mcell,1);
                [dm,r,cg_iter(count)] = CGLSQ( dm, A , b );

                %% Step length, line search                
                ncg = ncg+cg_iter(count); % Record the number of CG iterations

                gamma = 2;

                % Reduce step length in order to reduce phid
                phi_temp = [1 2];   
                while (phi_temp(1) > phi(end) || gamma == 2) && (phi_temp(1)/phi_temp(2))<1

                    phi_temp(2) = phi_temp(1);
                    
                    gamma = 0.5 * gamma;

                    gdm = gamma * dm;

                    ddm(2) = norm(gdm);

                    m_temp = invmod + gdm;

                    phi_temp(1) = objfunc(m_temp,phim,lambda(count));

                end
              
                % Mesure relative change in dm
                if npcg == 1

                    rddm = 1;
                    ddm(1) = ddm(2);

                else

                    rddm = ddm(2)/ddm(1);

                end


                % Update model
                invmod = invmod + dm;
                
                fprintf('GN iter %i |g| rel:\t\t %8.5e\n',npcg,rddm);
                npcg = npcg + 1;

            end                

%% Save iteration details
               
            
            rdm(count) = norm(m_in - invmod) / norm(m_in);
            
            phi_d(count) = sum((G*(invmod)-d).^2);
            phi_m(count) = (invmod)'*(phim)*(invmod);
            phi(count) = objfunc(invmod,phim,lambda(count)); 
              
            if rdm(count) < 5e-3 && switcher >= 1 && phi_d(count) > target *0.98 && phi_d(count) < target *1.02

                switcher = 3;
                fprintf('---------->')
                fprintf(' misfit:\t %8.5e \n',phi_d(count))
                fprintf('Relative dm:\t %8.5e ', rdm(count));
                fprintf('<----------\n')
                
                continue

            end    
            
            % Check to see if overshooted the target misfit,
            % If yes, then compare with the previous iteration and keep
            % closest                
            if phi_d(count) < target *0.98
                
                phi_d(count) = sum((G*(invmod)-d).^2);                       
                
                fprintf('---------->')
                fprintf(' misfit:\t %8.5e ',phi_d(count))
                fprintf('<----------\n')
                fprintf('\n# BACKTRACKING STEP #\n');
                
                if switcher == 0
                    switcher = 1;
                end
                
                beta(count+1) = 1.02 * beta(count);
                
                            
            % Else reduce beta and continue inversion
            elseif phi_d(count) > target * 1.02
                
                
               if switcher == 0
                    
                    beta(count+1) = 0.75 * beta(count); 
                    
                else
                    
                    beta(count+1) = 0.98 * beta(count);
                    
                end  
                
                fprintf('---------->')
                fprintf(' misfit:\t %8.5e ',phi_d(count))
                fprintf('<----------\n')    
                fprintf('\n# REDUCE STEP #\n');
                
            else
                
                
                fprintf('---------->')
                fprintf(' misfit:\t %8.5e \n',phi_d(count))
                fprintf('Relative dm:\t %8.5e ', rdm(count));
                fprintf('<----------\n')
                fprintf('\n# STATIONNARY STEP #\n');
                
                if switcher == 0
                    switcher = 1;
                end
                
                beta(count+1) = beta(count);

                 
            end
           
            % OPTIONAL: Plot model
            set(figure(counter+3), 'Position', [50 200 750 750]);
            plot(x,model,'b','LineWidth',2); hold on
            plot(x, invmod,'r','LineWidth',3); hold off
            axis([0 1 -.025 0.425]);
            axis square

            
            end
            
            load('m_wr')
            hold on 
            plot(x,m_wr,'k--','LineWidth',3)
            L = legend('True model',...
                '$\phi_m \; = \alpha \| \mathbf{ W_r^{-1} \; R\; \;f}(m) \|_2^2$',...
                '$\phi_m \; = \alpha \| \mathbf{ \tilde{R}\; \; f} (\tilde{m}) \|_2^2$','Location','NorthWest');
            set(L, 'interpreter', 'latex','FontSize',15)
            xlabel('\bf Z');
            ylabel('\bf Model m(Z)')
            axis square
            tlabel=(['$\mathbf{\phi_m\;=\;\| m \|_' num2str(pvec(pp)) '\;+\; \| \nabla m \|_' num2str(qvec(qq)) '}$']);
            title(tlabel, 'interpreter', 'latex','Fonts',18)
            grid on
            % Plot final phid
%             figure(3); 
%             plot((G*(invmod)).*wd,'go');

            % Load other results and plot
%             load('m_IWR');
%             load('m_phim');
%             figure(3+counter); hold on
%             plot(m_phim,'--','LineWidth',3);
%             xlabel('\bfZ');
%             ylabel('\bfm(Z)');
%             legend('\bfTrue model','\bf Wr || m ||_0 + Wr || \nabla m ||_0','\bf || Wr m ||_0 + || \nabla (Wr m) ||_0','Location','NorthWest')

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
            misfit(pp,qq,ll) = phi_d(count);
            models_out(pp,qq,ll,:) = invmod(:);
            
%             figure; plotyy(1:count,phi_d,1:count,phi_m);

            fprintf('End of lp inversion. Number of iterations: %i\n',count)
            fprintf('Final model misfit: %8.3e\n\n',norm(model-invmod,2))
            
            counter = counter+1;
            
        end
    end
end

%% Plot result 
% sliceomatic(l1,pvec,qvec,lvec);
% load('l1');
[I,J]=find(l2==min(min(l2)));
set(figure, 'Position', [50 200 750 750]);
imagesc(pvec,qvec,l2);colorbar;colormap winter;hold on
contour(pvec,qvec,l2,'ShowText','on','LineColor',[0 0 0])
set(gca,'YDir','Normal')
ylabel('$\mathbf{\| m \|_p}$', 'interpreter', 'latex','Fonts',15)
xlabel('$\mathbf{\| \nabla m \|_q}$', 'interpreter', 'latex','Fonts',15)
plot(qvec(J),pvec(I),'ro','LineWidth',2);hold on
axis square
title('\textbf{Model error} ($\mathbf{\| m - m^{*} \|_2}$)', 'interpreter', 'latex','Fonts',15)

% load('misfit');
[I,J]=find(l2==min(min(l2)));
set(figure, 'Position', [50 200 750 750]);
imagesc(pvec,qvec,misfit);colorbar;colormap winter;hold on
% contour(pvec,qvec,misfit,'ShowText','on','LineColor',[0 0 0])
set(gca,'YDir','Normal')
ylabel('$\mathbf{\| m \|_p}$', 'interpreter', 'latex','Fonts',15)
xlabel('$\mathbf{\| \nabla m \|_q}$', 'interpreter', 'latex','Fonts',15)
plot(qvec(J),pvec(I),'ro','LineWidth',2);hold on
axis square
title('\textbf{Final data misifit} $\phi_d$', 'interpreter', 'latex','Fonts',15)
caxis([38 42])