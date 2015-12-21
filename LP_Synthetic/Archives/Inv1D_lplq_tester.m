% Inv1D_lplq
% Program creates a 1D toy problem and generates data using kernels of the
% form k(z) = exp^(-ajx) cos(b j 2pi x)
% 
% Data can then be inverted using various lp-norms applied on the model and
% its gradient. Depth weighting is applied to the model objective function,
% but outside the norm.
% 
% Written by: D. Fournier
% Last Update: November 11, 2014

close all
clear all

addpath data\
addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\Lp_norm_codes\Dual_norm_1D\functions
addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\Lp_norm_codes\Dual_norm_1D\sliceomatic\
% Set up 1D problem

%% INPUT VARIABLES

% Model space
nx = 200;

% Kernel
nk = 40; %number of frequencies

decay = -0.10;

basis = 0.20;

% Data noise
amp_pct = 0.02;
floor_pct = 0.02;

% Lp-norm parameters
pvec = 0%0:0.5:2;
qvec = 0%0:0.5:2;
lvec = [1.0]%0.1:0.4:1.9;

% Add an extra test to find best starting beta for lplq
% chi_vec = [2 3 4 5 7.5 10 20 30];
% multip = 100;

%% SCRIPT STARTS HERE
%% % Generate kernel functions and depth weighting
x = (0 : 1/(nx) : 1) * nx/10;
mcell=nx;
dx=ones(1,mcell) * abs((min(x)-max(x))/mcell);


set(figure, 'Position', [50 200 1000 500])
G = zeros(nk,mcell);
wr = zeros(1,mcell);
for ii = 1 : nk 
    
    b = basis * 2*pi* (ii) / (nx/10);
    a = decay * (ii) / (nx/10);
    
    G(ii,1:mcell) = exp( a * x(2:end) ) .* (a/(a^2+b^2)*cos (b * x(2:end)) + b/(a^2+b^2) *sin (b * x(2:end))) -...
        exp( a * x(1:end-1) ) .* (a/(a^2+b^2)*cos (b * x(1:end-1)) + b/(a^2+b^2) *sin (b * x(1:end-1)));

    G(ii,1:mcell) = G(ii,1:mcell) / max(G(ii,1:mcell));
    
    wr = wr + ((exp( a*x(2:end))-exp( a * x(1:end-1) ))/a) ;
    plot(G(ii,:));hold on

end

wr = (wr / max(wr)).^0.5;

IWr = spdiags(1./wr',0,mcell,mcell);
Wr = spdiags(wr',0,mcell,mcell);
Wrx = spdiags(wr(1:end-1)',0,mcell-1,mcell-1);

figure(1); plot(wr','r--','LineWidth',3);
axis square;grid on
temp = legend('\bfg(z) = e^{-jz} cos(j 2\pi z)','\bfDepth Weighting (Wr)');
ylabel('\bfAmplitude');
xlabel('\bfZ');
a=get(temp,'children');
set(a(2),'LineWidth',3,'Color',[1 0 0],'LineStyle','--'); 
%% Generate model

x = (1/nx : 1/(nx) : 1) * nx/10;
% Create mdoel: Square and tanh functions
model = 0.4* exp(-abs(((x-15)*10 / (nx/10)).^2)*2); 

model(30:70)=0.20;

model=model(:);

set(figure, 'Position', [50 200 1000 500]);
plot(x,model,'LineWidth',2);axis([x(1) x(end) -0.1 0.5]);hold on
axis square
grid on


%% FORWARD MODELING
% Generate data and corrupt with noise

data= G * model;

% Corrupt data with Gaussian noise
rand_noise = rand(nk,1);
save ('noise_1D','rand_noise')
load('noise_1D');

amp_noise = amp_pct * rand_noise .* abs(data);
floor_noise = floor_pct * rand_noise * max(abs(data));

noise = amp_noise + floor_noise;

d = data + noise;

wd =  abs(d)*0.05 + 0.05 * std(data);

figure; 
plot(data);
hold on
errorbar(d,wd,'r*')
legend('Data','Data+noise')
grid on

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
% RMS = zeros(length(pvec),length(qvec),length(lvec));
% linf = zeros(length(pvec),length(qvec),length(lvec));
% l2 = zeros(length(pvec),length(qvec),length(lvec));
% l1 = zeros(length(pvec),length(qvec),length(lvec));
% misfit = zeros(length(pvec),length(qvec),length(lvec));
% models_out = zeros(length(pvec),length(qvec),length(lvec),mcell);

counter = 1;  

objfunc = @(m,phi,b) sum( ( G * m - d ).^2 ) + ( m' * b * phi * m );

% Iterate over all the combination of norms as specified by the vectors
% pvec, qvec and lvec.

chi_start = zeros(length(pvec),length(qvec),length(lvec),length(chi_vec));
RMS = zeros(length(pvec),length(qvec),length(lvec),length(chi_vec));
linf = zeros(length(pvec),length(qvec),length(lvec),length(chi_vec));
l2 = zeros(length(pvec),length(qvec),length(lvec),length(chi_vec));
l1 = zeros(length(pvec),length(qvec),length(lvec),length(chi_vec));
misfit = zeros(length(pvec),length(qvec),length(lvec),length(chi_vec));
models_out = zeros(length(pvec),length(qvec),length(lvec),length(chi_vec),1,mcell);


for ll= 1:length(lvec)

    for pp = 1:length(pvec)
        
        for qq = 1:length(qvec)
            
            for bb = 1 : length(chi_vec)
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
            count_lp = 0;      
            while phi_d(end) > target*1.05
                
                    count=count+1;

                    if phi_d(end)/nk > (target/nk * 1.05)    %First iteration

                        WxtWx = ( Vx * Wx)' * ( Vx * Wx );
                        WstWs = ( V * Ws)' * ( V * Ws);
                        
                        Rs = speye(mcell);
                        Rx = speye(size(Wx,1));
                    
                        phim_start =  ax * WxtWx + as * WstWs;
                        
                        phim = phim_start; 
                        
                        
                        
                        % Initial beta trace(G'G) / trace(phim)
                        if count == 1
                            
                            beta = sum(sum(G.^2,1)) / sum(diag(phim,0).^2) *1e+5;  
                        
                        end
                        
                        lambda(count) = beta(count);
                        
                        phi(count) = objfunc(invmod,phim,beta(count));
                        
                    else

                        fprintf('\nLP-LQ ITER');
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
                        
                        count_lp = count_lp + 1;
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
                
            while npcg < 10 && rddm > 1e-3

                A= [ G  ;...
                    sqrt( lambda(count) * scale_s ) * ( Wr * V * Rs * Ws ) ;...
                    sqrt( lambda(count) * scale_x ) * ( Wrx * Vx * Rx * Wx  )] ;

                b = [- (G *invmod - d) ; ...
                    -sqrt( lambda(count) * scale_s ) * ( Wr * V * Rs * Ws * invmod) ;...
                    -sqrt( lambda(count) * scale_x ) * ( Wrx * Vx * Rx * Wx  * invmod )] ;

                %% Gauss-Newton step
                dm = zeros(mcell,1);
                [dm,r,iter] = CGLSQ( dm, A , b );

                %% Step length, line search                
                ncg = ncg+iter; % Record the number of CG iterations

                gamma = 2;

                % Reduce step length in order to reduce phid
                phi_temp = 0;   
                while (phi_temp > phi(end) || gamma == 2) && gamma>1e-4

                    gamma = 0.5 * gamma;

                    gdm = gamma * dm;

                    ddm(2) = norm(dm);

                    m_temp = invmod + gdm;

                    phi_temp = objfunc(m_temp,phim,lambda(count));

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
               
            

            phi_d(count) = sum((G*(invmod)-d).^2);
                
            % Check to see if overshooted the target misfit,
            % If yes, then compare with the previous iteration and keep
            % closest                
            if phi_d(count) < target *0.95 || (phi_d(end)/nk < (target/nk * chi_vec(bb) * 0.95) && count_lp==0)
                
                phi_d(count) = sum((G*(invmod)-d).^2);                       
                
                fprintf('---------->')
                fprintf(' misfit:\t %8.5e ',phi_d(count))
                fprintf('<----------\n')
                fprintf('\n# BACKTRACKING STEP #\n');
                
                beta(count) = beta(count)*1.1;
                invmod = m_in;
                phi_d(count) = phi_d(count-1);
                
                count = count-1;
                continue
            
            % Else reduce beta and continue inversion
            else
                
                phi_d(count) = sum((G*(invmod)-d).^2);
                
                if count_lp == 0
                    
                    chi_start(pp,qq,ll,bb) = phi_d(count)/nk;
                    
                end
                
                
                phi_m(count) = (invmod)'*(phim)*(invmod);
                phi(count) = objfunc(invmod,phim,lambda(count)); 
                
                fprintf('---------->')
                fprintf(' misfit:\t %8.5e ',phi_d(count))
                fprintf('<----------\n')
                
                
                beta(count+1) = 0.5 * beta(count);           
                
            end
           
            % OPTIONAL: Plot model
            set(figure(4), 'Position', [50 200 1000 500])
            plot(1:mcell,model); hold on
            plot(1:mcell, invmod,'r','LineWidth',3); hold off
            axis([0 mcell -.1 0.6]);
            axis square

            
            end
            
%             legend('\bfTrue model','\bf Wr || m ||_0 + Wr || \nabla m ||_0');
%             xlabel('\bfZ');
%             ylabel('\bfm(Z)');
            
            
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
            RMS(pp,qq,ll,bb)=(sum((invmod-model).^2)/length(data))^(0.5);
            linf(pp,qq,ll,bb) = norm(model-invmod,'inf');
            l2(pp,qq,ll,bb) = norm(model-invmod,2);
            l1(pp,qq,ll,bb) = norm(model-invmod,1);
            misfit(pp,qq,ll,bb) = phi_d(count);
            models_out(pp,qq,ll,bb,:) = invmod(:);
            
%             figure; plotyy(1:count,phi_d,1:count,phi_m);

            fprintf('End of lp inversion. Number of iterations: %i\n',count)
            fprintf('Final model misfit: %8.3e\n\n',norm(model-invmod,2))
            
            counter = counter+1;
            
            end
        end
    end
end



%% Plot result 
% sliceomatic(l1,pvec,qvec,lvec);
set(figure, 'Position', [50 200 1000 500]);
ltext = [];
count = 1;
for pp = 1:length(pvec)
        
        for qq = 1:length(qvec)
            
            for ll = 1 : length(lvec)
                xx = chi_start(pp,qq,ll,:); xx = xx(:);
                yy = l1(pp,qq,ll,:); yy = yy(:);

               
                % Find best fitting quadratic
                cc = [xx.^0 xx.^1 xx.^2 xx.^3]\yy;
                
%                 plot(xx,yy,'Color',[ll/length(lvec) pp/length(pvec) qq/length(qvec)].^1.5,'LineWidth',2); hold on
                plot(xx,yy,'Color',[pp/length(pvec) rand(1) qq/length(qvec)],'LineWidth',2); hold on
                
                ltext{count} = ['pvec ' num2str(pvec(pp)) ' qvec ' num2str(qvec(qq)) ' lvec ' num2str(lvec(ll))];
                
                count = count + 1;
            end
            
        end
        
end
legend(ltext)

% Fit polynomial
set(figure, 'Position', [50 200 1000 500]);
for pp = 1:length(pvec)
        
        for qq = 1:length(qvec)
            
            for ll = 1 : length(lvec)
                xx = chi_start(pp,qq,ll,:); xx = xx(:);
                yy = l1(pp,qq,ll,:); yy = yy(:);

                % Find best fitting quadratic
                cc = [xx.^0 xx.^1 xx.^2 xx.^3]\yy;
                
                func= cc(4)*xx.^3 + cc(3)*xx.^2 + cc(2).*xx + cc(1);
                plot(xx,func,'r:','LineWidth',2); hold on
                
                index = find(min(func)==func);

                plot(xx(index),func(index),'bo')
            end
            
        end
        
end
xlabel('\bf Starting chi factor');ylabel('\bf|m_{true} - m_{rec}|')