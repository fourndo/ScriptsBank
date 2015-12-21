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
addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB
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
pvec = [0 2]%0:0.25:2;
qvec = [0 2]%0:0.25:2;
lvec = ones(1,1)*1.0%[0.25:0.25:1.75]%0.1:0.4:1.9;

% Add an extra test to find best starting beta for lplq
chi_vec = 1%[2 3 4 5 7.5 10 20 30];
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

x = (1/nx : 1/(nx) : 1) ;

cntr = round(nx/4);

% Create mdoel: Square and tanh functions
model = 0.4* exp(-(((x-3*cntr/nx) ).^2 * nx)); 

model(cntr-round(cntr/4):cntr+round(cntr/4))=0.25;

model=model(:);

set(figure, 'Position', [50 200 1000 500]);
plot(x,model,'LineWidth',2);axis([x(1) x(end) -0.1 0.5]);hold on
axis square
grid on


%% FORWARD MODELING
% Generate data and corrupt with noise

data= G * model;

% Corrupt data with Gaussian noise
% rand_noise = rand(nk,1);
% save ('noise_1D','rand_noise')
load('noise_1D');

amp_noise = amp_pct .* abs(data);
floor_noise = min(abs(data));

wd = amp_noise + floor_noise;
% noise = amp_noise + floor_noise;

d = data + rand_noise.*wd;

% wd =  abs(d)*0.05 + 0.05 * std(data);

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

                       % Active cell
            count=0; % Initiate iteration count 
            switcher = 0;
            ncg(counter) = 0;
            lp_count = 0;
            cg_iter = 0;
            delta = 1;
            lambda = [];
            
        while switcher ~= 3

                    count=count+1;

                    if switcher == 0   %First iteration

                        delta(count) = 1;
                        
                        WxtWx = ( Vx * Wx)' * ( Vx * Wx );
                        WstWs = ( V * Ws)' * ( V * Ws);
                        
                        Rs = speye(mcell);
                        Rx = speye(size(Wx,1));
                    
                        phim_start =  ax * WxtWx + as * WstWs;
                        
                        MOF = phim_start; 
                        
                        aVRWs = sqrt(as) * Wr * V * Ws;
                        aVRWx = sqrt(ax) * Wrx * Vx * Wx;
                        
                        % Initial beta trace(G'G) / trace(phim)
                        if count == 1
                            
                            beta = sum(sum(G.^2,1)) / sum(diag(MOF,0).^2) *1e+4;  
                        
                        end
                        
                        lambda(count) = beta(count);
                        
                        phi(count) = objfunc(invmod,MOF,beta(count));
                        
                    else
                        
                        fprintf('# # LP-LQ ITER# #');
                        lp_count = lp_count+1;
                        
                        if cg_iter(end) < 0.5*mcell && switcher == 1
                            
                            delta(count) = 10^-(lp_count/2);
                            
                        else
                            
                            switcher = 2;
                            delta(count) = delta(count-1);
                            
                        end

                        Rs = spdiags( 1./ ( abs(invmod) .^( 2-pvec(pp) ) + delta(count) ) .^0.5 ,0,mcell,mcell);
                        Rx = spdiags( 1./ ( abs(Wx * invmod) .^( 2-qvec(qq) ) + delta(count) ) .^0.5 ,0,mcell-1,mcell-1);
                       
                        
                        WxtRxWx = ( Vx * Rx * Wx)' * ( Vx * Rx * Wx);
                        WstRsWs = ( V * Rs * Ws)' * ( V * Rs * Ws);

                        scale_s = 0.5 * ...
                            (  phi_m(end) ) /...
                            ( as * invmod' * WstRsWs * invmod );
                        
                        scale_x = 0.5 * ...
                            (  phi_m(end) ) /...
                            ( ax * invmod' * WxtRxWx * invmod );
                        
                        scales(1) = lvec(ll) * scale_s * as;
                        scales(2) = (2.0 - lvec(ll)) * scale_x * ax;
                        
                        MOF =  scales(1) * WstRsWs +...
                                scales(2) * WxtRxWx;
                        
                        aVRWs = sqrt(scales(1)) * Wr * V * Rs * Ws;
                        aVRWx = sqrt(scales(2)) * Wrx * Vx * Rx * Wx;                        
                        
                        mu = ( invmod'* phim_start * invmod ) / (invmod'*MOF*invmod ) ;
       
                        lambda(count) = beta(count) * mu;
                        

                    end

                
                fprintf('\n# # # # # #\n');
                fprintf('BETA ITER: \t %i  \nbeta: \t %8.5e \n',count,beta(count));

                m_in = invmod;  % Save initial model for backtracking steps
                
%                 eigval_GtG = eig(G'*G); 
%                 eigval_MOF = eig( lambda(count) * (MOF) );
%                 
%                 eigval = eig( G'*G + lambda(count) * (MOF));
%                 
%                 condnum = max(abs(eigval)) / min(abs(eigval));
%                 
%                 figure(4); 
% %                 plot(eigval_GtG,'*'); hold on
% %                 plot(eigval_MOF,'ro'); hold on
%                 plot(eigval,'*'); hold off
%                 fprintf('Condition number: %f\n',condnum);

                [invmod,cg_iter] = GNsolver( G, invmod, d, phi(end), MOF, lambda(count) , aVRWs, aVRWx ,[],[]);             
                ncg(counter) = ncg(counter) + cg_iter;
                
                rdm(count) = norm(m_in - invmod) / norm(m_in);
                
               
                phi_d(count) = sum((G*(invmod)-d).^2);
                phi_m(count) = (invmod)'*(MOF)*(invmod);
                phi(count) = objfunc(invmod,MOF,beta(count));

                if (rdm(count) < 1e-3 && switcher >= 1 && phi_d(count) > target *0.98 && phi_d(count) < target *1.02) || count > 50 

                    switcher = 3;
                    fprintf('---------->')
                    fprintf(' misfit:\t %8.5e ',phi_d(count))
                    fprintf(' ** Relative dm:\t %8.5e ', rdm(count));
                    fprintf('<----------\n')

                    continue

                end    
                
            % Check to see if overshooted the target misfit,
            % If yes, then compare with the previous iteration and keep
            % closest                
            if phi_d(count) < target * 0.98
                
                phi_d(count) = sum((G*(invmod)-d).^2);                       
                
                fprintf('---------->');
                fprintf(' misfit:\t %8.5e ',phi_d(count));
                fprintf('Relative dm:\t %8.5e ', rdm(count));
                fprintf('<----------\n');
                fprintf('\n# NEXT ITER: INCREASING BETA #\n');
                
                if switcher == 0
                    switcher = 1;
                end
%                 beta(count+1) = beta(count) * target / phi_d(count);
                beta(count+1) = beta(count) * 1.1;
                
%                 invmod = m_in;
%                 phi_d(count) = phi_d(count-1);
                
%                 count = count-1;
%                 continue
            
            % Else reduce beta and continue inversion
            elseif phi_d(count) > target * 1.02
                
                fprintf('---------->')
                fprintf(' misfit:\t %8.5e ',phi_d(count))
                fprintf(' ** Relative dm:\t %8.5e ', rdm(count));
                fprintf('<----------\n')
                fprintf('\n# NEXT ITER - REDUCING BETA #\n');
                if switcher == 0
                    
                    beta(count+1) = 0.75 * beta(count); 
                    
                else
                    
%                     beta(count+1) = beta(count) * target / phi_d(count);
                    beta(count+1) = beta(count) * 0.95;
                    
                end
                
            else
                
                
                fprintf('---------->')
                fprintf(' misfit:\t %8.5e ',phi_d(count))
                fprintf(' ** Relative dm:\t %8.5e ', rdm(count));
                fprintf('<----------\n')
                fprintf('\n# NEXT ITER - STATIONNARY STEP #\n');
                
                if switcher == 0
                    switcher = 1;
                end
                
                beta(count+1) = beta(count);

                 
            end
           
            
            % OPTIONAL: Plot model
            set(figure(4+counter), 'Position', [50 200 1000 500])
            plot(1:mcell,model); hold on
            plot(1:mcell, invmod,'r','LineWidth',3); hold off
            axis([0 mcell -.1 0.6]);
            axis square

            
            
        end
            tlabel=(['$\mathbf{\phi_m\;=\;\| m \|_{' num2str(pvec(pp)) '}\;+\; \| \nabla m \|_{' num2str(qvec(qq)) '}}$']);
            title(tlabel, 'interpreter', 'latex','Fonts',18)
            
            RMS(pp,qq,ll)=(sum((invmod-model).^2)/length(data))^(0.5);
            linf(pp,qq,ll) = norm(model-invmod,'inf');
            l2(pp,qq,ll) = norm(model-invmod,2);
            l1(pp,qq,ll) = norm(model-invmod,1);
            misfit(pp,qq,ll) = phi_d(count);
            models_out(pp,qq,ll,:) = invmod(:);
            
%             figure; plotyy(1:count,phi_d,1:count,phi_m);

            ldml(counter) = norm(model-invmod,2);
            fprintf('Relative dm:\t %8.5e \n', rdm(count));
            fprintf('End of lp inversion. Number of iterations: %i\n',count)
            fprintf('Final model misfit: %8.3e\n\n',ldml(counter))
            
            counter = counter+1;
            

        end
    end
end



%% Plot result 
% sliceomatic(l1,pvec,qvec,lvec);
% set(figure, 'Position', [50 200 1000 500]);
% ltext = [];
% count = 1;
% for pp = 1:length(pvec)
%         
%         for qq = 1:length(qvec)
%             
%             for ll = 1 : length(lvec)
%                 xx = chi_start(pp,qq,ll,:); xx = xx(:);
%                 yy = l1(pp,qq,ll,:); yy = yy(:);
% 
%                
%                 % Find best fitting quadratic
%                 cc = [xx.^0 xx.^1 xx.^2 xx.^3]\yy;
%                 
% %                 plot(xx,yy,'Color',[ll/length(lvec) pp/length(pvec) qq/length(qvec)].^1.5,'LineWidth',2); hold on
%                 plot(xx,yy,'Color',[pp/length(pvec) rand(1) qq/length(qvec)],'LineWidth',2); hold on
%                 
%                 ltext{count} = ['pvec ' num2str(pvec(pp)) ' qvec ' num2str(qvec(qq)) ' lvec ' num2str(lvec(ll))];
%                 
%                 count = count + 1;
%             end
%             
%         end
%         
% end
% legend(ltext)
% 
% % Fit polynomial
% set(figure, 'Position', [50 200 1000 500]);
% for pp = 1:length(pvec)
%         
%         for qq = 1:length(qvec)
%             
%             for ll = 1 : length(lvec)
%                 xx = chi_start(pp,qq,ll,:); xx = xx(:);
%                 yy = l1(pp,qq,ll,:); yy = yy(:);
% 
%                 % Find best fitting quadratic
%                 cc = [xx.^0 xx.^1 xx.^2 xx.^3]\yy;
%                 
%                 func= cc(4)*xx.^3 + cc(3)*xx.^2 + cc(2).*xx + cc(1);
%                 plot(xx,func,'r:','LineWidth',2); hold on
%                 
%                 index = find(min(func)==func);
% 
%                 plot(xx(index),func(index),'bo')
%             end
%             
%         end
%         
% end
% xlabel('\bf Starting chi factor');ylabel('\bf|m_{true} - m_{rec}|')

figure;
[hAx,hLine1,hLine2] = plotyy(delta,ncg,delta,ldml,'semilogx','semilogx');
ylabel(hAx(1),'\bfNumber of CG steps');
ylabel(hAx(2),'\bf|m_{true} - m_{k}|_1');
xlabel('\bf Factor (\epsilon)')
set(hAx(2),'YColor','k')
set(hLine2,'LineWidth',2,'Color','k')
set(hLine1,'LineWidth',2)
grid on