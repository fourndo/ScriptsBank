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
% addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\Lp_norm_codes\Dual_norm_1D\functions
% addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\Lp_norm_codes\Dual_norm_1D\sliceomatic\
addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB
% Set up 1D problem

%% INPUT VARIABLES

% Model space
nx = 500;

% Kernel
nk = 5; %number of frequencies

decay = -0.20;

basis = 0.5;

% Data noise
amp_pct = 0.05;
floor_pct = 0;0.005;

tol = 0.05;
% Lp-norm parameters
pQ = 1;%[0 1 2];%0:0.1:2;%0%
qQ = [0 1 2];%[0 1 2];%0:0.1:2;%0%
lQ = 0;%ones(1,1)*1.0%[0.25:0.25:1.75]%0.1:0.4:1.9;

% Add an extra test to find best starting beta for lplq
chi_vec = 1%[2 3 4 5 7.5 10 20 30];
% multip = 100;

% Percentile for cutoff
pct_cutoff = 75;

eps_p = 1e-5;
eps_q = 1e-5;


%% SCRIPT STARTS HERE
%% % Generate kernel functions and depth weighting
x = (0 : 1/(nx) : 1);
mcell=nx;
dx=ones(mcell,1) ;%* abs((min(x)-max(x))/mcell);


set(figure, 'Position', [50 0 775 1000]);
% plot_h = tight_subplot(2,2,[.05 .01],[.1 0.1]);

G = zeros(nk,mcell);
G1 = zeros(nk,mcell);
G2 = zeros(nk,mcell);
data = zeros(nk,1);
% wr = zeros(1,mcell);
subplot(2,2,2)

% Gerate functions of the form x cos(wx)
for ii = 1 : nk 
    
    b = basis * 2*pi* (ii);
%     a = decay * (ii);
    
%     G1(ii,1:mcell) = ( x(2:end).*sin(b*x(2:end)) + cos(b*x(2:end)) ) -...
%         ( x(1:end-1).*sin(b*x(1:end-1)) + cos(b*x(1:end-1)) );
% 
%     G(ii,1:mcell) = ( x(2:end).*cos(b*x(2:end)) - (sin(b*x(2:end)) + x(2:end).*sin(b*x(2:end))) -...
%         ( x(1:end-1).*sin(b*x(1:end-1)) );
%     G(ii,1:mcell) = ( x(2:end).*sin(b*x(2:end))/b + cos(b*x(2:end))/(b^2) ) -...
%         ( x(1:end-1).*sin(b*x(1:end-1))/b + cos(b*x(1:end-1))/(b^2) ) ;
    G(ii,1:mcell) = 1e+5*(sin(b*x(2:end)) - sin(b*x(1:end-1)))/b;

      
%     wr = wr + ((exp( a*x(2:end))-exp( a * x(1:end-1) ))/a) ;
    plot( x(2:end) , G(ii,:) );hold on

end

wr = (sum(G.^2,1)).^0.0;
wr = (wr / max(wr)).^0.5;

IWr = spdiags(1./wr',0,mcell,mcell);
Wr = spdiags(wr',0,mcell,mcell);
% Wrx = spdiags(wr(1:end-1)',0,mcell-1,mcell-1);

% figure(1); plot(wr','r--','LineWidth',3);
axis square;grid on
% temp = legend('\bfg(z) = e^{-jz} cos(j 2\pi z)','\bfDepth Weighting (Wr)');
ylabel('$\mathbf{g}(m)$', 'interpreter', 'latex','FontSize',16)
set(get(gca,'YLabel'),'Rotation',360);
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') - [0.05 0.0 0.00]);
text(.45,-1.2,'(b)', 'interpreter', 'latex','FontSize',14);
xlabel('$\mathbf{z}$', 'interpreter', 'latex','FontSize',16);
% a=get(temp,'children');
% set(a(2),'LineWidth',3,'Color',[1 0 0],'LineStyle','--'); 
%% Generate model

x = (1/nx : 1/(nx) : 1) ;

cntr = round(nx/2);

% Create mdoel: Square and tanh functions
% model = 0.4* exp(-(((x-3*cntr/nx) ).^2 * nx)); 
model = zeros(nx,1);
model(cntr-round(cntr/4):cntr+round(cntr/4))=1;

% model = model(2:end) - model(1:end-1);
% model = [model;model(1)];
% model([100 250 400]) = [-.25 .5 -0.5];
% model=model(:);

subplot(2,2,1)
plot(x,model,'LineWidth',2);axis([x(1) x(end) -0.6 0.6]);hold on
ylabel('$\mathbf{m}$', 'interpreter', 'latex','FontSize',16)
set(get(gca,'YLabel'),'Rotation',360);
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') - [0.05 0.0 0.00]);
text(.45,-0.23,'(a)', 'interpreter', 'latex','FontSize',14);
xlabel('$\mathbf{z}$', 'interpreter', 'latex','FontSize',16);
axis square
grid on


%% FORWARD MODELING
% Generate data and corrupt with noise

% data= G1 * model - G2 * model;
data= G * model;
% Corrupt data with Gaussian noise
rand_noise = floor_pct * randn(nk,1);
% save ('noise_1D','rand_noise')
% load('noise_1D');

% amp_noise = amp_pct .* abs(data);
% floor_noise = floor_pct * min(abs(data));

wd = ones(nk,1);
% noise = amp_noise + floor_noise;

d = data + rand_noise;

% wd =  abs(d)*0.05 + 0.05 * std(data);

axes('Position',[.3 .15 .4 .4])
plot(data);
hold on
errorbar(d,wd,'r*')
legend('Data','Data+noise')
ylabel('$\mathbf{d}$', 'interpreter', 'latex','FontSize',16)
set(get(gca,'YLabel'),'Rotation',360);
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') - [0.05 0.0 0.00]);
axis square
grid on
text(11.5,-9,'(c)', 'interpreter', 'latex','FontSize',14);

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
V = spdiags(sqrt(dx),0,nx,nx);
s= ones(mcell,1);
t= ones(mcell,1);

% Global constant
as = 1.0 / min(dx) ^2;  %Smallnest term
ax = 1.0;               %Smoothness term                     

% Store all the final models
% RMS = zeros(length(pQ),length(qQ),length(lQ));
% linf = zeros(length(pQ),length(qQ),length(lQ));
% l2 = zeros(length(pQ),length(qQ),length(lQ));
% l1 = zeros(length(pQ),length(qQ),length(lQ));
% misfit = zeros(length(pQ),length(qQ),length(lQ));
% models_out = zeros(length(pQ),length(qQ),length(lQ),mcell);
% iter_num = zeros(length(pQ),length(qQ),length(lQ));

Pac = speye(mcell);
            
dkdt = @(p,ep) ((ep).^(1./(4-2*p)) - ep);

counter = 1;  

objfunc = @(m,phi,b) sum( ( G * m - d ).^2 ) + ( m' * b * phi * m );
set(figure, 'Position', [50 0 775 1000]);
plot_h = tight_subplot(2,3,[.05 .01],[.1 0.1]);
% Iterate over all the combination of norms as specified by the vectors
% pvec, qvec and lvec.
for ll= 1:length(lQ)

    for pp = 1:length(pQ)
        
        for qq = 1:length(qQ)
            
            
            
            % Message prompt
            head = ['lp: ' num2str(pQ(pp)) ' lq: ' num2str(qQ(qq)) ' psi: ' num2str(lQ(ll))];
            fprintf('Starting lp inversion %s\n',...
                head)
            
            invmod      = ones(mcell,1)*1e-1;       % Initial model       
                        
            phi_d       = sum((G*invmod - d).^2);   % Initial misfit

                       % Active cell
            count=0; % Initiate iteration count 
            switcher = 0;
            ncg(counter) = 0;
            lp_count = 0;
            beta_count = 0;
            cg_iter = 0;
%             delta_p = 1;
            lambda = [];
            phi_m = [];
            phi = [];
                        
            delta_p = [];
            delta_q = [];        

            % Define two zones with different lp,lq
%             zone{1,1} = 1:75;  
%             zone{2,1} = 125:200; 
% 
%             % Transition zone with linear interp
%             zone{3,1} = 76:124;
%             linterp = cumsum(ones(length(zone{3,1}),1)/length(zone{3,1}));
            
            % Z{1} = ones(mcell,1);
            z{1} = ones(mcell,1);%zeros(nx,1); z{1}(zone{1,1}) = 1;
            trans{1} = ones(mcell,1);%z{1}; trans{1}(zone{3,1}) = 1-linterp;
            
%             z{2} = zeros(nx,1); z{2}(zone{2,1}) = 1;
%             trans{2} = z{2}; trans{2}(zone{3,1}) = 1-linterp(end:-1:1);
                        
%             pvec = ones(nx,1) * pQ(pp);%(0.*trans{1} + 2.*trans{2}) ./ (trans{1} + trans{2}) ;
%             
%             qvec = ones(nx,1) * qQ(qq);%(0.*trans{1} + 2.*trans{2}) ./ (trans{1} + trans{2}) ;
%             
%             lvec = ones(nx,1) * lQ(ll);
%             lvec(z{1}==1) = 0;
%             lvec(z{2}==1) = 2;
           LP = [pQ(pp) qQ(qq) lQ(ll)];
           dphi_m = 100;
            
        while switcher ~= 3

                    count=count+1;

                    if switcher == 0   %First iteration

                        delta_p(count) = 1e-1;%max(invmod);%prctile(abs(invmod(invmod ~= 0)),pct_cutoff);
                        dmdx = abs(Wx * invmod);
                        delta_q(count) = 1e-1;%max(dmdx);%prctile(gradm(gradm~=0),pct_cutoff);
%                         delta_q(count) = max(abs(invmod));                        
                        % Initial beta trace(G'G) / trace(phim)
 
%                         scale_p = 1 /max(abs( as * WstWs * invmod ));
%                         
%                         scale_x = 1/max(abs( ax * WxtWx * invmod ));
                        
                         
                        aVRWs = sqrt(as) * Wr * V * Ws;
                        aVRWx = sqrt(ax * (2.0 - lQ(ll))) * Wr * Vx * Wx;
                        
                        lvec = ones(mcell,1);
                        
                        MOF = aVRWs'*aVRWs + aVRWx'*aVRWx;
                        
                        if count==1
                            beta = 1e+4;%sum(sum(G.^2,1)) / sum(diag(MOF,0).^2) *1e+6;
                        end
                        
%                         lambda(count) = beta(count);
                        
                        phi(count) = objfunc(invmod,MOF,beta(count));
                        gamma = 1;
%                         x_tresh = dkdt(delta(count));

                        p_tresh = 0.01;
                        q_tresh = 0.01;
                        
                    else
                        
%                         fprintf('# # LP-LQ ITER# #');
                        
                            
                        lp_count = lp_count+1;

%                         if lp_count == 1
%                             
%                             [eps_p,eps_q] = get_eps(invmod,5,Wx,[],[]);
%                             
%                         end
                        
                        % Only matters for the third case where epsilon
                        % starts high and cooled
                        if delta_p(end)> eps_p %&& switcher == 1
                            
                            delta_p(count) = delta_p(count-1)*.5;
                            
                        else 
                            
                            delta_p(count) = delta_p(count-1);%delta_p(count-1);
                            
                        end
                        
                        if delta_q(end)> eps_q %&& switcher == 1
                            
                            delta_q(count) = delta_q(count-1)*.5;
                            
                        else

                            delta_q(count) = delta_q(count-1);%delta_q(count-1);
                            
                        end
                        
                        % Only for first scenario
                        % Check if the algo has converged and need to
                        % adjust beta
                        if dphi_m(end) < 1  %&& traffic_p(end)*100 < 2 && traffic_q(end)*100 < 2 
                            
%                             delta_q(count) = eps_q;
%                             delta_p(count) = eps_p;
                            switcher = 2;
                            
                        end
                            
       
%                         x_tresh = dkdt(pvec,delta(count));
                       
                       
                        
                        tRs = spalloc(mcell,mcell,mcell);
                        tRx = spalloc(mcell,mcell,mcell);
                        
%                         scale_p = zeros(mcell,1);
%                         scale_x = zeros(mcell,1);
                        
                        lvec = zeros(mcell,1);
                            
       
%                         x_tresh = dkdt(pvec,delta(count));
                        
                       
                        
                        tRs = spalloc(mcell,mcell,mcell);
                        tRx = spalloc(mcell,mcell,mcell);
                        
%                         scale_p = zeros(mcell,1);
%                         scale_x = zeros(mcell,1);
                        lvec = zeros(mcell,1);
                        for ii = 1 : size(s,2)
                            
                            T = spdiags(t(:,ii),0,mcell,mcell);
                            S = spdiags(s(:,ii),0,mcell,mcell);
                            
                            rs = 1./ ( abs(invmod) .^2 + delta_p(count).^2 ) .^( (1-LP(ii,1)/2) );
                            rx = 1./ ( abs(Wx * invmod) .^2 + delta_q(count).^2 ) .^( (1-LP(ii,2)/2) );

                            Rs = spdiags( rs.^0.5 ,0,mcell,mcell);
                            Rx = spdiags( rx.^0.5 ,0,mcell,mcell);
                        
                            avrws =  Rs *  Ws ;
                            avrwx =  Rx *  Wx;
                            
                            
                                
                            grad_s = delta_p(count).^ (1-LP(ii,1)/2);                        
                            grad_x = delta_q(count).^ (1-LP(ii,2)/2);

                            
                            tRs = tRs + sqrt( grad_s )  * T * Rs;
                            tRx = tRx + sqrt( grad_x )  * T * Rx;
                            
                            lvec = T * ones(mcell,1) * LP(ii,3);
                            
                        end
                                                                       
                        aVRWs =  spdiags( sqrt( as *  lvec ),0,mcell,mcell)  * Wr * tRs * V * Ws ;
                        aVRWx =  spdiags( sqrt( ax * (2.0 - lvec)  ),0,mcell,mcell) * Wr * tRx * Vx * Wx;
%                          figure(100);plot((aVRWx)'*(aVRWx) * invmod );
                           
                        gamma = phi_m(count-1)*gamma /...
                            ( invmod' * (aVRWs' * aVRWs + aVRWx' * aVRWx) * invmod );
                                                                                                 
                        aVRWs = sqrt(gamma) * aVRWs;
                        aVRWx = sqrt(gamma) * aVRWx;
                        
                        MOF = ( aVRWs' * aVRWs + aVRWx' * aVRWx );
                        
%                         mu = ( phi_m(count-1) ) / (invmod'*MOF*invmod ) ;
       
%                         beta(count) = beta(count);
                        
%                         figure(100)
%                         [mm,idx] = sort(invmod);
%                         gradphi = abs(mm)./(mm.^2 + delta_p(end).^2).^(1-pvec/2);
%                         hist(mm,15); hold on
%                         plot(abs(mm),gradphi/max(gradphi)*50,'r');
%                         plot([p_tresh p_tresh],[0 200],'--');
%                         ylim([0 55])
%                         xlim([0 0.5])
%                         hold off
                        
                    end

                
                
%                 fprintf('\n# # # # # #\n');
%                 fprintf('BETA ITER: \t %i  \nbeta: \t %8.5e \n',count,beta(count));

                m_in = invmod;  % Save initial model for backtracking steps
                

                diagA = sum(G.^2,1) + beta(count)*spdiags(MOF,0)';
                PreC     = Pac * spdiags(1./diagA(:),0,mcell,mcell);
                
                A= [ G  ;...
                sqrt( beta(count) ) * ( aVRWs ) ;...
                sqrt( beta(count) ) * ( aVRWx )] ;

                b = [(d) ; ...
                sqrt( beta(count) ) * ( aVRWs * zeros(mcell,1)) ;...
                sqrt( beta(count) ) * ( aVRWx * zeros(mcell,1))] ;

                [invmod,~,cg_iter] = PCGLSQ( invmod, A , b, PreC, Pac );
%                 tncg = tncg + ncg;
%                 [invmod,cg_iter] = GNsolver( G, invmod, d, phi(end), MOF, PreC, Pac, lambda(count) , aVRWs, aVRWx ,[],[]);             
                ncg(counter) = ncg(counter) + cg_iter;
                
                
                % Measure the change in model update
                if count==1 

                    rdm(count) =  1;
                    gdm(1) = norm(m_in - invmod);

                else

                    gdm(2) = norm(m_in - invmod);
                    rdm(count) = abs( gdm(2) - gdm(1) ) / norm(invmod);

                    gdm(1) = gdm(2);

                end
                
               
                phi_d(count) = sum((G*(invmod)-d).^2);
                phi_m(count) = (invmod)'*(MOF)*(invmod)/gamma; 
                phi_p(count) = sum((aVRWs*invmod).^2);
                phi_x(count) = sum((aVRWx*invmod).^2);
                phi(count) = objfunc(invmod,MOF,beta(count));

                if count > 1
                dphi_m(count) = abs(phi_m(count) - phi_m(count-1))/phi_m(count) *100;
                
                end

                
            
                if switcher ~= 0;
                    
                    % Plot phi_m_out
%                      aVRWs =  spdiags( sqrt( as * lvec  ),0,mcell,mcell)  * Wr * tRs * V * Ws ;
%                         aVRWx =  spdiags( sqrt((2.0 - lvec) * ax  ),0,mcell,mcell) * Wr * tRx * Vx * Wx;
            
                else
                    
                    aVRWx =  spdiags( sqrt((2.0 - lvec) * ax ),0,mcell,mcell) * Wr * Vx * Wx;

                    aVRWs =  spdiags( sqrt((lvec) * as ),0,mcell,mcell) * Wr * V * Ws;
                    
                end
                
                 phim_out(count) = invmod'*(aVRWx'*aVRWx + aVRWs'*aVRWs) * invmod;
                 
                
                
            % Check to see if overshooted the target misfit,               
            if phi_d(count) < target * (1-tol) || phi_d(count) > target * (1+tol)%&& switcher~=1
  
                if switcher == 0 && phi_d(count) < target * (1-tol)
                    
                    % Come back by a fraction and start IRLS
                    beta(count+1) = beta(count);
                    switcher = 1;                    

                elseif switcher == 0 
                    
                    beta(count+1) = 0.9 * beta(count);
                    
                elseif phi_d(count) > target * (1+tol)
                    
                    beta(count+1) = beta(count)*(1-tol);
                    
                else
                    
                   beta(count+1) = beta(count) * target / phi_d(count);
                   
                end
          
            % Else keep beta and iterate               
            else

                if switcher == 0
                    switcher = 1;
                   beta_count = 0;
                end
                
                beta_count = beta_count+1;
                beta(count+1) = beta(count);

                 
            end
           
            if lp_count > 0  && dphi_m(end) < tol && switcher==2 && phi_d(count) > target *(1-tol) && phi_d(count) < target *(1+tol)

                    switcher = 3;

                    continue

            end      
            
            % OPTIONAL: Plot model
%             set(figure(4), 'Position', [50 200 1000 500])
% axes(plot_h((ll)*2+1));
%             axes(plot_h(counter));
%             plot(1:mcell,model); hold on
%             plot(1:mcell, invmod,'r','LineWidth',3); hold off
%             axis([0 mcell -.1 0.6]);
%             axis square

            
            
        end
                    % OPTIONAL: Plot model

            axes(plot_h(counter));
            plot(1:mcell,model); hold on
            plot(1:mcell, invmod,'r','LineWidth',1); hold off
            axis([0 mcell min([min(model) min(invmod)]) max([max(model) max(invmod)]) ]);
            arg = ['\phi_d:' num2str(phi_d(end))];
            title(arg);
            text(25,0.9*max([max(model) max(invmod)]),['l2= ' num2str(round(sum(invmod.^2)))]);
            text(25,0.75*max([max(model) max(invmod)]),['l1= ' num2str(round(sum((invmod.^2)./((invmod.^2 + eps_p.^2).^0.5))))]);
            text(25,0.6*max([max(model) max(invmod)]),['l0= ' num2str(round(sum((invmod.^2)./(invmod.^2 + eps_p.^2))))]);
            text(400,0.9*max([max(model) max(invmod)]),['$\mathbf{l_' num2str(pQ(pp)) '}$'],'interpreter','latex','FontSize',12);

            axis square

            axes(plot_h(counter+3));
%             plot(1:mcell,model); hold on
            plot(1:mcell, cumsum(invmod).*dx,'r','LineWidth',3); hold off
            axis([0 mcell -.2 1]);
            axis square

            ldml(counter) = norm(model-invmod,2);
%             fprintf('Relative dm:\t %8.5e \n', rdm(count));
%             fprintf('End of lp inversion. Number of iterations: %i\n',count)
%             fprintf('Final model misfit: %8.3e\n\n',ldml(counter))
            
            counter = counter+1;
            

        end
    end
end


