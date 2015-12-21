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
nx = 200;

% Kernel
nk = 20; %number of frequencies

decay = -0.20;

basis = 0.20;

% Data noise
amp_pct = 0.05;
floor_pct = 0.02;

% Lp-norm parameters
LP = [2 2 1
	 0 0 2
	0 0 2];

% Add an extra test to find best starting beta for lplq
% chi_vec = 1%[2 3 4 5 7.5 10 20 30];
% multip = 100;

% Percentile for cutoff
pct_cutoff = 75;
eps_p = [1e-8 1e-3 1e-3];
eps_q = [1e-8 1e-8 1e-8];

%% SCRIPT STARTS HERE
%% % Generate kernel functions and depth weighting
z = (0 : 1/(nx) : 1);
mcell=nx;
dx=ones(1,mcell) * abs((min(z)-max(z))/mcell);


set(figure, 'Position', [50 0 775 1000]);
% plot_h = tight_pubplot(2,2,[.05 .01],[.1 0.1]);

G = zeros(nk,mcell);
wr = zeros(1,mcell);
subplot(2,2,2)
for ii = 1 : nk 
    
    b = basis * 2*pi* (ii);
    a = decay * (ii);
    
    G(ii,1:mcell) = exp( a * z(2:end) ) .* (a/(a^2+b^2)*cos (b * z(2:end)) + b/(a^2+b^2) *sin (b * z(2:end))) -...
        exp( a * z(1:end-1) ) .* (a/(a^2+b^2)*cos (b * z(1:end-1)) + b/(a^2+b^2) *sin (b * z(1:end-1)));

    G(ii,1:mcell) = G(ii,1:mcell) / max(G(ii,1:mcell));
    
    wr = wr + ((exp( a*z(2:end))-exp( a * z(1:end-1) ))/a) ;
    plot( z(2:end) , G(ii,:) );hold on

end

% wr = (sum(G.^2,1)).^0.5;
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

z = (1/nx : 1/(nx) : 1) ;

cntr = round(nx/4);

% Create mdoel: Square and tanh functions
model = 0.4* exp(-(((z-3*cntr/nx) ).^2 * nx)); 

model(cntr-round(cntr/4):cntr+round(cntr/4))=0.25;

model=model(:);

subplot(2,2,1)
plot(z,model,'LineWidth',2);axis([z(1) z(end) -0.1 0.5]);hold on
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

data= G * model;

% Corrupt data with Gaussian noise
% rand_noise = rand(nk,1);
% save ('noise_1D','rand_noise')
load('noise_1D');

amp_noise = amp_pct .* abs(data);
floor_noise = floor_pct * min(abs(data));

wd = amp_noise + floor_noise;
% noise = amp_noise + floor_noise;

d = data + rand_noise.*wd;

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
V = spdiags(sqrt(dx'),0,nx,nx);


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

Pac = speye(mcell);
            
dkdt = @(p,ep) ((ep).^(1./(4-2*p)) - ep);

counter = 1;  

objfunc = @(m,phi,b) sum( ( G * m - d ).^2 ) + ( m' * b * phi * m );
% set(figure, 'Position', [50 0 775 1000]);
% plot_h = tight_pubplot(3,3,[.05 .01],[.1 0.1]);
% Iterate over all the combination of norms as specified by the vectors
% pvec, qvec and lvec.
set(figure, 'Position', [50 0 775 1000]);
plot_h = tight_subplot(2,3,[.01 .01],[.25 0.05]);
for ll= 1:3           
            
            
            % Message prompt
%             head = ['lp: ' num2str(pQ(ll)) ' lq: ' num2str(qQ(ll)) ' psi: ' num2str(lQ(ll))];
%             fprintf('Starting lp inversion %s\n',...
%                 head)
            
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
            gamma = [];
            phi_m = [];
            grad_MOF= [];
            phi = [];
            
            % Paremeters for epsilon cooling
            traffic_p = 1;
            traffic_q = 1;
            group_p = 0;
            p_tresh = 1;
            q_tresh = 1;
            
            delta_p = [];
            delta_q = [];
                        
            WxtWx = ( Wr * Vx * Wx)' * ( Wr * Vx * Wx );
            WstWs = ( Wr * V * Ws)' * ( Wr * V * Ws);

            % Define two zones with different lp,lq
%             zone{1,1} = 1:75;  
%             zone{2,1} = 125:200; 
% 
%             % Transition zone with linear interp
%             zone{3,1} = 76:124;
%             linterp = cumsum(ones(length(zone{3,1}),1)/length(zone{3,1}));
%             
%             % Z{1} = ones(mcell,1);
%             Z{1} = ones(mcell,1);%zeros(nx,1); z{1}(zone{1,1}) = 1;
%             trans{1} = ones(mcell,1);%z{1}; trans{1}(zone{3,1}) = 1-linterp;
%             
%             z{2} = zeros(nx,1); z{2}(zone{2,1}) = 1;
%             trans{2} = z{2}; trans{2}(zone{3,1}) = 1-linterp(end:-1:1);
                        
%             pvec = ones(nx,1) * pQ(ll);%(0.*trans{1} + 2.*trans{2}) ./ (trans{1} + trans{2}) ;
%             
%             qvec = ones(nx,1) * qQ(ll);%(0.*trans{1} + 2.*trans{2}) ./ (trans{1} + trans{2}) ;
%             
%             lvec = ones(nx,1) * lQ(ll);
%             lvec(z{1}==1) = 0;
%             lvec(z{2}==1) = 2;
           
            traffic_p = 1;
            traffic_q = 1;
            group_p = mcell;
            dphi_m = 1;
        while switcher ~= 3

                    count=count+1;

                    if switcher == 0   %First iteration

                        if ll == 2
                        delta_p(count) = 1e-1;%prctile(abs(invmod(invmod ~= 0)),pct_cutoff)*2;
                        delta_q(count) = 1e-1;%abs(prctile(gradm,pct_cutoff)*2);
                 
                        else
                        delta_p(count) = eps_p(ll);%prctile(abs(invmod(invmod ~= 0)),pct_cutoff)*2;
                        delta_q(count) = eps_q(ll);%abs(prctile(gradm,pct_cutoff)*2);
                        
                        end
%                         delta_q(count) = max(abs(invmod));                        
                        % Initial beta trace(G'G) / trace(phim)

%                         scale_p = 1 /max(abs( as * WstWs * invmod ));
%                         
%                         scale_x = 1/max(abs( ax * WxtWx * invmod ));
                        
                         
                        aVRWs = sqrt(as) * Wr * V * Ws;
                        aVRWx = sqrt(ax) * Wr * Vx * Wx;
                        
                        MOF = aVRWs'*aVRWs + aVRWx'*aVRWx;
                        
                        if count==1
                            MOF_l2 =MOF;
                            beta = sum(sum(G.^2,1)) / sum(diag(MOF,0).^2) *1e+1;
                        end
                        
%                         lambda(count) = beta(count);
                        
                        phi(count) = objfunc(invmod,MOF,beta(count));
                        gamma(count) = 1;
%                         x_tresh = dkdt(delta(count));
                        lvec = ones(mcell,1);
                        p_tresh = 0.01;
                        q_tresh = 0.01;
                        
                    else
                        
%                         fprintf('# # LP-LQ ITER# #');
                        if switcher ~=2
                            lp_count = lp_count+1;
                            
                        else
                            
                            beta_count = beta_count+1;
                            
                        end
                        
                       
                        if lp_count == 1 && switcher == 1
                             
%                              [eps_p,eps_q] = get_eps_v2(invmod,LP(ll,1),LP(ll,2),Wx,[],[]);
%                             [eps_p,eps_q] = get_eps(invmod,10,Wx,[],[]);
%                             eps_q= 1e-4;
                            p_tresh = eps_p(ll);%0.01;%prctile(abs(invmod(invmod > 0)),pct_cutoff);
%                             delta_p(count) = eps_p;

        %                     gradm = Wx * invmod;
        %                     eps_x = std(gradm);
                            q_tresh = eps_q(ll);%prctile(abs(gradm(gradm > 0)),pct_cutoff);
%                             delta_q(count) = eps_q;%eps_q;%x_tresh*2;
                            
%                             eps_q = prctile(abs(gradm),mpct(idx));
                            
                            %% TEST ZHADNOV ALGO
%                             s_test = zeros(1,10);
%                             eps_temp = 10.^(-(1:10)/2);
%                             for ee = 1 : 10
%                                 
%                                 
%                                 r = spdiags(invmod.^2 + eps_temp(ee).^2,0,mcell,mcell);
%                                 s_test(ee) = invmod'*r*invmod;
%                             end
%                             
%                             figure;plot(eps_temp,s_test);
                            
                        end
                        
                        if switcher == 1 && delta_p(end)> eps_p(ll) %&& switcher == 1 && lp_count>1;%traffic_p(end)*100 > 1 && switcher == 1
                            
                            delta_p(count) = delta_p(count-1)/2;
                            
                        else 
                            
                            delta_p(count) = eps_p(ll);
                            
                        end
                        
                        if switcher == 1 && delta_q(end)> eps_q(ll) %&& switcher == 1 && lp_count>1;%traffic_q(end)*100 > 1 
                            
                            delta_q(count) = delta_q(count-1)/2;
                            
                        else

                            delta_q(count) = eps_q(ll);
                            
                        end
                        
                        if dphi_m(end) < 1 && switcher ~= 2 %&& delta_q(end)< q_tresh% delta_p(end)<= eps_p && delta_q(end)<= eps_q && switcher ~= 2;%traffic_p(end)*100 < 1 && traffic_q(end)*100 < 1 && switcher == 1
                            
%                             delta_q(count) = eps_q;
%                             delta_p(count) = eps_p;
                            switcher = 2;
                            
                        end
                            
                       
                        
                        tRs = zeros(mcell,1);
                        tRx = zeros(mcell,1);

                        
                        lvec = zeros(mcell,1);
                        
                        t = ones(mcell,1);
                        for ii = 1 : size(t,2)
                            
                            T = spdiags(t(:,ii),0,mcell,mcell);


                            rs = sqrt(1./ ( abs(invmod) .^2 + delta_p(count).^2 ) .^( (1-LP(ll,1)/2) ));
                            rx = sqrt(1./ ( abs(Wx * invmod) .^2 + delta_q(count).^2 ) .^( (1-LP(ll,2)/2) ));

                            Rs = spdiags( rs ,0,mcell,mcell);
                            Rx = spdiags( rx ,0,mcell,mcell);

                            avrws =  Rs ;
                            avrwx =  Rx * Wx;

                            grad_s =(avrws)'*(avrws) * invmod ;                        
                            grad_x =(avrwx)'*(avrwx) * invmod ;

                            eta_s =1;%eps_p^(1-LP(ll,1)/2);
                            eta_x =1;%eps_q^(1-LP(ll,2)/2);

                            tRs = tRs + sqrt( eta_s )  * T * rs;
                            tRx = tRx + sqrt( eta_x )  * T * rx;

                            lvec = lvec + T * ones(mcell,1) * LP(ll,3);
                        end
                        
                        tRs = spdiags( tRs , 0 , mcell , mcell);
                        tRx = spdiags( tRx , 0 , mcell , mcell);
                        
                        aVRWs =  spdiags( sqrt( as * lvec  ),0,mcell,mcell)  * Wr * tRs * V * Ws ;
                        aVRWx =  spdiags( sqrt((2.0 - lvec) * ax  ),0,mcell,mcell) * Wr * tRx * Vx * Wx;
                        
                        MOF =  ( aVRWs' * aVRWs + aVRWx' * aVRWx );
                        
%                         if ll == 2 || ll ==1
                            
                            gamma(count) = 1;
                            
%                         else
%                             
%                             gamma(count) = phi_m(count-1) /...
%                                 ( invmod' * (aVRWs' * aVRWs + aVRWx' * aVRWx) * invmod );
% 
% %                             gamma(count) = grad_MOF(count-1) /...
% %                                 norm( (aVRWs' * aVRWs + aVRWx' * aVRWx) * invmod );
% 
% %                             gamma(count) = (max(MOF_l2*invmod) /(max(MOF*invmod)));
%                             
%                         end                                              
                        
                        aVRWs = sqrt(gamma(count)) * aVRWs;
                        aVRWx = sqrt(gamma(count)) * aVRWx;
                        
                        MOF =  ( aVRWs' * aVRWs + aVRWx' * aVRWx );
                        
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

                group_p(count) = sum(abs(invmod) <= p_tresh);
                group_q(count) = sum(abs(Vx*Wx*invmod) <= q_tresh);
                
%                 fprintf('\n# # # # # #\n');
%                 fprintf('BETA ITER: \t %i  \nbeta: \t %8.5e \n',count,beta(count));

                m_in = invmod;  % Save initial model for backtracking steps
                grad_MOF(count) = norm((MOF)*(invmod));

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
                
                
            % Check to see if overshooted the target misfit,               
            if phi_d(count) < target * 0.95 %&& switcher~=1
                
                phi_d(count) = sum((G*(invmod)-d).^2);                       
                
                
                if switcher == 0
                    
                    % Come back almost a full beta step
                    beta(count+1) = beta(count);
                
                    switcher = 1;                    

                else
                   
                    beta(count) = beta(count) * target / phi_d(count);
                    invmod = m_in;
                    count = count -1;
                    beta_count = beta_count + 1;
                    continue
                    
                end
          
            % Else reduce beta and continue inversion
            elseif phi_d(count) > target * 1.05 %&& switcher == 0

                if switcher == 0
                    
                    beta(count+1) = 0.9 * beta(count); 
                    
                else
                                        
                    beta(count) = beta(count) * target / phi_d(count);
                    invmod = m_in;
                    count = count -1;
                    beta_count = beta_count + 1;
                    continue
                    
                end
                
            else

                if switcher == 0
                    switcher = 1;
                end
                
                beta(count+1) = beta(count);
                beta_count = beta_count + 1;
                 
            end
           
            phi_m(count) = (invmod)'*(MOF)*(invmod);                
            phi_p(count) = sum(aVRWs*invmod).^2;
            phi_x(count) = sum(aVRWx*invmod).^2;
            phi(count) = objfunc(invmod,MOF,beta(count));

            if count > 1
                temp = sum(abs(invmod) <= p_tresh);                     
                traffic_p(count) = abs(group_p(count) - temp) / group_p(count);

                temp = sum(abs(Vx*Wx*invmod) <= q_tresh);
                traffic_q(count) = abs(group_q(count) - temp) / group_q(count);

                dphi_m(count) = abs(phi_m(count) - phi_m(count-1))/phi_m(count) *100;
            end

            if switcher > 0  && delta_p(count) <= eps_p(ll) && dphi_m(end) < 1 && phi_d(count) > target *0.95 && phi_d(count) < target *1.05

                switcher = 3;

                continue

            end    
                
            % OPTIONAL: Plot model
%             set(figure(4), 'Position', [50 200 1000 500])
            axes(plot_h(counter));
            plot(1:mcell,model); hold on
            plot(1:mcell, invmod,'r','LineWidth',3); hold off
            axis([0 mcell -.1 0.6]);
            axis square

            
            
        end
                    % OPTIONAL: Plot model
%             set(figure, 'Position', [50 200 1000 500])

%             RMS(pp,qq,ll)=(sum((invmod-model).^2)/length(data))^(0.5);
%             linf(pp,qq,ll) = norm(model-invmod,'inf');
%             l2(pp,qq,ll) = norm(model-invmod,2);
%             l1(pp,qq,ll) = norm(model-invmod,1);
%             misfit(pp,qq,ll) = phi_d(count);
%             models_out(pp,qq,ll,:) = invmod(:);
%%           
            axes(plot_h(ll));
            plot(z,model); hold on
            plot(z, invmod,'r','LineWidth',2); 
            text(0.30,0.55,['\underline{Final Iteration}'],'interpreter', 'latex','FontSize',10)
            text(0.30,0.50,['$\mathbf{\phi_d} = ' num2str(round(phi_d(count)*100)/100) '$'],'interpreter', 'latex','FontSize',10)
            text(0.30,0.45,['$\epsilon = 10^{' num2str(round(log10(delta_p(end))*10)/10) '}$'],'interpreter', 'latex','FontSize',10)
            text(0.30,0.40,['$\|\delta m\|_1$ = ' num2str(round(norm(model-invmod,1)*100)/100) ],'interpreter', 'latex','FontSize',10)
                            xlabel('$z$','interpreter', 'latex','FontSize',16);

            text(0.30,0.35,['$k^{th} iter=\;' num2str(beta_count) '$'],'interpreter', 'latex','FontSize',10)
            
            if ll == 1
            set(gca, 'YAxisLocation', 'left')
            ylabel('$m$','interpreter', 'latex','FontSize',16);
            set(get(gca,'YLabel'),'Rotation',360);
            else
                set(gca,'YTickLabel',[]);
            end
            if ll ==1 
                text(-0.05,0.7,'$\phi_m=\|\mathbf{W}_s \mathbf{m} \|_2^2 + \| \mathbf{W}_x \mathbf{G}_x \mathbf{m} \|_2^2$','interpreter', 'latex','FontSize',12)
%                 set(gca,'XTickLabel',[]);
                text(0.45,-.35,'$(a)$','interpreter', 'latex','FontSize',16)
                
            elseif ll == 2
                text(0.15,0.7,['$\phi_m= \| \mathbf{W}_x \; \mathbf{R}_x \; \mathbf{G}_x \; \mathbf{m} \|_2^2$'],'interpreter', 'latex','FontSize',12)
%                 set(gca,'XTickLabel',[]);
                text(0.45,-.35,'$(b)$','interpreter', 'latex','FontSize',16)
            else
                text(0.1,0.7,['$\phi_m= \gamma^{(k)} [\| \mathbf{W}_x \; \mathbf{R}_x \; \mathbf{G}_x \; \mathbf{m} \|_2^2]$'],'interpreter', 'latex','FontSize',12)
                xlabel('$z$','interpreter', 'latex','FontSize',16);
                text(0.45,-.35,'$(c)$','interpreter', 'latex','FontSize',16)
                
            end
            
            


            hold off
            grid on
            axis([z(1) z(end) -.1 0.6]);
%             set(gca,'YTickLabel',[]);
            axis square
%             tlabel=(['$\mathbf{\phi_m\;=\;\| m \|_{' num2str(pQ(pp)) '}\;+\; \| \nabla m \|_{' num2str(qQ(qq)) '}}$']);
%             title(tlabel, 'interpreter', 'latex','FontSize',18)
            
            n_iter = count - (lp_count+beta_count);
            axes(plot_h(ll+3));           
            [h,line1,line2] = plotyy(1:count,phi_d,1:count,phi_m,'semilogy'); hold on
            set(line1,'LineWidth',1.5);
            set(line2,'LineWidth',1.5);
            
            grid on
            xlim(h(1),[0 count])
            xlim(h(2),[0 count])
            ylim(h(1),[1e-1 5e+5])
            ylim(h(2),[1e-1 5e+5])
            axis(h(1),'square')
            axis(h(2),'square')
%             set(gca, 'YAxisLocation', 'right')
% if ll == 2
            xlabel('Iteration')
            
% end
if ll ==1
%             ylabel(h(1),'$\phi_d$' ,'interpreter', 'latex','FontSize',14);
%             set(get(h(1),'YLabel'),'Rotation',360);
%             set(ylabh,'Position',get(ylabh,'Position') - [4 0 0.00]);
            text(-5,1e+5,'$\phi_d$' ,'interpreter', 'latex','FontSize',14,'Color','b');
else
    
            set(h(1),'YTickLabel',[]);
end
if ll ==3
%             ylabh = get(h(2),'YLabel');
%             set(ylabh,'Position',get(ylabh,'Position') + [3 0 0.00]);
%             ylabel(h(2),'$\phi_m$' ,'interpreter', 'latex','FontSize',14);
%             set(get(h(2),'YLabel'),'Rotation',360);
            text(count+2,1e+5,'$\phi_m$' ,'interpreter', 'latex','FontSize',14,'Color','r');
            set(h(2),'YTickLabel',[])
else
    
    set(h(2),'YTickLabel',[]);
end          
            
            n_iter = count - lp_count - beta_count;
            xx = [0 0 n_iter n_iter 0];
            yy = [1e-1 5e+5 5e+5 1e-1 1e-1];
            h = fill(xx,yy,'b');
            set(h,'FaceAlpha',0.1,'LineStyle','none');
            
                       
%             n_iter = n_iter + beta_count;
            xx = [n_iter n_iter n_iter+lp_count n_iter+lp_count n_iter];
            yy = [1e-1 5e+5 5e+5 1e-1 1e-1];
            h = fill(xx,yy,'r');
            set(h,'FaceAlpha',0.1,'LineStyle','none');
            
            plot([0 count],[target target],'b--');
            text(0,target,'$\phi_d^*$','interpreter', 'latex','FontSize',12,'Color','b','BackgroundColor','w','EdgeColor','b')
            
            % Plot phi_m*
            rx = 1./ ( abs(Wx * model) .^2 + delta_q(count).^2 ) .^( (1-LP(ll,2)/2)/2 );
            Rx = spdiags( rx ,0,mcell,mcell);
            aVRWx =   sqrt((2.0 - LP(ll,3)) * ax )* Wr * Vx * Rx * Wx;
            
            rs = 1./ ( abs(Ws * model) .^2 + delta_p(count).^2 ) .^( (1-LP(ll,1)/2)/2 );
            Rs = spdiags( rs ,0,mcell,mcell);
            aVRWs =   sqrt((LP(ll,3)) * as )* Wr * V * Rs * Ws;
            
            phim_true = model'*(aVRWx'*aVRWx + aVRWs'*aVRWs) * model;
            
            plot([0 count],[phim_true phim_true],'r--');
            text(1,phim_true,'$\phi_m^*$','interpreter', 'latex','FontSize',12,'Color','r','BackgroundColor','w','EdgeColor','r','HorizontalAlignment','left')
            

            if ll ==1 
                
%                 set(gca,'XTickLabel',[]);
%             title('$\phi_d$' ,'interpreter', 'latex','FontSize',16);
            elseif ll == 2
                
%                 set(gca,'XTickLabel',[]);
                
            else
                % Create legend for color scheme
                axes('Position',[.36 .21 .28 .03])
                axis([0 3 0 1])

               h = fill([0 0 1.0 1.0 0],[0 1 1 0 0],'b');set(h,'FaceAlpha',0.1); hold on
                h = fill([1.0 1.0 2 2 1.0],[0 1 1 0 0],'r');set(h,'FaceAlpha',0.1);
                h = fill([2.0 2.0 3 3 2.0],[0 1 1 0 0],'w');set(h,'FaceAlpha',0.1);
                set(gca,'XTickLabel',[])
                set(gca,'YTickLabel',[])
%                 set(gca,'Visible','off');  
                
                text(0.5,0.4,'$l_2-norm$','interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
                text(1.5,0.4,'$IRLS$','interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
                text(2.5,0.4,'$\phi_d^{(k)} \rightarrow \phi_d^*$','interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')

            end
%             figure; plotyy(1:count,phi_d,1:count,phi_m);

            ldml(counter) = norm(model-invmod,2);
%             fprintf('Relative dm:\t %8.5e \n', rdm(count));
%             fprintf('End of lp inversion. Number of iterations: %i\n',count)
%             fprintf('Final model misfit: %8.3e\n\n',ldml(counter))
            
            counter = counter+1;
            

end

