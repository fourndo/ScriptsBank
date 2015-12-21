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
addpath 'C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\Thesis_Plots'
% addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\Lp_norm_codes\Dual_norm_1D\functions
% addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\Lp_norm_codes\Dual_norm_1D\sliceomatic\
addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB
% Set up 1D problem

%% INPUT VARIABLES

% Model space
nx = 200;

% Kernel
nk = 30; %number of frequencies

decay = -0.20;

basis = 0.20;

% Data noise
amp_pct = 0.05;
floor_pct = 0.02;

% Lp-norm parameters
LP = [2 2 1
     0 0 0
	 0 0 0
     0 0 0
     0 0 0
	0 0 0];

% Add an extra test to find best starting beta for lplq
% chi_vec = 1%[2 3 4 5 7.5 10 20 30];
% multip = 100;

% Percentile for cutoff
pct_cutoff = 75;
eps_p = (1e-1);
eps_q = (1e-8);

tol = 0.1;
%% SCRIPT STARTS HERE
%% % Generate kernel functions and depth weighting
z = (0 : 1/(nx) : 1);
mcell=nx;
dx=ones(1,mcell) * abs((min(z)-max(z))/mcell);


set(figure, 'Position', [50 0 775 775]);
% plot_h = tight_pubplot(2,2,[.05 .01],[.1 0.1]);

G = zeros(nk,mcell);
wr = zeros(1,mcell);
axes('Position',[0.6 .575 .38 .38]);
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
ylabel('$\mathbf{f}_j(x)$', 'interpreter', 'latex','FontSize',14)
set(get(gca,'YLabel'),'Rotation',360);
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') - [0.05 0.0 0.00]);
text(.9,0.9,'(b)', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle');
xlabel('$x$', 'interpreter', 'latex','FontSize',16);
% a=get(temp,'children');
% set(a(2),'LineWidth',3,'Color',[1 0 0],'LineStyle','--'); 
%% Generate model

z = (1/nx : 1/(nx) : 1) ;

cntr = round(nx/4);

% Create mdoel: Square and tanh functions
model = 0.4* exp(-(((z-3*cntr/nx) ).^2 * nx)); 

model(cntr-round(cntr/4):cntr+round(cntr/4))=0.25;

model=model(:);

axes('Position',[0.1 .575 .38 .38]);
plot(z,model,'LineWidth',2);axis([z(1) z(end) -0.1 0.5]);hold on
ylabel('$\mathbf{m}$', 'interpreter', 'latex','FontSize',14)
set(get(gca,'YLabel'),'Rotation',360);
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') - [0.05 0.0 0.00]);
text(.9,0.45,'(a)', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle');
xlabel('$x$', 'interpreter', 'latex','FontSize',16);
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

axes('Position',[.3 .05 .4 .4])
plot(data);
hold on
errorbar(d,wd,'r*')
legend('Data','Data+noise','Location','SouthEast')
ylabel('$\mathbf{d}$', 'interpreter', 'latex','FontSize',16)
set(get(gca,'YLabel'),'Rotation',360);
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') - [0.05 0.0 0.00]);
xlim([0 length(data)])
axis square
grid on
xlabel('$i^{th} datum$','interpreter','latex')
text(27.5,13,'(c)', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle');

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
set(figure, 'Position', [50 0 775 775]);
set(figure, 'Position', [50 0 775 775]);
set(figure, 'Position', [50 0 775 775]);
% plot_h = tight_subplot(2,3,[.01 .01],[.25 0.05]);
for ll= 1:size(LP,1)           
            
            
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
            l2_count = 0;
            beta_count = 0;
            cg_iter = 0;
%             delta_p = 1;
            gamma = [];
            phi_m = [];
            phi_d = [];
            grad_MOF= [];
            phi = [];
            phi_lp=[];
            
            % Paremeters for epsilon cooling
            traffic_p = 1;
            traffic_q = 1;
%             group_p = 0;
%             p_tresh = 1;
%             q_tresh = 1;
            
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
           
%             traffic_p = 1;
%             traffic_q = 1;
%             group_p = mcell;
            dphi_m = 1;
            gamma = 1;
            l2_count = 0;
        while switcher ~= 3

                    count=count+1;

                    if switcher == 0   %First iteration

                        l2_count = l2_count+1;
                        delta_p(count) = eps_q;%prctile(abs(invmod(invmod ~= 0)),pct_cutoff)*2;
                        gradm = abs(Wx * invmod);
                        
                        if ll ~= 5 && ll ~= 6
                            
                            delta_q(count) = 1e-8;
                            
                        else
%                             
                            delta_q(count) = 0.01;%abs(prctile(gradm,pct_cutoff)*2);
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
                        
%                         x_tresh = dkdt(delta(count));
                        lvec = ones(mcell,1);
                        p_tresh = 0.01;
                        q_tresh = 0.01;
                        
                    else
                        
%                         fprintf('# # LP-LQ ITER# #');
%                         if switcher ~= 2
                            lp_count = lp_count+1;
                                                     
%                         end
                        
                        
                        if lp_count == 1 && ll == 6
                                
                            
%                             [eps_p,eps_q] = get_eps(invmod,5,Wx,[],[]);
                            [eps_p,eps_q] = get_eps_v2(invmod,LP(ll,1),LP(ll,2),Wx,[],[]);
                            [~,~] = get_eps_Lcurve(invmod,10,Wx,[],[]);
                                                    
                        end
                        
                        % Only matters for the third case where epsilon
                        % starts high and cooled
                        if delta_p(end)> eps_p && switcher == 1
                            
                            delta_p(count) = delta_p(count-1)*.5;
                            
                        else 
                            
                            delta_p(count) = delta_p(count-1);%delta_p(count-1);
                            
                        end
                        
                        if delta_q(end)> eps_q && switcher == 1
                            
                            delta_q(count) = delta_q(count-1)*.5;
                            
                        else

                            delta_q(count) = delta_q(count-1);%delta_q(count-1);
                            
                        end
                        
                        % Only for first scenario
                        % Check if the algo has converged and need to
                        % adjust beta
                        if dphi_m(end) < 2  && ((lp_count > 1 && delta_q(end)<= eps_q) || ll==4 || ll==5 ) %&& traffic_p(end)*100 < 2 && traffic_q(end)*100 < 2 
                            
%                             delta_q(count) = eps_q;
%                             delta_p(count) = eps_p;
                            switcher = 2;
                            
                        end
                            
       
% %                         x_tresh = dkdt(pvec,delta(count));
%                         rs = 1./ ( abs(invmod) .^2 + delta_p(count).^2 ) .^( (1-pvec/2)/2 );
%                         rx = 1./ ( abs(Wx * invmod) .^2 + delta_q(count).^2 ) .^( (1-qvec/2)/2 );
%                         
%                         Rs = spdiags( rs ,0,mcell,mcell);
%                         Rx = spdiags( rx ,0,mcell,mcell);
                       
                        
                        tRs = zeros(mcell,1);
                        tRx = zeros(mcell,1);
                        
%                         scale_p = zeros(mcell,1);
%                         scale_x = zeros(mcell,1);
                        
                        lvec = zeros(mcell,1);
                        
                        t = ones(mcell,1);
                        for ii = 1 : size(t,2)
                            
                            T = spdiags(t(:,ii),0,mcell,mcell);
%                             S = spdiags(s(:,ii),0,mcell,mcell);
                        
                         rs = sqrt(1./ ( abs(invmod) .^2 + delta_p(count).^2 ) .^( (1-LP(ll,1)/2) ));
                        rx = sqrt(1./ ( abs(Wx * invmod) .^2 + delta_q(count).^2 ) .^( (1-LP(ll,2)/2) ));
                        
                        Rs = spdiags( rs ,0,mcell,mcell);
                        Rx = spdiags( rx ,0,mcell,mcell);
                        
                            avrws =  Rs ;
                            avrwx =  Rx * Wx;
                            
                            grad_s =(avrws)'*(avrws) * invmod ;                        
                            grad_x =(avrwx)'*(avrwx) * invmod ;
                                                        
                            tRs = tRs + T * rs;
                            tRx = tRx + T * rx;
 
                            lvec = lvec + T * ones(mcell,1) * LP(ll,3);
                            
                        end
                        
                        tRs = spdiags( tRs , 0 , mcell , mcell);
                        tRx = spdiags( tRx , 0 , mcell , mcell);
                        
                        aVRWs =  spdiags( sqrt( as * lvec  ),0,mcell,mcell)  * Wr * tRs * V * Ws ;
                        aVRWx =  spdiags( sqrt((2.0 - lvec) * ax  ),0,mcell,mcell) * Wr * tRx * Vx * Wx;
                        
                        MOF =  ( aVRWs' * aVRWs + aVRWx' * aVRWx );
                        
                        % The third scenario scales beta between each IRLS
                        if ll == 4 || ll == 5 || ll == 6 
                                                       
                            gamma = (phi_m(count-1)*gamma) /...
                                ( invmod' * (aVRWs' * aVRWs + aVRWx' * aVRWx) * invmod );
                            
                        end                                              
                        
                        aVRWs = sqrt(gamma) * aVRWs;
                        aVRWx = sqrt(gamma) * aVRWx;
                        
                        MOF =  ( aVRWs' * aVRWs + aVRWx' * aVRWx );
                                                
                    end

%                 group_p(count) = sum(abs(invmod) <= p_tresh);
%                 group_q(count) = sum(abs(Wx*invmod) <= q_tresh);
                
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
                
                
%                 % Measure the change in model update
%                 if count==1 
% 
%                     rdm(count) =  1;
%                     gdm(1) = norm(m_in - invmod);
% 
%                 else
% 
%                     gdm(2) = norm(m_in - invmod);
%                     rdm(count) = abs( gdm(2) - gdm(1) ) / norm(invmod);
% 
%                     gdm(1) = gdm(2);
% 
%                 end
                
               
             phi_d(count) = sum((G*(invmod)-d).^2);
             fprintf('%8.3e\n',phi_d(count));   
                
            % Check to see if overshooted the target misfit,               
            if phi_d(count) < target * (1-tol) || phi_d(count) > target * (1+tol)%&& switcher~=1
  
                beta_count = beta_count+1;
                
                if switcher == 0 && phi_d(count) < target * (1-tol)
                    
                    % Come back by a fraction and start IRLS
                    beta(count+1) = beta(count) * target / phi_d(count);
                    switcher = 1;                    

                elseif switcher == 0 
                    
                    beta(count+1) = 0.9 * beta(count);
                    
                else
                    
%                    if ll == 3 && switcher==2%Scenario #3, adjust beta and iterate 
%                         
%                         beta(count+1) = beta(count) * target / phi_d(count);
%                         
                   if (ll == 2 && lp_count~=1)  %Scenario #1, let it be
                       
                       beta(count+1) = beta(count);
                                              
                                          
                   else

                        beta(count+1) =beta(count) * target / phi_d(count);

                       
                       if ll == 3 || (ll == 2 && lp_count==1)
                        beta(count)=beta(count+1);
                        invmod = m_in;
                        count = count -1;
                        lp_count = lp_count-1;
%                         ncg(counter) = ncg(counter) - cg_iter;
                        continue
                       end
                        
                   end
                   
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
           
            phi_m(count) = (invmod)'*(MOF)*(invmod)/gamma; 
            phi_lp(count) = sum((Wx*invmod)>0);
            phi_p(count) = sum(aVRWs*invmod).^2;
            phi_x(count) = sum(aVRWx*invmod).^2;
            phi(count) = objfunc(invmod,MOF,beta(count));

            if count > 1
                dphi_m(count) = abs(phi_m(count) - phi_m(count-1))/phi_m(count) *100;
                dphi_lp(count) = abs(phi_lp(count) - phi_lp(count-1))/phi_lp(count) *100;
            end

            if lp_count > 0  && dphi_m(end) < 2 && ((switcher==2 && phi_d(count) > target *(1-tol) && phi_d(count) < target *(1+tol)) || ll == 2)

                switcher = 3;

                continue

            end    
           
            
            % OPTIONAL: Plot model
%             set(figure(4), 'Position', [50 200 1000 500])
%             axes(plot_h(counter));
%             plot(1:mcell,model); hold on
%             plot(1:mcell, invmod,'r','LineWidth',3); hold off
%             axis([0 mcell -.1 0.6]);
%             axis square

            
            
        end
                    % OPTIONAL: Plot model
            if ll == 1
                figure(2)
                axes('Position',[0.1 .55 .38 .38]);
                set(gca, 'YAxisLocation', 'left')
                ylabel('$m$','interpreter', 'latex','FontSize',16);
                set(get(gca,'YLabel'),'Rotation',360);
%                 text(0.5,0.68,['$\phi_m=\gamma^{(k)} \Big [  \|$ \boldmath$ W_s  \; R_s \; m \|_2^2$ $+ \|$ \boldmath$ W_x \;  R_x \; G_x m \|_2^2 \Big ]$'],'interpreter', 'latex','FontSize',11,'HorizontalAlignment','center')
                hold on
                text(0.45,-.25,'$(a)$','interpreter', 'latex','FontSize',16)
                text(0.5,.65,['$\mathbf{L_2-norm}$'],'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
                
            elseif ll==2
                figure(2)
                axes('Position',[0.55 .55 .38 .38]);
%                 text(0.5,0.68,['$\phi_m=\gamma^{(k)} \Big [  \|$ \boldmath$ W_s  \; \hat R_s \; m \|_2^2$ $+ \|$ \boldmath$ W_x \;  \hat R_x \; G_x m \|_2^2 \Big ]$'],'interpreter', 'latex','FontSize',11,'HorizontalAlignment','center')
                text(0.45,-.25,'$(b)$','interpreter', 'latex','FontSize',16)
                hold on
                text(0.5,.65,['$\mathbf{IRLS: \; Fix\;\{\epsilon_q,\;\beta\}}$'],'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
                set(gca,'YTickLabel',[]);
                
            elseif ll==3
                figure(3)
                axes('Position',[0.1 .55 .38 .38]);
                set(gca, 'YAxisLocation', 'left')
                ylabel('$m$','interpreter', 'latex','FontSize',16);
                set(get(gca,'YLabel'),'Rotation',360);
%                 text(0.5,0.68,['$\phi_m=\gamma^{(k)} \Big [  \|$ \boldmath$ W_s  \; R_s \; m \|_2^2$ $+ \|$ \boldmath$ W_x \;  R_x \; G_x m \|_2^2 \Big ]$'],'interpreter', 'latex','FontSize',11,'HorizontalAlignment','center')
                hold on
                text(0.45,-.25,'$(a)$','interpreter', 'latex','FontSize',16)
                text(0.5,.65,['$\mathbf{IRLS:\; Fix\;\{\epsilon_q\},\; \beta-search}$'],'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
                
            elseif ll == 4 
                figure(3)
                axes('Position',[0.55 .55 .38 .38]);
%                 text(0.5,0.68,['$\phi_m=\gamma^{(k)} \Big [  \|$ \boldmath$ W_s  \; \hat R_s \; m \|_2^2$ $+ \|$ \boldmath$ W_x \;  \hat R_x \; G_x m \|_2^2 \Big ]$'],'interpreter', 'latex','FontSize',11,'HorizontalAlignment','center')
                text(0.45,-.25,'$(b)$','interpreter', 'latex','FontSize',16)
                hold on
                text(0.5,.65,['$\mathbf{IRLS: \; Fix\;\{\epsilon_q\},\; \beta^{(k)},\; \gamma^{(k)}}$'],'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
                set(gca,'YTickLabel',[]);
                
            elseif ll == 5
                
                figure(4)
                axes('Position',[0.1 .55 .38 .38]);
                set(gca, 'YAxisLocation', 'left')
                ylabel('$m$','interpreter', 'latex','FontSize',16);
                set(get(gca,'YLabel'),'Rotation',360);
%                 text(0.5,0.68,['$\phi_m=\gamma^{(k)} \Big [  \|$ \boldmath$ W_s  \; R_s \; m \|_2^2$ $+ \|$ \boldmath$ W_x \;  R_x \; G_x m \|_2^2 \Big ]$'],'interpreter', 'latex','FontSize',11,'HorizontalAlignment','center')
                hold on
                text(0.45,-.25,'$(a)$','interpreter', 'latex','FontSize',16)
                text(0.5,.65,['$\mathbf{IRLS: \; Cool\;\{\epsilon_q^{(k)}\rightarrow 0\},\; \beta^{(k)},\; \gamma^{(k)}}$'],'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
                
            else
                
                 figure(4)
                axes('Position',[0.55 .55 .38 .38]);
%                 text(0.5,0.68,['$\phi_m=\gamma^{(k)} \Big [  \|$ \boldmath$ W_s  \; \hat R_s \; m \|_2^2$ $+ \|$ \boldmath$ W_x \;  \hat R_x \; G_x m \|_2^2 \Big ]$'],'interpreter', 'latex','FontSize',11,'HorizontalAlignment','center')
                text(0.45,-.25,'$(b)$','interpreter', 'latex','FontSize',16)
                hold on
                text(0.5,.65,['$\mathbf{IRLS: \; Cool\;\{\epsilon_q^{(k)}\rightarrow \epsilon^*\},\; \beta^{(k)},\; \gamma^{(k)}}$'],'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
                set(gca,'YTickLabel',[]);
                
                
            end
            
            
            plot(z,model,'b'); hold on
            plot(z, invmod,'r','LineWidth',2); 
            
            axis([z(1) z(end) -.1 0.6]);
            xlabel('$z$','interpreter', 'latex','FontSize',16);
                set(get(gca,'YLabel'),'Rotation',360);
            axis square
            text(0.10,0.55,['\underline{q=' num2str(LP(ll,2)) ', $\epsilon^{(k)}=\;' num2str(round(delta_q(end)/10^floor(log10(delta_q(end))))) 'e' num2str(floor(log10(delta_q(end)))) '$}'],'interpreter', 'latex','FontSize',12)
            text(0.10,0.50,['$\mathbf{\phi_d^{(k)}} = ' num2str(round(phi_d(count)*100)/100) '$'],'interpreter', 'latex','FontSize',10)
            text(0.10,0.40,['$\|\delta m\|_1$ = ' num2str(round(norm(model-invmod,1)*100)/100) ],'interpreter', 'latex','FontSize',10)
            text(0.10,0.45,['$\phi_m^{(k)}$ = ' num2str(round(phi_m(count)*100)/100) ],'interpreter', 'latex','FontSize',10)

            text(0.90,0.55,['\underline{Iterations}'],'interpreter', 'latex','FontSize',10,'HorizontalAlignment','right')
            text(0.90,0.50,['$\beta = ' num2str(beta_count) '$'],'interpreter', 'latex','FontSize',10,'HorizontalAlignment','right')
            text(0.90,0.45,['$CG$ = ' num2str(ncg(counter)) ],'interpreter', 'latex','FontSize',10,'HorizontalAlignment','right')
            

            hold off
            grid on

%             set(gca,'YTickLabel',[]);
            axis square
%             tlabel=(['$\mathbf{\phi_m\;=\;\| m \|_{' num2str(pQ(pp)) '}\;+\; \| \nabla m \|_{' num2str(qQ(qq)) '}}$']);
%             title(tlabel, 'interpreter', 'latex','FontSize',18)
            
            
            if ll == 1
                
                figure(2)
                axes('Position',[0.1 .075 .38 .38]);
                
            elseif ll == 2
                
                figure(2)
                axes('Position',[0.55 .075 .38 .38]);
                
            elseif ll == 3
                figure(3)
                axes('Position',[0.1 .075 .38 .38]);
                
            elseif ll ==4
                
                figure(3)
                axes('Position',[0.55 .075 .38 .38]);
                
            elseif ll == 5
                
                figure(4)
                axes('Position',[0.1 .075 .38 .38]);
                
            else 
                
                figure(4)
                axes('Position',[0.55 .075 .38 .38]);
                
            end
%             axes(plot_h((ll)+3));           
            [h,line1,line2] = plotyy(1:count,phi_d,1:count,phi_m,'semilogy'); hold on
            set(line1,'LineWidth',1.5);
            set(line2,'LineWidth',1.5);
            
%             plot(1:count,phi_x,'g--','LineWidth',1.5)

                        
%             plot(1:count,phim_p,'c-.','LineWidth',1.5)
%             
%             plot(1:count,phim_Ekb_p,'co','MarkerSize',5)
%             plot(1:count,phim_Ekb_q,'gx','MarkerSize',5)
            
%             plot(1:count,grad_s,'k-.','LineWidth',1.5)
%             plot(1:count,grad_x,'k--','LineWidth',1.5)
            
            grid on
            xlim(h(1),[0 count])
            xlim(h(2),[0 count])
            ylim(h(1),[1 5e+4])
            ylim(h(2),[1 5e+4])
            
            axis(h(1),'square')
            axis(h(2),'square')
%             set(gca, 'YAxisLocation', 'right')
                
            xlabel('IRLS Steps')
            
            if ll ==1 || ll ==3
            %             ylabel(h(1),'$\phi_d$' ,'interpreter', 'latex','FontSize',14);
            %             set(get(h(1),'YLabel'),'Rotation',360);
            %             set(ylabh,'Position',get(ylabh,'Position') - [4 0 0.00]);
            text(-2,5e+4,'$\phi_d,\;\phi_m$' ,'interpreter', 'latex','FontSize',14,'Color','k','HorizontalAlignment','left','VerticalAlignment','bottom');
            else

            set(h(1),'YTickLabel',[]);
            end
%             if ll ==3
%             %             ylabh = get(h(2),'YLabel');
%             %             set(ylabh,'Position',get(ylabh,'Position') + [3 0 0.00]);
%             %             ylabel(h(2),'$\phi_m$' ,'interpreter', 'latex','FontSize',14);
%             %             set(get(h(2),'YLabel'),'Rotation',360);
%             text(count+2,1e+4,'$\phi_m$' ,'interpreter', 'latex','FontSize',14,'Color','r','VerticalAlignment','bottom');
%             set(h(2),'YTickLabel',[])
%             else

            set(h(2),'YTickLabel',[]);
%             end          

%             n_iter = count - lp_count - beta_count;
            xx = [0 0 l2_count l2_count 0];
            yy = [0.1 5e+4 5e+4 0.1 0.1];
            h = fill(xx,yy,'b');
            set(h,'FaceAlpha',0.1,'LineStyle','none');
            
            
%             n_iter = n_iter + lp_count;
            xx = [l2_count l2_count count count l2_count];
            yy = [0.1 5e+4 5e+4 0.1 0.1];
            h = fill(xx,yy,'r');
            set(h,'FaceAlpha',0.1,'LineStyle','none');
            
            if ll ~=1
            rx = 1./ ( abs(Wx * model) .^2 + delta_q(count).^2 ) .^( (1-LP(ll,2)/2)/2 );
            Rx = spdiags( rx ,0,mcell,mcell);
            aVRWx =   sqrt((2.0 - LP(ll,3)) * ax )* Wr * Vx * Rx * Wx;
            phim_true = model'*aVRWx'*aVRWx*model;
            
            plot([l2_count count],[phim_true phim_true],'r--');
            text(l2_count,phim_true,'$\phi_{m^*}$','interpreter', 'latex','FontSize',12,'Color','r','BackgroundColor','w','EdgeColor','r','HorizontalAlignment','right')
            end
            
            plot([0 count],[target target],'b--');
            text(20,target,'$\phi_d^*$','interpreter', 'latex','FontSize',12,'Color','b','BackgroundColor','w','EdgeColor','b','HorizontalAlignment','center','VerticalAlignment','middle')
%             text(2,phim_q(1),'$\phi_x^{(k)}$','interpreter', 'latex','FontSize',12,'Color','g','BackgroundColor','w','EdgeColor','g','HorizontalAlignment','left','VerticalAlignment','bottom')
%             text(2,phim_p(1),'$\phi_s^{(k)}$','interpreter', 'latex','FontSize',12,'Color','r','BackgroundColor','w','EdgeColor','r','HorizontalAlignment','left','VerticalAlignment','bottom')
            text(2,phi_d(1),'$\phi_d^{(k)}$','interpreter', 'latex','FontSize',12,'Color','b','BackgroundColor','w','EdgeColor','b','HorizontalAlignment','left','VerticalAlignment','bottom')
            text(2,phi_m(1),'$\phi_m^{(k)}$','interpreter', 'latex','FontSize',12,'Color','r','BackgroundColor','w','EdgeColor','r','HorizontalAlignment','left','VerticalAlignment','top')

            % Plot phi_m*
%             rx = 1./ ( abs(Wx * model) .^2 + delta_q(count).^2 ) .^( (1-qvec/2)/2 );
% %             eta_x = mcell/(rx'*rx);
%             Rx = spdiags( sqrt(eta_x) * rx ,0,mcell,mcell);
%             aVRWx =  spdiags( sqrt((2.0 - lvec) * ax ),0,mcell,mcell) * Wr * Vx * Rx * Wx;
%             
%             rs = 1./ ( abs(model) .^2 + delta_p(count).^2 ) .^( (1-pvec/2)/2 );
% %             eta_s = mcell/(rs'*rs);
%             Rs = spdiags( sqrt(eta_s) * rs ,0,mcell,mcell);
%             aVRWs =  spdiags( sqrt((lvec) * as ),0,mcell,mcell) * Wr * V * Rs * Ws;
%             
%             phim_true = model'*(aVRWx'*aVRWx + aVRWs'*aVRWs) * model;
%             
%             plot([0 count],[phim_true phim_true],'r--');
%             text(2,phim_true,'$\phi_m^*$','interpreter', 'latex','FontSize',12,'Color','r','BackgroundColor','w','EdgeColor','r','HorizontalAlignment','left','VerticalAlignment','bottom')
            
            
            if ll == 2 || ll ==4
                
               
                   
                   axes('Position',[.38 .01 .28 .03])

                axis([0 3 0 1])

                h = fill([0 0 1.5 1.5 0],[0 1 1 0 0],'b');set(h,'FaceAlpha',0.1); hold on
                h = fill([1.5 1.5 3 3 1.5],[0 1 1 0 0],'r');set(h,'FaceAlpha',0.1);
                set(gca,'XTickLabel',[])
                set(gca,'YTickLabel',[])
%                 set(gca,'Visible','off');  
                
                text(0.75,0.4,'$l_2-norm$','interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
                text(2.25,0.4,'$l_p-norm$','interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
            end
%             figure; plotyy(1:count,phi_d,1:count,phi_m);

            ldml(counter) = norm(model-invmod,2);
%             fprintf('Relative dm:\t %8.5e \n', rdm(count));
%             fprintf('End of lp inversion. Number of iterations: %i\n',count)
%             fprintf('Final model misfit: %8.3e\n\n',ldml(counter))
            
            counter = counter+1;
            

end

