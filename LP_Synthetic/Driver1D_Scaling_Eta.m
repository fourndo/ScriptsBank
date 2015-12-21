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
nk = 30; %number of frequencies

decay = -0.20;

basis = 0.20;

% Data noise
amp_pct = 0.05;
floor_pct = 0.02;

% multip = 100;

% Percentile for cutoff
pct_cutoff = 75;
eps_p = (1e-3);
eps_q = (1e-4);

tol = 0.05
%% SCRIPT STARTS HERE
%% % Generate kernel functions and depth weighting
z = (0 : 1/(nx) : 1);
mcell=nx;
dx=ones(1,mcell) * abs((min(z)-max(z))/mcell);


set(figure, 'Position', [50 0 775 1000]);
% plot_h = tight_subplot(2,2,[.05 .01],[.1 0.1]);

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


%% Create inversion parameters
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
set(figure, 'Position', [50 0 775 400]);
% plot_h = tight_subplot(2,3,[.01 .01],[.25 0.05]);
for ll= 1:3           

    %% Create zones and transition for lp,lq,lu

% Create smoothing operator for transitions
av = @(n) spdiags (kron(ones(n,1),[0.25,0.5,0.25]),[-1,0,1],n,n);

A = av(mcell); avcx(1,1:2) = 0.5; avcx(end,end-1:end) = 0.5;
% avcz = av(nz); avcz(1,1:2) = 0.5; avcz(end,end-1:end) = 0.5;

% Avcx = kron ( avcx, speye (nz));
% Avcz = kron ( speye (nx), avcz);

% A = (Avcx * Avcz);
A = A^3;
A = spdiags(1./sum(A,2),0,mcell,mcell) *A;


if ll~=3
LP = [0 2 1];  
% Define zones {x,y,lp,lqx,lqz,ll}
zone{1,1} = 1:mcell; 
% zone{2,1} = 101:200; 
else
LP = [0 0 1;
      1 2 1];
% Define zones {x,y,lp,lqx,lqz,ll}
zone{1,1} = 1:100; 
zone{2,1} = 101:200;  
end    
    
% Build background tile
s = zeros(mcell,size(zone,1)+1);
s(:,end) = ones(mcell,1);

for ii = 1 : size(zone,1)
       
    temp = zeros(mcell,1); 
    temp(zone{ii,1}) = 1;
    s(:,ii) = temp(:);
    
    % Zero out tiles from the background
    s(s(:,ii)==1,end) = 0;
    
end
%Remove background zone if empty
if sum(s(:,end)) ==0
    
    s = s(:,1:end-1);
    
end
% Build transition matrices and lp,lq,ll vector
t = A * s;

% set(figure, 'Position', [50 0 775 1000]);
% lvec = zeros(mcell,1);
% qvec = zeros(mcell,1);
% pvec = zeros(mcell,1);

% for jj = 1 : size(s,2)
%     
%     lvec =  lvec + t(:,jj)*LP(jj,3);
%     qvec =  qvec + t(:,jj)*LP(jj,2);
%     pvec =  pvec + t(:,jj)*LP(jj,1);
%     
% end

            
 %%         % START - Message prompt
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
            phi = [];
            phim_p = [];
            phim_q = [];
            
            phim_Ekb_p = [];
            phim_Ekb_q = [];
            
            
            % Paremeters for epsilon cooling
            traffic_p = 1;
            traffic_q = 1;
            group_p = 0;
            p_tresh = 1;
            q_tresh = 1;
            grad_s = [];
            grad_x = [];
            delta_p = [];
            delta_q = [];
                        
            WxtWx = ( Wr * Vx * Wx)' * ( Wr * Vx * Wx );
            WstWs = ( Wr * V * Ws)' * ( Wr * V * Ws);

%             lvec(z{1}==1) = 0;
%             lvec(z{2}==1) = 2;
           
            traffic_p = 1;
            traffic_q = 1;
            group_p = mcell;
            l2_count = 0;
        while switcher ~= 3

                    count=count+1;

                    if switcher == 0   %First iteration

                        l2_count = l2_count+1;
                        delta_p(count) = 1e-1;%prctile(abs(invmod(invmod ~= 0)),pct_cutoff);
                        gradm = abs(Wx * invmod);
                        delta_q(count) = 1e-1;%prctile(abs(gradm(gradm ~= 0)),pct_cutoff);
%                         delta_q(count) = max(abs(invmod));                        
                        % Initial beta trace(G'G) / trace(phim)

%                         scale_p = 1 /max(abs( as * WstWs * invmod ));
%                         
%                         scale_x = 1/max(abs( ax * WxtWx * invmod ));
                        
%                         avrws =  Ws;
%                         avrwx =  Wx;                         
% 
%                         grad_s(count) = max(abs( (avrws)'*(avrws) * invmod ));                        
%                         grad_x(count) = max(abs( (avrwx)'*(avrwx) * invmod ));
                        
                        aVRWs = sqrt(as) * Wr * V * Ws;
                        aVRWx = sqrt(ax) * Wr * Vx * Wx;
                        
                        MOF = aVRWs'*aVRWs + aVRWx'*aVRWx;
                        
                        if count==1
                            beta = sum(sum(G.^2,1)) / sum(diag(MOF,0).^2) *1e+1;
                        end
                        
%                         lambda(count) = beta(count);
                        
                        phi(count) = objfunc(invmod,MOF,beta(count));
                        gamma = 1;
%                         x_tresh = dkdt(delta(count));
                        lvec = ones(mcell,1);
                        
                        p_tresh = 0.01;
                        q_tresh = 0.01;
                        
                    else
                        
%                         fprintf('# # LP-LQ ITER# #');
                        fprintf('# # LP-LQ ITER# #\n');
%                         if switcher ~= 2
                            lp_count = lp_count+1;
                                                     
%                         end
                        
                        
                        if lp_count == 1 
                                
                            
%                             [eps_p,eps_q] = get_eps(invmod,5,Wx,[],[]);
                            [eps_p,eps_q] = get_eps_v2(invmod,LP(1,1),LP(1,2),Wx,[],[]);
%                             [~,~] = get_eps_Lcurve(invmod,10,Wx,[],[]);
                                                    
%                             eps_p = 1e-4;
%                             eps_q = 1e-4;
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
                        if dphi_m(end) < 1  && beta_count > 1 && delta_q(end)<= eps_q && delta_p(end)<= eps_p%&& traffic_p(end)*100 < 2 && traffic_q(end)*100 < 2 
                            
%                             delta_q(count) = eps_q;
%                             delta_p(count) = eps_p;
                            switcher = 2;
                            
                        end
                            
       
%                         x_tresh = dkdt(pvec,delta(count));
                       
                       
                        
                        tRs = zeros(mcell,1);
                        tRx = zeros(mcell,1);
                        
%                         scale_p = zeros(mcell,1);
%                         scale_x = zeros(mcell,1);
                        
                        lvec = zeros(mcell,1);
                        for ii = 1 : size(s,2)
                            
                            T = spdiags(t(:,ii),0,mcell,mcell);
                            S = spdiags(s(:,ii),0,mcell,mcell);
                        
                         rs = sqrt(1./ ( abs(invmod) .^2 + delta_p(count).^2 ) .^( (1-LP(ii,1)/2) ));
                        rx = sqrt(1./ ( abs(Wx * invmod) .^2 + delta_q(count).^2 ) .^( (1-LP(ii,2)/2) ));
%                          rs = sqrt(1./ (( abs(invmod).^ (2-LP(ii,1))  + delta_p(count) )  ));
%                         rx = sqrt(1./ (( abs(Wx * invmod).^ (2-LP(ii,2))  + delta_q(count) )   ));                        
                        Rs = spdiags( rs ,0,mcell,mcell);
                        Rx = spdiags( rx ,0,mcell,mcell);
                        
                            avrws =  Rs ;
                            avrwx =  Rx * Wx;
                            
                            grad_s =(avrws)'*(avrws) * invmod ;                        
                            grad_x =(avrwx)'*(avrwx) * invmod ;
                            
                            if ll == 1
                                
                                eta_s = 1;
                                eta_x = 1;
                                
                            else
                                % WORKING
%                                 eta_s = 1./max(abs(grad_s));
%                                 eta_x = 1./max(abs(grad_x));
%                                     eta_s = 1/norm(grad_s);
%                                     eta_x = 1/norm(grad_x);
%                                     eta_s = sqrt(mcell)/norm(rs);
%                                     eta_x = sqrt(mcell)/norm(rx);
                                    
                                % NOT WORKING
%                                     eta_s = 1/max(rs);
%                                     eta_x = 1/max(rx);
                                    eta_s = delta_p(count).^( (1-LP(ii,1)/2) );
                                    eta_x = delta_q(count).^( (1-LP(ii,2)/2) );
%                                     eta_s = mcell/sum(rs.^2);
%                                     eta_x = mcell/sum(rx.^2);
%                                     eta_s =norm(invmod)/norm(grad_s);
%                                     eta_x =norm(Wx'*Wx*invmod)/norm(grad_x);


                            end
                            
                            tRs = tRs + sqrt( eta_s )  * T * rs;
                            tRx = tRx + sqrt( eta_x )  * T * rx;
 
                            lvec = lvec + T * ones(mcell,1) * LP(ii,3);
                        end
                        
                        temp=sum(tRs)+sum(tRx);
                        tRs = spdiags( tRs , 0 , mcell , mcell);
                        tRx = spdiags( tRx , 0 , mcell , mcell);
                        
                        aVRWs =  spdiags( sqrt( as * lvec  ),0,mcell,mcell)  * Wr * tRs * V * Ws ;
                        aVRWx =  spdiags( sqrt((2.0 - lvec) * ax  ),0,mcell,mcell) * Wr * tRx * Vx * Wx;
%                          %                          figure(100);plot((aVRWx)'*(aVRWx) * invmod );
                           
                        gamma = phi_m(count-1) * gamma /...
                            ( invmod'*(aVRWs' * aVRWs + aVRWx' * aVRWx) * invmod );
%                                          gamma(count)=2*mcell/temp;                                                        
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
                
                
                phi_d(count) = sum((G*(invmod)-d).^2);
%              fprintf('%8.3e\n',phi_d(count));   
                
            % Check to see if overshooted the target misfit,               
            if phi_d(count) < target * (1-tol) || phi_d(count) > target * (1+tol)%&& switcher~=1
  
                if switcher == 0 && phi_d(count) < target * (1-tol)
                    
                    % Come back by a fraction and start IRLS
                    beta(count+1) = beta(count) * target / phi_d(count);
                    switcher = 1;                    

                elseif switcher == 0 
                    
                    beta(count+1) = 0.9 * beta(count);
                    
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
           
            phi_m(count) = (invmod)'*(MOF)*(invmod)/gamma; 
            phi_lp(count) = sum((Wx*invmod)>0);
            phi_p(count) = sum(aVRWs*invmod).^2;
            phi_x(count) = sum(aVRWx*invmod).^2;
            phi(count) = objfunc(invmod,MOF,beta(count));

            if count > 1
                dphi_m(count) = abs(phi_m(count) - phi_m(count-1))/phi_m(count) *100;
                dphi_lp(count) = abs(phi_lp(count) - phi_lp(count-1))/phi_lp(count) *100;
            end

            if lp_count > 0  && dphi_m(end) < 5 && switcher==2 && phi_d(count) > target *(1-tol) && phi_d(count) < target *(1+tol)

                switcher = 3;

                continue

            end  
            
        end
            
            if ll == 1
                figure(2)
                axes('Position',[0.1 .55 .38 .38]);
                set(gca, 'YAxisLocation', 'left')
                ylabel('$m$','interpreter', 'latex','FontSize',16);
                set(get(gca,'YLabel'),'Rotation',360);
                text(0.5,0.68,['$\phi_m=\gamma^{(k)} \Big [  \|$ \boldmath$ W_s  \; R_s \; m \|_2^2$ $+ \|$ \boldmath$ W_x \;  R_x \; G_x m \|_2^2 \Big ]$'],'interpreter', 'latex','FontSize',11,'HorizontalAlignment','center')
                hold on
                text(0.45,-.25,'$(a)$','interpreter', 'latex','FontSize',16)
                text(0.35,-0.05,['$\mathbf{p:' num2str(LP(1,1)) ' ,\; q:' num2str(LP(1,2)) '}$'],'interpreter', 'latex','FontSize',12)
                
            elseif ll==2
                figure(2)
                axes('Position',[0.55 .55 .38 .38]);
                text(0.5,0.68,['$\phi_m=\gamma^{(k)} \Big [  \|$ \boldmath$ W_s  \; \hat R_s \; m \|_2^2$ $+ \|$ \boldmath$ W_x \;  \hat R_x \; G_x m \|_2^2 \Big ]$'],'interpreter', 'latex','FontSize',11,'HorizontalAlignment','center')
                text(0.45,-.25,'$(b)$','interpreter', 'latex','FontSize',16)
                hold on
                text(0.35,-0.05,['$\mathbf{p:' num2str(LP(1,1)) ',\; q:' num2str(LP(1,2)) '}$'],'interpreter', 'latex','FontSize',12)
                set(gca,'YTickLabel',[]);
                
            else
                figure(3)
                axes('Position',[-.1 .2 .75 .75]);
%                 text(-.3,0.7,['$\phi_m\;=\mu^{(k)} \Big [ \alpha_p \|$ \boldmath$\eta_p \hat W_p \; R_p \; m \|_2^2$ $+ \alpha_q\|$ \boldmath$\eta_q \hat W_q \; R_q \; G m \|_2^2 \Big ]$'],'interpreter', 'latex','FontSize',12)
                set(gca, 'YAxisLocation', 'left')
                ylabel('$m$','interpreter', 'latex','FontSize',16);

                hold on
                text(0.1,-0.05,['$\mathbf{p:' num2str(LP(1,1)) ' ,\; q:' num2str(LP(1,2)) '}$'],'interpreter', 'latex','FontSize',12)
                text(0.6,-0.05,['$\mathbf{p:' num2str(LP(2,1)) ' ,\; q:' num2str(LP(2,2)) '}$'],'interpreter', 'latex','FontSize',12)
                hold on
                xx = [0 0 0.5 0.5 0];
                yy = [-.1 0.6 0.6 -.1 -.1];
                h = fill(xx,yy,[0.25 0.25 0.25]);
                set(h,'FaceAlpha',0.1,'LineStyle','none');
                
                xx = [0.5 0.5 1 1 0.5];
                yy = [-.1 0.6 0.6 -.1 -.1];
                h = fill(xx,yy,'w');
                set(h,'FaceAlpha',0.1,'LineStyle','none');
                text(0.9,0.55,'$(a)$' ,'interpreter', 'latex','FontSize',14,'Color','k','HorizontalAlignment','center','VerticalAlignment','middle')
            end
            
            
            plot(z,model,'b'); hold on
            plot(z, invmod,'r','LineWidth',2); 
            
            axis([z(1) z(end) -.1 0.6]);
            xlabel('$x$','interpreter', 'latex','FontSize',16);
                set(get(gca,'YLabel'),'Rotation',360);
            axis square
            text(0.30,0.55,['\underline{Final Iteration}'],'interpreter', 'latex','FontSize',10)
            text(0.30,0.50,['$\mathbf{\phi_d} = ' num2str(round(phi_d(count)*100)/100) '$'],'interpreter', 'latex','FontSize',10)
            text(0.30,0.45,['$\epsilon = 10^{' num2str(round(log10(delta_p(end))*10)/10) '}$'],'interpreter', 'latex','FontSize',10)
            text(0.30,0.40,['$\|\delta m\|_1$ = ' num2str(round(norm(model-invmod,1)*100)/100) ],'interpreter', 'latex','FontSize',10)


            hold off
            grid on

%             set(gca,'YTickLabel',[]);
            axis square
%             tlabel=(['$\mathbf{\phi_m\;=\;\| m \|_{' num2str(pQ(pp)) '}\;+\; \| \nabla m \|_{' num2str(qQ(qq)) '}}$']);
%             title(tlabel, 'interpreter', 'latex','FontSize',18)
            
            n_iter = count - (lp_count+beta_count);
            
            if ll == 1
                figure(2)
                axes('Position',[0.1 .075 .38 .38]);
                
            elseif ll == 2
                
                figure(2)
                axes('Position',[0.55 .075 .38 .38]);
                
            else
               
                figure(3)
                axes('Position',[0.4 .2 .75 .75]);
                
            end
%             axes(plot_h((ll)+3));           
            [h,line1,line2] = plotyy(1:count,phi_d,1:count,phi_m,'semilogy'); hold on
            set(line1,'LineWidth',1.5);
            set(line2,'LineWidth',1.5);
            
%             plot(1:count,phim_q,'g--','LineWidth',1.5)

                        
%             plot(1:count,phim_p,'c-.','LineWidth',1.5)
%             
%             plot(1:count,phim_Ekb_p,'co','MarkerSize',5)
%             plot(1:count,phim_Ekb_q,'gx','MarkerSize',5)
            
%             plot(1:count,grad_s,'k-.','LineWidth',1.5)
%             plot(1:count,grad_x,'k--','LineWidth',1.5)
            
            grid on
            xlim(h(1),[0 count])
            xlim(h(2),[0 count])
            ylim(h(1),[0.1 5e+4])
            ylim(h(2),[0.1 5e+4])
            
            axis(h(1),'square')
            axis(h(2),'square')
%             set(gca, 'YAxisLocation', 'right')
                
            xlabel('Iteration')
            
            if ll ==1
            %             ylabel(h(1),'$\phi_d$' ,'interpreter', 'latex','FontSize',14);
            %             set(get(h(1),'YLabel'),'Rotation',360);
            %             set(ylabh,'Position',get(ylabh,'Position') - [4 0 0.00]);
            text(-2,1e+4,'$\phi_d$' ,'interpreter', 'latex','FontSize',14,'Color','b','HorizontalAlignment','right','VerticalAlignment','bottom');
            else

%             set(h(1),'YTickLabel',[]);
            end
            if ll ==3
            %             ylabh = get(h(2),'YLabel');
            %             set(ylabh,'Position',get(ylabh,'Position') + [3 0 0.00]);
            %             ylabel(h(2),'$\phi_m$' ,'interpreter', 'latex','FontSize',14);
            %             set(get(h(2),'YLabel'),'Rotation',360);
            text(0,2e+4,'$\phi_m,\; \phi_d$' ,'interpreter', 'latex','FontSize',14,'Color','k','VerticalAlignment','bottom','HorizontalAlignment','right');
            set(h(2),'YTickLabel',[])
            text(count-2,2e+4,'$(b)$' ,'interpreter', 'latex','FontSize',14,'Color','k','HorizontalAlignment','right')
            else

            set(h(2),'YTickLabel',[]);
            end          

            n_iter = count - lp_count;
            xx = [0 0 l2_count l2_count 0];
            yy = [0.1 5e+4 5e+4 0.1 0.1];
            h = fill(xx,yy,'b');
            set(h,'FaceAlpha',0.1,'LineStyle','none');
            
            
%             n_iter = n_iter + lp_count;
            xx = [l2_count l2_count count count l2_count];
            yy = [0.1 5e+4 5e+4 0.1 0.1];
            h = fill(xx,yy,'r');
            set(h,'FaceAlpha',0.1,'LineStyle','none');
            
            plot([0 count],[target target],'b--');
            text(2,target,'$\phi_d^*$','interpreter', 'latex','FontSize',12,'Color','b','BackgroundColor','w','EdgeColor','b','HorizontalAlignment','left','VerticalAlignment','top')
%             text(2,phim_q(1),'$\phi_x^{(k)}$','interpreter', 'latex','FontSize',12,'Color','g','BackgroundColor','w','EdgeColor','g','HorizontalAlignment','left','VerticalAlignment','bottom')
            text(2,phi_m(1),'$\phi_m^{(k)}$','interpreter', 'latex','FontSize',12,'Color','r','BackgroundColor','w','EdgeColor','r','HorizontalAlignment','left','VerticalAlignment','bottom')
            text(2,phi_d(1),'$\phi_d^{(k)}$','interpreter', 'latex','FontSize',12,'Color','b','BackgroundColor','w','EdgeColor','b','HorizontalAlignment','left','VerticalAlignment','bottom')

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
            
            
            if ll == 2 || ll ==3
                
               if ll == 2
                   
                   axes('Position',[.375 .01 .28 .03])
                   
               else
                   
                   axes('Position',[.68 .2 .2 .075])
                   
               end
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

%% Export table for true model
fprintf('\nWHOLE MODEL\n')
fprintf('Norm\tphi_lp\\t\t\tphi_lq\t\t\tmax(GRADp)\t\tmin(GRADq)\n')
for pp = 0:0.5:2
    
    rs = 1./ ( abs(model) .^2 + delta_p(count).^2 ) .^( (1-pp/2) );
    rx = 1./ ( abs(Wx * model) .^2 + delta_q(count).^2 ) .^( (1-pp/2) );

    Rs = spdiags( rs.^0.5 ,0,mcell,mcell);
    Rx = spdiags( rx.^0.5 ,0,mcell,mcell);
                        
    aVRWs =  spdiags( sqrt( as * lvec  ),0,mcell,mcell)  * Wr * Rs * V * Ws ;
    aVRWx =  spdiags( sqrt((2.0 - lvec) * ax  ),0,mcell,mcell) * Wr * Rx * Vx * Wx;
                    
    nrm_lp = model' * aVRWs' * aVRWs * model;
    nrm_lq = model' * aVRWx' * aVRWx * model;
    
    Ekbl_p = sum( abs((sqrt(as) * Wr * V * model).^2 + delta_p(count)^2).^pp/2);
    Ekbl_q = sum( abs((sqrt(ax) * Wr * Vx * Wx * model).^2 + delta_q(count)^2).^pp/2);
    
    grad_p = abs((aVRWs' * aVRWs) * model);
    grad_q = abs((aVRWx' * aVRWx) * model);

    fprintf('%4.2f\t%12.8e\t%12.8e\t%12.8e\t%12.8e\t\n',pp,nrm_lp,nrm_lq,max(grad_p),max(grad_q));
end

% fprintf('\nSQUARE ANOMALY\n')
% fprintf('Norm\tphi_p\t\t\tphi_q\n')
% for pp = 0:0.5:2
%     
%     rs = 1./ ( abs(model) .^2 + delta_p(count).^2 ) .^( (1-pp/2) );
%     rx = 1./ ( abs(Wx * model) .^2 + delta_q(count).^2 ) .^( (1-pp/2) );
% 
%     Rs = spdiags( rs.^0.5 ,0,mcell,mcell);
%     Rx = spdiags( rx.^0.5 ,0,mcell,mcell);
%              
%     temp = [zeros(100,1);ones(100,1)];
%     S = spdiags(temp,0,mcell,mcell);
%     aVRWs =  spdiags( sqrt( as * lvec  ),0,mcell,mcell)  * Wr * Rs * V * S * Ws ;
%     aVRWx =  spdiags( sqrt((2.0 - lvec) * ax  ),0,mcell,mcell) * Wr * Rx * Vx * S * Wx;
%                     
%     nrm_p = model' * aVRWs' * aVRWs * model;
%     nrm_q = model' * aVRWx' * aVRWx * model;
%     
%     fprintf('%4.2f\t%12.8e\t%12.8e\t\n',pp,nrm_p,nrm_q);
% end
% 
% fprintf('\nSMOOTH ANOMALY\n')
% fprintf('Norm\tphi_p\t\t\tphi_q\n')
% for pp = 0:0.5:2
%     
%     rs = 1./ ( abs(model) .^2 + delta_p(count).^2 ) .^( (1-pp/2) );
%     rx = 1./ ( abs(Wx * model) .^2 + delta_q(count).^2 ) .^( (1-pp/2) );
% 
%     Rs = spdiags( rs.^0.5 ,0,mcell,mcell);
%     Rx = spdiags( rx.^0.5 ,0,mcell,mcell);
%              
%     temp = [ones(100,1);zeros(100,1)];
%     S = spdiags(temp,0,mcell,mcell);
%     aVRWs =  spdiags( sqrt( as * lvec  ),0,mcell,mcell)  * Wr * Rs * V * S * Ws ;
%     aVRWx =  spdiags( sqrt((2.0 - lvec) * ax  ),0,mcell,mcell) * Wr * Rx * Vx * S * Wx;
%                     
%     nrm_p = model' * aVRWs' * aVRWs * model;
%     nrm_q = model' * aVRWx' * aVRWx * model;
%     
%     fprintf('%4.2f\t%12.8e\t%12.8e\t\n',pp,nrm_p,nrm_q);
% end