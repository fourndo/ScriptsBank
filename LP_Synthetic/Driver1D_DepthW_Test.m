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

% Lp-norm parameters
pQ = [2 0]%0:.1:2;%0%
qQ = [0 2]%0:.1:2;%0%
lQ = [0 1]%ones(1,1)*1.0%[0.25:0.25:1.75]%0.1:0.4:1.9;

% Add an extra test to find best starting beta for lplq
chi_vec = 1%[2 3 4 5 7.5 10 20 30];
% multip = 100;

% Percentile for cutoff
pct_cutoff = 50;
eps_p = 1e-3;
eps_q =1e-3;
%% SCRIPT STARTS HERE
%% % Generate kernel functions and depth weighting
x = (0 : 1/(nx) : 1);
mcell=nx;
dx=ones(1,mcell) * abs((min(x)-max(x))/mcell);


set(figure, 'Position', [50 0 775 1000]);
% plot_h = tight_subplot(2,2,[.05 .01],[.1 0.1]);

G = zeros(nk,mcell);
wr = zeros(1,mcell);
subplot(2,2,2)
for ii = 1 : nk 
    
    b = basis * 2*pi* (ii);
    a = decay * (ii);
    
    G(ii,1:mcell) = exp( a * x(2:end) ) .* (a/(a^2+b^2)*cos (b * x(2:end)) + b/(a^2+b^2) *sin (b * x(2:end))) -...
        exp( a * x(1:end-1) ) .* (a/(a^2+b^2)*cos (b * x(1:end-1)) + b/(a^2+b^2) *sin (b * x(1:end-1)));

    G(ii,1:mcell) = G(ii,1:mcell) / max(G(ii,1:mcell));
    
    wr = wr + ((exp( a*x(2:end))-exp( a * x(1:end-1) ))/a) ;
    plot( x(2:end) , G(ii,:) );hold on

end

% wr = (sum(G.^2,1)).^0.5;
wr = (wr / max(wr)).^0.5;

IWr = spdiags(1./wr',0,mcell,mcell);
Wr = spdiags(wr',0,mcell,mcell);
Wrx = spdiags(wr(1:end-1)',0,mcell-1,mcell-1);

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

cntr = round(nx/4);

% Create mdoel: Square and tanh functions
model = 0.4* exp(-(((x-3*cntr/nx) ).^2 * nx)); 

model(cntr-round(cntr/4):cntr+round(cntr/4))=0.25;

model=model(:);

subplot(2,2,1)
plot(x,model,'LineWidth',2);axis([x(1) x(end) -0.1 0.5]);hold on
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
Gin = Wd * G ;

target = nk;

% Weight data with uncertainty
d = Wd * d;

%% INVERSION
% Create inversion parameters
[Wx,Vx] = getWx1D(nx,dx);
[Ws] = speye(mcell);

Wx = Wx(1:end-1,:);
Vx = Vx(1:end-1,1:end-1);

V = spdiags(sqrt(dx'),0,nx,nx);
s= ones(mcell,1);
t= ones(mcell,1);

% Global constant
as = 1.0 / min(dx) ^2;  %Smallnest term
ax = 1.0;               %Smoothness term                     

% Store all the final models


Pac = speye(mcell);
            
dkdt = @(p,ep) ((ep).^(1./(4-2*p)) - ep);



objfunc = @(m,phi,b) sum( ( Gin * m - d ).^2 ) + ( m' * b * phi * m );
set(figure, 'Position', [50 0 775 500]);
axs{1} = axes('Position',[0-.025 .3 .65 .65]);
axs{2} = axes('Position',[0.44 .3 .65 .65]);
% plot_h = tight_subplot(1,2,[.1 .1],[.2 0.2]);
% Iterate over all the combination of norms as specified by the vectors
% pvec, qvec and lvec.
for rr = [2 1]
 counter = 1;  
if rr == 2

G = Gin*IWr;

else
    
    G = Gin;
    
end
for ll= 1:length(lQ)

            
            
            
            % Message prompt
            head = ['lp: ' num2str(pQ(ll)) ' lq: ' num2str(qQ(ll)) ' psi: ' num2str(lQ(ll))];
            fprintf('Starting lp inversion %s\n',...
                head)
            
            invmod      = ones(mcell,1)*1e-1;       % Initial model
            
            if rr == 2
            invmod      = Wr*invmod;       % Initial model 
            end
            
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
            
            % Paremeters for epsilon cooling
            traffic_p = 1;
            traffic_q = 1;
            group_p = 0;
            p_tresh = 1;
            q_tresh = 1;
            
            delta_p = [];
            delta_q = [];
                        
            % Define two zones with different lp,lq
            zone{1,1} = 1:75;  
            zone{2,1} = 125:200; 

            % Transition zone with linear interp
            zone{3,1} = 76:124;
            linterp = cumsum(ones(length(zone{3,1}),1)/length(zone{3,1}));
            
            % Z{1} = ones(mcell,1);
            z{1} = ones(mcell,1);%zeros(nx,1); z{1}(zone{1,1}) = 1;
            trans{1} = ones(mcell,1);%z{1}; trans{1}(zone{3,1}) = 1-linterp;
            
%             z{2} = zeros(nx,1); z{2}(zone{2,1}) = 1;
%             trans{2} = z{2}; trans{2}(zone{3,1}) = 1-linterp(end:-1:1);
                        
            pvec = ones(nx,1) * pQ(ll);%(0.*trans{1} + 2.*trans{2}) ./ (trans{1} + trans{2}) ;
            
            qvec = ones(size(Wx,1),1) * qQ(ll);%(0.*trans{1} + 2.*trans{2}) ./ (trans{1} + trans{2}) ;
            
            lvec = ones(nx,1) * lQ(ll);
%             lvec(z{1}==1) = 0;
%             lvec(z{2}==1) = 2;
           
            traffic_p = 1;
            traffic_q = 1;
            group_p = mcell;
            
        while switcher ~= 3

                    count=count+1;

                    if switcher == 0   %First iteration

                        delta_p(count) = prctile(abs(invmod(invmod ~= 0)),pct_cutoff)*5;
                        
                        if rr == 2
                            gradm = abs(Wx * IWr * invmod);
                        
                        else
                            gradm = abs(Wx * invmod);
                        end
                        
                        delta_q(count) = prctile(gradm(gradm~=0),pct_cutoff)*5;
%                         delta_q(count) = max(abs(invmod));                        
                        % Initial beta trace(G'G) / trace(phim)

%                         scale_p = 1 /max(abs( as * WstWs * invmod ));
%                         
%                         scale_x = 1/max(abs( ax * WxtWx * invmod ));
                        
                         
                        aVRWs = sqrt(as) *  V * Ws;
                        aVRWx = sqrt(ax) *  Vx * Wx;
                        
                        if rr==1
                            
                            aVRWs = Wr * aVRWs ;
                            aVRWx = Wrx * aVRWx ;
                        
                        end
                        
                        MOF = aVRWs'*aVRWs + aVRWx'*aVRWx;
                        
                        if count==1
                            beta = sum(sum(G.^2,1)) / sum(diag(MOF,0).^2) *1e+1;
                        end
                        
%                         lambda(count) = beta(count);
                        
                        phi(count) = objfunc(invmod,MOF,beta(count));
                        
%                         x_tresh = dkdt(delta(count));

                        p_tresh = 0.01;
                        q_tresh = 0.01;
                        
                    else
                        
%                         fprintf('# # LP-LQ ITER# #');
                        if switcher ~=2
                            
                            lp_count = lp_count+1;
                            
                        else
                            
                            beta_count = beta_count+1;
                            
                        end
                        
                        
                        if lp_count == 1
                            
                            p_tresh = eps_p;%prctile(abs(invmod(invmod ~= 0)),pct_cutoff);
%                             gradm = abs(Vx * Wx * invmod);
                            q_tresh = eps_q;%prctile(gradm(gradm~=0),pct_cutoff);

                        
                            delta_p(count) = 1e-1;%delta_p(count-1)/2;
                            delta_q(count) = 1e-1;%delta_q(count-1)/2;
                            %x_tresh.^2;
%                             figure; hist(abs(invmod),20);
%                             ylim([0 100])
%                             set(figure(100), 'Position', [50 200 750 750])
%                             plot(x,model,'LineWidth',2);axis([x(1) x(end) -0.1 0.5]);hold on
%                             axis square
%                             grid on
            
                        end
                        
                        if delta_p(end)> eps_p;%traffic_p(end)*100 > 1 && switcher == 1
                            
                            delta_p(count) = delta_p(end)/2;
                            
                        else 
                            
                            delta_p(count) = delta_p(count-1);
                            
                        end
                        
                        if delta_q(end)> eps_q;%traffic_q(end)*100 > 1 && switcher == 1
                            
                            delta_q(count) = delta_q(end)/2;
                            
                        else

                            delta_q(count) = delta_q(count-1);
                            
                        end
                        
                        if delta_p(end)< eps_p && delta_q(end)< eps_q%traffic_p(end)*100 < 1 && traffic_q(end)*100 < 1 && switcher == 1
                            
                            delta_q(count) = eps_q;
                            delta_p(count) = eps_p;
                            switcher = 2;
                            
                        end
                            
       
%                         x_tresh = dkdt(pvec,delta(count));

                        if rr == 2
                           
                            rs = 1./ ( abs(invmod) .^2 + delta_p(count).^2 ) .^( (1-pvec/2) );
                            rx = 1./ ( abs(Wx * invmod) .^2 + delta_q(count).^2 ) .^( (1-qvec/2) );
                            
                        else
                            
                            rs = 1./ ( abs(invmod) .^2 + delta_p(count).^2 ) .^( (1-pvec/2) );
                            rx = 1./ ( abs(Wx * invmod) .^2 + delta_q(count).^2 ) .^( (1-qvec/2) );
                            
                        end
                        
                        Rs = spdiags( rs.^0.5 ,0,mcell,mcell);
                        Rx = spdiags( rx.^0.5 ,0,size(Wx,1),size(Wx,1));
                       
                                                
%                         scale_p = zeros(mcell,1);
%                         scale_x = zeros(mcell,1);
                        
                            avrws =  Rs * Ws ;
                            avrwx =  Rx * Wx ;
                            
                                                                
                            eta_s = delta_p(count)^(1-pQ(ll)/2);%1 / max(abs( (avrws)'*(avrws) * invmod ));                        
                            eta_x = delta_q(count)^(1-qQ(ll)/2);%1 / max(abs( (avrwx)'*(avrwx) * invmod ));

                            
%                             theta_s = theta_s + sqrt( max_qs )  * T  ;
%                             theta_x = theta_x + sqrt( max_qx )  * T  ;
 
                                                                       
                        aVRWs =  spdiags( sqrt( lvec * as *eta_s ),0,mcell,mcell)  *  Rs * V * Ws ;
                        aVRWx =  spdiags( sqrt((2.0 - lvec) * ax *eta_x ),0,size(Wx,1),size(Wx,1)) *  Rx * Vx * Wx;
%                          figure(100);plot((aVRWx)'*(aVRWx) * invmod );
                        if rr ==1
                            
                            aVRWs = Wr * aVRWs ;
                            aVRWx = Wrx * aVRWx ;
                            
                        end
                        
                        gamma = phi_m(count-1) /...
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

                group_p(count) = sum(abs(invmod) <= p_tresh);
                group_q(count) = sum(abs(Vx*Wx*invmod) <= q_tresh);
                
%                 fprintf('\n# # # # # #\n');
%                 fprintf('BETA ITER: \t %i  \nbeta: \t %8.5e \n',count,beta(count));

                m_in = invmod;  % Save initial model for backtracking steps
                
%                 eigval_qtG = eig(G'*G); 
%                 eigval_MOF = eig( lambda(count) * (MOF) );
%                 
%                 eigval = eig( G'*G + lambda(count) * (MOF));
%                 
%                 condnum = max(abs(eigval)) / min(abs(eigval));
%                 
%                 figure(4); 
% %                 plot(eigval_qtG,'*'); hold on
% %                 plot(eigval_MOF,'ro'); hold on
%                 plot(eigval,'*'); hold off
%                 fprintf('Condition number: %f\n',condnum);

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
                
                if count > 1
                    temp = sum(abs(invmod) <= p_tresh);                     
                    traffic_p(count) = abs(group_p(count) - temp) / group_p(count);
                    
                    temp = sum(abs(Vx*Wx*invmod) <= q_tresh);
                    traffic_q(count) = abs(group_q(count) - temp) / group_q(count);
                    
                end
                
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
                phi_m(count) = (invmod)'*(MOF)*(invmod);
                phi_p(count) = sum(aVRWs*invmod).^2;
                phi_x(count) = sum(aVRWx*invmod).^2;
                phi(count) = objfunc(invmod,MOF,beta(count));

                if (rdm(count) < 5e-3 && switcher == 2 && phi_d(count) > target *0.99 && phi_d(count) < target *1.01) 

                    switcher = 3;
%                     fprintf('---------->')
%                     fprintf(' misfit:\t %8.5e ',phi_d(count))
%                     fprintf(' ** Relative dm:\t %8.5e ', rdm(count));
%                     fprintf('<----------\n')

                    continue

                end    
                
            % Check to see if overshooted the target misfit,
            % If yes, then compare with the previous iteration and keep
            % closest                
            if phi_d(count) < target * 0.99 && switcher~=1
                
                phi_d(count) = sum((G*(invmod)-d).^2);                       
                
%                 fprintf('---------->');
%                 fprintf(' misfit:\t %8.5e ',phi_d(count));
%                 fprintf('Relative dm:\t %8.5e ', rdm(count));
%                 fprintf('<----------\n');
%                 fprintf('\n# NEXT ITER: INCREASING BETA #\n');
                
                if switcher == 0
                    
                    % Come back almost a full beta step
                    beta(count+1) = beta(count) * 1.1;
                
                    switcher = 1;
                    
%                 elseif rdm(count) < 5e-3
%                     
%                     beta(count+1) = beta(count) * target / phi_d(count);
% 
%                     
                else
                   
                    % Slowly backtrack
%                     beta(count+1) = beta(count) * 1.05;
                    beta(count+1) = beta(count) * target / phi_d(count);
                end
                
                
                
%                 invmod = m_in;
%                 phi_d(count) = phi_d(count-1);
                
%                 count = count-1;
%                 continue
            
            % Else reduce beta and continue inversion
            elseif phi_d(count) > target * 1.01 && switcher ~= 1
                
%                 fprintf('---------->')
%                 fprintf(' misfit:\t %8.5e ',phi_d(count))
%                 fprintf(' ** Relative dm:\t %8.5e ', rdm(count));
%                 fprintf('<----------\n')
%                 fprintf('\n# NEXT ITER - REDUCING BETA #\n');
                if switcher == 0
                    
                    beta(count+1) = 0.9 * beta(count); 
                    
                else
                    
                    
%                     beta(count+1) = beta(count) * 0.95;
                    beta(count+1) = beta(count) * target / phi_d(count);
                    
                    
                end
                
            else
                
                
%                 fprintf('---------->')
%                 fprintf(' misfit:\t %8.5e ',phi_d(count))
%                 fprintf(' ** Relative dm:\t %8.5e ', rdm(count));
%                 fprintf('<----------\n')
%                 fprintf('\n# NEXT ITER - STATIONNARY STEP #\n');
                
                if switcher == 0
                    switcher = 1;
                end
                
                beta(count+1) = beta(count);

                 
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
            if rr == 2
                invmod = IWr*invmod;
            end
            
            model_error = norm(model-invmod,1);

            
            if rr == 1
                
            axes(axs{counter});
            plot(x,model); hold on
            plot(x, invmod,'k--','LineWidth',2); 
            
            else
            axes(axs{counter});
            plot(x, invmod,'r','LineWidth',2); hold on
            
            if ll == 1
                
                text(.35,0.43,['$\mathbf{q:' num2str(qQ(ll)) '}$'],'interpreter', 'latex','FontSize',12)

            else
                text(.35,0.43,['$\mathbf{p: ' num2str(pQ(ll)) ',\; q:' num2str(qQ(ll)) '}$'],'interpreter', 'latex','FontSize',12)
            end
            text(.35,0.39,['$\mathbf{\phi_d} = ' num2str(round(phi_d(count)*100)/100) '$'],'interpreter', 'latex','FontSize',10)
            text(.35,0.35,['$\epsilon = 10^{' num2str(round(log10(delta_p(end))*10)/10) '}$'],'interpreter', 'latex','FontSize',10)
            text(.35,0.31,['$\|\mathbf{m - m^*}\|$ = ' num2str(round(model_error*100)/100)],'interpreter', 'latex','FontSize',10)

            grid on
            axis([0 1 -.025 0.5]);
            
            if ll == 1
            ylabel('\bf{m}', 'interpreter', 'latex','FontSize',14);
            
            end
                xlabel('z', 'interpreter', 'latex','FontSize',14);
                set(get(gca,'YLabel'),'Rotation',360);
                ylabh = get(gca,'YLabel');
                set(ylabh,'Position',get(ylabh,'Position') - [0.1 .025 0.00]);
                
            if ll==1
                
                
                text(0.9,0.425,'(a)', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle');
            else
                
                set(gca,'YTickLabel',[]);
                text(0.9,0.425,'(b)', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle');

                
            
            end
            axis square
                
            end
            
            
%             figure; plotyy(1:count,phi_d,1:count,phi_m);

            ldml(counter) = norm(model-invmod,2);
%             fprintf('Relative dm:\t %8.5e \n', rdm(count));
%             fprintf('End of lp inversion. Number of iterations: %i\n',count)
%             fprintf('Final model misfit: %8.3e\n\n',ldml(counter))
            
            counter = counter+1;
            

        end

end

%%
axes('Position',[.6 .35 .01 .01])
plot(1,1,'b');hold on
plot(1,1,'r','LineWidth',2);hold on
plot(1,1,'k--','LineWidth',2);hold on
set(gca,'Visible','off');

h = legend('$m^*$','$\phi({\hat m})$',...
'$\phi({m})$','Location','North');
% h = legend('$m^*$','$\phi({\hat m}) = \| \;\mathbf{\hat F} \; \mathbf{ \hat m } - \mathbf{d}\|_2^2 + \| \;\mathbf{\hat W}_{ m} \; \mathbf{ \hat m }\|_2^2$',...
% '$\phi({m}) = \| \;\mathbf{F} \; \mathbf{m } - \mathbf{d}\|_2^2 + \| \mathbf{W}_{r}\;\mathbf{W}_{m} \; \mathbf{ m }\|_2^2$','Location','North');
set(h,'Interpreter', 'latex','FontSize',12)

% save('l1','l1');
% save('l2','l2');
% save('misfit','misfit');
% save('models_out','models_out');

%% Plot result 
% figure; hist(abs(invmod),200);text(0.04,90,['\epsilon = ' num2str(delta(end))])
% axis([0 0.1 0 100])
% xx = 1e-4:1e-4:max(abs(invmod));
% mm = xx./(xx.^(2-pQ(pp)) + delta(end));
% hold on; plot(xx,mm/max(mm)*100,'LineWidth',2,'Color','r')
% yy = get(gca,'Ylim');
% plot([x_tresh x_tresh],yy,'--')
        
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

% figure;
% [hAx,hLine1,hLine2] = plotyy(delta,ncg,delta,ldml,'semilogx','semilogx');
% ylabel(hAx(1),'\bfNumber of CG steps');
% ylabel(hAx(2),'\bf|m_{true} - m_{k}|_1');
% xlabel('\bf Factor (\epsilon)')
% set(hAx(2),'YColor','k')
% set(hLine2,'LineWidth',2,'Color','k')
% set(hLine1,'LineWidth',2)
% grid on

% [Ilow,Jlow]=find(l1==min(min(l1)));
% [Iupp,Jupp]=find(l1==max(max(l1)));
% 
% set(figure, 'Position', [50 0 775 1000]);
% plot_h = tight_subplot(2,2,[.15 .05],[.1 0.05]);
% 
% axes(plot_h(1));
% imagesc(qQ,pQ,l1);colorbar('SouthOutside');
% hold on
% contour(qQ,pQ,l1,'ShowText','on','LineColor',[1 1 1])
% set(gca,'YDir','Normal')
% ylabel('$\mathbf{p}$', 'interpreter', 'latex','FontSize',12)
% set(get(gca,'YLabel'),'Rotation',360);
% ylabh = get(gca,'YLabel');
% set(ylabh,'Position',get(ylabh,'Position') - [0.1 0.05 0.00]);
% set(gca, 'XAxisLocation', 'top')
% xlabel('$\mathbf{q}$', 'interpreter', 'latex','FontSize',12)
% plot(qQ(Jlow),pQ(Ilow),'rs','LineWidth',2,'MarkerSize',14);hold on
% plot(qQ(Jupp),pQ(Iupp),'ks','LineWidth',2,'MarkerSize',14);hold on
% axis square
% text(0.65, -0.6 ,'$\mathbf{\| m^{(k)} - m^{*} \|_1}$', 'interpreter', 'latex','FontSize',12)
% text(0.9,-.8,'(a)', 'interpreter', 'latex','FontSize',14);
% % caxis([0 10])
% % for ii = 1 : 3
% %     
% %     for jj = 1 : 3
% %         
% %         iq = [1 11 21];
% %         
% %         text(qQ(iq(jj)),pQ(iq(ii)),['$\|m\|_{' num2str(qQ(iq(jj))) '} + \|m\|_{' num2str(pQ(iq(ii))) '}$'],'Color','w','FontSize',8,'BackgroundColor','k', 'interpreter', 'latex','FontSize',15)
% %         
% %     end
% %     
% % end
% % load('misfit');
% axes(plot_h(2));
% imagesc(qQ,pQ,misfit);%caxis([39 41]);
% colorbar('SouthOutside');
% hold on
% % contour(pvec,qvec,misfit,'ShowText','on','LineColor',[0 0 0])
% set(gca,'YDir','Normal')
% set(gca, 'XAxisLocation', 'top')
% set(gca,'YAxisLocation', 'right');
% set(get(gca,'YLabel'),'Rotation',360);
% ylabh = get(gca,'YLabel');
% set(ylabh,'Position',get(ylabh,'Position') + [0.1 .05 0.00]);
% ylabel('$\mathbf{p}$', 'interpreter', 'latex','FontSize',12)
% xlabel('$\mathbf{q}$', 'interpreter', 'latex','FontSize',12)
% plot(qQ(Jlow),pQ(Ilow),'rs','LineWidth',2,'MarkerSize',14);hold on
% plot(qQ(Jupp),pQ(Iupp),'ks','LineWidth',2,'MarkerSize',14);hold on
% % plot(qQ(J),pQ(I),'ro','LineWidth',2);hold on
% axis square
% caxis([19.75 20.25])
% text(0.7, -0.6 ,'$\mathbf{\| Fm - d} \|_2$', 'interpreter', 'latex','FontSize',12)
% text(0.9,-.8,'(b)', 'interpreter', 'latex','FontSize',14);
% % caxis([38 42])
% 
% % for ii = 1 : 3
% %     
% %     for jj = 1 : 3
% %         
% %         iq = [1 11 21];
% %         
% %         text(qQ(iq(jj)),pQ(iq(ii)),num2str(3*(ii-1)+jj),'Color','w','FontSize',15,'BackgroundColor','k')
% %         
% %     end
% %     
% % end
% % [I,J]=find(l1==min(min(l1)));
% % set(figure, 'Position', [50 200 750 750]);
% % imagesc(qQ,pQ,l1);colorbar('SouthOutside');colormap winter;hold on
% % contour(pQ,qQ,l1,'ShowText','on','LineColor',[0 0 0])
% % set(gca,'YDir','Normal')
% % ylabel('$\mathbf{\| m \|_p}$', 'interpreter', 'latex','Fonts',15)
% % xlabel('$\mathbf{\| \nabla m \|_q}$', 'interpreter', 'latex','Fonts',15)
% % plot(qQ(J),pQ(I),'ro','LineWidth',2);hold on
% % axis square
% % title('\textbf{Model error} ($\mathbf{\| m - m^{*} \|_1}$)', 'interpreter', 'latex','Fonts',15)
% 
% %% Plot models_out, one at a time
% mm = models_out(Ilow,Jlow,1,:);
% 
% axes(plot_h(4));
% plot(x,model,'LineWidth',2); hold on
% plot(0.1,0.425,'rs','LineWidth',2,'MarkerSize',14);hold on
% plot(x, mm(:),'r','LineWidth',3); hold off
% ylabel('$\mathbf{m}$', 'interpreter', 'latex','FontSize',14)
% set(get(gca,'YLabel'),'Rotation',360);
% ylabh = get(gca,'YLabel');
% set(ylabh,'Position',get(ylabh,'Position') - [0.05 0.0 0.00]);
% xlabel('$\mathbf{z}$', 'interpreter', 'latex','FontSize',14);
% axis square
% grid on
% axis square
% text(0.3,0.425,['$\mathbf{p: ' num2str(pQ(Ilow)) ',\; q:' num2str(pQ(Jlow)) '}$'], 'interpreter', 'latex','FontSize',12);
% text(0.3,0.395,['\underline{Final Iteration}'],'interpreter', 'latex','FontSize',12)
% text(0.3,0.36,['$\mathbf{\phi_d} = ' num2str(round(misfit(Ilow,Jlow,1)*1000)/1000) '$'],'interpreter', 'latex','FontSize',12)
% text(0.3,0.325,['$\|\delta m\|_1$ = ' num2str(round(norm(model-mm(:),1)*100)/100) ],'interpreter', 'latex','FontSize',12)
% axis([0 1 -0.05 0.45])
% text(0.47,-0.15,'(d)', 'interpreter', 'latex','FontSize',14);
% set(gca,'YTickLabel',[]);
% tlabel = ['$\mathbf{m^{(k)}}$'];
% title(tlabel, 'interpreter', 'latex','FontSize',14)
% mm = models_out(Iupp,Jupp,1,:);
% 
% axes(plot_h(3));
% plot(x,model,'LineWidth',2); hold on
% plot(0.1,0.425,'ks','LineWidth',2,'MarkerSize',14);hold on
% plot(x, mm(:),'r','LineWidth',3); hold off
% ylabel('$\mathbf{m}$', 'interpreter', 'latex','FontSize',14)
% set(get(gca,'YLabel'),'Rotation',360);
% ylabh = get(gca,'YLabel');
% set(ylabh,'Position',get(ylabh,'Position') - [0.05 0.0 0.00]);
% xlabel('$\mathbf{z}$', 'interpreter', 'latex','FontSize',14);
% axis square
% grid on
% axis square
% text(0.3,0.425,['$\mathbf{p: ' num2str(pQ(Iupp)) ',\; q:' num2str(pQ(Jupp)) '}$'], 'interpreter', 'latex','FontSize',12);
% text(0.3,0.395,['\underline{Final Iteration}'],'interpreter', 'latex','FontSize',12)
% text(0.3,0.360,['$\mathbf{\phi_d} = ' num2str(round(misfit(Iupp,Jupp,1)*1000)/1000) '$'],'interpreter', 'latex','FontSize',12)
% text(0.3,0.325,['$\|\delta m\|_1$ = ' num2str(round(norm(model-mm(:),1)*100)/100) ],'interpreter', 'latex','FontSize',12)
% axis([0 1 -0.05 0.45])
% tlabel = ['$\mathbf{m^{(k)}}$'];
% title(tlabel, 'interpreter', 'latex','FontSize',14)
% text(0.47,-0.15,'(c)', 'interpreter', 'latex','FontSize',14);
