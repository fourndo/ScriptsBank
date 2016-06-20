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
nk = 5; %number of frequencies

decay = -3;

basis = 1;

% Data noise
amp_pct = 0.00;
floor_pct = 0.1;

tol = 0.4;
% Lp-norm parameters
pQ = [2 0 0 0];%0:0.1:2;%[0 1 2];%0%
qQ = [2 0 0 2];%0:0.1:2;%[0 1 2];%0%
lQ = [1 2 0 1];%ones(1,1)*1.0%[0.25:0.25:1.75]%0.1:0.4:1.9;


% Percentile for cutoff
pct_cutoff = 75;

eps_p = 1e-8;
eps_q = 1e-8;
%% SCRIPT STARTS HERE
%% % Generate kernel functions and depth weighting
x = (0 : 1/(nx) : 1);
mcell=nx;
dx=ones(1,mcell) * abs((min(x)-max(x))/mcell);


set(figure, 'Position', [50 0 775 400]);
axs{1} = axes('Position',[0-.1 .175 .75 .75]);
axs{2} = axes('Position',[0.425 .175 .75 .75]);
% plot_h = tight_subplot(2,2,[.05 .01],[.1 0.1]);

G = zeros(nk,mcell);
wr = zeros(1,mcell);
axes(axs{2});
for ii = 1 : nk 
    
    b = basis *pi* (ii-1);
    a = decay ;
    
    G(ii,1:mcell) = exp( a * x(2:end) ) .* (a/(a^2+b^2)*cos (b * x(2:end)) + b/(a^2+b^2) *sin (b * x(2:end))) -...
        exp( a * x(1:end-1) ) .* (a/(a^2+b^2)*cos (b * x(1:end-1)) + b/(a^2+b^2) *sin (b * x(1:end-1)));

%     G(ii,1:mcell) = G(ii,1:mcell) / max(G(ii,1:mcell));
    
    wr = wr + ((exp( a*x(2:end))-exp( a * x(1:end-1) ))/a) ;
    plot( x(2:end) , G(ii,:), 'LineWidth',2,'Color',[0 0 0]+(ii)/(nk+2) );hold on

end

wr = (sum(G.^2,1)).^0.5;
wr = (wr / max(wr)).^0.5;
Wr = spdiags(wr',0,mcell,mcell);

% Wrx = spdiags(wr(1:end-1)',0,mcell-1,mcell-1);

% figure(1); plot(wr','r--','LineWidth',3);
axis square;grid on
% temp = legend('\bfg(z) = e^{-jz} cos(j 2\pi z)','\bfDepth Weighting (Wr)');
ylabel('$\mathbf{g}(m)$', 'interpreter', 'latex','FontSize',16)
% set(get(gca,'YLabel'),'Rotation',360);
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') - [0.05 0.0 0.00]);
text(0.9,max(G(:))*0.9,'(b)', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle');
xlabel('$\mathbf{x}$', 'interpreter', 'latex','FontSize',16);
% a=get(temp,'children');
% set(a(2),'LineWidth',3,'Color',[1 0 0],'LineStyle','--'); 
%% Generate model

x = (1/nx : 1/(nx) : 1) ;

cntr = round(nx/4);

% Create mdoel: Square and tanh functions
model = 0.4* exp(-(((x-3*cntr/nx) ).^2 * nx)); 

model(cntr-round(cntr/4):cntr+round(cntr/4))=0.25;

model=model(:);

axes(axs{1});
plot(x,model,'k','LineWidth',2);
axis([x(1) x(end) -0.075 0.425]);hold on
ylabel('$\mathbf{m}$', 'interpreter', 'latex','FontSize',16)
% set(get(gca,'YLabel'),'Rotation',360);
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') - [0.05 0.0 0.00]);
text(0.9,0.38,'(a)', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle');
xlabel('$\mathbf{x}$', 'interpreter', 'latex','FontSize',16);
axis square
grid on


%% FORWARD MODELING
% Generate data and corrupt with noise

data= G * model;

% Corrupt data with Gaussian noise
rand_noise = rand(nk,1);
% save ('noise_1D','rand_noise')
% load('noise_10nk');

amp_noise = amp_pct .* abs(data);
floor_noise = floor_pct * min(abs(data));

wd = amp_noise + floor_noise;
% noise = amp_noise + floor_noise;

d = data + rand_noise.*wd;
% wd =  abs(d)*0.05 + 0.05 * std(data);

figure
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
s= ones(mcell,1);
t= ones(mcell,1);

% Global constant
as = 1.0;  %Smarrnest term
ax = 1.0;               %Smoothness term                     

% Store arr the final models
RMS = zeros(length(pQ),length(qQ),length(lQ));
linf = zeros(length(pQ),length(qQ),length(lQ));
l2 = zeros(length(pQ),length(qQ),length(lQ));
l1 = zeros(length(pQ),length(qQ),length(lQ));
misfit = zeros(length(pQ),length(qQ),length(lQ));
models_out = zeros(length(pQ),length(qQ),length(lQ),mcell);
iter_num = zeros(length(pQ),length(qQ),length(lQ));

Pac = speye(mcell);
            
dkdt = @(p,ep) ((ep).^(1./(4-2*p)) - ep);

counter = 1;  

objfunc = @(m,phi,b) sum( ( G * m - d ).^2 ) + ( m' * b * phi * m );
%%
% set(figure, 'Position', [50 0 775 800]);
set(figure, 'Position', [50 0 775 400]);
axs{1} = axes('Position',[-.1 .175 .8 .8]);
p = [0 0.5 1 1.5 2];
epsl = 1e-2;
xx=-1:1e-5:1;
for pp = 1 : length(p)
    
    r = 1./(xx.^(2) + epsl.^2).^(1-p(pp)/2);
    
    dphi_dm = (xx).*r;
%     scale = 1./(max(sqrt(r)));

    scale = epsl^((1-p(pp)/2));
    dphi_dm = dphi_dm * scale;
    
%     dphi_dm = (x*(1e-2).^(1-p(pp)/2))./(x.^(2) + 1e-2.^2).^(1-p(pp)/2);
    
%     phim = phim./ ( max(x.*phim));
%     phim = phim.^0.6;


     plot(xx, dphi_dm,'k','LineWidth',p(pp)+0.25);hold on
     
     if p(pp) == 2 || p(pp)==1.5
         
              text(0.22,(0.22*(1e-2).^(1-p(pp)/2)) ./ (0.22.^(2) + (1e-2).^2).^(1 - p(pp)/2),['p =' num2str(p(pp))],...
        'BackgroundColor','w','EdgeColor','k')
         
     else
         
     text(0.02,(0.02*(1e-2).^(1-p(pp)/2)) ./ (0.02.^(2) + (1e-2).^2).^(1 - p(pp)/2),['p =' num2str(p(pp))],...
        'BackgroundColor','w','EdgeColor','k')
    
     end
end
text(sqrt(epsl),sqrt(epsl),['$\sqrt{\epsilon}$'],'BackgroundColor','w','EdgeColor','k','interpreter', 'latex','FontSize',14,'HorizontalAlignment','left','VerticalAlignment','bottom')
plot(sqrt(epsl),sqrt(epsl),'ko')

xlim([-0.01 0.25])
ylim([-0.1 0.6])
grid on
yy = get(gca,'Ylim');
plot([10^(-2) 10^(-2)],yy,'r-.')
ylabel('$\frac{ \partial \hat\phi_m}{\partial m}$', 'interpreter', 'latex','FontSize',16)
set(get(gca,'YLabel'),'Rotation',360);
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') - [0.025 0 0]);
xlabel('x', 'interpreter', 'latex','FontSize',14)
text(.215,0.53,'(a)', 'interpreter', 'latex','FontSize',14);
axis square
% axs{3} = axes('Position',[0.1 .55 .4 .4]);
% axs{4} = axes('Position',[0.55 .55 .4 .4]);
%%
% Iterate over arr the combination of norms as specified by the vectors
% pvec, qvec and lvec.
axs{2} = axes('Position',[0.35 .175 .8 .8]);
for rr = 4

        % Message prompt
        head = ['lp: ' num2str(pQ(rr)) ' lq: ' num2str(qQ(rr)) ' psi: ' num2str(lQ(rr))];
        fprintf('Starting lp inversion %s\n',...
            head)

        invmod      = ones(mcell,1)*1e-1;       % Initial model       

        phi_d       = sum((G*invmod - d).^2);   % Initial misfit

                   % Active cerr
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

        if rr ~=4
            
            lp    = ones(mcell,1)*pQ(rr);
            lq    = ones(mcell,1)*qQ(rr);
            ll    = ones(mcell,1)*lQ(rr);
        
        else

            LP = [0 0 1; 0.5 2 1];
            % Define zones {x,y,lp,lqx,lqz,ll}
            zone{1,1} = 1:100; 
            zone{2,1} = 101:200;
            zone{3,1} = 76:124;
            t = zeros(mcell,2);
            t(zone{1},1) = 1; t(zone{3},1)= 1:-1/length(zone{3}):1/length(zone{3});
            t(zone{2},2) = 1; t(zone{3},2)= 1/length(zone{3}):1/length(zone{3}):1;


            lp  = zeros(mcell,1);
            lq  = zeros(mcell,1);
            ll  = ones(mcell,1);
            
            % Generate lp vectors
            for jj = 1 : size(LP,1)

                lp =  lp + t(:,jj)*LP(jj,1);
                lq =  lq + t(:,jj)*LP(jj,2);
                %ll =  ll + t(:,jj)*LP(jj,3);
            end
            
        end

        dphi_m = 100;
        l2_count = 0;

        while switcher ~= 3

                count=count+1;

                if switcher == 0   %First iteration
                    l2_count = l2_count+1;
                    delta_p(count) = 1e-1;%max(invmod);%prctile(abs(invmod(invmod ~= 0)),pct_cutoff);
                    dmdx = abs(Wx * invmod);
                    delta_q(count) = 1e-1;%max(dmdx);%prctile(gradm(gradm~=0),pct_cutoff);

                    aVRWs = sqrt(as) * Wr * V * Ws;
                    aVRWx = sqrt(ax) * Wr * Vx * Wx;

                    lvec = ones(mcell,1);

                    MOF = aVRWs'*aVRWs + aVRWx'*aVRWx;

                    if count==1
                        beta = sum(sum(G.^2,1)) / sum(diag(MOF,0).^2) *1e+2;
                    end


                    phi(count) = objfunc(invmod,MOF,beta(count));
                    gamma = 1;


                else

                    lp_count = lp_count+1;
                    
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

%                         scale_p = zeros(mcell,1);
%                         scale_x = zeros(mcell,1);


                    rs = 1./ ( abs(invmod) .^2 + delta_p(count).^2 ) .^( (1-lp/2) );
                    rx = 1./ ( abs(Wx * invmod) .^2 + delta_q(count).^2 ) .^( (1-lq/2) );


                    eta_s = delta_p(count).^ (1-lp/2);                        
                    eta_x = delta_q(count).^ (1-lq/2);


                    sRs = spdiags( ( eta_s .* rs ).^0.5 ,0,mcell,mcell);
                    sRx = spdiags( ( eta_x .* rx ).^0.5 ,0,mcell,mcell);


                    aVRWs =  spdiags( sqrt( as *  ll ),0,mcell,mcell)  * Wr * sRs * Ws ;
                    aVRWx =  spdiags( sqrt( ax * (2.0 - ll)  ),0,mcell,mcell) * Wr * sRx * Wx ;
                    
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
                dphi_m(count) = abs(phi_m(count) - phi_m(count-1))/phi_m(1) *100;

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

            beta_count = beta_count+1;

            if switcher == 0 && phi_d(count) < target * (1-tol)

                % Come back by a fraction and start IRLS
                beta(count+1) = beta(count);
                switcher = 1;                    

            elseif switcher == 0 

                beta(count+1) = 0.9 * beta(count);

            else


                beta(count+1) =beta(count) * target / phi_d(count);



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

            if lp_count > 0  && dphi_m(end) < 1 && switcher==2 && phi_d(count) > target *(1-tol) && phi_d(count) < target *(1+tol)

                switcher = 3;

                continue

            end      
    
        % OPTIONAL: Plot model
%             set(figure(4), 'Position', [50 200 1000 500])
% axes(plot_h((rr)*2+1));
%             axes(axs{rr});
%             plot(1:mcell,model); hold on
%             plot(1:mcell, invmod,'r','LineWidth',3); hold off
%             axis([0 mcell -.1 0.6]);
%             axis square



        end
                % OPTIONAL: Plot model

        

        model_error = norm(model-invmod,1);




        
        if rr < 4
            axes(axs{1});
            
            if rr == 1
            plot(x,model,'k--','LineWidth',1); hold on
            plot(x, invmod,'k','LineWidth',2); hold on
            
            elseif rr == 2
                
                plot(x, invmod,'k','LineWidth',1); hold on
            else
            
            	plot(x, invmod,'k:','LineWidth',2); hold on
                legend('True','p=2, q=2','p=0','q=2','Location','North')
            end
            
        else
            axes(axs{2});
            
            plot(x,model,'k--','LineWidth',1); hold on
            plot(x, invmod,'k','LineWidth',2); hold on
            legend('True','Mixed norm','Location','North')
            plot([0.5,0.5],[-1,1],'k:')
            text(0.5,0.2,'Transition','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom')
            text(.1,0.-0.035,['$\mathbf{p: ' num2str(LP(1,1)) ',\; q: ' num2str(LP(1,2)) '}$'],'interpreter', 'latex','FontSize',12)
            text(.65,0.-0.035,['$\mathbf{p: ' num2str(LP(2,1)) ',\; q: ' num2str(LP(2,2)) '}$'],'interpreter', 'latex','FontSize',12)
            quiver(0.4,-0.035,0.2,0, 'LineWidth', 2, 'MaxHeadSize',0.5, 'Color','k');
            
            grid on
        end
       

        
%         plot(x, invmod,'r','LineWidth',2); hold on

%         if rr == 1
%             
%             text(.3,0.40,['$\mathbf{p: ' num2str(pQ(rr)) '}$'],'interpreter', 'latex','FontSize',12)
%         else
%             
%             text(.3,0.40,['$\mathbf{q:' num2str(qQ(rr)) '}$'],'interpreter', 'latex','FontSize',12)
%  
%         end
%         text(.3,0.36,['$\mathbf{\phi_d} = ' num2str(round(phi_d(count)*100)/100) '$'],'interpreter', 'latex','FontSize',10)
%         text(.35,0.32,['$\epsilon = 10^{' num2str(round(log10(delta_p(end))*10)/10) '}$'],'interpreter', 'latex','FontSize',10)
%         text(.3,0.32,['$\|\mathbf{m - m^*}\|$ = ' num2str(round(model_error*100)/100)],'interpreter', 'latex','FontSize',10)

        
        axis([0 1 -.075 0.425]);

        if rr == 1
            ylabel('\bf{m}', 'interpreter', 'latex','FontSize',14);

        
            
            set(get(gca,'YLabel'),'Rotation',360);
            ylabh = get(gca,'YLabel');
            set(ylabh,'Position',get(ylabh,'Position') - [0.075 0 0.00]);
              
            
        end
        
        xlabel('x', 'interpreter', 'latex','FontSize',14);
        if rr==1


            text(0.9,0.38,'(a)', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle');
        elseif rr ==4

            set(gca,'YTickLabel',[]);
            text(0.9,0.38,'(b)', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle');



        end
        axis square

        


%             figure; plotyy(1:count,phi_d,1:count,phi_m);

        ldml(counter) = norm(model-invmod,2);
%             fprintf('Relative dm:\t %8.5e \n', rdm(count));
%             fprintf('End of lp inversion. Number of iterations: %i\n',count)
%             fprintf('Final model misfit: %8.3e\n\n',ldml(counter))

        counter = counter+1;
            
end


