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

% addpath data\
addpath C:\Users\DominiqueFournier\Documents\GIT\CodeBank\FUNC_LIB
% addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\
% addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\Lp_norm_codes\Dual_norm_1D\functions
% addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\Lp_norm_codes\Dual_norm_1D\sliceomatic\
% Set up 1D problem

%% INPUT VARIABLES
% Lp-norm parameters
% lp = 2%0:0.1:2;
% lq = 0%0:0.1:2;
% mu = ones(1,1)*1.0%[0.25 0.5 1 1.5 1.75]%0.1:0.4:1.9;

% Inversion parameters
FLAG2 = 'SMOOTH_MOD'; % Choose between SMOOTH_MOD | SMOOTH_MOD_DIF

% Model space
nx = 60;
nz = 60;

mcell=nx*nz;

dx = 1 / nx;
dz = 1 / nz;

% Reference model
mref = zeros(nx*nz,1);

% Create observation locations
obsx = [zeros(1,9) ones(1,9)];
nobs = length(obsx);
obsz = [-0.1:-0.1:-0.9 -0.1:-0.1:-0.9];


dx = ones(nx,1) * dx;
dz = ones(nz,1) * dz;

x0 = 0;
z0 = 0;

xc = [x0 + cumsum(dx) - dx/2];
zc = [z0 - cumsum(dz) + dz/2];

% Topography (all ones if flat)
nullcell = ones(nz,nx);

% Kernel
nk = 10; %number of frequencies

decay = -0.25;

basis = 0.5;

% Data noise
amp_pct = 0.05;
floor_pct = 0.05;

% Percentile for cutoff
pct_cutoff = 50;

%% Create zones and transition for lp,lq,lu

% Create smoothing operator for transitions
av = @(n) spdiags (kron(ones(n,1),[0.25,0.5,0.25]),[-1,0,1],n,n);

avcx = av(nx); avcx(1,1:2) = 0.5; avcx(end,end-1:end) = 0.5;
avcz = av(nz); avcz(1,1:2) = 0.5; avcz(end,end-1:end) = 0.5;

Avcx = kron ( avcx, speye (nz));
Avcz = kron ( speye (nx), avcz);

A = (Avcx * Avcz);
A = A^3;
A = spdiags(1./sum(A,2),0,mcell,mcell) *A;

% Define zones {x,y,lp,lqx,lqz,ll}
zone{1,1} = 1:nx; zone{1,2} = 1:nz;
LP = [1 1 1 1];

% Mulitple zones
% zone{1,1} = 10:29; zone{1,2}=10:50; 
% zone{2,1} = 32:50; zone{2,2}=10:50; 
% zone{3,1} = []; zone{3,2}=[];
% LP = [2 0 0 1;
%       0 2 2 1;
%       0 1 1 1];

% zone{1,1} = 1:nx; zone{1,2}=1:nz; 
% LP = [1 2 2 1];
    
% Build background tile
s = zeros(mcell,size(zone,1));
s(:,end) = ones(mcell,1);

for ii = 1 : size(zone,1)-1
       
    temp = zeros(nz,nx); 
    temp(zone{ii,1},zone{ii,2}) = 1;
    s(:,ii) = temp(:);
    
    % Zero out tiles from the background
    s(s(:,ii)==1,end) = 0;
    
end

% Build transition matrices and lp,lq,ll vector
t = A * s;

% t = spdiags(1.,0,mcell,mcell) * t/sum(t,2);

% Get lp,lq and mu everywhere
% pvec    = T * pps';
% qxvec   = T * qqx';
% qzvec   = T * qqz';
% uvec    = T * uu';

set(figure, 'Position', [50 0 775 1000]);
lp_mod = zeros(mcell,1);
lqx_mod = zeros(mcell,1);
lqz_mod = zeros(mcell,1);

for jj = 1 : size(zone,1)
    
    lp_mod =  lp_mod + t(:,jj)*LP(jj,1);
    lqx_mod =  lqx_mod + t(:,jj)*LP(jj,2);
    lqz_mod =  lqz_mod + t(:,jj)*LP(jj,2);
    
    if jj==1
        axes('Position',[-0.25+jj/3.25 .2 .28 .28])    
        imagesc(xc,zc,reshape(t(:,jj),nz,nx));
        
    elseif jj == 2
    axes('Position',[-0.251+jj/3.25 .165 .282 .282])    
    imagesc(xc,zc,reshape(t(:,jj),nz,nx));
    colorbar('SouthOutside')
    set(gca,'YTickLabel',[]);
    xlabel('$x$', 'interpreter', 'latex','FontSize',16);
    text(.45,-1.35,'$\mathbf{T}_{jk}$', 'interpreter', 'latex','FontSize',14);
    text(.45,-1.5,'(c)', 'interpreter', 'latex','FontSize',14);

    else
        
        axes('Position',[-0.25+jj/3.25 .2 .28 .28])    
    imagesc(xc,zc,reshape(t(:,jj),nz,nx));
%         set(gca,'YTickLabel',[]);
        set(gca, 'YAxisLocation', 'right')
    end
    set(gca, 'XAxisLocation', 'top')
    set(gca,'YDir','normal')
    axis square
end
          
            
axes('Position',[-0.255+1/3.25 .6 .28 .28])
imagesc(xc,zc,reshape(lp_mod,nz,nx));
axis square
colorbar('SouthOutside')
set(gca, 'XAxisLocation', 'top')
set(gca,'YDir','normal')
text(.43,-1.4,'$\|m\|_p$', 'interpreter', 'latex','FontSize',14);
text(.45,-1.5,'(a)', 'interpreter', 'latex','FontSize',14);


axes('Position',[-0.23+2/3.25 .6 .28 .28])
imagesc(xc,zc,reshape(lqx_mod,nz,nx));
axis square
% set(gca,'YTickLabel',[]);
colorbar('SouthOutside')
set(gca, 'XAxisLocation', 'top')
set(gca,'YDir','normal')
text(.37,-1.4,'$\|\nabla_x m\|_p$', 'interpreter', 'latex','FontSize',14);
text(0.98,-1.5,'(b)', 'interpreter', 'latex','FontSize',14);
xlabel('$x$', 'interpreter', 'latex','FontSize',16);

axes('Position',[-0.255+3/3.25 .6 .28 .28])
imagesc(xc,zc,reshape(lqz_mod,nz,nx));
axis square
set(gca, 'YAxisLocation', 'right')
colorbar('SouthOutside')
set(gca, 'XAxisLocation', 'top')
set(gca,'YDir','normal')
text(.37,-1.4,'$\|\nabla_z m\|_p$', 'interpreter', 'latex','FontSize',14);

%% Generate model

model = zeros(nz,nx);

% Center the anomaly
dmx = floor(0.1 / dx(1));
dmz = floor(0.1 / dz(1));

centerx = floor(nx/2);
centerz = floor(nz/3);

% Square bloc
model(centerz-5:centerz+5,centerx-5:centerx+5)= -0.5;

% Gaussian anomaly
zz = 2*centerz-12 :2*centerz+12;
xx = centerx-12:centerx+12;

[ZZ,XX] = ndgrid(cumsum(dz(zz)) - sum([dz(zz);dz(zz(1))])/2,cumsum(dx(xx)) - sum([dx(xx);dx(xx(1))])/2);
r = sqrt(XX.^2 + ZZ.^2);
dr = sqrt(dx(1)^2 + dz(1)^2);

model(zz,xx)= 0.75 * exp(-abs(((r)*10).^2));


set(figure, 'Position', [50 0 775 1000]);

subplot(2,2,1)
imagesc(xc,zc,model);hold on
contour(xc,zc,model,[-0.1 -0.05 0.05:0.1:0.5],'Color','w');
ylabel('$z$', 'interpreter', 'latex','FontSize',16)
set(get(gca,'YLabel'),'Rotation',360);
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') - [0.05 0.0 0.00]);
set(gca, 'XAxisLocation', 'top')
xlabel('$x$', 'interpreter', 'latex','FontSize',16);
text(.40,-1.3,'Model', 'interpreter', 'latex','FontSize',14);
text(.45,-1.4,'(a)', 'interpreter', 'latex','FontSize',14);
axis square
grid on
            

caxis([-.6 0.6]);
colorbar('SouthOutside')
set(gca,'YDir','normal')
% axis([-1.1 1.1 -1.1 0.1])
axis square

hold on
plot(obsx,obsz,'r^','MarkerFaceColor','y','LineWidth',2)

%% SCRIPT STARTS HERE
%% % Generate kernel functions and depth weighting
xn = [x0;x0 + cumsum(dx)];
zn = [z0;z0 - cumsum(dz)];

[Zn,Xn] = ndgrid(zn,xn);

[Zc,Xc] = ndgrid(zc,xc);

[dZ , dX] = ndgrid(dz,dx);

ndata = nk * nobs;

%% Generate forward operator and depth weighting
G = zeros(ndata,mcell);
aa = zeros(ndata,nx);
wr = zeros(mcell,1);
count = 1;

for jj = 1 : nobs;
     
    for ii = 1 : nk 
 
        dXn = Xn - obsx(jj);
        dZn = Zn - obsz(jj);

        a = decay * (ii) ;
        b = basis * 2*pi* (ii) ;

        expx = @(m) exp( a * abs(m) );

        trig1 = @(m) cos(b*(m));       

        xupb = dXn(2:end,1:end-1);
        xlob = dXn(2:end,2:end);
        
        zlob = dZn(1:end-1,2:end);
        zupb = dZn(2:end,2:end);
        
        fupxupz = expx( sqrt( xupb.^2 + zupb.^2 ) ) .* trig1(sqrt( xupb.^2 + zupb.^2 ));

        fupxloz = expx( sqrt( xupb.^2 + zlob.^2 ) ) .* trig1(sqrt( xupb.^2 + zlob.^2 ));

        floxupz = expx( sqrt( xlob.^2 + zupb.^2 ) ) .* trig1(sqrt( xlob.^2 + zupb.^2 ));

        floxloz = expx( sqrt( xlob.^2 + zlob.^2 ) ) .* trig1(sqrt( xlob.^2 + zlob.^2 ));

        temp = (fupxupz + floxupz+ fupxloz + floxloz)/4 .* (xupb-xlob).*(zupb-zlob);
            
            
        G(count,:) = temp(:);


        temp = (expx( sqrt( xlob.^2 + zlob.^2 ) ) + expx( sqrt( xlob.^2 + zupb.^2 ) ) +...
            expx( sqrt( xupb.^2 + zlob.^2 ) ) + expx( sqrt( xupb.^2 + zupb.^2 ) )) / 4;
       

        wr = wr + temp(:);
 
        if jj == 5 && ii == 10
            
            subplot(2,2,2)
            imagesc(xc,zc,reshape(G(count,:),nx,nz)); hold on
            contour(xc,zc,model,[-0.1 -0.05 0.05:0.1:0.5],'Color','w');
            caxis([-2e-4 2e-4])
            axis square;grid on
            % temp = legend('\bfg(z) = e^{-jz} cos(j 2\pi z)','\bfDepth Weighting (Wr)');
            ylabel('$z$', 'interpreter', 'latex','FontSize',16)
            set(get(gca,'YLabel'),'Rotation',360);
            ylabh = get(gca,'YLabel');
            set(ylabh,'Position',get(ylabh,'Position') - [0.05 0.0 0.00]);
            text(.40,-1.3,'Kernel', 'interpreter', 'latex','FontSize',14);
            text(.45,-1.4,'(b)', 'interpreter', 'latex','FontSize',14);
            set(gca, 'XAxisLocation', 'top')
            xlabel('$x$', 'interpreter', 'latex','FontSize',16);
            colorbar('SouthOutside')
            set(gca,'YDir','normal')
            % axis([-1.1 1.1 -1.1 0.1])
            axis square
            plot(obsx,obsz,'r^','MarkerFaceColor','y','LineWidth',2)
        end

        count = count + 1;

        
    end
    
end
% wr = sqrt(sum(G.^2,1));
wr = wr / max(wr);
wr = (wr(:)).^0.5;


%% FORWARD MODELING
% Generate data and corrupt with noise
model=model(:);
data= G * model;

% Corrupt data with Gaussian noise
rand_noise = rand(ndata,1);
save ('noise','rand_noise')
load('noise');

amp_noise = amp_pct .* abs(data);
floor_noise = floor_pct * min(abs(data));

wd = amp_noise + floor_noise;
% noise = amp_noise + floor_noise;

d = data + rand_noise.*wd;

axes('Position',[.2 .2 .6 .3])
for ii = 1:nobs
    
    plot(data((1+(ii-1)*nk):(nk+(ii-1)*nk)));
    hold on
    errorbar(d((1+(ii-1)*nk):(nk+(ii-1)*nk)),wd((1+(ii-1)*nk):(nk+(ii-1)*nk)),'r*')
    legend('Data','Data+noise')
    grid on
end

grid on
yy = get(gca,'Ylim');
plot([10^(-2) 10^(-2)],yy,'r-.')
% ylabel('$Amplitude$', 'interpreter', 'latex','FontSize',16)
set(get(gca,'YLabel'),'Rotation',360);
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') - [0.01 0 0]);
xlabel('$Data$', 'interpreter', 'latex','FontSize',14)
text(5.7,-0.028,'(c)', 'interpreter', 'latex','FontSize',14);
% axis square
% Create uncertainty weighting matrix
Wd = spdiags(1./wd,0,ndata,ndata);

% Normalize d and G by standard deviation:
G = Wd * G ;

target = ndata;

% Weight data with uncertainty
d = Wd * d;

%% INVERSION
% Create inversion parameters
[ Wx, Wz, Vx, Vz ] = get_GRAD_op2D_SQUARE(dx,dz,nullcell);
Ws = speye(mcell);
V = spdiags(sqrt(dX(:).*dZ(:)),0,mcell,mcell);


% Global constant
alpha(1) = 1.0 / min(dx) ^2;  %Smallnest term
alpha(2) = 1.0;               %Smoothness term                     
alpha(3) = 1.0;

set(figure(5), 'Position', [50 0 775 400]);
% plot_h = tight_subplot(1,2,[.1 .025],[.1 0.1]);
for rr = 1 : 2

    if rr == 1
        
        FLAG1 = 'GRADm'; % Choose between GRADm | lGRADm|
        
    else
        
        FLAG1 = 'lGRADml'; % Choose between GRADm | lGRADm|
        
    end
    
phi_m = [];
phi = [];
counter = 1;  
dkdt = @(p,ep) ep .^(1/(2*(2-p)));
% dkdt = @(p,ep)(ep).^(1/(2*(2-p))) - p*ep;
polyn = @(a,ep) [1 0 2*ep-a 0 (ep^2 + a *ep)];
l0norm = @(x,ep) x./(x.^(2) + ep);
objfunc = @(m,phi,b) sum( ( G * m - d ).^2 ) + ( m' * b * phi * m );

% Iterate over all the combination of norms as specified by the vectors
% pvec, qvec and lvec.

% for ll= 1:length(mu)
% 
%     for pp = 1:length(lp)
%         
%         for qq = 1:length(lq)
            
            
%             delta(counter) = 1e-5;
            % Message prompt
%             head = ['lp: ' num2str(lp(pp)) ' lq: ' num2str(lq(qq)) ' psi: ' num2str(mu(ll))];
%             fprintf('Starting lp inversion %s\n',...
%                 head)
            
invmod      = ones(mcell,1)*1e-2;       % Initial model       

phi_d       = sum((G*invmod - d).^2);   % Initial misfit

count=0; % Initiate iteration count 
switcher = 0;
ncg(counter) = 0;
lp_count = 0;
cg_iter = 0;
delta_p = 1;
traffic = 1;
group = 0;
x_tresh = 1;
Pac = speye(mcell);
set(figure, 'Position', [50 0 775 1000]);   

while switcher ~= 3

            count=count+1;

            if switcher == 0   %First iteration
                delta_p(count) = 0.1;
                delta_x(count) = 0.1;

                [MOF,aVRWs,aVRWx,aVRWz] = get_lp_MOF_2D(invmod,mref,1,s,t,V,Ws,Vx,Wx,Vz,Wz,wr,alpha,LP,FLAG1,FLAG2,switcher,delta_p(count),delta_x(count));                                      
%                         [MOF,aVRWs,aVRWx,aVRWz] = get_lp_MOF_2D(invmod,mref,nx,nz,V,Ws,Vx,Wx,Vz,Wz,wr,alpha,2,2,1,FLAG1,FLAG2,delta(count));

                if count==1
                    % Initial beta trace(G'G) / trace(phim)
                    beta = sum(sum(G.^2,1)) / sum(diag(MOF,0).^2) * 1e-1; 
                    phi(count) = objfunc(invmod,MOF,beta(count));

                end

                x_tresh = 0.01;
                p_tresh = 0.01;

            else

                fprintf('# # LP-LQ ITER# #');
                lp_count = lp_count+1;



                if lp_count == 1
                    x_tresh = prctile(abs(invmod(invmod > 0)),pct_cutoff)*2;
                    
                    
                    delta_p(count) = 0.01;
                    delta_x(count) = 0.01;

                    %x_tresh.^2;
%                     figure; hist(abs(invmod),20);
%                     ylim([0 100])

                    subplot(2,2,1)
                    imagesc(xc,zc,reshape(invmod,nz,nx));
                    hold on;
                    contour(xc,zc,reshape(invmod,nz,nx),[-0.1 -0.05 0.01:0.1:0.5],'Color','w');
                    caxis([-.6 0.6]);
                    colorbar('SouthOutside')
                    set(gca,'YDir','normal')
                    axis square
                    set(gca, 'XAxisLocation', 'top')
%                     xlabel('$x$', 'interpreter', 'latex','FontSize',16);
                    text(.40,-1.3,'Model', 'interpreter', 'latex','FontSize',14);
                    text(.45,-1.4,'(a)', 'interpreter', 'latex','FontSize',14);
                    text(0.05,-0.05,'$\phi_m = \|\mathbf{m}\|_2 + \|\mathbf{\nabla m}\|_2$','interpreter', 'latex','FontSize',10,'EdgeColor','k','BackgroundColor','w');
                    ylabel('$z$', 'interpreter', 'latex','FontSize',16)
                    set(get(gca,'YLabel'),'Rotation',360);
                    xlabel('$x$', 'interpreter', 'latex','FontSize',14)
                    axis square
                    
                elseif (traffic_s(end)*100 > 1 || traffic_x(end)*100 > 1) && switcher == 1 

                    delta_p(count) = delta_p(end);
                    delta_x(count) = delta_x(end);
                    
                else

%                             switcher = 2;
                    delta_p(count) = delta_p(end);
                    delta_x(count) = delta_x(end);
                    switcher = 2;
                    
                    fprintf('\nADJUSTING BETA\n')
                end

                [MOF,aVRWs,aVRWx,aVRWz] = get_lp_MOF_2D(invmod,mref,phi_m(count-1),s,t,V,Ws,Vx,Wx,Vz,Wz,wr,alpha,LP,FLAG1,FLAG2,switcher,delta_p(count),delta_x(count));
                         
%                         lambda(count) = beta(count) * scale(count);
%                         scale = 1 / max(l0norm(invmod,delta(count)));
%                         x_tresh = dkdt(LP(:,1),delta(count));

%                         polynom = polyn(scale,delta(count));
%                         x_tresh = max(roots(polynom));

            end


         group_s(count) = sum(abs(invmod) <= p_tresh);
         
        dmdx = sqrt( (Wx * invmod).^2 + (Wz * invmod).^2 );
        group_x(count) = sum(abs(dmdx) <= x_tresh);
%                 eigval = eig(G'*G + lambda(count) * (phim) ); 
%                 condnum = max(abs(eigval)) / min(abs(eigval));
%                 
%                 figure(4); plot(eigval,'*');
%                 fprintf('Condition number: %f\n',condnum);
        fprintf('\n# # # # # # # # #\n');
        fprintf('BETA ITER: \t %i  \nbeta: \t %8.5e \n',count,beta(count));


        m_in = invmod;  % Save initial model for backtracking steps

        diagA = sum(G.^2,1) + beta(count)*spdiags(MOF,0)';
        PreC     = Pac * spdiags(1./diagA(:),0,mcell,mcell);

    [invmod,cg_iter(count),Pac] = GNsolver( G, invmod, d, phi(end), MOF, PreC, Pac, beta(count) , aVRWs, aVRWx, aVRWz, [] );                
    ncg(counter) = ncg(counter) + cg_iter(count);
%% Save iteration details

    if count > 1
        temp = sum(abs(invmod) <= p_tresh);                     
        traffic_s(count) = abs(group_s(count) - temp) / group_s(count);
        
        dmdx = sqrt( (Wx * invmod).^2 + (Wz * invmod).^2 );
        temp = sum(abs(dmdx) <= x_tresh);
        traffic_x(count) = abs(group_x(count) - temp) / group_x(count);
        
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
    phi(count) = objfunc(invmod,MOF,beta(count));

    if (rdm(count) < 1e-2 && switcher == 2 && phi_d(count) > target *0.98 && phi_d(count) < target *1.02) || count > 50 

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
    if phi_d(count) < target * 0.98 && switcher~=1


        fprintf('---------->')
        fprintf(' misfit:\t %8.5e ',phi_d(count))
        fprintf(' ** Relative dm:\t %8.5e ', rdm(count));
        fprintf('<----------\n')
        fprintf('\n# NEXT ITER - INCREASING BETA #\n');

%                 beta(count) = beta(count-1);
        if switcher == 0

            beta(count+1) = beta(count) * 1.5;
            switcher = 1;

% %                         x_tresh = dkdt(LP(:,1),delta(count));
%                     polynom = polyn(scale,delta(count));
%                     delta(count) = max(roots(polynom));%max(abs(invmod)).^((2-LP(:,1)));
        else

             beta(count+1) = beta(count) * 1.05;

        end
%                 beta(count+1) = beta(count) * target / phi_d(count);

%                 invmod = m_in;
%                 phi_d(count) = phi_d(count-1);
%                 phi_m(count) = phi_m(count-1);
%                 phi(count) = phi(count-1);
%                 
%                 count = count-1;
%                 lp_count = lp_count-1;

    elseif phi_d(count) > target * 1.02 && switcher ~= 1


        fprintf('---------->')
        fprintf(' misfit:\t %8.5e ',phi_d(count))
        fprintf(' ** Relative dm:\t %8.5e ', rdm(count));
        fprintf('<----------\n')
        fprintf('\n# NEXT ITER - REDUCING BETA #\n');
        if switcher == 0

            beta(count+1) = 0.5 * beta(count); 

%         elseif rdm(count) < 1e-2
% 
%             beta(count+1) = beta(count) * target / phi_d(count);

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

%     set(gca, 'XAxisLocation', 'top')
%     set(gca,'YDir','normal')
%     axis square



end
        
        %%    %% OPTIONAL: Plot model
        figure(5)
if rr == 1
              
    axes('Position',[-.1 .2 .75 .75]);
    imagesc(xc,zc,reshape(invmod,nz,nx)); hold on

else
    
    axes('Position',[0.4 .2 .75 .75]);
    imagesc(xc,zc,reshape(invmod,nz,nx)); hold on
    
end

caxis([-.6 0.6]);
% colorbar('SouthOutside')
set(gca,'YDir','normal')
axis square
contour(xc,zc,reshape(invmod,nz,nx),[-0.4:0.2:-.1 0.2:0.2:0.4],'Color','w');

%         xlabel('$x$', 'interpreter', 'latex','FontSize',16);
% text(.40,-1.3,'Model', 'interpreter', 'latex','FontSize',14);

if rr == 1
text(0.05,-0.1,'$\mathbf{R}_x^{(k)}  ={diag \Big( {(\mathbf{G}_x\;\mathbf{m}^{(k-1)})}^{2} + \epsilon_q^2 \Big)}^{\mathbf{q}/2 - 1} $','interpreter', 'latex','FontSize',10,'EdgeColor','k','BackgroundColor','w');
text(.45,-1.2,'(a)', 'interpreter', 'latex','FontSize',14);
% set(gca, 'XAxisLocation', 'top')
% set(gca, 'YAxisLocation', 'right')
ylabel('$z$', 'interpreter', 'latex','FontSize',16)
set(get(gca,'YLabel'),'Rotation',360);
else
   text(0.05,-0.1,'$\mathbf{R}_x^{(k)}  ={diag \Big( {(\mathbf{\nabla m}^{(k-1)})}^{2} + \epsilon_q^2 \Big)}^{\mathbf{q}/2 - 1} $','interpreter', 'latex','FontSize',10,'EdgeColor','k','BackgroundColor','w');
text(.45,-1.2,'(b)', 'interpreter', 'latex','FontSize',14); 
set(gca,'YTickLabel',[])
end


xlabel('$x$', 'interpreter', 'latex','FontSize',14)
axis square
% set(gca, 'XAxisLocation', 'top')
% set(gca, 'YAxisLocation', 'right')
end        

%% Add colorbar
ax = axes('Position',[0.375 0.3 .2 .5]);
cbar = colorbar('EastOutside');
set(cbar,'Ticks',[0 0.5 1])
set(cbar,'TickLabels',([-0.6 0 0.6]))
set(gca,'Visible','off');
text(1.25,-0.1,'$model$', 'interpreter', 'latex','FontSize',12)
            