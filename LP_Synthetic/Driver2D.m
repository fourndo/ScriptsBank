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
addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\
% addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\Lp_norm_codes\Dual_norm_1D\functions
% addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\Lp_norm_codes\Dual_norm_1D\sliceomatic\
% Set up 1D problem

%% INPUT VARIABLES
% Lp-norm parameters
% lp = 2%0:0.1:2;
% lq = 0%0:0.1:2;
% mu = ones(1,1)*1.0%[0.25 0.5 1 1.5 1.75]%0.1:0.4:1.9;

% Inversion parameters
FLAG1 = 'GRADm'; % Choose between GRADm | lGRADml
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
pct_cutoff = 75;

% Rotation angle 
theta=0;
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
% zone{1,1} = 10:29; zone{1,2}=10:50; 
% zone{2,1} = 32:50; zone{2,2}=10:50; 
% zone{3,1} = []; zone{3,2}=[];
% LP = [1 0 2 1;
%       1 2 0 1;
%       2 2 2 1];

zone{1,1} = 1:nx; zone{1,2}=1:nz; 
LP = [0 0.5 0.5 1];
    
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
    lqz_mod =  lqz_mod + t(:,jj)*LP(jj,3);
    
    if jj==1
        axes('Position',[-0.185+1/3.25 .165 .25 .25])    
        imagesc(xc,zc,reshape(t(:,jj),nz,nx));
axis square
text(-.4,-0.5,'(b)', 'interpreter', 'latex','FontSize',14);
text(-.25,-0.5,'$z$', 'interpreter', 'latex','FontSize',14);
text(.45,-1.1,'$\mathbf{T}_{1}$', 'interpreter', 'latex','FontSize',14);
        
    elseif jj == 2
    axes('Position',[-0.225+2/3.25 .165 .25 .25])    
    imagesc(xc,zc,reshape(t(:,jj),nz,nx));
    axis square
%     colorbar('SouthOutside')
    set(gca,'YTickLabel',[]);
    xlabel('$x$', 'interpreter', 'latex','FontSize',16);
    text(.45,-1.1,'$\mathbf{T}_{2}$', 'interpreter', 'latex','FontSize',14);
    

    else
        
        axes('Position',[-0.265+3/3.25 .120 .34 .34])    
    imagesc(xc,zc,reshape(t(:,jj),nz,nx));
    axis square
    colorbar('EastOutside')
        set(gca,'YTickLabel',[]);
        set(gca, 'YAxisLocation', 'right')
        text(.45,-1.1,'$\mathbf{T}_{3}$', 'interpreter', 'latex','FontSize',14);
    end
    set(gca, 'XAxisLocation', 'top')
    set(gca,'YDir','normal')
    
end
          
            
axes('Position',[-0.185+1/3.25 .6 .25 .25])
imagesc(xc,zc,reshape(lp_mod,nz,nx));
axis square
% colorbar('SouthOutside')
caxis([0 2])
set(gca, 'XAxisLocation', 'top')
set(gca,'YDir','normal')
text(.43,-1.1,'$p$', 'interpreter', 'latex','FontSize',14);
text(-.4,-0.5,'(a)', 'interpreter', 'latex','FontSize',14);
text(-.25,-0.5,'$z$', 'interpreter', 'latex','FontSize',14);

axes('Position',[-0.225+2/3.25 .6 .25 .25])
imagesc(xc,zc,reshape(lqx_mod,nz,nx));
axis square
% set(gca,'YTickLabel',[]);
% colorbar('SouthOutside')
set(gca, 'XAxisLocation', 'top')
set(gca,'YTickLabel',[]);
set(gca,'YDir','normal')
text(0.45,-1.1,'$q_x$', 'interpreter', 'latex','FontSize',14);
% text(0.98,-1.5,'(b)', 'interpreter', 'latex','FontSize',14);
xlabel('$x$', 'interpreter', 'latex','FontSize',16);

axes('Position',[-0.265+3/3.25 .555 .34 .34])
imagesc(xc,zc,reshape(lqz_mod,nz,nx));
axis square
set(gca,'YTickLabel',[]);
colorbar('EastOutside')
set(gca, 'XAxisLocation', 'top')
set(gca,'YDir','normal')
text(0.45,-1.1,'$q_z$', 'interpreter', 'latex','FontSize',14);

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
            
colormap(jet)
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
% rand_noise = rand(ndata,1);
% save ('noise','rand_noise')
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
[Ws, gx, gz, V, Vx, Vz ] = get_GRAD_op2D_SQUARE_kron(dx,dz,nullcell,'FWR');
[~, gx_b, gz_b, ~,~,~ ] = get_GRAD_op2D_SQUARE_kron(dx,dz,nullcell,'BACK');
% Ws = speye(mcell);
% V = spdiags(sqrt(dX(:).*dZ(:)),0,mcell,mcell);

%% Rotate gradient operators

Rz = @(x)   [cosd(x) -sind(x);
            sind(x) cosd(x)];
               
 
Rot = Rz(theta);


rxx = spdiags(ones(mcell,1)*Rot(1,1),0,mcell,mcell);
rxy = spdiags(ones(mcell,1)*Rot(1,2),0,mcell,mcell);
ryx = spdiags(ones(mcell,1)*Rot(2,1),0,mcell,mcell);
ryy = spdiags(ones(mcell,1)*Rot(2,2),0,mcell,mcell);


Wx = rxx * gx + rxy * gz_b;
Wz = ryx * gx + ryy * gz;

% Show gradient
% mm = zeros(nz,nx);
% mm(round(nz/2),round(nx/2)) = 1;
% 
% figure; 
% axes('Position',[0.05 .25 .5 .5])
% imagesc(xc,zc,reshape(Wx*mm(:),nz,nx)); hold on
% scatter(xc(round(nx/2)),zc(round(nz/2)),'ro')
% caxis([-1 1]);
% colorbar('EastOutside')
% set(gca,'YDir','normal')
% axis square
% 
% axes('Position',[0.5 .25 .5 .5])
% imagesc(xc,zc,reshape(Wz*mm(:),nz,nx)); hold on
% scatter(xc(round(nx/2)),zc(round(nz/2)),'ro')
% caxis([-1 1]);
% % colorbar('EastOutside')
% set(gca,'YDir','normal')
% axis square

%% Global constant
alpha(1) = 0.5 / min(dx) ^2;  %Smallnest term
alpha(2) = 1.0;               %Smoothness term                     
alpha(3) = 1.0;



% Store all the final models
% RMS = zeros(length(lp),length(lq),length(mu));
% linf = zeros(length(lp),length(lq),length(mu));
% l2 = zeros(length(lp),length(lq),length(mu));
% l1 = zeros(length(lp),length(lq),length(mu));
% misfit = zeros(length(lp),length(lq),length(mu));
% models_out = zeros(length(lp),length(lq),length(mu),mcell);
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
beta_count = 0;
cg_iter = 0;
delta_p = 1;
traffic_s = 1;
traffic_x = 1;
group_s = 0;
p_tresh = 1;
x_tresh = 1;
switch_p = 0;
switch_q = 0;
dphi_m = 100;
gamma = 1;
Pac = speye(mcell);
set(figure(3), 'Position', [50 0 775 1000]);   
while switcher ~= 3

            count=count+1;

            if switcher == 0   %First iteration
                delta_p(count) = 0.5;
                delta_q(count) = 0.5;
                
                [MOF,aVRWs,aVRWx,aVRWz] = get_lp_MOF_2D(invmod,mref,1,s,t,V,Ws,Vx,Wx,Vz,Wz,wr,alpha,[2 2 2 1],FLAG1,FLAG2,switcher,delta_p(count));                                      
%                         [MOF,aVRWs,aVRWx,aVRWz] = get_lp_MOF_2D(invmod,mref,nx,nz,V,Ws,Vx,Wx,Vz,Wz,wr,alpha,2,2,1,FLAG1,FLAG2,delta(count));

                if count==1
                    % Initial beta trace(G'G) / trace(phim)
                    beta = sum(sum(G.^2,1)) / sum(diag(MOF,0).^2) * 1e+1; 
                    phi(count) = objfunc(invmod,MOF,beta(count));

                end

                p_tresh = 0.01;
                x_tresh = 0.01;

            else

                fprintf('# # LP-LQ ITER# #');
                if switcher ~=2
                            
                    lp_count = lp_count+1;

                else

                    beta_count = beta_count+1;

                end


                
                if lp_count == 1
                    
%                     [eps_p,eps_q] = get_eps_v2(invmod,lp_mod,lqx_mod,Wx,Wz,[]);
                    [eps_p,eps_q] = get_eps(invmod,10,Wx,Wz,[]);
                    eps_p = 1e-4;
                    eps_q = 1e-4;
                    p_tresh = eps_p;%0.01;%prctile(abs(invmod(invmod > 0)),pct_cutoff);
                    
%                     gradm = Wx * invmod;
%                     eps_x = std(gradm);
                    x_tresh = eps_q;%prctile(abs(gradm(gradm > 0)),pct_cutoff);
                    
%                     figure; hist(abs(invmod),20);
%                     ylim([0 100])
%%
                    figure(3)
                    axes('Position',[0.075 .57 .405 .405])
                    imagesc(xc,zc,reshape(invmod,nz,nx));
                    hold on;
%                     pp = prctile(invmod,[10 25 50 75 90]);
%                     pp_cut = prctile(invmod,pct_cutoff);
%                     contour(xc,zc,reshape(invmod,nz,nx),[-x_tresh/2 -x_tresh/5 x_tresh/5 x_tresh/2],'Color','w');
                    contour(xc,zc,reshape(invmod,nz,nx),[-p_tresh p_tresh],'Color','r','LineStyle','--');
                    caxis([-.6 0.6]);
%                     colorbar('EastOutside')
                    set(gca,'YDir','normal')
                    axis square
                    set(gca, 'XAxisLocation', 'top')
                    colormap(jet)
%                     xlabel('$x$', 'interpreter', 'latex','FontSize',16);
                    text(.40,-1.1,'Model', 'interpreter', 'latex','FontSize',14);
                    text(.45,-1.2,'(a)', 'interpreter', 'latex','FontSize',14);
%                     text(0.025,-0.125,['\underline{$L_2$-norm Solution}' char(10)...
%                         '$\mathbf{\phi_d} = ' num2str(round(phi_d(end)*100)/100) '$' char(10)...
%                         '$\|\delta m\|_1$ = ' num2str(round(norm(model-invmod,1)*100)/100) ],...
%                         'interpreter', 'latex','FontSize',10,'EdgeColor','k','BackgroundColor','w')
%                     text(0.05,-0.10,['$\mathbf{\phi_d} = ' num2str(round(phi_d(end)*100)/100) '$'],'interpreter', 'latex','FontSize',10)
                    ylabel('$z$', 'interpreter', 'latex','FontSize',16)
                    set(get(gca,'YLabel'),'Rotation',360);
                    xlabel('$x$', 'interpreter', 'latex','FontSize',14)
                    axis square
                            
                end
                
%                 if delta_p(end)> 1e-8;%traffic_p(end)*100 > 1 && switcher == 1
%                             
%                             delta_p(count) = delta_p(end)/2;
%                             
%                         else 
%                             
%                             delta_p(count) = delta_p(count-1);
%                             
%                 end
                        
                if delta_p(end)> eps_p %&& switcher == 1%(dphi_p(end) > 2 || lp_count ==1) && switch_p == 0
                            
                    switch_p = 1;
                    delta_p(count) = delta_p(count-1)*.5;

                    if delta_p(count) < eps_p
                        
                        delta_p(count) = eps_p;
                        
                    end
                else 

                    delta_p(count) = delta_p(count-1);%delta_p(count-1);

                end

                if delta_q(end)> eps_q %&& switcher == 1%(dphi_q(end) > 2 || lp_count ==1) && switch_q == 0

                    delta_q(count) = delta_q(count-1)*.5;
                    
                    if delta_q(count) < eps_q
                        
                        delta_q(count) = eps_q;
                        
                    end

                else

                    switch_q = 1;
                    delta_q(count) = delta_q(count-1);%delta_q(count-1);

                end

                if dphi_m(end) < 1  && delta_p(count) == eps_p && delta_q(count) == eps_q 

%                             delta_q(count) = eps_q;
%                             delta_p(count) = eps_p;
                    switcher = 2;

                end

                [MOF,aVRWs,aVRWx,aVRWz,gamma] = get_lp_MOF_2D(invmod,mref,phi_m(count-1),s,t,V,Ws,Vx,Wx,Vz,Wz,wr,alpha,LP,FLAG1,FLAG2,switcher,delta_p(count),delta_q(count));                         
%                         [MOF,aVRWs,aVRWx,aVRWz] = get_lp_MOF_2D(invmod,mref,nx,nz,V,Ws,Vx,Wx,Vz,Wz,wr,alpha,pvec(pp),qvec(qq),lvec(ll),FLAG1,FLAG2,delta(counter));

%                         lambda(count) = beta(count) * scale(count);
%                         scale = 1 / max(l0norm(invmod,delta(count)));
%                         x_tresh = dkdt(LP(:,1),delta(count));

%                         polynom = polyn(scale,delta(count));
%                         x_tresh = max(roots(polynom));

            end


            %%
    if lp_count==1
        
        mm = sort(abs(invmod));
        gradmm = mm./(mm.^2 + delta_p(end).^2).^(1-LP(1,1)/2);
        axes('Position',[-0.21+1/3.4 .275 .27 .27])
        [n, xout] =hist(abs(invmod),100); hold off
        [h_plot,h1,h2] = plotyy(xout,n,mm,gradmm/max(gradmm),'bar','plot');  
        hold(h_plot(1),'on');
        hold(h_plot(2),'on');
        set(h2,'LineWidth',2)
                
%         ylim(h_plot(1),[0 100]);
        xlim(h_plot(1),[0 0.2])
        xlim(h_plot(2),[0 0.2])
        
        axis(h_plot(1),'square')
        axis(h_plot(2),'square')
%                 [n, xout] =hist(h_plot(1),abs(invmod),100);
        set(h1,'barwidth', 1, 'basevalue', 1,'FaceColor',[0.7 0.7 0.7],'LineWidth',0.5);

        plot(h_plot(2),[p_tresh p_tresh],[0 max(n)],'r--','LineWidth',2)

        % Change scale
        set(h_plot(1),'yscale','log')
        ylim(h_plot(1),[1 3000]);
        set(h_plot(1),'YTick',[1;10;100;1000]);
        
%         text(0.1,90,['\epsilon = ' num2str(delta(end))],'EdgeColor','k')
        set(h_plot(2),'YTickLabel',[]);
        xlabel(h_plot(1),'$|m|$', 'interpreter', 'latex','FontSize',12)
%         temp = title(['$\epsilon_p^{(1)}$=' num2str(delta_p(end)) '$\;:\;\epsilon_q^{(1)}$=' num2str(delta_q(end)) '$\;:\;\delta\phi_m^{(1)}$=' num2str(round(dphi_m(end))) '$\%$']);
temp = title(['IRLS Iteration: 0']);        
set(temp,'interpreter', 'latex','FontSize',13);
        

        
        ylabel('$Hist(m)$', 'interpreter', 'latex','FontSize',16)
        hold on
            
        
    elseif lp_count == 3
        mm = sort(abs(invmod));
        gradmm = mm./(mm.^2 + delta_p(end).^2).^(1-LP(1,1)/2);
        axes('Position',[-0.21+2/3.4 .275 .27 .27])
        [n, xout] =hist(abs(invmod),100); hold off
        [h_plot,h1,h2] = plotyy(xout,n,mm,gradmm/max(gradmm),'bar','plot');  
        hold(h_plot(1),'on');
        hold(h_plot(2),'on');
        set(h2,'LineWidth',2)
                
%         ylim(h_plot(1),[0 100]);
        xlim(h_plot(1),[0 0.2])
        xlim(h_plot(2),[0 0.2])
        
        axis(h_plot(1),'square')
        axis(h_plot(2),'square')
%                 [n, xout] =hist(h_plot(1),abs(invmod),100);
        set(h1,'barwidth', 1, 'basevalue', 1,'FaceColor',[0.7 0.7 0.7],'LineWidth',0.5);

        plot(h_plot(2),[p_tresh p_tresh],[0 max(n)],'r--','LineWidth',2)

        % Change scale
        set(h_plot(1),'yscale','log')
        ylim(h_plot(1),[1 3000]);
%         set(h_plot(1),'YTick',[1 10 100 1000]);
        
        xlabel(h_plot(1),'$|m|$', 'interpreter', 'latex','FontSize',12)

        text(0.075,1e-1,'$(c)$', 'interpreter', 'latex','FontSize',14)
        
        set(h_plot(2),'YTickLabel',[]);
        set(h_plot(1),'YTickLabel',[]);
%         temp = title(['$\epsilon_p^{(3)}$=' num2str(delta_p(end)) '$\;:\;\epsilon_q^{(3)}$=' num2str(delta_q(end)) '$\;:\;\delta\phi_m^{(3)}$=' num2str(round(dphi_m(end))) '$\%$']);
temp = title(['IRLS Iteration: 3']);        
set(temp,'interpreter', 'latex','FontSize',13);
        temp = legend('Hist','$\mathbf{\hat g_p}(m)$',['$\epsilon^*:\;$' num2str(round(eps_p(end)/10^floor(log10(eps_p(end))))) 'e' num2str(floor(log10(eps_p(end))))]) ;
        set(temp,'interpreter', 'latex','FontSize',10)
        hold on

    elseif lp_count == 6
        
        mm = sort(abs(invmod));
        gradmm = mm./(mm.^2 + delta_p(end).^2).^(1-LP(1,1)/2);
        axes('Position',[-0.21+3/3.4 .275 .27 .27])
        [n, xout] =hist(abs(invmod),100); hold off
        [h_plot,h1,h2] = plotyy(xout,n,mm,gradmm/max(gradmm),'bar','plot');  
        hold(h_plot(1),'on');
        hold(h_plot(2),'on');
        set(h2,'LineWidth',2)
                
%         ylim(h_plot(1),[0 100]);
        xlim(h_plot(1),[0 0.2])
        xlim(h_plot(2),[0 0.2])
        
        axis(h_plot(1),'square')
        axis(h_plot(2),'square')
%                 [n, xout] =hist(h_plot(1),abs(invmod),100);
        set(h1,'barwidth', 1, 'basevalue', 1,'FaceColor',[0.7 0.7 0.7],'LineWidth',0.5);

        plot(h_plot(2),[p_tresh p_tresh],[0 max(n)],'r--','LineWidth',2)
        xlabel(h_plot(1),'$|m|$', 'interpreter', 'latex','FontSize',12)
        % Change scale
        set(h_plot(1),'yscale','log')
        ylim(h_plot(1),[1 3000]);
%         set(h_plot(1),'YTick',[1 10 100 1000]);
        set(h_plot(1),'YTickLabel',[])
%         text(0.2,0.2/(0.2^2 + delta(end).^2).^(1 - LP(1,1)/2),'$\mathbf{R^TRm}$', 'interpreter', 'latex','FontSize',10,'EdgeColor','k','BackgroundColor','w')
%         temp = title(['$\epsilon_p^{(5)}$=' num2str(delta_p(end)) '$\;:\;\epsilon_q^{(5)}$=' num2str(delta_q(end)) '$\;:\;\delta\phi_m^{(5)}$=' num2str(round(dphi_m(end))) '$\%$']);
        temp = title(['IRLS Iteration: 6']);        
        set(temp,'interpreter', 'latex','FontSize',13);
        ylabel(h_plot(2),'$\mathbf{\hat g_p}(m)$', 'interpreter', 'latex','FontSize',12)
%         set(get(h_plot(2),'YLabel'),'Rotation',360);
%         ylabh = get(h_plot(2),'YLabel');
%         set(ylabh,'Position',get(ylabh,'Position') + [0.025 0.05 0.00]);
%         set(gca, 'YLabelLocation', 'top')

        hold on
        
        
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

%     if switcher==1
%     [inv_temp,~,~] = GNsolver( G, invmod, d, phi(end), MOF, PreC, Pac, beta(count) , aVRWs, aVRWx, aVRWz, [] );                
%     
%     %% TEST Re-adjust beta on the fly
%     gamma =  norm(G'*(G*inv_temp-d))/norm(G'*(G*m_in-d));
%     beta(count) = gamma * beta(count);
%     
%     end
    
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
    phi_p(count) = (aVRWs*invmod)'*(aVRWs*invmod);
    phi_q(count) = (aVRWx*invmod)'*(aVRWx*invmod) + (aVRWz*invmod)'*(aVRWz*invmod);
    
    if count > 1
                dphi_m(count) = abs(phi_m(count) - phi_m(count-1))/phi_m(count) *100;
                dphi_p(count) = abs(phi_p(count) - phi_p(count-1))/phi_p(count) *100;
                dphi_q(count) = abs(phi_q(count) - phi_q(count-1))/phi_q(count) *100;
    end

    if switcher == 0
    grad_MOF = norm(MOF*invmod);
    end
    
    phi(count) = objfunc(invmod,MOF,beta(count));

    if (dphi_m(end) < 1 && switcher == 2 && phi_d(count) > target *0.95 && phi_d(count) < target *1.05) || count > 50 

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
    if phi_d(count) < target * 0.95 %&& switcher~=1


        fprintf('---------->')
        fprintf(' misfit:\t %8.5e ',phi_d(count))
        fprintf(' ** Relative dm:\t %8.5e ', rdm(count));
        fprintf('<----------\n')
        fprintf('\n# NEXT ITER - INCREASING BETA #\n');

%                 beta(count) = beta(count-1);
        if switcher == 0

            beta(count+1) = beta(count);
            switcher = 1;

% %                         x_tresh = dkdt(LP(:,1),delta(count));
%                     polynom = polyn(scale,delta(count));
%                     delta(count) = max(roots(polynom));%max(abs(invmod)).^((2-LP(:,1)));
        else

%              beta(count+1) = beta(count) * 1.1;
            beta(count+1) = beta(count) * target / phi_d(count);
        end
%                 

%                 invmod = m_in;
%                 phi_d(count) = phi_d(count-1);
%                 phi_m(count) = phi_m(count-1);
%                 phi(count) = phi(count-1);
%                 
%                 count = count-1;
%                 lp_count = lp_count-1;

    elseif phi_d(count) > target * 1.05 %&& switcher ~= 1


        fprintf('---------->')
        fprintf(' misfit:\t %8.5e ',phi_d(count))
        fprintf(' ** Relative dm:\t %8.5e ', rdm(count));
        fprintf('<----------\n')
        fprintf('\n# NEXT ITER - REDUCING BETA #\n');
        if switcher == 0

            beta(count+1) = 0.75 * beta(count); 

%         elseif rdm(count) < 1e-2
% 
%             beta(count+1) = beta(count) * target / phi_d(count);

        else

                    beta(count+1) = beta(count) * target / phi_d(count);
%             beta(count+1) = beta(count) * 0.95;

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
        figure(3)
        axes('Position',[0.5 .524 .495 .495])
        imagesc(xc,zc,reshape(invmod,nz,nx)); hold on
        caxis([-.6 0.6]);
        colorbar('EastOutside')
        set(gca,'YDir','normal')
        axis square
        
%         contour(xc,zc,reshape(invmod,nz,nx),pp,'Color','w');
%                     contour(xc,zc,reshape(invmod,nz,nx),[-x_tresh/2 -x_tresh/5 x_tresh/5 x_tresh/2],'Color','w');
                    contour(xc,zc,reshape(invmod,nz,nx),[-p_tresh p_tresh],'Color','r','LineStyle','--');
        set(gca, 'XAxisLocation', 'top')
        xlabel('$x$', 'interpreter', 'latex','FontSize',14)
        set(gca,'YTickLabel',[])
%         xlabel('$x$', 'interpreter', 'latex','FontSize',16);
        text(.40,-1.1,'Model', 'interpreter', 'latex','FontSize',14);
        text(.45,-1.2,'(b)', 'interpreter', 'latex','FontSize',14);
%         text(0.025,-0.15,['\underline{Mixed $L_p-norm$}' char(10) ...
%             '$\epsilon^{(k)}=\;' num2str(round(delta_p(end)/10^floor(log10(delta_p(end))))) 'e' num2str(floor(log10(delta_p(end)))) '$' char(10)...
%             '$\mathbf{\phi_d} = ' num2str(round(phi_d(end)*100)/100) '$' char(10)...
%             '$\|\delta m\|_1$ = ' num2str(round(norm(model-invmod,1)*1000)/1000) ],...
%             'interpreter', 'latex','FontSize',10,'EdgeColor','k','BackgroundColor','w')
%         text(0.05,-0.10,['$\mathbf{\phi_d} = ' num2str(round(phi_d(count)*100)/100) '$'],'interpreter', 'latex','FontSize',10)
%         text(0.05,-0.25,['$\epsilon = 10^{' num2str(round(log10(delta(end))*10)/10) '}$'],'interpreter', 'latex','FontSize',10)
%         text(0.05,-0.3,['$\|\delta m\|_1$ = ' num2str(round(norm(model-invmod,1)*100)/100) ],'interpreter', 'latex','FontSize',10)
            
        annotation('arrow',[0.2 0.2],[0.58 0.55],'LineWidth',5,'HeadWidth',10);
        annotation('arrow',[0.82 0.82],[0.55 0.58],'LineWidth',5,'HeadWidth',10);
        annotation('arrow',[0.33 0.37],[0.53 0.53],'LineWidth',5,'HeadWidth',10);
        annotation('arrow',[0.63 0.67],[0.53 0.53],'LineWidth',5,'HeadWidth',10);
        
        n_iter = count - (lp_count+beta_count);
        
        set(figure, 'Position', [50 0 775 1000]);
        plot_h = tight_subplot(2,1,[.1 .01],[.1 0.1]);
        
        axes(plot_h(1))
        h = plotyy(1:count,phi_d,1:count,phi_m,'semilogy'); hold on
        set(h(1),'LineWidth',1);
        set(h(2),'LineWidth',1);

        grid on
        xlim(h(1),[0 count])
        xlim(h(2),[0 count])
        ylim(h(1),[min(phi_d)*0.9 max(phi_d)*1.1])
        ylim(h(2),[min(phi_m)*0.9 max(phi_m)*1.1])
        axis(h(1),'square')
        axis(h(2),'square')
%             set(gca, 'YAxisLocation', 'right')

        xlabel('Iteration')

        ylabel(h(1),'$\phi_d$' ,'interpreter', 'latex','FontSize',16);
        set(get(h(1),'YLabel'),'Rotation',360);
        ylabh = get(h(1),'YLabel');
        set(ylabh,'Position',get(ylabh,'Position') - [2 0 0.00]);

        ylabh = get(h(2),'YLabel');
        set(ylabh,'Position',get(ylabh,'Position') + [2 0 0.00]);
        ylabel(h(2),'$\phi_m$' ,'interpreter', 'latex','FontSize',14);
        set(get(h(2),'YLabel'),'Rotation',360);


        xx = [0 0 n_iter n_iter 0];
        yy = [1 1e+4 1e+4 1 1];
        h = fill(xx,yy,'g');
        set(h,'FaceAlpha',0.1,'LineStyle','none');

        n_iter = n_iter + lp_count;
        xx = [0 0 n_iter n_iter 0];
        yy = [1 1e+4 1e+4 1 1];
        h = fill(xx,yy,'r');
        set(h,'FaceAlpha',0.1,'LineStyle','none');

        n_iter = n_iter + beta_count;
        xx = [0 0 n_iter n_iter 0];
        yy = [1 1e+4 1e+4 1 1];
        h = fill(xx,yy,'b');
        set(h,'FaceAlpha',0.1,'LineStyle','none');

        plot([0 count],[target target],'b--');
        text(0,target,'$\phi_d^*$','interpreter', 'latex','FontSize',12,'Color','b','BackgroundColor','w','EdgeColor','b')
        text(-7,180,'(d)', 'interpreter', 'latex','FontSize',14);
        
        %% ADD LEGEND
        axes('Position',[.375 .025 .28 .03])
        axis([0 3 0 1])

        h = fill([0 0 1 1 0],[0 1 1 0 0],'g');set(h,'FaceAlpha',0.1); hold on
        h = fill([0 0 2 2 0],[0 1 1 0 0],'r');set(h,'FaceAlpha',0.1);
        h = fill([0 0 3 3 0],[0 1 1 0 0],'b');set(h,'FaceAlpha',0.1);
        set(gca,'XTickLabel',[])
        set(gca,'YTickLabel',[])
%                 set(gca,'Visible','off');  

        text(0.05,0.4,'$$\phi_d^{(1)} \rightarrow \phi_d^*$','interpreter', 'latex','FontSize',12)
        text(1.05,0.4,'$\epsilon^{(1)} \rightarrow \epsilon^*$','interpreter', 'latex','FontSize',12)
        text(2.05,0.4,'$\phi_d^{(k)} \rightarrow \phi_d^*$','interpreter', 'latex','FontSize',12)
            
        %%
        
%         axes(plot_h(2))
%         plot(traffic);
        %% 
        figure(6);
        axis([0 0.1 0 100])
        xx = 1e-4:1e-4:max(abs(invmod));
        mm = xx./(xx.^(2-LP(1)) + delta_p(end));
        figure(6); hold on; plot(xx,mm/max(mm)*100,'LineWidth',2,'Color','r')
        var = ['LP: ' num2str(LP(1))];
        yy = get(gca,'Ylim');
        plot([p_tresh p_tresh],yy,'--')
        title(var)
        set(gca, 'XAxisLocation', 'top')

            xlabel('\bfX');
            ylabel('\bfZ');


            ldml(counter) = norm(model-invmod,1);
            fprintf('End of lp inversion. Number of iterations: %i\n',count)
            fprintf('Final model misfit: %8.3e\n',ldml(counter))
            fprintf('TOTAL NUMBER OF CG STEPS: %i\n\n',ncg(counter))
            
            
            counter = counter+1;
            



