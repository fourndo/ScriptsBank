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
% addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\
% addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\Lp_norm_codes\Dual_norm_1D\functions
% addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\Lp_norm_codes\Dual_norm_1D\sliceomatic\
% Set up 1D problem

%% INPUT VARIABLES
% Lp-norm parameters
pvec = [0 0.1]%0:0.5:2;
qvec = [0:0.1:2]%0:0.5:2;
lvec = ones(1,1)*1.0%[0.25 0.5 1 1.5 1.75]%0.1:0.4:1.9;

% Inversion parameters
FLAG1 = 'lGRADml'; % Choose between GRADm | lGRADm|
FLAG2 = 'SMOOTH_MOD'; % Choose between SMOOTH_MOD | SMOOTH_MOD_DIF

% Model space
nx = 60;
nz = 60;

dx = 1 / nx;
dz = 1 / nz;

% Reference model
mref = zeros(nx*nz,1);

% Create observation locations
obsx = [zeros(1,11) ones(1,11)];
nobs = length(obsx);
obsz = [0:-0.1:-1 0:-0.1:-1];


dx = ones(nx,1) * dx;
dz = ones(nz,1) * dz;

x0 = 0;
z0 = 0;

xc = [x0 + cumsum(dx) - dx/2];
zc = [z0 - cumsum(dz) + dz/2];

% Topography (all ones if flat)
nullcell = ones(nz,nx);

% Kernel
nk = 15; %number of frequencies

decay = -0.25;

basis = 0.5;

% Data noise
amp_pct = 0.00;
floor_pct = 0.02;

%% Generate model

model = zeros(nz,nx);

% Center the anomaly
dmx = floor(0.1 / dx(1));
dmz = floor(0.1 / dz(1));

centerx = floor(nx/2);
centerz = floor(nz/3);

% Square bloc
model(centerz-5:centerz+5,centerx-5:centerx+5)= 0.5;

% Gaussian anomaly
zz = 2*centerz-12:2*centerz+12;
xx = centerx-12:centerx+12;

[ZZ,XX] = ndgrid(cumsum(dz(zz)) - sum([dz(zz);dz(zz(1))])/2,cumsum(dx(xx)) - sum([dx(xx);dx(xx(1))])/2);
r = sqrt(XX.^2 + ZZ.^2);
dr = sqrt(dx(1)^2 + dz(1)^2);

model(zz,xx)= 0.75 * exp(-abs(((r)*10).^2));


set(figure, 'Position', [50 200 750 750])
imagesc(xc,zc,model);
caxis([0 0.75]);
colorbar
set(gca,'YDir','normal')
% axis([-1.1 1.1 -1.1 0.1])
axis square

hold on
plot(obsx,obsz,'r^','MarkerFaceColor','y','LineWidth',2)
model=model(:);
%% SCRIPT STARTS HERE
%% % Generate kernel functions and depth weighting
xn = [x0;x0 + cumsum(dx)];
zn = [z0;z0 - cumsum(dz)];

[Zn,Xn] = ndgrid(zn,xn);



[Zc,Xc] = ndgrid(zc,xc);



[dZ , dX] = ndgrid(dz,dx);
% dZ = dZ(:);
% dX = dX(:);

mcell=nx*nz;

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
       

        R = sqrt( (Zc-obsz(jj)).^2 + (obsx(jj)-Xc).^2 );
        dR = sqrt( dZ.^2 + dX.^2 );
%         x = Xn - obsx(jj);
%         z = Zn - obsz(jj);

        a = decay * (ii) ;
        b = basis * 2*pi* (ii) ;
%         c = decay * (ii) ;
%         d = basis * 2*pi* (ii) ;

        const1 = a / (a^2 + b^2);
        const2 = b / (a^2 + b^2);
%         const3 = c / (c^2 + d^2);
%         const4 = d / (c^2 + d^2);

        expx = @(m) exp( a * abs(m) );
%         expz = @(m) exp( c * abs(m) );

        trig1 = @(m) cos(b*(m));       
        trig2 = @(m) sin(b*(m));


%         xupb = (const1 * expx(R+dR)) .* (const2*trig1(b,R+dR));
%         xlob = (const1 * expx(R)) .* (const2*trig1(b,R));
%         temp = (xupb - xlob);

        xupb = dXn(2:end,1:end-1);
        xlob = dXn(2:end,2:end);
        
        zlob = dZn(1:end-1,2:end);
        zupb = dZn(2:end,2:end);
        
        fupxupz = const1 * expx( sqrt( xupb.^2 + zupb.^2 ) ) .* trig1(sqrt( xupb.^2 + zupb.^2 ))+...
                const2 * expx( sqrt( xupb.^2 + zupb.^2 ) ) .* trig2(sqrt( xupb.^2 + zupb.^2 ));

        fupxloz = const1 * expx( sqrt( xupb.^2 + zlob.^2 ) ) .* trig1(sqrt( xupb.^2 + zlob.^2 ))+...
                const2 * expx( sqrt( xupb.^2 + zlob.^2 ) ) .* trig2(sqrt( xupb.^2 + zlob.^2 ));

        floxupz = const1 * expx( sqrt( xlob.^2 + zupb.^2 ) ) .* trig1(sqrt( xlob.^2 + zupb.^2 ))+...
                const2 * expx( sqrt( xlob.^2 + zupb.^2 ) ) .* trig2(sqrt( xlob.^2 + zupb.^2 ));

        floxloz = const1 * expx( sqrt( xlob.^2 + zlob.^2 ) ) .* trig1(sqrt( xlob.^2 + zlob.^2 ))+...
                const2 * expx( sqrt( xlob.^2 + zlob.^2 ) ) .* trig2(sqrt( xlob.^2 + zlob.^2 ));
% 
        temp = (fupxupz - fupxloz) - (-floxupz + floxloz);
            
%             if sign(obsx(jj))<1
%                 
%                 temp= fliplr(temp);
%                 
%             end
            
        G(count,:) = temp(:);

    %     G(ii,1:mcell) = G(ii,1:mcell) / max(G(ii,1:mcell));

%         temp = expx(R) - expx(R+dR);
       
%         temp = ( expx( sqrt( xupb.^2 + zupb.^2 ) ) -...
%             expx( sqrt( xupb.^2 + zlob.^2 ) ) )/a -...
%             ( expx( sqrt( xlob.^2 + zlob.^2 ) ) -...
%             expx( sqrt( xlob.^2 + zupb.^2 ) ) )/a;
        
%          temp = (( expx(x(2:end,1:end-1))/a .* expz(z(1:end-1,2:end))/c )-...
%            ( expx(x(2:end,2:end))/a .* expz(z(2:end,2:end))/c ) -...
%            ( expx(x(2:end,1:end-1))/a .* expz(z(1:end-1,2:end))/c ) +...
%            ( expx(x(2:end,2:end))/a .* expz(z(2:end,2:end))/c )).^2;

%        if sign(obsx(jj))<1
%                 
%                 temp= fliplr(temp);
%                 
%        end
            
%         wr = wr + temp(:);
        
%         set(figure(3), 'Position', [10 250 1000 500]);
%         imagesc(xc,zc,reshape(G(count,:),nz,nx));
%         colorbar
%         set(gca,'YDir','normal')


        count = count + 1;

        
    end
    
end

wr = sqrt(sum(G.^2,1));

wr = (wr(:)).^0.5;


%% FORWARD MODELING
% Generate data and corrupt with noise

data= G * model;

% Corrupt data with Gaussian noise
% rand_noise = rand(ndata,1);
% save ('noise','rand_noise')
load('noise');

amp_noise = amp_pct * rand_noise .* abs(data);
floor_noise = floor_pct * rand_noise * max(abs(data));

noise = amp_noise + floor_noise;

d = data + noise;

% wd =  abs(d)*0.02 + 0.05 * std(d);
wd =  ones(ndata,1) * 0.1*std(d);

set(figure, 'Position', [10 250 1000 500]);
for ii = 1:nobs
    
    plot(data((1+(ii-1)*nk):(nk+(ii-1)*nk)));
    hold on
    errorbar(d((1+(ii-1)*nk):(nk+(ii-1)*nk)),wd((1+(ii-1)*nk):(nk+(ii-1)*nk)),'r*')
    legend('Data','Data+noise')
    grid on
end
% Create uncertainty weighting matrix
Wd = spdiags(1./wd,0,ndata,ndata);

% Normalize d and G by standard deviation:
G = Wd * G ;

target = ndata;

% Weight data with uncertainty
d = Wd * d;

%% INVERSION
% Create inversion parameters
[ Wx, Wz, Vx, Vz ] = get_GRAD_op2D(dx,dz,nullcell);
Ws = speye(mcell);
V = spdiags(sqrt(dX(:).*dZ(:)),0,mcell,mcell);


% Global constant
alpha(1) = 1.0 / min(dx) ^2;  %Smallnest term
alpha(2) = 1.0;               %Smoothness term                     
alpha(3) = 1.0;



% Store all the final models
RMS = zeros(length(pvec),length(qvec),length(lvec));
linf = zeros(length(pvec),length(qvec),length(lvec));
l2 = zeros(length(pvec),length(qvec),length(lvec));
l1 = zeros(length(pvec),length(qvec),length(lvec));
misfit = zeros(length(pvec),length(qvec),length(lvec));
models_out = zeros(length(pvec),length(qvec),length(lvec),mcell);
phi_m = [];
phi = [];
counter = 1;  

objfunc = @(m,phi,b) sum( ( G * m - d ).^2 ) + ( m' * b * phi * m );

% Iterate over all the combination of norms as specified by the vectors
% pvec, qvec and lvec.

for ll= 1:length(lvec)

    for pp = 1:length(pvec)
        
        for qq = 1:length(qvec)
            
            
%             delta(counter) = 1e-5;
            % Message prompt
            head = ['lp: ' num2str(pvec(pp)) ' lq: ' num2str(qvec(qq)) ' psi: ' num2str(lvec(ll))];
            fprintf('Starting lp inversion %s\n',...
                head)
            
            invmod      = ones(mcell,1)*1e-2;       % Initial model       
                        
            phi_d       = sum((G*invmod - d).^2);   % Initial misfit
           
            count=0; % Initiate iteration count 
            switcher = 0;
            ncg(counter) = 0;
            lp_count = 0;
            cg_iter = 0;
            delta = 1;
            
        while switcher ~= 3
                
                    count=count+1;

                    if switcher == 0   %First iteration
                        delta(count) = 1;
                        
                        [MOF,aVRWs,aVRWx,aVRWz] = get_lp_MOF_2D_v2(invmod,mref,1,nx,nz,V,Ws,Vx,Wx,Vz,Wz,wr,alpha,2,2,1,FLAG1,FLAG2,delta(count));                                      
%                         [MOF,aVRWs,aVRWx,aVRWz] = get_lp_MOF_2D(invmod,mref,nx,nz,V,Ws,Vx,Wx,Vz,Wz,wr,alpha,2,2,1,FLAG1,FLAG2,delta(count));

                        if count==1
                            % Initial beta trace(G'G) / trace(phim)
                            beta = sum(sum(G.^2,1)) / sum(diag(MOF,0).^2) * 1e+2; 
                            phi(count) = objfunc(invmod,MOF,beta(count));
                            
                        end
                       
                        
                    else

                        fprintf('# # LP-LQ ITER# #');
                        lp_count = lp_count+1;
                        
                        if cg_iter(end) < 0.5*mcell && switcher == 1
                            
                            delta(count) = 10^-(lp_count/2);
                            
                        else
                            
                            switcher = 2;
                            delta(count) = delta(count-1);
                            
                        end
                        
                        [MOF,aVRWs,aVRWx,aVRWz] = get_lp_MOF_2D_v2(invmod,mref,phi_m(end),nx,nz,V,Ws,Vx,Wx,Vz,Wz,wr,alpha,pvec(pp),qvec(qq),lvec(ll),FLAG1,FLAG2,delta(count));                         
%                         [MOF,aVRWs,aVRWx,aVRWz] = get_lp_MOF_2D(invmod,mref,nx,nz,V,Ws,Vx,Wx,Vz,Wz,wr,alpha,pvec(pp),qvec(qq),lvec(ll),FLAG1,FLAG2,delta(counter));
       
%                         lambda(count) = beta(count) * scale(count);
                        
                        
                    end

                
 
%                 eigval = eig(G'*G + lambda(count) * (phim) ); 
%                 condnum = max(abs(eigval)) / min(abs(eigval));
%                 
%                 figure(4); plot(eigval,'*');
%                 fprintf('Condition number: %f\n',condnum);
                fprintf('\n# # # # # # # # #\n');
                fprintf('BETA ITER: \t %i  \nbeta: \t %8.5e \n',count,beta(count));

               
                m_in = invmod;  % Save initial model for backtracking steps
                
            [invmod,cg_iter(count)] = GNsolver( G, invmod, d, phi(end), MOF, beta(count) , aVRWs, aVRWx, aVRWz, [] );                
            ncg(counter) = ncg(counter) + cg_iter(count);
%% Save iteration details
               
            rdm(count) = norm(m_in - invmod) / norm(m_in);
                
               
            phi_d(count) = sum((G*(invmod)-d).^2);
            phi_m(count) = (invmod)'*(MOF)*(invmod);
            phi(count) = objfunc(invmod,MOF,beta(count));
            
            if (rdm(count) < 1e-2 && switcher >= 1 && phi_d(count) > target *0.95 && phi_d(count) < target *1.05) || count > 50 

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
            if phi_d(count) < target * 0.95
                                      
                
                fprintf('---------->')
                fprintf(' misfit:\t %8.5e ',phi_d(count))
                fprintf(' ** Relative dm:\t %8.5e ', rdm(count));
                fprintf('<----------\n')
                fprintf('\n# NEXT ITER - INCREASING BETA #\n');
                
%                 beta(count) = beta(count-1);
                if switcher == 0
                    switcher = 1;
                end
%                 beta(count+1) = beta(count) * target / phi_d(count);
                beta(count+1) = beta(count) * 1.25;
%                 invmod = m_in;
%                 phi_d(count) = phi_d(count-1);
%                 phi_m(count) = phi_m(count-1);
%                 phi(count) = phi(count-1);
%                 
%                 count = count-1;
%                 lp_count = lp_count-1;

            elseif phi_d(count) > target * 1.05
                 
                
                fprintf('---------->')
                fprintf(' misfit:\t %8.5e ',phi_d(count))
                fprintf(' ** Relative dm:\t %8.5e ', rdm(count));
                fprintf('<----------\n')
                fprintf('\n# NEXT ITER - REDUCING BETA #\n');
                if switcher == 0
                    
                    beta(count+1) = 0.75 * beta(count); 
                    
                else
                    
%                     beta(count+1) = beta(count) * target / phi_d(count);
                    beta(count+1) = beta(count) * 0.9;
                    
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
            set(figure(3+counter), 'Position', [50 200 750 750])
            imagesc(xc,zc,reshape(invmod,nz,nx));
            caxis([0 0.75]);
            colorbar
            set(gca,'YDir','normal')
            axis square

            
        end
%             set(figure(mod(counter,10)+2), 'Position', [50 200 750 750])
%             imagesc(xc,zc,reshape(invmod,nz,nx));
%             caxis([0 0.75]);
%             colorbar
%             set(gca,'YDir','normal')
%             axis square
%             
% %             legend('\bfTrue model','\bf Wr || m ||_0 + Wr || \nabla m ||_0');
            xlabel('\bfX');
            ylabel('\bfZ');
            tlabel=(['$\mathbf{\phi_m\;=\;\| m \|_{' num2str(pvec(pp)) '}\;+\; \| \nabla m \|_{' num2str(qvec(qq)) '}}$']);
            title(tlabel, 'interpreter', 'latex','Fonts',18)
            
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
            RMS(pp,qq,ll)=(sum((invmod-model).^2)/length(data))^(0.5);
            linf(pp,qq,ll) = norm(model-invmod,'inf');
            l2(pp,qq,ll) = norm(model-invmod,2);
            l1(pp,qq,ll) = norm(model-invmod,1);
            misfit(pp,qq,ll) = phi_d(count);
            models_out(pp,qq,ll,:) = invmod(:);
            
%             figure; plotyy(1:count,phi_d,1:count,phi_m);

            ldml(counter) = norm(model-invmod,1);
            fprintf('End of lp inversion. Number of iterations: %i\n',count)
            fprintf('Final model misfit: %8.3e\n',ldml(counter))
            fprintf('TOTAL NUMBER OF CG STEPS: %i\n\n',ncg(counter))
            
            
            counter = counter+1;
            
        end
    end
end

save('l2_A','l2');
save('l1_A','l1');
save('misfit_A','misfit');
save('models_out_A','models_out');

%% Plot result 
% sliceomatic(l1,pvec,qvec,lvec);
% Plot tradeoff between model error and CG iterations

% figure;
% [hAx,hLine1,hLine2] = plotyy(delta,rdm,delta,cg_iter,'loglog','loglog');
% ylabel(hAx(2),'\bfCG iterations');
% ylabel(hAx(1),'\bf|m_{k-1} - m_{k}|/m_{k}');
% xlabel('\bf Factor log(\epsilon)')


% figure;
% [hAx,hLine1,hLine2] = plotyy(delta,ncg,delta,ldml,'semilogx','semilogx');
% ylabel(hAx(1),'\bfNumber of CG steps');
% ylabel(hAx(2),'\bf|m_{true} - m_{k}|_1');
% xlabel('\bf Factor (\epsilon)')
% set(hAx(2),'YColor','k')
% set(hLine2,'LineWidth',2,'Color','k')
% set(hLine1,'LineWidth',2)

% grid on

% [I,J]=find(l2==min(min(l2)));
% 
% set(figure, 'Position', [50 200 750 750]);
% imagesc(pvec,qvec,l2);colorbar;colormap winter;hold on
% contour(pvec,qvec,l2,'ShowText','on','LineColor',[0 0 0])
% set(gca,'YDir','Normal')
% ylabel('$\mathbf{\| m \|_p}$', 'interpreter', 'latex','Fonts',15)
% xlabel('$\mathbf{\| \nabla m \|_q}$', 'interpreter', 'latex','Fonts',15)
% plot(qvec(J),pvec(I),'ro','LineWidth',2);hold on
% axis square
% title('\textbf{Model error} ($\mathbf{\| m - m^{*} \|_2}$)', 'interpreter', 'latex','Fonts',15)
% 
% % load('misfit');
% set(figure, 'Position', [50 200 750 750]);
% imagesc(pvec,qvec,misfit);colorbar;colormap winter;hold on
% % contour(pvec,qvec,misfit,'ShowText','on','LineColor',[0 0 0])
% set(gca,'YDir','Normal')
% ylabel('$\mathbf{\| m \|_p}$', 'interpreter', 'latex','Fonts',15)
% xlabel('$\mathbf{\| \nabla m \|_q}$', 'interpreter', 'latex','Fonts',15)
% plot(qvec(J),pvec(I),'ro','LineWidth',2);hold on
% axis square
% title('\textbf{Final data misifit} $\phi_d$', 'interpreter', 'latex','Fonts',15)
% % caxis([38 42])
% 
% [I,J]=find(l1==min(min(l1)));
% set(figure, 'Position', [50 200 750 750]);
% imagesc(pvec,qvec,l1);colorbar;colormap winter;hold on
% % contour(pvec,qvec,misfit,'ShowText','on','LineColor',[0 0 0])
% set(gca,'YDir','Normal')
% ylabel('$\mathbf{\| m \|_p}$', 'interpreter', 'latex','Fonts',15)
% xlabel('$\mathbf{\| \nabla m \|_q}$', 'interpreter', 'latex','Fonts',15)
% plot(qvec(J),pvec(I),'ro','LineWidth',2);hold on
% axis square
% title('\textbf{Final data misifit} $\phi_d$', 'interpreter', 'latex','Fonts',15)