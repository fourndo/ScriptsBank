% Dual lp-norm algorithm in 1D
% This program is a demo version produced for the course 
% UBC-CPSC 546: Numerical Optimization
% 
% The program generates a simple 1D toy model and a forward model operator.
% The goal is to solve an ill-pose inverse problem of the type Ax = b. A
% regularization term is added to stabilize the problem and introduce some
% soft constraint on the solution. An approximate lp-norm penality on the
% model and the gradients are applied for a range of lp-norm. Solutions for
% a range of combinaison of lp-norms are computed and dsiplayed.
%
% Sub-function required: 
% CGLS.m
%
% Professor: Michael Friedlander
% Author: Dominique Fournier
% Last update: April 15h, 2014

close all
clear all

% Welcome message
fprintf('Starting lp-norm algorithm\n')
fprintf('Author: Dominique Fournier\n')
fprintf('Last update: April15th, 2014\n')

% Set up 1D problem
x = 3.14/1000:3.14/1000:3.14;
x = x(:) * 1000;

%% Set up model size and plot
mcell=length(x);
dx=ones(1,mcell) * abs((min(x)-max(x))/mcell);

% Create 1D model: Square and sin waves
model = 0.4 * exp(-abs((x*4/1000-pi*2.1).^6)).*tanh((x*8/1000-pi*2.1)); 

% Square wave
model(780:860)=0.6;

% Triangular wave
model(175:225)=0.6*[0:0.02:1];model(225:275)=0.6*[1:-0.02:0];

% Vector
model=model(:);


mref = zeros(mcell,1) ;
figure;plot(x,model,'LineWidth',2);axis([x(1) x(end) -0.1 1]);hold on
plot(x,mref,'r','LineWidth',2);
xlabel('Distance')
ylabel('Amplitude')
title('\bfTrue model')

%% Generate kernel
nfreq = 50; %number of frequencies
kernel = zeros(nfreq*2,mcell); % Pre-allocate space

figure;
for ii = 1 : nfreq 
    
    kernel(ii,1:mcell) = exp(- x*2 / ( pi  * 2000)  * (ii-1)/10) .*cos (x /2000 * (ii-1));
%     + ...
%         exp(-0.25 * x2/(pi)  * (ii-1)) .*cos (0.25 * x2 *2 * (ii-1));
    plot(kernel(ii,:));hold on

end

for ii = 1 : nfreq 
    
    kernel(ii+nfreq,1:mcell) = exp(- (x(end) - x)*2 / (pi * 2000)  * (ii-1)/10) .*cos ((x(end) - x) /2000 * (ii-1));
%     + ...
%         exp(-0.25 * x2/(pi)  * (ii-1)) .*cos (0.25 * x2 *2 * (ii-1));
    plot(kernel(ii+nfreq,:),'r');hold on

end
xlabel('\bfDistance');
ylabel('\bfAmplitude');
title('\bf Kernel functions e^{(\omega x)} cos(\omega x)')

% Total number of data
ndata = size(kernel,1);

%% Generate raw data and corrupt with noise
data= kernel * model;

% Corrupt data with Gaussian noise
rand_noise = randn(ndata,1);
pct_noise=0.05; % Percent noise
floor = ones(ndata,1) * pct_noise * max(abs(data));
noise = floor .* rand_noise;

d = data + noise; % Add noise

figure; plot(data); hold on;
plot(d,'r*')
xlabel('Frequency');
ylabel('Amplitude');
legend('Data','Data + noise');
title('Data from Ax = B and noise')

% Create uncertainty weighting matrix
Wd = spdiags(1./floor,0,ndata,ndata);

% Normalize d and G by standard deviation:
G = Wd * kernel;

% Target data misfit \phi_d^*
target = ndata;

% Weight data with uncertainty
d = Wd * d;

% Derivative operator for gradient norm
Wx = spdiags([-dx(:).^-0.5 dx(:).^-0.5],[0 1],mcell-1,mcell);

% Mass matrix for model norm
Ws = speye(mcell);

% Global constantclc
as = 1.0 / min(dx) ^2;  %Smallnest term
ax = 1.0;               %Smoothness term                     

% Lp-norm parameters chosen by the user
% We do a nice spread for the demo

fprintf('A solution will be found for all the combination of:\n')
pnorm = [0 1 2]
qnorm = [0 1 2]
lvec = [0.5 1.5]

% Stabilizing term
epsilon = 10^-(10);

% Store all the final models
models_out = zeros(length(pnorm),length(qnorm),length(lvec),mcell);

  

objfunc = @(m,dm,a,phi,l) norm(G*(m + a * dm) - d) +...
    (m + a * dm)' * l * phi * (m + a * dm);

% Iterate over all the combination of norms as specified by the vectors
% pvec, qvec and lvec.
counter = 1;
set(figure(4), 'Position', [25 100 1800 900])

for ll= 1:length(lvec)

    for pp = 1:length(pnorm)
        
        for qq = 1:length(qnorm)
            
            % Message prompt
            head = ['lp: ' num2str(pnorm(pp)) ' lq: ' num2str(qnorm(qq)) ' psi: ' num2str(lvec(ll))];
            fprintf('Starting Inversion %i of %i\n',counter, (length(pnorm)*length(qnorm)*length(lvec)))
            fprintf('%s\n',head)
            
            invmod      = ones(mcell,1)*1e-4;       % Initial model       
                        
            phi_d       = sum((G*invmod - d).^2);   % Initial misfit
            phi_m       = [];
            phim        = [];
            phi         = [];
            D = speye(mcell);           % Active cell
            count=0; % Initiate iteration count
            
            beta = [];
            eta = [];
            CGiter(counter) = 0;
            tic
            
            while phi_d(end) > target
                
                    count=count+1;
                    
                    %First iteration, do a two norm on model and gradient
                    if count==1   

                        WxtWx = (Wx)' * (Wx);
                        WstWs = (Ws)' * (Ws);
                        
                        Rs = speye(mcell);
                        Rx = speye(size(Wx,1));
                        
                        % Scale for gradient norm
                        scale_x     = ax * abs( 2.0-lvec(ll) ); 
                        scale_s     = as * lvec(ll) ;
                        
                        % Starting Model objective function 
                        phim_start =  ax * WxtWx + as * WstWs;
                        
                        phim = phim_start; 
                        
                        % Initial trade-off parameter
                        beta = sum(sum(G.^2,1)) / sum(diag(phim)) * 1e+2;                
                        eta = beta(count);
                        
                    else

                        rs = 1./ ( abs(Ws * invmod) .^( 2-qnorm(qq) ) + epsilon );
                        rs = rs / max(rs) + 10^(-8);
                        
                        Rs = spdiags( rs.^ 0.5 ,0,mcell,mcell);
                        
                        rx = 1./ ( abs(Wx * invmod) .^( 2-pnorm(pp) ) + epsilon ) ;
                        rx = rx / max(rx) +  10^(-8);
                        
                        Rx = spdiags( rx .^0.5,0,mcell-1,mcell-1);

                        WxtRxWx = (Rx * Wx)' * (Rx * Wx);
                        WstRsWs = (Rs * Ws)' * (Rs * Ws);

                        scale_s = as * (invmod' * WstWs * invmod ) /...
                            ( invmod' * WstRsWs * invmod );

                        scale_x = ax * (invmod' * scale_s * WstRsWs * invmod ) /...
                            ( invmod' * WxtRxWx * invmod );
                        
                        scale_s = lvec(ll) * scale_s ;
                        scale_x = (2.0 - lvec(ll)) * scale_x ;
                        
                        phim =  scale_s * WstRsWs + scale_x * WxtRxWx;
                        
                        eta = (invmod' * phim_start * invmod) / (invmod' * phim * invmod) ;
       
                        eta = beta(count) * eta;
                        
                    end

                phi(count) = norm(G*invmod - d) +...
                    invmod' * eta * phim * invmod;
                
                % Generate approximated Hessian
                A= [ G ;...
                    sqrt( eta * scale_s ) * ( Rs * Ws )  ;...
                    sqrt( eta * scale_x ) * ( Rx * Wx ) ];

                % Generate gradient (RHS)
                b = [- (G *invmod - d) ; ...
                    -sqrt( eta * scale_s ) * ( Rs * Ws * (invmod - mref)) ;...
                    -sqrt( eta * scale_x ) * ( Rx * Wx * (invmod - mref))];

                %% Gauss-Newton step
                dm = zeros(mcell,1);
                [dm,r,iter] = CGLSQ( dm, A , b );
                
                % Count the total number of CG iterations per inversion
                CGiter(counter) = CGiter(counter) + iter;
                
                % Starting step length
                alpha = 1; 
                
                % Evaluate current objective function 
                temp = objfunc(invmod,dm,alpha,phim,eta);
               
                % Reduce step length until \phi decreases
                while temp > phi(count)
                    
                    temp = objfunc(invmod,dm,alpha,phim,eta);
                    alpha = 0.75 * alpha;
                    
                end
                      
                % Update model
                invmod = invmod + alpha * dm;
                

            %% Save iteration details
            phi_d(count) = sum((G*(invmod)-d).^2);
            phi_m(count) = (invmod)'*( eta * phim)*(invmod);
                      
            % Reduce trade-off parameter
            if phi_d(count) < target*10 && count~=1
                
              beta(count+1) = 0.9*beta(count);

            else
                
              beta(count+1) = 0.5*beta(count);
              
            end

            
            end

            % Plot Results
            subplot((length(pnorm)*length(qnorm)*length(lvec))/3,3,counter)
            plot(1:mcell,model,1:mcell,invmod,'r','LineWidth',2);axis([0 mcell -.1 1])
            name_it = ['\bf|dm|p = ' num2str(pnorm(pp)) ' , |m|q = ' num2str(qnorm(qq)) ' , scale = ' num2str(lvec(ll))];
            title(name_it)
            
            % Compute model error and post on graph
            errl1 = norm(model - invmod,1);
            text(50,0.75,['|err|_{1}: ' num2str(errl1)]);
            
            errl2 = norm(model - invmod,2);
            text(500,0.75,['|err|_{2}: ' num2str(errl2)]);
            
            steps(counter) = count;

                
            fprintf('End of lp inversion in: %f sec\n',toc)
            fprintf('Number of beta iterations: %i\n',count)
            fprintf('Number of CG iterations: %i of max: %i\n',CGiter(counter),count*400);
            fprintf('Final data misfit: %8.3e. \nFinal model misfit: %8.3e\n\n',phi_d(count),norm(model-invmod,2))
            
            counter = counter+1;
            
        end
    end
end
 