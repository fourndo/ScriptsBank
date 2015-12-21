% Test for eigenvalue
% Solve a problem of the form Ax = b + e
% Analyse eigenvalues and convergence of CG

clear all
close all
addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB

x = [0.1;1];

mcell = length(x);

% Forward operator
G = [1 1.2;
     1e-1 5e-1];

d = G*x;

ndata = length(d);

noise = rand(ndata,1) * 0.02 * min(d);

d = d + noise;

wd = 1 ./ (abs(d) * 0.1);
Wd = spdiags(wd,0,ndata ,ndata );

G = Wd * G;
d = Wd * d;

% Normal equation
GtG = G'* G;

 % Regularization
B = [1 0
     0 1];

% Map function
x1 = -0.5:0.2:0.5;
x2 = 0:0.2:1.2;

nx = length(x1);
ny = length(x2);

% Objective function calculator
objfunc = @(m,A) m' * (A) * m - 2*m'*(G'*d) + d'*d;

% Create objective function to solve underdetermine system
% Solve optimization problem ||Ax - b|| + I

% Solve the system for x
phid = 10;
betaB = 1000;
invxB = ones(mcell,1) * 1e-3;
count = 1;
% figure
while phid > ndata 
    

    
    H = GtG + betaB(count)*B;
    
    [invxB(:,count+1),~,~]=CGiter(invxB(:,end), H , G'*d);
    phid = norm( d - G*invxB(:,end) ) / length(x);
    
    invzB(count+1) = objfunc(invxB(:,end),H);
    
    tempx = invxB(1,count)-0.2:0.02:invxB(1,count)+0.2;
    tempy = invxB(2,count)-0.2:0.02:invxB(2,count)+0.2;
    
    [X1,X2] = ndgrid(tempx,tempy);
    X1  = X1(:);
    X2  = X2(:);
    X = [X1 X2]; X = X';

    func = zeros(size(X,2),1);
    func2 = zeros(size(X,2),1);
    for ii = 1 : size(X,2)
        
        func(ii) = objfunc(X(:,ii),H) ;

        func2(ii) = objfunc(X(:,ii),GtG) ;

    end
    
    func = reshape(func,length(tempx),length(tempy));
    func2 = reshape(func,length(tempx),length(tempy));
    
    if mod(count,1)==0
    figure(1)
    surf(tempy,tempx,func); hold on
    surf(tempy,tempx,func2); hold on
    view(65,20)
    xlim([0 1.0]);
    ylim([-.5 0.5]);
    axis square 
    plot3(invxB(2,count),invxB(1,count),(invzB(count)),'bo','MarkerFaceColor','w','MarkerSize',10,'LineWidth',2);hold on
    plot3(invxB(2,count+1),invxB(1,count+1),(invzB(count+1)),'bo','MarkerFaceColor','r','MarkerSize',10,'LineWidth',2);hold on
   
    hold off
    end
    
    betaB(count+1) = betaB(count) * 0.75;
    
    count = count + 1;
%     surf(tempy,tempx,func);hold on
%     view(-90,20)
%     axis equal
    
end
phim = (invxB(:,end)'*B*invxB(:,end));

fprintf('Final model B --> x1: %f , x2 %f\n',invxB(1,end),invxB(2,end))
fprintf(' Phi_d: %f \n',phid)
eps = 5e-3;
iter_start = 4;
 
% Solve the system for x lp_norm
phid = 1;
betaC = betaB(iter_start);
invxC = invxB(:,iter_start);
invzC = invzB(iter_start);
count = 1;
dm = ones(mcell,1);
ddm = [1 1];
rddm = 1;
while rddm > 1e-3 || phid > ndata * 1.25 || phid < ndata * 0.8;
    

    
    % Re-compute regularization
    r = 1 ./ (invxC(:,end).^2 + eps^2);
    Rc = spdiags(r,0,mcell,mcell);
    
    H = GtG + betaC(count)*Rc;
    
    [dm,~,~]=CGiter(invxC(:,end), H , -(H)*invxC(:,end) + G'*d);
    
    ddm(2) = norm(dm);
    
    tempx = invxC(1,end)-0.5:0.01:invxC(1,end)+0.5;
    tempy = invxC(2,end)-0.5:0.01:invxC(2,end)+0.5;
    
    invxC(:,count+1) = invxC(:,end) + dm;
    
    phid = norm( d - G * invxC(:,end) ) / length(x);
    
    invzC(count+1) = objfunc(invxC(:,end),H);
    
    if count==1%mod(count,2)==0   
    [X1,X2] = ndgrid(tempx,tempy);
    X1  = X1(:);
    X2  = X2(:);
    X = [X1 X2]; X = X';

    func = zeros(size(X,2),1);
    func2 = zeros(size(X,2),1);
    func3 = zeros(size(X,2),1);
    for ii = 1 : size(X,2)
        
        % Re-compute regularization
        r = 1 ./ (X(:,ii).^2 + eps^2);
        rc = spdiags(r,0,mcell,mcell);
    
        
        func(ii) = objfunc(X(:,ii),GtG + betaC(count)*rc) ;

        func2(ii) = objfunc(X(:,ii),GtG + betaC(count)*Rc ) ;

        func3(ii) = objfunc(X(:,ii),GtG) ;
        
    end
    
    func = reshape(func,length(tempx),length(tempy));
    func2 = reshape(func2,length(tempx),length(tempy));
    func3 = reshape(func3,length(tempx),length(tempy));
    
    limx = 1500;
    func2(func2>limx) = NaN;
%     func3(func3>200) = NaN;
%     func = X' * ((H) * X) ;
%     func = (diag(func) - 2*X'*(G'*d) + d'*d);
%     func = reshape(func,length(tempx),length(tempy));
%     
%     func2 = X' * ((GtG) * X) ;
%     func2 = (diag(func2) - 2*X'*(G'*d) + d'*d);
%     func2 = reshape(func2,length(tempx),length(tempy));
%     
%     func3 = X' * ((betaC(count) * Rc) * X) ;
%     func3 = (diag(func3) - 2*X'*(G'*d) + d'*d);
%     func3 = reshape(func3,length(tempx),length(tempy));
    
    
    figure(2)
    surf(tempy,tempx,func/limx); hold on; lighting phong
    surf(tempy,tempx,func2/limx); hold on; lighting phong
    surf(tempy,tempx,func3/limx); hold on; lighting phong
    view(100,5)
    xlim([invxC(1,count)-1 invxC(1,count)+1]);
    ylim([invxC(2,count)-1 invxC(2,count)+1]);
    axis equal 
    
    plot3(invxC(2,count),invxC(1,count),objfunc(invxC(:,count),H)/limx,'bo','MarkerFaceColor','w','MarkerSize',10,'LineWidth',2);hold on
    plot3(invxC(2,count+1),invxC(1,count+1),invzC(count+1)/limx,'bo','MarkerFaceColor','r','MarkerSize',10,'LineWidth',2);hold on
    hold off
    end
    
    if phid < ndata*0.75
        
        betaC(count+1) = betaC(count) * 1.25;
     
    elseif phid > ndata*1.25
        
        betaC(count+1) = betaC(count) * 0.8;
        
    else
        
        betaC(count+1) = betaC(count);
        
    end
    
    % Mesure relative change in dm
    if count == 1

        rddm = 1;
        ddm(1) = ddm(2);

    else

        rddm = ddm(2)/ddm(1);

    end
    
    count = count + 1;
        
end

fprintf('Final model C --> x1: %f , x2 %f\n',invxC(1,end),invxC(2,end))
fprintf(' Phi_d: %f \n',phid)
% Solve the system for x using normalized lp_norm
betaD = betaB(iter_start);
invxD = invxB(:,iter_start);
invzD = invzB(iter_start);
count = 1;
dm = ones(mcell,1);
ddm = [1 1];
rddm = 1;
while rddm > 1e-3 || phid > ndata * 1.25 || phid < ndata * 0.8;
    
    
    
    % Re-compute regularization
    r = 1 ./ (invxD(:,end).^2 + eps^2);
    Rd = spdiags(r,0,mcell,mcell);
     
    % Scale regularization on phim2
    scale = ( invxD(:,end)'*B*invxD(:,end) ) / ( invxD(:,end)'*Rd*invxD(:,end) );
    
    H = GtG + scale * betaD(count) * Rd;
    
    [dm,~,~]=CGiter(invxD(:,end), H , -(H)*invxD(:,end) + G'*d);
    
    ddm(2) = norm(dm);
    
    tempx = invxD(1,count)-0.5:0.01:invxD(1,count)+0.5;
    tempy = invxD(2,count)-0.5:0.01:invxD(2,count)+0.5;
    
    invxD(:,count+1) = invxD(:,end) + dm;
    
    phid = norm( d - G * invxD(:,end) ) / length(x);
    
    invzD(count+1) =  objfunc(invxD(:,end),H);
    
   
    if count==2%mod(count,1)==0
    [X1,X2] = ndgrid(tempx,tempy);
    X1  = X1(:);
    X2  = X2(:);
    X = [X1 X2]; X = X';

    func = zeros(size(X,2),1);
    func2 = zeros(size(X,2),1);
    func3 = zeros(size(X,2),1);
    for ii = 1 : size(X,2)
        
        % Re-compute regularization
        r = 1 ./ (X(:,ii).^2 + eps^2);
        rd = spdiags(r,0,mcell,mcell);
%         ss = ( X(:,ii)'*B*X(:,ii) ) / ( X(:,ii)'*rd*X(:,ii) );

        func(ii) = objfunc(X(:,ii),GtG) ;

        func2(ii) = objfunc(X(:,ii),GtG + scale*betaD(count)*rd) ;

        func3(ii) = objfunc(X(:,ii),GtG + scale*betaD(count)*Rd) ;
        
    end
    
    func = reshape(func,length(tempx),length(tempy));
    func2 = reshape(func2,length(tempx),length(tempy));
    func3 = reshape(func3,length(tempx),length(tempy));
    limx = 300;
    func(func>limx) = NaN;
    func2(func2>limx) = NaN;
    func3(func3>limx) = NaN;
%     func = X' * ((H) * X) ;
%     func = (diag(func) - 2*X'*(G'*d) + d'*d);
%     func = reshape(func,length(tempx),length(tempy));
%     
%     func2 = X' * ((GtG) * X) ;
%     func2 = (diag(func2) - 2*X'*(G'*d) + d'*d);
%     func2 = reshape(func2,length(tempx),length(tempy));
%     
%         
%     func3 = X' * ((GtG + betaD(count) * B) * X) ;
%     func3 = (diag(func3) - 2*X'*(G'*d) + d'*d);
%     func3 = reshape(func3,length(tempx),length(tempy));
    
    
    figure(3)
    
    surf(tempy,tempx,func/limx); hold on;
    surf(tempy,tempx,func2/limx); hold on;
    surf(tempy,tempx,func3/limx); hold on; 
    view(100,5)
    xlim([invxD(1,count)-1 invxD(1,count)+1]);
    ylim([invxD(2,count)-1 invxD(2,count)+1]);
    axis equal 
    
    plot3(invxD(2,count),invxD(1,count),objfunc(invxD(:,count),H)/limx,'bo','MarkerFaceColor','w','MarkerSize',10,'LineWidth',2);hold on
    plot3(invxD(2,count+1),invxD(1,count+1),invzD(count+1)/limx,'bo','MarkerFaceColor','r','MarkerSize',10,'LineWidth',2);hold on
   
    hold off
    end
    
    if phid < ndata * 0.8
        
        betaD(count+1) = betaD(count) * 1.25;
     
    elseif phid > ndata * 1.25
        
        betaD(count+1) = betaD(count) * 0.8;
        
    else
        
        betaD(count+1) = betaD(count);
        
    end
    
    % Mesure relative change in dm
    if count == 1

        rddm = 1;
        ddm(1) = ddm(2);

    else

        rddm = ddm(2)/ddm(1);

    end
    
    phim = (invxD(:,end)'*Rd*invxD(:,end));
    count = count + 1;
    
end
fprintf('Final model B --> x1: %f , x2 %f\n',invxD(1,end),invxD(2,end))
fprintf(' Phi_d: %f \n',phid)
% 
% cmap1 = zeros(nx,ny,3);
% cmap1(:,:,1) = abs(phi3)/max(abs(phi3(:))); 
% cmap1(:,:,2) = abs(phi3)/max(abs(phi3(:))); 
% cmap1(:,:,3) = 1; 
% 
% cmap2 = zeros(nx,ny,3);
% cmap2(:,:,1) = 1; 
% cmap2(:,:,2) = abs(phi4)/max(abs(phi4(:))); 
% cmap2(:,:,3) = abs(phi4)/max(abs(phi4(:))); 
% 
% figure;
% % surf(x1,x2,func1);
% surf(x2,x1,phi3,cmap1);hold on
% surf(x2,x1,phi4,cmap2);
% view(-200,40)
% quiver3(x(2),x(1),-1,0,0,1,'LineWidth',2)
% title('C vs norm(C)');
% 
