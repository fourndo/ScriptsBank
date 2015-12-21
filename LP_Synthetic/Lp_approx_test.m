% Test approximation of l0 norm and plot objective function

close all
clear all

addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB


%% Set up problem
lob = -0.2;
hib = 1;

xx = lob : 1e-2 : hib;
yy = lob : 1e-2 : hib;

[Y,X] = ndgrid(xx,yy);


G = [1 2];
d = [1];

GtG = G'*G;
[eivec,eival] = eig(GtG);

epsl = 0.001;%10.^(-5:0);
p = [1 0];
beta = [0.01];

% Vector gradient location
xvec = lob:0.2:hib;
yvec = lob:0.2:hib;

[Xvec,Yvec] = ndgrid(xvec,yvec);

mref = [0. 0.];
%% Function
l2_func = @(r,x) sum((sqrt(r) * x').^2);

mis_func = @(x) ( (G * x').^2 )' - 2 * x * G' * d + d'*d;

lp_func = @(x,pp,eps) sum(x'.^2./(x'.^2 + eps.^2).^(1-pp/2));

wt_func = @(x,pp,eps) (1./(x'.^2 + eps.^2).^(1-pp/2));

grad_l2 = @(x,W) W*x';

grad_regfunc = @(x,pp,eps) (x'./(x'.^2 + eps.^2).^(1-pp/2));

grad_objfunc = @(x,reg) ((G'*G + reg)* x') - kron(G' * d,ones(1,size(x,1))) ;

%% Solve least-norm problem
lnorm = G' * (d/(G*G'));
x_0 = [lnorm(1) 0.6];
y_0 = [lnorm(2) 0.2];


%% Misfit function
misfit = mis_func([X(:) Y(:)]);
misfit = reshape(misfit,size(X,1),size(X,2));

% Compute gradient of function
UV = -grad_objfunc([Xvec(:) Yvec(:)],0);
U = reshape(UV(1,:),size(Xvec,1),size(Xvec,2));
V = reshape(UV(2,:),size(Xvec,1),size(Xvec,2));

% figure; imagesc(xx,yy,misfit); colorbar

set(figure, 'Position', [50 0 775 1000]);
plot_h = tight_subplot(3,3,[.01 .01],[.1 0.1]);

axes(plot_h(1));
% plot the origin for reference
plot([lob hib],[0 0],'k:');
plot([0 0],[lob hib],'k:');
contour(xx,yy,misfit,10,'LineColor','k'); hold on
plot(0,0,'k*','LineWidth',1)
quiver(Xvec,Yvec,U,V,0.5,'LineWidth',1);
set(gca,'Ydir','normal')

% Plot the null space of G 
plot( [lob hib]  ,[(y_0(1) - (x_0(1) - lob)*eivec(1,2)/eivec(1,1)) ...
   (y_0(1) - (x_0(1) - hib)*eivec(1,2)/eivec(1,1))] ,'--');

% Plot the least-norm axis
plot( [0 x_0(1)]  ,[0 y_0(1)],'k--','LineWidth',2);
plot(x_0(1),y_0(1),'ko','MarkerFaceColor','k');

set(gca,'XTickLabel',[]);
axis square
grid on
axis([lob hib lob hib])
xlabel('(a)')
text(-.15,0.9,'$\phi_d = \|\mathbf{F \; m - d}\|_2^2$', 'interpreter', 'latex','FontSize',12,'BackgroundColor','w','EdgeColor','k');
% quiver([0;0],[0;0],eivec(:,1),eivec(:,2),1,'LineWidth',2);

%% L2 regularization

smallest = l2_func([1 0;0 1],[X(:)-mref(1) Y(:)-mref(2)]);
smallest = reshape(smallest,size(X,1),size(X,2));
% set(figure, 'Position', [50 200 500 500]); axis square
axes(plot_h(2));
% plot the origin for reference
plot([lob hib],[0 0],'k:');
plot([0 0],[lob hib],'k:');
UV = -grad_l2([Xvec(:)-mref(1) Yvec(:)-mref(2)],[1 0;0 1]);
U = reshape(UV(1,:),size(Xvec,1),size(Xvec,2));
V = reshape(UV(2,:),size(Xvec,1),size(Xvec,2));
quiver(Xvec,Yvec,U,V,0.5,'LineWidth',1); hold on

contour(yy,xx,smallest,10,'LineColor','k'); hold on
plot(x_0(1),y_0(1),'ko','MarkerFaceColor','k');hold on;
text(x_0(1)-0.1,y_0(1)+0.1,['[' num2str(x_0(1)) ' ; ' num2str(y_0(1)) ']'],...
        'BackgroundColor','w','EdgeColor','k')
        
plot(x_0(2),y_0(2),'ko','MarkerFaceColor','w');
text(x_0(2)-0.1,y_0(2)+0.1,['[' num2str(x_0(2)) ' ; ' num2str(y_0(2)) ']'],...
        'BackgroundColor','w','EdgeColor','k')
plot(0,0,'k*','LineWidth',1)

% Plot vector gradient at location xval,yval
uv = -grad_l2([x_0(1)-mref(1),y_0(1)-mref(2)],[1 0;0 1]); %sqrt(uv(1).^2 + uv(2).^2)
uv = uv/norm(uv);
quiver(x_0(1),y_0(1),uv(1),uv(2),0.1,'LineWidth',1,'Color','k','MaxHeadSize',1);

uv = -grad_l2([x_0(2)-mref(1),y_0(2)-mref(2)],[1 0;0 1]); %sqrt(uv(1).^2 + uv(2).^2)
uv = uv/norm(uv);
quiver(x_0(2),y_0(2),uv(1),uv(2),0.1,'LineWidth',1,'Color','k','MaxHeadSize',1);


axis square
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
grid on
xlabel('(b)')
title('\boldmath$\phi_m^{(1)}$', 'interpreter', 'latex','FontSize',16);

%% Regularized misfit
misfit_reg = misfit + beta*smallest;
[yval2,xval2] = find(misfit_reg==min(min(misfit_reg))); mis_func([xx(xval2) yy(yval2)]);
% set(figure, 'Position', [50 200 500 500]); axis square
axes(plot_h(3));
% plot the origin for reference
plot([lob hib],[0 0],'k:');
plot([0 0],[lob hib],'k:');
plot(x_0(1),y_0(1),'ko','MarkerFaceColor','k');hold on;
plot(x_0(2),y_0(2),'ko','MarkerFaceColor','w');

% Plot null space of F
plot( [lob hib]  ,[(y_0(1) - (x_0(1) - lob)*eivec(1,2)/eivec(1,1)) ...
   (y_0(1) - (x_0(1) - hib)*eivec(1,2)/eivec(1,1))] ,'r--');
plot( [0 x_0(1)]  ,[0 y_0(1)],'k--','LineWidth',1);

for jj = 1 : 2

    x0 = x_0(jj);
    y0 = y_0(jj);
    
    phi_out = l2_func([1 0;0 1],[x0 y0]);

    beta_in = beta;

    

    R = spdiags(wt_func([x0,y0],2,epsl),0,2,2);
    phi_in = l2_func(R,[x0,y0]);

    scale = phi_out/phi_in;

    IRLS = l2_func(R,[X(:) Y(:)]);
    IRLS = reshape(IRLS,size(X,1),size(X,2))*scale;

    reg_objfunc = misfit + beta_in*IRLS;

%     [yval,xval] = find(reg_objfunc==min(min(reg_objfunc))); 
    m_out = (G'*G + beta*R)\(G'*d);
    x_iter = m_out(1);
    y_iter = m_out(2);

%     mis_func([x_iter y_iter])

%     UV = -grad_objfunc([Xvec(:) Yvec(:)],beta_in*scale*R);
%     U = reshape(UV(1,:),size(Xvec,1),size(Xvec,2));
%     V = reshape(UV(2,:),size(Xvec,1),size(Xvec,2));
    hold on
    plot(0,0,'k*','LineWidth',1)
%     quiver(Xvec,Yvec,U,V,0.5,'LineWidth',1);

    % plot the origin for reference
    plot([lob hib],[0 0],'k:');
    plot([0 0],[lob hib],'k:');



    phi_out = scale * phi_in;

    R = spdiags(wt_func([x_iter,y_iter],2,epsl*100),0,2,2);

    IRLS = l2_func(R,[X(:) Y(:)]);
    IRLS = reshape(IRLS,size(X,1),size(X,2))*scale;

    if jj == 1
        
        contour(xx,yy,misfit + beta_in*50*IRLS,10,'LineColor','k','LineStyle','-'); hold on


    end
    
    if jj == 2
    quiver(x0,y0,[x_iter - x0],[y_iter - y0],0.9,'LineWidth',2,'Color','r','MaxHeadSize',0.5);
    
    end
    plot(x_iter,y_iter,'k^','MarkerSize',10);
    text(x_iter-0.15,y_iter+0.1,['[' num2str(round(100*x_iter)/100) ' ; ' num2str(round(100*y_iter)/100) ']'],...
        'BackgroundColor','w','EdgeColor','k')
%     text(x0+0.05,y0+(-1)^jj*0.1,['\beta=' num2str(beta_in)],...
%         'BackgroundColor','w','EdgeColor','k')
end

phid_l2 = norm(G * [x_iter;y_iter] - d);

set(gca,'Ydir','normal')
grid on
axis square
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
axis([lob hib lob hib])
xlabel('(c)')
text(1.05,0.4,'\boldmath$l_2$', 'interpreter', 'latex','FontSize',14);
title('\boldmath$\phi_d + \phi_m^{(n)}$', 'interpreter', 'latex','FontSize',16);

%% L1

for pp = 1 : 2
    
reg = lp_func([X(:) Y(:)], p(pp) , epsl*100);
reg = reshape(reg,size(X,1),size(X,2));
% set(figure, 'Position', [50 200 500 500]); axis square
axes(plot_h(pp*3+1));
% plot the origin for reference
plot([lob hib],[0 0],'k:');hold on;
plot([0 0],[lob hib],'k:'); 
contour(xx,yy,reg,10,'LineColor','k');hold on; 

UV = -grad_regfunc([Xvec(:) Yvec(:)],p(pp),epsl);
U = reshape(UV(1,:),size(Xvec,1),size(Xvec,2));
V = reshape(UV(2,:),size(Xvec,1),size(Xvec,2));

% plot location of least-norm
plot(x_0(1),y_0(1),'ko','MarkerFaceColor','k');hold on;
plot(x_0(2),y_0(2),'ko','MarkerFaceColor','w');
quiver(Xvec,Yvec,U,V,0.5,'LineWidth',1);

if  pp == 1
    set(gca,'XTickLabel',[]);
    
end

grid on
axis square
if pp == 1

    xlabel('(d)')
    text(-.15,0.9,'$\sum_{i=1}^2 {(x_i^2 + \epsilon^2)}^{1/2}$', 'interpreter', 'latex','FontSize',12,'BackgroundColor','w','EdgeColor','k');

else

    xlabel('(g)')
    text(-.15,0.9,'$\sum_{i=1}^2 {(x_i^2 + \epsilon^2)}$', 'interpreter', 'latex','FontSize',12,'BackgroundColor','w','EdgeColor','k');

end
% figure; imagesc(xx,yy,misfit + reg); colorbar
% figure; contour(xx,yy,misfit + reg,20,'LineColor','k');
% set(gca,'Ydir','normal')

%% Plot Base map
axes(plot_h(pp*3+2));
% plot the origin for reference
plot([lob hib],[0 0],'k:');hold on;
plot([0 0],[lob hib],'k:');
plot(x_0(1),y_0(1),'ko','MarkerFaceColor','k');hold on;
% text(x_0(1)+0.05,y_0(1)+0.1,['[' num2str(x_0(1)) ' ; ' num2str(y_0(1)) ']'],...
%         'BackgroundColor','w','EdgeColor','k')
        
plot(x_0(2),y_0(2),'ko','MarkerFaceColor','w');
% text(x_0(2)+0.05,y_0(2)+0.1,['[' num2str(x_0(2)) ' ; ' num2str(y_0(2)) ']'],...
%         'BackgroundColor','w','EdgeColor','k')
    
plot(0,0,'k*','LineWidth',1)
if  pp == 1
    set(gca,'XTickLabel',[]);
    
end
axes(plot_h(pp*3+3));
% plot the origin for reference
plot([lob hib],[0 0],'k:');hold on;
plot([0 0],[lob hib],'k:');
plot(x_0(1),y_0(1),'ko','MarkerFaceColor','k');hold on;
plot(x_0(2),y_0(2),'ko','MarkerFaceColor','w');
plot(0,0,'k*','LineWidth',1)
plot( [lob hib]  ,[(y_0(1) - (x_0(1) - lob)*eivec(1,2)/eivec(1,1)) ...
(y_0(1) - (x_0(1) - hib)*eivec(1,2)/eivec(1,1))] ,'r--');
plot( [0 x_0(1)]  ,[0 y_0(1)],'k--','LineWidth',1);
    if  pp == 1
    set(gca,'XTickLabel',[]);
    
end
for jj = 1 : 2
        
    x0 = x_0(jj);
    y0 = y_0(jj);
    
    R = spdiags(wt_func([x0,y0],p(pp),epsl),0,2,2);

    scale = l2_func([1 0;0 1],[x0 y0])/l2_func(R,[x0 y0]);

    IRLS = l2_func(R,[X(:) Y(:)]);
    IRLS = reshape(IRLS,size(X,1),size(X,2))*scale;
      
    % set(figure, 'Position', [50 200 500 500]); axis square
    axes(plot_h(pp*3+2));
               
    UV = -grad_l2([Xvec(:) Yvec(:)],R);
    U = reshape(UV(1,:),size(Xvec,1),size(Xvec,2));
    V = reshape(UV(2,:),size(Xvec,1),size(Xvec,2));

    % Plot unit gradient at the starting model
    uv = -grad_l2([x0,y0],R); %uv*scale(1)
    uv = uv/norm(uv);
    quiver(x0,y0,uv(1),uv(2),0.1,'LineWidth',1,'Color','k','MaxHeadSize',1);
    
    if jj == 1
        
        contour(xx,yy,IRLS,10,'LineColor','k','LineStyle','-'); hold on
%         quiver(Xvec,Yvec,U,V,0.5,'Color','b','LineWidth',1,'LineStyle','-');
        
    else

        contour(xx,yy,IRLS,10,'LineColor','k','LineStyle','--'); hold on
%         quiver(Xvec,Yvec,U,V,0.5,'Color','b','LineWidth',1,'LineStyle','--');

    end
    

    set(gca,'Ydir','normal')
    grid on
    axis square
    set(gca,'YTickLabel',[]);
    
    if pp == 1
        
        xlabel('(e)')
        
    else
        
        xlabel('(h)')
        
    end

    % Iterate until convergence
    axes(plot_h(pp*3+3));
    
    
    plot([lob hib],[0 0],'k:');
    plot([0 0],[lob hib],'k:');

    beta_in = beta;

    x_iter = x0;
    y_iter = y0;
    
    phi_out = l2_func([1 0;0 1],[x_iter y_iter]); 
    
    for ii = 1 : 15

        if ii ~=1
        beta_in = beta_in * (phid_l2 / norm(G*[x_iter;y_iter]-d+1e-4)+1e-4) ;
        end
        
        R = spdiags(wt_func([x_iter,y_iter],p(pp),epsl),0,2,2);
        phi_in = l2_func(R,[x_iter,y_iter]);

        scale = phi_out/phi_in;

        IRLS = l2_func(R,[X(:) Y(:)]);
        IRLS = reshape(IRLS,size(X,1),size(X,2))*scale;

        reg_objfunc = misfit + beta_in*IRLS;

%         [yval,xval] = find(reg_objfunc==min(min(reg_objfunc))); 
        m_out = (G'*G + beta_in*scale*R)\(G'*d);
        x_iter = m_out(1);
        y_iter = m_out(2);

%         mis_func([x_iter y_iter])



        phi_out = scale * phi_in;

    end

     mis_func([x_iter y_iter])
     
%     plot(x_iter,y_iter,'k^','MarkerSize',8);
%     text(x0+0.05,y0+(-1)^jj*0.1,['\beta=' num2str(beta_in)],...
%             'BackgroundColor','w','EdgeColor','k')

%% FOR PLOTTING ONLY
    R = spdiags(wt_func([x_iter,y_iter],p(pp),epsl*100),0,2,2);

    IRLS = l2_func(R,[X(:) Y(:)]);
    IRLS = reshape(IRLS,size(X,1),size(X,2))*scale;

    if jj == 1
        
        contour(xx,yy,misfit + beta_in*10*IRLS,10,'LineColor','k','LineStyle','-'); hold on

    else

        contour(xx,yy,misfit + beta_in*10*IRLS,10,'LineColor','k','LineStyle','--'); hold on

    end
            
    if jj == 1
        quiver(x0,y0,[x_iter - x0],[y_iter - y0],0.9,'LineWidth',2,'Color','r','MaxHeadSize',0.5);
        plot(x_iter,y_iter,'k^','MarkerSize',8,'MarkerFaceColor','k');
         text(x_iter+0.05,y_iter+0.1,['[' num2str(round(100*x_iter)/100) ' ; ' num2str(round(100*y_iter)/100) ']'],...
        'BackgroundColor','w','EdgeColor','k')
    else
        quiver(x0,y0,[x_iter - x0],[y_iter - y0],0.9,'LineWidth',2,'Color','r','MaxHeadSize',0.25);
        
        if pp~=1
        plot(x_iter,y_iter,'k^','MarkerSize',8);
         text(x_iter-0.2,y_iter-0.1,['[' num2str(num2str(round(100*x_iter)/100)) ' ; ' num2str(num2str(round(100*y_iter)/100)) ']'],...
        'BackgroundColor','w','EdgeColor','k')
    
        end
    end
    
end

grid on
axis square
set(gca,'YTickLabel',[]);
axis([lob hib lob hib])
if pp == 1

    xlabel('(f)')
    text(1.05,0.4,'\boldmath$l_1$', 'interpreter', 'latex','FontSize',14);
    
else

    xlabel('(i)')
    text(1.05,0.4,'\boldmath$l_0$', 'interpreter', 'latex','FontSize',14);

end

end


