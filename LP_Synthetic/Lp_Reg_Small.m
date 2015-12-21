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
p = [2 1 0];
beta = [0.01];

% Vector gradient location
xvec = lob:0.2:hib;
yvec = lob:0.2:hib;

[Xvec,Yvec] = ndgrid(xvec,yvec);

mref = [0. 0.];

axs  = 0.28;
%% Function
l2_func = @(r,x) sum((sqrt(r) * x').^2);

mis_func = @(x) ( (G * x').^2 )' - 2 * x * G' * d + d'*d;

lp_func = @(x,pp,eps) sum(x'.^2./(x'.^2 + eps.^2).^(1-pp/2),2);

wt_func = @(x,pp,eps) (1./(x'.^2 + eps.^2).^(1-pp/2));

grad_l2 = @(x,W) W*x';

grad_regfunc = @(x,pp,eps) (x'./(x'.^2 + eps.^2).^(1-pp/2));

grad_objfunc = @(x,reg) ((G'*G + reg)* x') - kron(G' * d,ones(1,size(x,1))) ;

%% Solve least-norm problem
lnorm = G' * (d/(G*G'));
x_0 = [lnorm(1) 0.6];
y_0 = [lnorm(2) 0.2];

%% Trace lines linking the different plots
set(figure, 'Position', [50 0 800 1000]);

%% Misfit function
misfit = mis_func([X(:) Y(:)]);
misfit = reshape(misfit,size(X,1),size(X,2));

% Compute gradient of function
UV = -grad_objfunc([Xvec(:) Yvec(:)],0);
U = reshape(UV(1,:),size(Xvec,1),size(Xvec,2));
V = reshape(UV(2,:),size(Xvec,1),size(Xvec,2));

% figure; imagesc(xx,yy,misfit); colorbar


% plot_h = tight_subplot(3,3,[.01 .01],[.1 0.1]);
axes('Position',[0.38 0.72 axs axs]);
% axes(plot_h(1));
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
xlabel('$m_1$','interpreter','latex','FontSize',12) 
ylabel('$m_2$','interpreter','latex','FontSize',12,'rot',360)
text(-.15,0.9,'$\phi_d = \|\mathbf{F \; m - d}\|_2^2$', 'interpreter', 'latex','FontSize',12,'BackgroundColor','w','EdgeColor','k');
text(0.9,-0.1,'(a)','HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor','w','Color','k','EdgeColor','k')


%% L2 regularization
for pp = 1:3
% smallest = l2_func([1 0;0 1],[X(:)-mref(1) Y(:)-mref(2)]);
% smallest = reshape(smallest,size(X,1),size(X,2));
reg = lp_func([X(:)';Y(:)'], p(pp) , epsl*100);
reg = reshape(reg,size(X,1),size(X,2));
% set(figure, 'Position', [50 200 500 500]); axis square

if pp == 1
    
    axes('Position',[0.075 0.4 axs axs]);
    title('$\phi_m = \sum_{i=1}^2 {m_i^2}$', 'interpreter', 'latex','FontSize',14,'BackgroundColor','w','EdgeColor','k');
    hold on
    xlabel('$m_1$','interpreter','latex','FontSize',12) 
    ylabel('$m_2$','interpreter','latex','FontSize',12,'rot',360)
elseif pp == 2
    
    axes('Position',[0.38 0.4 axs axs]);
    title('$\phi_m = \sum_{i=1}^2 \frac{m_i^2}{{(m_i^2 + \epsilon^2)}^{1/2}}$', 'interpreter', 'latex','FontSize',16,'BackgroundColor','w','EdgeColor','k');
    hold on
    xlabel('$m_1$','interpreter','latex','FontSize',12) 
    set(gca,'YTickLabel',[]);
else
    
    axes('Position',[0.68 0.4 axs axs]);
    title('$\phi_m = \sum_{i=1}^2 \frac{m_i^2}{{(m_i^2 + \epsilon^2)}}$', 'interpreter', 'latex','FontSize',16,'BackgroundColor','w','EdgeColor','k');
    hold on
    xlabel('$m_1$','interpreter','latex','FontSize',12)
    set(gca,'yaxislocation','right');
end


% plot the origin for reference
plot([lob hib],[0 0],'k:');
plot([0 0],[lob hib],'k:');
UV = -grad_regfunc([Xvec(:) Yvec(:)],p(pp),epsl);
U = reshape(UV(1,:),size(Xvec,1),size(Xvec,2));
V = reshape(UV(2,:),size(Xvec,1),size(Xvec,2));

% plot location of least-norm
plot(x_0(1),y_0(1),'ko','MarkerFaceColor','k');hold on;
plot(x_0(2),y_0(2),'ko','MarkerFaceColor','w');
quiver(Xvec,Yvec,U,V,0.5,'LineWidth',1); hold on

contour(yy,xx,reg,10,'LineColor','k'); hold on
plot(x_0(1),y_0(1),'ko','MarkerFaceColor','k');hold on;
plot(x_0(2),y_0(2),'ko','MarkerFaceColor','w');

plot(0,0,'k*','LineWidth',1)

if pp == 1
    text(x_0(1)-0.1,y_0(1)+0.1,['[' num2str(x_0(1)) ' ; ' num2str(y_0(1)) ']'],...
        'BackgroundColor','w','EdgeColor','k')

    text(x_0(2)-0.1,y_0(2)+0.1,['[' num2str(x_0(2)) ' ; ' num2str(y_0(2)) ']'],...
        'BackgroundColor','w','EdgeColor','k')
end

% Plot vector gradient at location xval,yval
R = spdiags(wt_func([x_0(1),y_0(1)],p(pp),epsl),0,2,2);
uv = -grad_l2([x_0(1)-mref(1),y_0(1)-mref(2)],R); %sqrt(uv(1).^2 + uv(2).^2)
uv = uv/norm(uv);
quiver(x_0(1),y_0(1),uv(1),uv(2),0.1,'LineWidth',1,'Color','k','MaxHeadSize',1);

R = spdiags(wt_func([x_0(2),y_0(2)],p(pp),epsl),0,2,2);
uv = -grad_l2([x_0(2)-mref(1),y_0(2)-mref(2)],R); %sqrt(uv(1).^2 + uv(2).^2)
uv = uv/norm(uv);
quiver(x_0(2),y_0(2),uv(1),uv(2),0.1,'LineWidth',1,'Color','k','MaxHeadSize',1);


axis square
% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);
grid on
box on
% xlabel('(b)')
% title('\boldmath$\phi_m^{(1)}$', 'interpreter', 'latex','FontSize',16);
if pp == 1
    text(0.9,-0.1,'(b)','HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor','w','Color','k','EdgeColor','k')

elseif pp == 2
    
    text(0.9,-0.1,'(d)','HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor','w','Color','k','EdgeColor','k')

else
    text(0.9,-0.1,'(f)','HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor','w','Color','k','EdgeColor','k')

end
    
    
%% Regularized misfit
misfit_reg = misfit + beta*reg;
if pp == 1
    axes('Position',[0.075 0.075 axs axs]);
    title('$\phi_{\;l2} =\phi_d + \phi_m$', 'interpreter', 'latex','FontSize',14,'BackgroundColor','w','EdgeColor','k');
    hold on
    xlabel('$m_1$','interpreter','latex','FontSize',12) 
    ylabel('$m_2$','interpreter','latex','FontSize',12,'rot',360)

elseif pp == 2
    
    axes('Position',[0.38 0.075 axs axs]);
    title('$\phi_{\;l1} =\phi_d + \phi_m$', 'interpreter', 'latex','FontSize',14,'BackgroundColor','w','EdgeColor','k');
    hold on
    xlabel('$m_1$','interpreter','latex','FontSize',12) 
    set(gca,'YTickLabel',[]);
else
    
    axes('Position',[0.68 0.075 axs axs]);
    title('$\phi_{\;l0} = \phi_d + \phi_m$', 'interpreter', 'latex','FontSize',14,'BackgroundColor','w','EdgeColor','k');
    hold on
    xlabel('$m_1$','interpreter','latex','FontSize',12)
    set(gca,'yaxislocation','right');
end
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

    beta_in = 1;

    

    R = spdiags(wt_func([x0,y0],p(pp),epsl),0,2,2);
    phi_in = l2_func(R,[x0,y0]);

    scale = phi_out/phi_in;

    IRLS = l2_func(R,[X(:) Y(:)]);
    IRLS = reshape(IRLS,size(X,1),size(X,2))*scale;

    reg_objfunc = misfit + beta_in*IRLS;

%     [yval,xval] = find(reg_objfunc==min(min(reg_objfunc))); 
%     m_out = (G'*G + beta*R)\(G'*d);
%     x_iter = m_out(1);
%     y_iter = m_out(2);

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

%     R = spdiags(wt_func([x_iter,y_iter],2,epsl*100),0,2,2);

%     IRLS = l2_func(R,[X(:) Y(:)]);
%     IRLS = reshape(IRLS,size(X,1),size(X,2))*scale;

    if jj == 1
        
        contour(xx,yy,reg_objfunc,10,'LineColor','k','LineStyle','-'); hold on
%         quiver(Xvec,Yvec,U,V,0.5,'Color','b','LineWidth',1,'LineStyle','-');
        
    else

        contour(xx,yy,reg_objfunc,10,'LineColor','k','LineStyle','--'); hold on
%         quiver(Xvec,Yvec,U,V,0.5,'Color','b','LineWidth',1,'LineStyle','--');

    end
%     plot(x_iter,y_iter,'k^','MarkerSize',10);
%     text(x_iter-0.15,y_iter+0.1,['[' num2str(round(100*x_iter)/100) ' ; ' num2str(round(100*y_iter)/100) ']'],...
%         'BackgroundColor','w','EdgeColor','k')
%     text(x0+0.05,y0+(-1)^jj*0.1,['\beta=' num2str(beta_in)],...
%         'BackgroundColor','w','EdgeColor','k')
end

% phid_l2 = norm(G * [x_iter;y_iter] - d);

set(gca,'Ydir','normal')
grid on
box on
axis square

axis([lob hib lob hib])
% text(1.05,0.4,'\boldmath$l_2$', 'interpreter', 'latex','FontSize',14);
% title('\boldmath$\phi_d + \phi_m^{(n)}$', 'interpreter', 'latex','FontSize',16);

if pp == 1
    text(0.9,-0.1,'(c)','HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor','w','Color','k','EdgeColor','k')

elseif pp == 2
    
    text(0.9,-0.1,'(e)','HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor','w','Color','k','EdgeColor','k')

else
    text(0.9,-0.1,'(g)','HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor','w','Color','k','EdgeColor','k')

end
end

%% FIGURE 2 - Solve the systems and plot
set(figure, 'Position', [50 0 800 400]);
axs = axs*2;
for pp = 1 : 3
    
reg = lp_func([X(:)';Y(:)'], p(pp) , epsl*100);
reg = reshape(reg,size(X,1),size(X,2));

if pp == 1
    axes('Position',[-.075 0.2 axs axs]);
    title('$\phi_{\;l2} =\phi_d + \beta\;\phi_m^{(k)}$', 'interpreter', 'latex','FontSize',14,'BackgroundColor','w','EdgeColor','k');
    hold on
    xlabel('$m_1$','interpreter','latex','FontSize',12) 
    ylabel('$m_2$','interpreter','latex','FontSize',12,'rot',360)

elseif pp == 2
    
    axes('Position',[0.225 0.2 axs axs]);
    title('$\phi_{\;l1} =\phi_d + \beta\;\phi_m^{(k)}$', 'interpreter', 'latex','FontSize',14,'BackgroundColor','w','EdgeColor','k');
    hold on
    xlabel('$m_1$','interpreter','latex','FontSize',12) 
    set(gca,'YTickLabel',[]);
else
    
    axes('Position',[0.53 0.2 axs axs]);
    title('$\phi_{\;l0} =\phi_d + \beta\;\phi_m^{(k)}$', 'interpreter', 'latex','FontSize',14,'BackgroundColor','w','EdgeColor','k');
    hold on
    xlabel('$m_1$','interpreter','latex','FontSize',12)
    set(gca,'yaxislocation','right');
end

% % plot the origin for reference
% plot([lob hib],[0 0],'k:');hold on;
% plot([0 0],[lob hib],'k:'); 
% contour(xx,yy,reg,10,'LineColor','k');hold on; 
% 
% UV = -grad_regfunc([Xvec(:) Yvec(:)],p(pp),epsl);
% U = reshape(UV(1,:),size(Xvec,1),size(Xvec,2));
% V = reshape(UV(2,:),size(Xvec,1),size(Xvec,2));
% 
% % plot location of least-norm
% plot(x_0(1),y_0(1),'ko','MarkerFaceColor','k');hold on;
% plot(x_0(2),y_0(2),'ko','MarkerFaceColor','w');
% quiver(Xvec,Yvec,U,V,0.5,'LineWidth',1);
% 
% 
% 
% grid on
% axis square
% if pp == 1
% 
%     xlabel('(d)')
%     text(-.15,0.9,'$\sum_{i=1}^2 {(x_i^2 + \epsilon^2)}^{1/2}$', 'interpreter', 'latex','FontSize',12,'BackgroundColor','w','EdgeColor','k');
% 
% else
% 
%     xlabel('(g)')
%     text(-.15,0.9,'$\sum_{i=1}^2 {(x_i^2 + \epsilon^2)}$', 'interpreter', 'latex','FontSize',12,'BackgroundColor','w','EdgeColor','k');
% 
% end
% figure; imagesc(xx,yy,misfit + reg); colorbar
% figure; contour(xx,yy,misfit + reg,20,'LineColor','k');
% set(gca,'Ydir','normal')

%% Plot Base map
plot([lob hib],[0 0],'k:');hold on;
plot([0 0],[lob hib],'k:');
plot(x_0(1),y_0(1),'ko','MarkerFaceColor','k');hold on;
plot(x_0(2),y_0(2),'ko','MarkerFaceColor','w');
plot(0,0,'k*','LineWidth',1)
plot( [lob hib]  ,[(y_0(1) - (x_0(1) - lob)*eivec(1,2)/eivec(1,1)) ...
(y_0(1) - (x_0(1) - hib)*eivec(1,2)/eivec(1,1))] ,'r--');
plot( [0 x_0(1)]  ,[0 y_0(1)],'k--','LineWidth',1);

for jj = 1 : 2
        
    x0 = x_0(jj);
    y0 = y_0(jj);
    
    R = spdiags(wt_func([x0,y0],p(pp),epsl),0,2,2);

    scale = l2_func([1 0;0 1],[x0 y0])/l2_func(R,[x0 y0]);


    beta_in = 0.005;

    x_iter = x0;
    y_iter = y0;
    
    phi_out = l2_func([1 0;0 1],[x_iter y_iter]); 
    
%     phid = sum((G * [x_iter;y_iter] - d)^2);
    for ii = 1 : 20

%         if ii ~=1
%         beta_in = beta_in * (phid_l2 / norm(G*[x_iter;y_iter]-d+1e-4)+1e-4) ;
%         end
%         beta_in = beta_in/5;
        
        R = spdiags(wt_func([x_iter,y_iter],p(pp),epsl),0,2,2);
        phi_in = l2_func(R,[x_iter,y_iter]);

        scale = 1;%phi_out/phi_in;

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
         text(x_iter+0.05,y_iter+0.1,['[' num2str(round(10*x_iter)/10) ' ; ' num2str(round(10*y_iter)/10) ']'],...
        'BackgroundColor','w','EdgeColor','k')
    else
        quiver(x0,y0,[x_iter - x0],[y_iter - y0],0.9,'LineWidth',2,'Color','r','MaxHeadSize',0.25);
        
        if pp~=1
        plot(x_iter,y_iter,'k^','MarkerSize',8);
         text(x_iter-0.2,y_iter-0.1,['[' num2str(num2str(round(10*x_iter)/10)) ' ; ' num2str(num2str(round(10*y_iter)/10)) ']'],...
        'BackgroundColor','w','EdgeColor','k')
    
        end
    end
    
end

grid on
box on
set(gca,'Ydir','normal')
axis square
axis([lob hib lob hib])

if pp == 1
    text(0.9,0.9,'(a)','HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor','w','Color','k','EdgeColor','k')

elseif pp == 2
    
    text(0.9,0.9,'(b)','HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor','w','Color','k','EdgeColor','k')

else
    text(0.9,0.9,'(c)','HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor','w','Color','k','EdgeColor','k')

end

end


%% Plot for Doug - Value of norms with l2 and l0 along null(G)

% x1 = 0 : 0.0001 : 1;
% x2 = 0.5 - 0.5 * x1;
% 
% phil0 = lp_func([x1;x2],0,epsl);
% phil2 = lp_func([x1;x2],2,epsl);
% 
% phi_tot = phil0 + phil2;
% 
% gradl0 = max(grad_regfunc([x1;x2],0,epsl),2);
% gradl2 = max(grad_regfunc([x1;x2],2,epsl),2);
% 
% grad_tot = max(grad_regfunc([x1;x2],0,epsl) + grad_regfunc([x1;x2],2,epsl),2);
% figure; plot(x1,phil0);
% hold on
% plot(x1,phil2)
% plot(x1,phi_tot)
% 
% figure; semilogy(x1,gradl0);
% hold on
% semilogy(x1,gradl2)
% semilogy(x1,grad_tot)