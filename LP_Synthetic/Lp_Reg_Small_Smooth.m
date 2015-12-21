% Test approximation of l0 norm and plot objective function

close all
clear all

addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB


%% Set up problem
lob = -0.5;
hib = 1.5;

xx = lob : 1e-2 : hib;
yy = lob : 1e-2 : hib;

[Y,X] = ndgrid(xx,yy);


G = [1e-5 1];
d = [1];

% G = [1 2];
% d = [1];

GtG = G'*G;
[eivec,eival] = eig(GtG);

epsl = 0.1;%10.^(-5:0);
p = [0];
beta = [0.005];

% Vector gradient location
xvec = lob:0.2:hib;
yvec = lob:0.2:hib;

[Xvec,Yvec] = ndgrid(xvec,yvec);

mref = [0. 0.];

axs  = 0.36;
%% Function
l2_func = @(r,x) sum((r * x').^2,1);

mis_func = @(x) ( (G * x').^2 )' - 2 * x * G' * d + d'*d;

lp_func = @(x,pp,eps) sum(x'.^2./(x'.^2 + eps.^2).^(1-pp/2));

wt_func = @(x,pp,eps) sqrt(1./(x'.^2 + eps.^2).^(1-pp/2));

grad_l2 = @(x,W) (W'*W)*x';

grad_regfunc = @(x,pp,eps) (x'./(x'.^2 + eps.^2).^(1-pp/2));

grad_objfunc = @(x,reg) ((G'*G + reg)* x') - kron(G' * d,ones(1,size(x,1))) ;

%% Solve least-norm problem
lnorm = G' * (d/(G*G'));
% x_0 = [lnorm(1) 1];
% y_0 = [lnorm(2) epsl];
x_0 = [0 1];
y_0 = [1 1];

%% Misfit function
misfit = mis_func([X(:) Y(:)]);
misfit = reshape(misfit,size(X,1),size(X,2));

% Compute gradient of function
UV = -grad_objfunc([Xvec(:) Yvec(:)],0);
U = reshape(UV(1,:),size(Xvec,1),size(Xvec,2));
V = reshape(UV(2,:),size(Xvec,1),size(Xvec,2));

% figure; imagesc(xx,yy,misfit); colorbar

set(figure, 'Position', [50 0 775 600]);

axes('Position',[0.025 0.56 axs axs]);

% plot the origin for reference
plot([lob hib],[0 0],'k:');
plot([0 0],[lob hib],'k:');
contour(xx,yy,misfit,15,'LineColor','k'); hold on
plot(0,0,'k*','LineWidth',1)
quiver(Xvec,Yvec,U,V,0.5,'LineWidth',1);
set(gca,'Ydir','normal')

% Plot the null space of G 
plot( [lob hib]  ,[(lnorm(2) - (lnorm(1) - lob)*eivec(1,2)/eivec(1,1)) ...
   (lnorm(2) - (lnorm(1) - hib)*eivec(1,2)/eivec(1,1))] ,'r--');

% Plot the least-norm axis
plot( [0 lnorm(1)]  ,[0 lnorm(2)],'k--','LineWidth',2);
plot(lnorm(1),lnorm(2),'ko','MarkerFaceColor','k');

set(gca,'XTickLabel',[]);
axis square
grid on
axis([lob hib lob hib])
xlabel('$m_1$','interpreter','latex','FontSize',12)
ylabel('$m_2$','interpreter','latex','FontSize',12)
title('$\phi_d = \|\mathbf{F \; m - d}\|_2^2$', 'interpreter', 'latex','FontSize',14);
% quiver([0;0],[0;0],eivec(:,1),eivec(:,2),1,'LineWidth',2);
text(1.25,-.25,'(a)','HorizontalAlignment','center','BackgroundColor','w','Color','k','EdgeColor','w','FontSize',12,'HorizontalAlignment','center','VerticalAlignment','middle')
%% L2 regularization
Wx = [1 -1];
% Wx= Wx'*Wx;

smoothness = l2_func(Wx,[X(:) Y(:)]);
smoothness = reshape(smoothness,size(X,1),size(X,2));
% set(figure, 'Position', [50 200 500 500]); axis square
axes('Position',[0.65 0.56 axs axs]);
% plot the origin for reference
plot([lob hib],[0 0],'k:');
plot([0 0],[lob hib],'k:');
UV = -grad_l2([Xvec(:)-mref(1) Yvec(:)-mref(2)],Wx);
U = reshape(UV(1,:),size(Xvec,1),size(Xvec,2));
V = reshape(UV(2,:),size(Xvec,1),size(Xvec,2));
quiver(Xvec,Yvec,U,V,0.5,'LineWidth',1); hold on

contour(yy,xx,smoothness,10,'LineColor','k'); hold on
% plot( [0 lnorm(1)]  ,[0 lnorm(2)],'k--','LineWidth',2);
plot(x_0(1),y_0(1),'ko','MarkerFaceColor','k');hold on;
      
plot(x_0(2),y_0(2),'ko','MarkerFaceColor','w');
plot(0,0,'k*','LineWidth',1)

% Plot vector gradient at location xval,yval
uv = -grad_l2([x_0(1)-mref(1),y_0(1)-mref(2)],Wx); %sqrt(uv(1).^2 + uv(2).^2)
uv = uv/norm(uv);
quiver(x_0(1),y_0(1),uv(1),uv(2),0.2,'LineWidth',1,'Color','k','MaxHeadSize',1);

uv = -grad_l2([x_0(2)-mref(1),y_0(2)-mref(2)],Wx); %sqrt(uv(1).^2 + uv(2).^2)
uv = uv/norm(uv);
quiver(x_0(2),y_0(2),uv(1),uv(2),0.25,'LineWidth',1,'Color','k','MaxHeadSize',1);



plot( [lob hib]  ,[(lnorm(2) - (lnorm(1) - lob)*eivec(1,2)/eivec(1,1)) ...
   (lnorm(2) - (lnorm(1) - hib)*eivec(1,2)/eivec(1,1))] ,'r--');

axis square
set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);
grid on
xlabel('$m_1$','interpreter','latex','FontSize',12)
title('\boldmath${\phi_x = \|G_x m\|}_2^2$', 'interpreter', 'latex','FontSize',14);
text(1.25,-.25,'(c)','HorizontalAlignment','center','BackgroundColor','w','Color','k','EdgeColor','w','FontSize',12,'HorizontalAlignment','center','VerticalAlignment','middle')


%% l2 Model norm
smallness = l2_func([1 0;0 1],[X(:)-mref(1) Y(:)-mref(2)]);
smallness = reshape(smallness,size(X,1),size(X,2));
% set(figure, 'Position', [50 200 500 500]); axis square
axes('Position',[0.335 0.56 axs axs]);
% plot the origin for reference
plot([lob hib],[0 0],'k:');
plot([0 0],[lob hib],'k:');
UV = -grad_l2([Xvec(:)-mref(1) Yvec(:)-mref(2)],[1 0;0 1]);
U = reshape(UV(1,:),size(Xvec,1),size(Xvec,2));
V = reshape(UV(2,:),size(Xvec,1),size(Xvec,2));
quiver(Xvec,Yvec,U,V,0.5,'LineWidth',1); hold on

contour(yy,xx,smallness,10,'LineColor','k'); hold on
plot(x_0(1),y_0(1),'ko','MarkerFaceColor','k');hold on;

% plot( [0 lnorm(1)]  ,[0 lnorm(2)],'k--','LineWidth',2);        
plot(x_0(2),y_0(2),'ko','MarkerFaceColor','w');

plot(0,0,'k*','LineWidth',1)

% Plot vector gradient at location xval,yval
uv = -grad_l2([x_0(1)-mref(1),y_0(1)-mref(2)],[1 0;0 1]); %sqrt(uv(1).^2 + uv(2).^2)
uv = uv/norm(uv);
quiver(x_0(1),y_0(1),uv(1),uv(2),0.2,'LineWidth',1,'Color','k','MaxHeadSize',1);

uv = -grad_l2([x_0(2)-mref(1),y_0(2)-mref(2)],[1 0;0 1]); %sqrt(uv(1).^2 + uv(2).^2)
uv = uv/norm(uv);
quiver(x_0(2),y_0(2),uv(1),uv(2),0.25,'LineWidth',1,'Color','k','MaxHeadSize',1);

% text(hib,1.1,'B"','HorizontalAlignment','right','FontSize',12,'BackgroundColor','w')
% text(lob,1.1,'B','FontSize',12,'BackgroundColor','w')

plot( [lob hib]  ,[(lnorm(2) - (lnorm(1) - lob)*eivec(1,2)/eivec(1,1)) ...
   (lnorm(2) - (lnorm(1) - hib)*eivec(1,2)/eivec(1,1))] ,'r--');

axis square
set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);
grid on
xlabel('$m_1$','interpreter','latex','FontSize',12)
title('\boldmath${\phi_s = \|m\|}_2^2$', 'interpreter', 'latex','FontSize',14);
text(1.25,-.25,'(b)','HorizontalAlignment','center','BackgroundColor','w','Color','k','EdgeColor','w','FontSize',12,'HorizontalAlignment','center','VerticalAlignment','middle')
%% Plot total objective function

obj_func_tot = 0.01*(smallness + smoothness) + misfit;
% set(figure, 'Position', [50 200 500 500]); axis square
axes('Position',[0.025 0.075 axs axs]);
% plot the origin for reference
plot([lob hib],[0 0],'k:');
plot([0 0],[lob hib],'k:');
% UV = -grad_l2([Xvec(:)-mref(1) Yvec(:)-mref(2)],[1 0;0 1]);
% U = reshape(UV(1,:),size(Xvec,1),size(Xvec,2));
% V = reshape(UV(2,:),size(Xvec,1),size(Xvec,2));
% quiver(Xvec,Yvec,U,V,0.5,'LineWidth',1); hold on

[~,temp] = min(obj_func_tot(:));
[i,j] = ind2sub(size(X),temp);

plot(xx(j),yy(i),'^k');hold on;


contour(yy,xx,obj_func_tot,15,'LineColor','k'); hold on
plot(x_0(1),y_0(1),'ko','MarkerFaceColor','k');hold on;
quiver(x_0(1),y_0(1),0.5,0,'LineWidth',2,'Color','r','MaxHeadSize',0.75)
% plot( [0 lnorm(1)]  ,[0 lnorm(2)],'k--','LineWidth',2);        
plot(x_0(2),y_0(2),'ko','MarkerFaceColor','w');
quiver(x_0(2),y_0(2),-0.5,0,'LineWidth',2,'Color','r','MaxHeadSize',0.75)
plot(0,0,'k*','LineWidth',1)

plot( [lob hib]  ,[(lnorm(2) - (lnorm(1) - lob)*eivec(1,2)/eivec(1,1)) ...
   (lnorm(2) - (lnorm(1) - hib)*eivec(1,2)/eivec(1,1))] ,'r--');

text(hib,1.0,'$A^\prime$','interpreter','latex','HorizontalAlignment','right','FontSize',12,'BackgroundColor','w','EdgeColor','k')
text(lob,1.0,'$A$','interpreter','latex','FontSize',12,'BackgroundColor','w','EdgeColor','k')
axis square
% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);
grid on
xlabel('$m_1$','interpreter','latex','FontSize',12)
ylabel('$m_2$','interpreter','latex','FontSize',12)
title('\boldmath$\phi(m) = \phi_d + \beta(\;\phi_s + \phi_x)$', 'interpreter', 'latex','FontSize',13);
text(1.25,-.25,'(d)','HorizontalAlignment','center','BackgroundColor','w','Color','k','EdgeColor','w','FontSize',12,'HorizontalAlignment','center','VerticalAlignment','middle')

%% Plot for Doug - Value of norms with l2 and l0 along null(G)
Wx = [1 -1];
% Wx = Wx'*Wx;
x2 = lob : 0.001 : hib;
x1 = ones(size(x2));

phiWx = l2_func(Wx,[x1;x2]') ;
phiWs = l2_func(speye(2),[x1;x2]') ;
% phiG = l2_func(G,[x1;x2]') ;
% phiWs = lp_func([x1;x2]',0,0.2);

phi_tot = phiWx + phiWs ;

axes('Position',[0.43 0.075 1.5*axs axs]);

plot(x2,phiWx,'k-','LineWidth',1);hold on
plot(x2,phiWs,'k:','LineWidth',2)
plot(x2,phi_tot,'k','LineWidth',2)
plot(x2,x2*0,'r--','LineWidth',2);
plot(0.5,1.5,'^k','MarkerFaceColor','w','MarkerSize',10);hold on;
% Do the same for l0

% axis([lob hib 0 4])
rs = wt_func([x1;x2]',p,epsl) ;
scale = 1./(max(rs,[],1))';
% scale = ones(size(x2'))*epsl^(1-p/2);
rs = rs * spdiags(sqrt(scale),0,length(x2),length(x2));
% rs = rs * spdiags(ones(size(x2'))*sqrt(epsl^(1-p/2)),0,length(x2),length(x2));

grid on

text(-.25,1.75,'\boldmath${\phi_x}$',...
    'BackgroundColor','w','EdgeColor','k','LineStyle','-','LineWidth',2,...
    'interpreter','latex','FontSize',10,...
    'HorizontalAlignment','center','VerticalAlignment','middle');

text(-0.25,1.0,'\boldmath${\phi_s}$',...
    'BackgroundColor','w','EdgeColor','k','LineStyle',':','LineWidth',2,...
    'interpreter','latex','FontSize',12,...
    'HorizontalAlignment','center','VerticalAlignment','middle');


text(0,2,0,'${\phi(m)}$',...
    'BackgroundColor','w','EdgeColor','k','LineWidth',2,...
    'interpreter','latex','FontSize',12,...
    'HorizontalAlignment','center','VerticalAlignment','middle');

% quiver(0.5,1.75,0,-0.25,1,'LineWidth',1,'Color','k','MaxHeadSize',0.75)
% text(0.5,1.75,'${min\;\phi(m)}$',...
%     'BackgroundColor','w','EdgeColor','k','LineWidth',1,...
%     'interpreter','latex','FontSize',12,...
%     'HorizontalAlignment','center','VerticalAlignment','bottom');

quiver(0.75,1,-.25,0.25,1,'LineWidth',1,'Color','k','MaxHeadSize',0.75)
text(0.75,1,'\boldmath${\frac{\partial\phi_s}{\partial\;m_1}}= 1$',...
    'BackgroundColor','w','EdgeColor','k','LineWidth',1,...
    'interpreter','latex','FontSize',12,...
    'HorizontalAlignment','center','VerticalAlignment','middle');

quiver(0.75,0.25,-0.25,0,1,'LineWidth',1,'Color','k','MaxHeadSize',0.75)
text(0.75,0.25,'\boldmath${\frac{\partial\phi_x}{\partial\;m_1}}= -1$',...
    'BackgroundColor','w','EdgeColor','k','LineWidth',1,...
    'interpreter','latex','FontSize',12,...
    'HorizontalAlignment','left','VerticalAlignment','middle');

title('Section $A-A^\prime$','interpreter','latex')

% axis square
% text(hib,0,'$A^\prime$','interpreter','latex','HorizontalAlignment','right','FontSize',12,'BackgroundColor','w','EdgeColor','k')
% text(lob,0,'A','FontSize',12,'BackgroundColor','w','EdgeColor','k')

xlabel('$m_1$','interpreter','latex','FontSize',12)
ylabel('$\phi$','interpreter','latex','FontSize',12,'Rotation',360)

axis([lob hib -0.5 2.5])

text(-.25,0,'$\boldmath{\phi_d}$','interpreter','latex',...
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'BackgroundColor','w','LineStyle','--',...
    'Color','k','EdgeColor','r','FontSize',12,'LineWidth',2)

text(1.4,2.25,'(e)','HorizontalAlignment','center','BackgroundColor','w','Color','k','EdgeColor','w','FontSize',12,'VerticalAlignment','middle')

% phiWs_l0 = sum((rs.*[x1;x2]).^2) ;
% 
% phi_tot = phiWs_l0 + phiWx;
% figure; plot(x2,phiWx);
% hold on
% plot(x2,phiWs_l0)
% plot(x2,phi_tot)
% leg = legend('$\phi_x:l_2$','$\phi_s:l_0$','$\phi_{m}$');
% set(leg,'interpreter', 'latex');
% [xval,idx] = min(phi_tot);
% plot(x2(idx),xval,'kx')
% text(x2(idx),xval,num2str(x2(idx)),'VerticalAlignment','bottom')
% 
% figure; plot(x2(2:end),abs(phiWx(2:end) - phiWx(1:end-1)));
% hold on
% plot(x2(2:end),abs(phiWs_l0(2:end) - phiWs_l0(1:end-1)))
% plot(x2(2:end),abs(phi_tot(2:end) - phi_tot(1:end-1)))
% leg = legend('$\mathbf{\hat g_x}$','$\mathbf{\hat g_s}$','$\mathbf{\hat g}$');
% set(leg,'interpreter', 'latex');