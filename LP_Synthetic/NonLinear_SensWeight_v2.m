% Test approximation of l0 norm and plot objective function

% close all
% clear all

addpath '..\FUNC_LIB'


%% Set up problem
lowb = -pi;
uppb = pi;

subfact = 8;

xx = lowb : (uppb-lowb)/64 : uppb;
yy = lowb : (uppb-lowb)/64 : uppb;

[X,Y] = ndgrid(xx,yy);

nC = length(X(:));

indSub = false(size(X));
indSub(1:subfact:end,1:subfact:end) = true;


nvx = sum(indSub(1,:));
nvy = sum(indSub(:,1));

Xvec = reshape(X(indSub), nvx, nvy);
Yvec = reshape(Y(indSub), nvx, nvy);

indSub = indSub(:);

d = 1;

% G = [1 2];
% d = [1];

% GtG = G'*G;
% [eivec,eival] = eig(GtG);

epsl = 1e-10;%10.^(-5:0);
beta = 1e-1;

% % Vector gradient location
% xvec = lowb:pi/4:uppb;
% yvec = lowb:pi/4:uppb;

% [Xvec,Yvec] = ndgrid(xvec,yvec);

mref = [0. 0.];

x_0 = [49 33];
%        9 41];

axs  = 0.36;
offs = 0.4;
tol = 1e-1;
mark = {'ks','k^'};
arrow = [-0.4 0.5; 0 -1;0.4 0.5]/4;
%% Function
% F = @(x1,x2)(exp(-x1/pi).*(x2));
% J = @(x1,x2)([-1/pi*exp(-x1/pi).*(x2) exp(-x1/pi)]);
% F = @(x1,x2)(cos(x1*2).*sin(x2*2));
% J = @(x1,x2)([-2*sin(x1*2).*sin(x2*2) 2*cos(x1*2).*cos(x2*2)]);
F = @(x1,x2)(sin(x1+pi/3).*sin(x2+pi/4));
J = @(x1,x2)([cos(x1+pi/4).*sin(x2+pi/4) sin(x1+pi/4).*cos(x2+pi/4)]);
% a = 1; b = 8; F = @(x)((x(:,1)*a)+b*(x(:,2)));
% J = @(x)([a b]);
% F = @(x1,x2)(x1.^2+x2.^2)^0.5;
% J = @(x1,x2)([2*x1.*(x1.^2+x2.^2)^-0.5 2*x2.*(x1.^2+x2.^2)^-0.5]);
% F = @(x1,x2)(tanh(x1+epsl).*tanh(x2+epsl));
% J = @(x1,x2)([sech(x1+epsl).^2.*tanh(x2+epsl) tanh(x1+epsl).*sech(x2+epsl).^2]);
% F = @(x1,x2)(x1.*x2);
% J = @(x1,x2)([x2 x1]);

reg_func = @(W,x) (W * (x)')'*(W * (x)');

mis_func = @(x) (F(x(:,1),x(:,2))-d).^2;

lp_func = @(x,pp,eps) sum(x'.^2./(x'.^2 + eps.^2).^(1-pp/2));

wt_func = @(x,pp,eps) sqrt(1./(x'.^2 + eps.^2).^(1-pp/2));

grad_l2 = @(x) J(x(:,1),x(:,2))'*x;

grad_regfunc = @(x,pp,eps) (x'./(x'.^2 + eps.^2).^(1-pp/2));

grad_misfitfunc = @(x) J(x(:,1),x(:,2))'*spdiags(F(x(:,1),x(:,2))-d,0,size(x,1),size(x,1));
% grad_objfunc = @(x,reg) ((J(x).*kron((F(x)-d),ones(1,size(x,2))) + beta*reg * x)) ;

% Wr = @(x) ()
%% Solve least-norm problem
% lnorm = J(0)' * (d/(J(0)*J(0)'));
% x_0 = [lnorm(1) 1];
% y_0 = [lnorm(2) epsl];
W = [1 0;0 1];
phid = zeros(size(X));
d_phid = zeros(length(X(:)),2);
phim_wr = zeros(size(X));
phim_l2 = zeros(size(X));
d_phim_wr = zeros(length(X(:)),2);
d_phim_l2 = zeros(length(X(:)),2);
for ii = 1:length(X(:))
    
    m = [X(ii) Y(ii)];
    j = J(m(1),m(2));
    fm = F(m(1),m(2));
    wr = abs(j)';
    wr = wr + max(wr)*epsl;
    Wr = spdiags((wr/max(wr)).^0.5,0,2,2);
    phid(ii) = (fm - d)^2;
    d_phid(ii,:) = -j * (fm - d);
    phim_wr(ii) = reg_func( Wr, m);%( Wr * W * m')' * ( Wr * W * m');
    phim_l2(ii) = reg_func( W, m);%' * ( W * m');
    d_phim_wr(ii,:) = -( Wr * W )'*( Wr * W) * m';
    d_phim_l2(ii,:) = -( W )'*( W) * m';

end
%% Misfit function
U = reshape(d_phid(indSub,1),nvx,nvy);
V = reshape(d_phid(indSub,2),nvx,nvy);

set(figure, 'Position', [50 0 775 600]);

axes('Position',[0.035 0.56 axs axs]);

% plot the origin for reference
plot([lowb uppb],[0 0],'k:');
plot([0 0],[lowb uppb],'k:');
contour(X,Y,phid,10); hold on

[val,~] = min(phid(:));
contour(X,Y,phid,[0 val+tol],'k--', 'LineWidth',1); hold on

% plot(0,0,'k*','LineWidth',1)
quiver(X(indSub), Y(indSub), U(:), V(:), 0.6, 'LineWidth', 1.,'Color','k');
set(gca,'Ydir','normal')
% plot(x_0(1),y_0(1),'ko','MarkerFaceColor','k', 'MarkerSize',8);hold on;
% plot(x_0(2),y_0(2),'ko','MarkerFaceColor','w');

for ii = 1 : size(x_0,1)
     
    x_in = x_0(ii,:);
    plot(X(x_in(1),1),Y(1,x_in(2)),mark{ii},'MarkerFaceColor','k', 'MarkerSize',8);hold on;
end

set(gca,'Ytick',[-pi -pi/2 0 pi/2 pi])
set(gca,'YTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'})
set(gca,'Xtick',[ -pi/2 0 pi/2])
set(gca,'XTickLabel',{'-\pi/2','0','\pi/2'})
% set(gca,'YTickLabel',{'\pi'},pi);
axis square
grid on
axis([lowb uppb lowb uppb])
xlabel('$m_1$','interpreter','latex','FontSize',12)
ylabel('$m_2$','interpreter','latex','FontSize',12,'Rotation',0)
title('$\phi_d = \|\mathbf{F \;[ m] - d}\|_2^2$', 'interpreter', 'latex','FontSize',14);
% quiver([0;0],[0;0],eivec(:,1),eivec(:,2),1,'LineWidth',2);
text(lowb+offs,lowb+offs,'(a)','HorizontalAlignment','center','BackgroundColor','w','Color','k','EdgeColor','w','FontSize',12,'HorizontalAlignment','center','VerticalAlignment','middle')
% 
%% l2 Model norm

axes('Position',[0.34 0.56 axs axs]);
% plot the origin for reference
plot([lowb uppb],[0 0],'k:');
plot([0 0],[lowb uppb],'k:');

contour(X,Y,phim_l2,10); hold on

[val,~] = min(phim_l2(:));
contour(X,Y,phim_l2,[0 val+tol],'k--', 'LineWidth',1); hold on

indSub = false(size(X));
indSub(1:8:end,1:8:end) = true;
nvx = sum(indSub(1,:));
nvy = sum(indSub(:,1));

U = reshape(d_phim_l2(indSub,1),nvx,nvy);
V = reshape(d_phim_l2(indSub,2),nvx,nvy);
quiver(X(indSub), Y(indSub), U(:), V(:) ,0.75, 'LineWidth', 1.5,'Color','k');

for ii = 1 : size(x_0,1)
     
    x_in = x_0(ii,:);
    plot(X(x_in(1),1),Y(1,x_in(2)),mark{ii},'MarkerFaceColor','k', 'MarkerSize',8);hold on;
end
axis square
% set(gca,'Ytick',[-pi -pi/2 0 pi/2 pi])
set(gca,'YTickLabel',[])
set(gca,'Xtick',[ -pi/2 0 pi/2])
set(gca,'XTickLabel',{'-\pi/2','0','\pi/2'})
grid on
% xlabel('$m_1$','interpreter','latex','FontSize',12)
title('\boldmath${\phi_m = \|\mathbf{m}\|}_2^2$', 'interpreter', 'latex','FontSize',14);
text(lowb+offs,lowb+offs,'(b)','HorizontalAlignment','center','BackgroundColor','w','Color','k','EdgeColor','w','FontSize',12,'HorizontalAlignment','center','VerticalAlignment','middle')
%% Plot total objective function
% 
dphi_tot = beta*d_phim_l2 + d_phid;
U = reshape(dphi_tot(:,1),size(X));
V = reshape(dphi_tot(:,2),size(X));

obj_func_tot_l2 = beta*(phim_l2) + phid;%abs(U)+abs(V);%
% % set(figure, 'Position', [50 200 500 500]); axis square
axes('Position',[0.64 0.56 axs axs]);

[val,temp] = min(obj_func_tot_l2(:));
[i,j] = ind2sub(size(X),temp);

contour(X,Y,obj_func_tot_l2,10); hold on
contour(X,Y,obj_func_tot_l2,[0 val+tol],'k--', 'LineWidth',1); hold on
% plot(xx(i),yy(j),'^k');hold on;
% plot(lnorm(1),lnorm(2),'ko','MarkerFaceColor','r');

cmap = colormap(gray);
colormap(cmap(end:-1:1,:))
% plot(x_0(1),y_0(1),'ko','MarkerFaceColor','k', 'MarkerSize',8);hold on;
% plot(x_0(2),y_0(2),'ko','MarkerFaceColor','w');
% 

% quiver(X(indSub), Y(indSub), U(:), V(:) ,0.5, 'LineWidth', 1,'Color','k');
for ii = 1 : size(x_0,1)
     
    x_in = x_0(ii,:);
    m = [X(x_in(1),1); Y(1,x_in(2))];
    plot(m(1), m(2) ,mark{ii},'MarkerFaceColor','k', 'MarkerSize',8);hold on;
%     
%     h = streamline(X',Y',U',V',X(x_in(1),1),Y(1,x_in(2)),[0.5,100]);
%     xx = h.XData;
%     yy = h.YData;
%     ang = atan2((xx(9) - xx(10)),(yy(9) - yy(10)));
%     
%     rot = [cos(ang) sin(ang); -sin(ang) cos(ang)];
%     
%     arrow = (rot*arrow')';
%     patch(arrow(:,1)+xx(10), arrow(:,2)+yy(10),'b')

    % Solve the system with GN
    
    
    
    
    dm = [1 0];
    count(ii) = 0;
    while norm(dm) > 1e-3
%         G = [m(1), m(2)];
        mref = [0; 0];
        phi_in = (F(m(1),m(2))-d)^2 + beta * (W*m)'*(W*m);
        PreC = speye(2);
        Pac = speye(2);
        lowBvec = [lowb; lowb];
        uppBvec = [uppb; uppb];
        [m_out, tncg, Pac] = GN_solver( F, J(m(1),m(2)), m, mref, d, phi_in, beta , PreC, Pac, lowBvec, uppBvec , W, W);
        
        dm = m_out-m;
        width = 4*((m(1)-m_out(1))^2 + (m(2)-m_out(2))^2)^0.25 / pi;
        plot([m(1) m_out(1)], [m(2) m_out(2)], 'k', 'LineWidth',width)
        count(ii) = count(ii) + 1;
        if count(ii) == 1
%             plot(m_out(1), m_out(2),'kx');
            ang = atan2(dm(1),dm(2));

            rot = [cos(ang) sin(ang); -sin(ang) cos(ang)];

            vec = -(rot*arrow')';
            patch(vec(:,1)+m(1)+dm(1)/3, vec(:,2)+m(2)+dm(2)/3,'k')
        end
        
        m = m_out;
    end
    
%     plot(m_out(1), m_out(2),'kx');

%     sprintf(['Number GN iter: ' num2str(count)]);
end

axis square
set(gca,'Ytick',[ -pi/2 0 pi/2])
set(gca,'YTickLabel',{'-\pi/2','0','\pi/2'})
set(gca,'YAxisLocation','right')
set(gca,'Xtick',[ -pi/2 0 pi/2])
set(gca,'XTickLabel',{'-\pi/2','0','\pi/2'})

grid on
% xlabel('$m_1$','interpreter','latex','FontSize',12)
% ylabel('$m_2$','interpreter','latex','FontSize',12)
title('\boldmath$\phi(m) = \phi_d + \beta\;\phi_m$', 'interpreter', 'latex','FontSize',13);
text(lowb+offs,lowb+offs,'(c)','HorizontalAlignment','center','BackgroundColor','w','Color','k','EdgeColor','w','FontSize',12,'HorizontalAlignment','center','VerticalAlignment','middle')

%% Weighted l2 Model norm

% axes('Position',[0.34 0.09 axs axs]);
% % plot the origin for reference
% plot([lowb uppb],[0 0],'k:');
% plot([0 0],[lowb uppb],'k:');
% 
% contour(X,Y,phim_wr,5); hold on
% 
% [val,~] = min(phim_wr(:));
% contour(X,Y,phim_wr,[0 val+tol],'r--', 'LineWidth',1); hold on
% 
% U = reshape(d_phim_wr(indSub,1),nvx,nvy);
% V = reshape(d_phim_wr(indSub,2),nvx,nvy);
% quiver(X(indSub), Y(indSub), U(:), V(:) ,1., 'LineWidth', 1.5,'Color','k');
% 
% for ii = 1 : size(x_0,1)
%      
%     x_in = x_0(ii,:);
%     plot(X(x_in(1),1),Y(1,x_in(2)),mark{ii},'MarkerFaceColor','k', 'MarkerSize',8);hold on;
% end
% axis square
% % set(gca,'Ytick',[-pi -pi/2 0 pi/2 pi])
% set(gca,'YTickLabel',[])
% set(gca,'Xtick',[ -pi/2 0 pi/2])
% set(gca,'XTickLabel',{'-\pi/2','0','\pi/2'})
% grid on
% xlabel('$m_1$','interpreter','latex','FontSize',12)
% % ylabel('$m_2$','interpreter','latex','FontSize',12,'Rotation',0)
% title('\boldmath${\phi_m = \|\mathbf{W}_r\;\mathbf{m}\|}_2^2$', 'interpreter', 'latex','FontSize',14);
% text(lowb+offs,lowb+offs,'(d)','HorizontalAlignment','center','BackgroundColor','w','Color','k','EdgeColor','w','FontSize',12,'HorizontalAlignment','center','VerticalAlignment','middle')
% %% Plot total objective function
% % 
% 
% % % set(figure, 'Position', [50 200 500 500]); axis square
% axes('Position',[0.64 0.09 axs axs]);
% dphi_tot = beta*d_phim_wr + d_phid;
% U = reshape(dphi_tot(:,1),size(X));
% V = reshape(dphi_tot(:,2),size(X));
% % [val,temp] = min(obj_func_tot(:));
% [i,j] = ind2sub(size(X),temp);
% obj_func_tot_wr =  beta*(phim_wr) + phid; %abs(U)+abs(V);%
% [val,ind] = min(obj_func_tot_wr(:));
% contour(X,Y,obj_func_tot_wr,10, 'LineWidth', 0.5); hold on
% contour(X,Y,obj_func_tot_wr,[0 val+tol],'r--', 'LineWidth',1); hold on
% % plot(xx(i),yy(j),'^k');hold on;
% % plot(lnorm(1),lnorm(2),mark{ii},'MarkerFaceColor','r');
% 
% cmap = colormap(gray);
% colormap(cmap(end:-1:1,:))

% plot(x_0(1),y_0(1),mark{ii},'MarkerFaceColor','k', 'MarkerSize',8);hold on;
% plot(x_0(2),y_0(2),mark{ii},'MarkerFaceColor','w');
% 
% quiver(X(indSub), Y(indSub), U(:), V(:) ,0.5, 'LineWidth', 1,'Color','k');

% for ii = 1 : size(x_0,1)
%      
%     x_in = x_0(ii,:);
%     plot(X(x_in(1),1),Y(1,x_in(2)),mark{ii},'MarkerFaceColor','k', 'MarkerSize',8);hold on;
%     
%         h = streamline(X',Y',U',V',X(x_in(1),1),Y(1,x_in(2)),[0.5,100]);
%     xx = h.XData;
%     yy = h.YData;
%     ang = atan2((xx(9) - xx(10)),(yy(9) - yy(10)));
%     arrow = [-0.5 0.5; 0 -1;0.5 0.5]/4;
%     rot = [cos(ang) sin(ang); -sin(ang) cos(ang)];
%     
%     arrow = (rot*arrow')';
%     patch(arrow(:,1)+xx(10), arrow(:,2)+yy(10),'b')
% end
for ii = 1 : size(x_0,1)
     
    x_in = x_0(ii,:);
    m = [X(x_in(1),1); Y(1,x_in(2))];
    plot(m(1), m(2) ,mark{ii},'MarkerFaceColor','k', 'MarkerSize',8);hold on;
%     
%     h = streamline(X',Y',U',V',X(x_in(1),1),Y(1,x_in(2)),[0.5,100]);
%     xx = h.XData;
%     yy = h.YData;
%     ang = atan2((xx(9) - xx(10)),(yy(9) - yy(10)));
%     arrow = [-0.5 0.5; 0 -1;0.5 0.5]/4;
%     rot = [cos(ang) sin(ang); -sin(ang) cos(ang)];
%     
%     arrow = (rot*arrow')';
%     patch(arrow(:,1)+xx(10), arrow(:,2)+yy(10),'b')

    % Solve the system with GN
    
    
    
    
    dm = [1 0];
    count_wr(ii) = 0;
    while norm(dm) > 1e-3
%         G = [m(1), m(2)];
        mref = [0; 0];
        
        
        wr = abs(J(m(1),m(2)))';
%         wr = wr + max(wr)*epsl;
        Wr = spdiags((wr/max(wr) + epsl).^0.5,0,2,2);
        phi_in = (F(m(1),m(2))-d)^2 + beta * (Wr*m)'*(Wr*m);
        PreC = speye(2);
        Pac = speye(2);
        lowBvec = [lowb; lowb];
        uppBvec = [uppb; uppb];
        [m_out, tncg, Pac] = GN_solver(F, J(m(1),m(2)), m, mref, d, phi_in, beta , PreC, Pac, lowBvec, uppBvec , Wr'*Wr, Wr);
        
        dm = m_out-m;
        

        width = 4*((m(1)-m_out(1))^2 + (m(2)-m_out(2))^2)^0.25 / pi;
        plot([m(1) m_out(1)], [m(2) m_out(2)], 'r', 'LineWidth',width)
        

        count_wr(ii) = count_wr(ii) + 1;        
        if count_wr(ii) == 1
            ang = atan2(dm(1),dm(2));

            rot = [cos(ang) sin(ang); -sin(ang) cos(ang)];

            vec = -(rot*arrow')';
            patch(vec(:,1)+m(1)+dm(1)/2, vec(:,2)+m(2)+dm(2)/2,'r')
        end
        m = m_out;
    end
%     plot(m_out(1), m_out(2),'kx');
%     sprintf(['Number GN iter: ' num2str(count)]);
end
% axis square
% set(gca,'Ytick',[ -pi/2 0 pi/2])
% set(gca,'YTickLabel',{'-\pi/2','0','\pi/2'})
% set(gca,'YAxisLocation','right')
% set(gca,'Xtick',[ -pi/2 0 pi/2])
% set(gca,'XTickLabel',{'-\pi/2','0','\pi/2'})
% grid on
% xlabel('$m_1$','interpreter','latex','FontSize',12)
% % ylabel('$m_2$','interpreter','latex','FontSize',12)
% % title('\boldmath$\phi(m) = \phi_d + \beta(\;\phi_s + \phi_x)$', 'interpreter', 'latex','FontSize',13);
% text(lowb+offs,lowb+offs,'(e)','HorizontalAlignment','center','BackgroundColor','w','Color','k','EdgeColor','w','FontSize',12,'HorizontalAlignment','center','VerticalAlignment','middle')

