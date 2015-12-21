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
xvec = 0.3:0.05:0.7;
yvec = 0.1:0.05:0.5;

[Xvec,Yvec] = ndgrid(xvec,yvec);

xvec = Xvec(:);
yvec = Yvec(:);

mref = [0. 0.];
%% Function
l2_func = @(r,x) sum((sqrt(r) * x').^2);

mis_func = @(x) ( (G * x').^2 )' - 2 * x * G' * d + d'*d;

lp_func = @(x,pp,eps) sum(x'.^2./(x'.^2 + eps.^2).^(1-pp/2));

wt_func = @(x,pp,eps) (1./(x'.^2 + eps.^2).^(1-pp/2));

grad_l2 = @(x,W) W*x';

grad_regfunc = @(x,pp,eps) (x'./(x'.^2 + eps.^2).^(1-pp/2));

grad_objfunc = @(x,reg) ((G'*G + reg)* x') - kron(G' * d,ones(1,size(x,1))) ;

%% Misfit function
misfit = mis_func([X(:) Y(:)]);
misfit = reshape(misfit,size(X,1),size(X,2));

%% L2 solution
m_out = (G'*G + beta*eye(2))\(G'*d);
phid_l2 = norm(G * m_out - d);


%% Plot Base map
subplot(2,1,1);
% plot the origin for reference
plot([lob hib],[0 0],'k:');hold on;
plot([0 0],[lob hib],'k:');
plot(xvec,yvec,'k*','MarkerFaceColor','k');hold on;
% text(x_0(1)+0.05,y_0(1)+0.1,['[' num2str(x_0(1)) ' ; ' num2str(y_0(1)) ']'],...
%         'BackgroundColor','w','EdgeColor','k')
        
% plot(x_0(2),y_0(2),'ko','MarkerFaceColor','w');
% text(x_0(2)+0.05,y_0(2)+0.1,['[' num2str(x_0(2)) ' ; ' num2str(y_0(2)) ']'],...
%         'BackgroundColor','w','EdgeColor','k')
    
plot(0,0,'k*','LineWidth',1)

subplot(2,1,2);
% plot the origin for reference
plot([lob hib],[0 0],'k:');hold on;
plot([0 0],[lob hib],'k:');
plot(xvec,yvec,'ko','MarkerFaceColor','k');hold on;
% plot(x_0(2),y_0(2),'ko','MarkerFaceColor','w');
plot(0,0,'k*','LineWidth',1)
% plot( [lob hib]  ,[(y_0(1) - (x_0(1) - lob)*eivec(1,2)/eivec(1,1)) ...
% (y_0(1) - (x_0(1) - hib)*eivec(1,2)/eivec(1,1))] ,'r--');
% plot( [0 x_0(1)]  ,[0 y_0(1)],'k--','LineWidth',1);
    
for jj = 1 : length(xvec)
        
    x0 = xvec(jj);
    y0 = yvec(jj);
    
    R = spdiags(wt_func([x0,y0],0,epsl),0,2,2);

    scale = l2_func([1 0;0 1],[x0 y0])/l2_func(R,[x0 y0]);

    IRLS = l2_func(R,[X(:) Y(:)]);
    IRLS = reshape(IRLS,size(X,1),size(X,2))*scale;
      
subplot(2,1,1)
               
    UV = -grad_l2([Xvec(:) Yvec(:)],R);
    U = reshape(UV(1,:),size(Xvec,1),size(Xvec,2));
    V = reshape(UV(2,:),size(Xvec,1),size(Xvec,2));

    % Plot unit gradient at the starting model
    uv = -grad_l2([x0,y0],R); %uv*scale(1)
    uv = uv/norm(uv);
    
    if x0~=0 & y0~=0
    quiver(x0,y0,uv(1),uv(2),0.1,'LineWidth',1,'Color','k','MaxHeadSize',1);
    end
    
%     if jj == 1
%         
%         contour(xx,yy,IRLS,10,'LineColor','k','LineStyle','-'); hold on
% %         quiver(Xvec,Yvec,U,V,0.5,'Color','b','LineWidth',1,'LineStyle','-');
%         
%     else
% 
%         contour(xx,yy,IRLS,10,'LineColor','k','LineStyle','--'); hold on
% %         quiver(Xvec,Yvec,U,V,0.5,'Color','b','LineWidth',1,'LineStyle','--');
% 
%     end
    

    set(gca,'Ydir','normal')
    grid on
    axis square
    set(gca,'YTickLabel',[]);
    xlabel('(e)')

    % Iterate until convergence
    subplot(2,1,2)
    
    
    plot([lob hib],[0 0],'k:');
    plot([0 0],[lob hib],'k:');

    beta_in = beta;

    x_iter = x0;
    y_iter = y0;
    
    phi_out = l2_func([1 0;0 1],[x_iter y_iter]); 
    
    dm = 1;
    m_in = [x0;y0];
    ii = 1;
    while dm > 1e-4

        if ii ~=1
        beta_in = beta_in * (phid_l2 / norm(G*[x_iter;y_iter]-d+1e-4)+1e-4) ;
        end
        
        R = spdiags(wt_func([x_iter,y_iter],0,epsl),0,2,2);
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

        dm = norm(m_in - m_out);

        phi_out = scale * phi_in;

        m_in = m_out;
        ii = ii + 1;
    end

     mis_func([x_iter y_iter]);
     
%     plot(x_iter,y_iter,'k^','MarkerSize',8);
%     text(x0+0.05,y0+(-1)^jj*0.1,['\beta=' num2str(beta_in)],...
%             'BackgroundColor','w','EdgeColor','k')

%% FOR PLOTTING ONLY
    R = spdiags(wt_func([x_iter,y_iter],0,epsl*100),0,2,2);

    IRLS = l2_func(R,[X(:) Y(:)]);
    IRLS = reshape(IRLS,size(X,1),size(X,2))*scale;

%     if jj == 1
%         
%         contour(xx,yy,misfit + beta_in*10*IRLS,10,'LineColor','k','LineStyle','-'); hold on
% 
%     else
% 
%         contour(xx,yy,misfit + beta_in*10*IRLS,10,'LineColor','k','LineStyle','--'); hold on
% 
%     end
            

        quiver(x0,y0,[x_iter - x0],[y_iter - y0],1,'LineWidth',1,'Color','r','MaxHeadSize',0.25);
        plot(x_iter,y_iter,'k^','MarkerSize',8,'MarkerFaceColor','k');
%          text(x_iter+0.05,y_iter+0.1,['[' num2str(round(100*x_iter)/100) ' ; ' num2str(round(100*y_iter)/100) ']'],...
%         'BackgroundColor','w','EdgeColor','k')


    
end

grid on
axis square
set(gca,'YTickLabel',[]);
axis([lob hib lob hib])
xlabel('(f)')



