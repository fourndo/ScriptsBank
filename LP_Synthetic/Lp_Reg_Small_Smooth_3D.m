% Test approximation of l0 norm and plot objective function

close all

addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB


%% Set up problem
lob = -0.3;
hib = 1.1;

xx = lob : 1e-3 : hib;
yy = lob : 1e-3 : hib;

[Y,X] = ndgrid(xx,yy);


G = [1 1 1];
d = 1;

GtG = G'*G;
[eivec,eival] = eig(GtG);

epsl_in = 0.01;%10.^(-5:0);
p = [1 0];
beta = 0.01;

Wx = [1 -1/2 -1/2;-1/2 1 -1/2;-1/2 -1/2 1]*1;

alpha_s = 100;
alpha_x = 1;

% Vector gradient location
xvec = lob:0.2:hib;
yvec = lob:0.2:hib;
zvec = lob:0.2:hib;

[Xvec,Yvec,Zvec] = ndgrid(xvec,yvec,zvec);

mref = [0. 0.];

setcam = [180 90];
%% Function
l2_func = @(r,x) sum((r * x).^2,1);

mis_func = @(x) ( (G * x').^2 )' - 2 * x * G' * d + d'*d;

lp_func = @(x,pp,eps) sum(x'.^2./(x'.^2 + eps.^2).^(1-pp/2));

wt_func = @(x,pp,eps) sqrt(1./(x.^2 + eps.^2).^(1-pp/2));

grad_l2 = @(x,W) (W'*W)*x;

grad_regfunc = @(x,pp,eps) (x'./(x'.^2 + eps.^2).^(1-pp/2));

grad_objfunc = @(x,reg) ((G'*G + reg)* x') - kron(G' * d,ones(1,size(x,1))) ;

%% Solve least-norm problem
lnorm = G' * (d/(G*G'));
% x_0 = [lnorm(1) 0.6];
% y_0 = [lnorm(2) 0.2];
% x_in = [-.5:0.1:0.5 -.5:0.1:0 0.5 0.5 0.25 -0.25];
% y_in = [ones(1,11)*0.5 .5:-0.1:0 0.25 -0.25 0.5 0.5];
x_in = -0.195:0.1:1.0;
y_in = -0.20:0.1:1.0;

[x_in,y_in] = ndgrid(x_in,y_in);
x_in = x_in(:);
y_in = y_in(:);
%% Generate surface of solutions

nspace = @(x,y) (eivec(1,3)*(x - lnorm(1)) + eivec(2,3)*(y - lnorm(2)))./(-eivec(3,3)) + (lnorm(3));

z_in = nspace(x_in,y_in)+0.01 ;

[X2D,Y2D] = ndgrid(xvec,yvec);

Z2D = nspace(X2D,Y2D);



%% Misfit function
for pp = 1 : 2
%% Plot Base map
set(figure, 'Position', [25 100 1800 900])
ss = surf(X2D,Y2D,Z2D);
set(ss,'LineStyle','none','FaceColor',[.25 .25 .25])
alpha(ss,.4)

view(setcam);
hold on


% plot least-norm solution
plot3(lnorm(1),lnorm(2),lnorm(3),'ro')

for jj = 1 : length(x_in)
    
plot3(x_in(jj),y_in(jj),z_in(jj),'ko','MarkerFaceColor','w');hold on;

end
% text(x_0(2)+0.05,y_0(2)+0.1,['[' num2str(x_0(2)) ' ; ' num2str(y_0(2)) ']'],...
%         'BackgroundColor','w','EdgeColor','k')
    
% plot(0,0,'k*','LineWidth',1)
% if  pp == 1
%     set(gca,'XTickLabel',[]);
%     
% end
% axes(plot_h(pp*3+3));
% plot the origin for reference
% plot([lob hib],[0 0],'k:');hold on;
% plot([0 0],[lob hib],'k:');
% plot(x_in(1),y_in(1),'ko','MarkerFaceColor','k');hold on;
% plot(x_in(2),y_in(2),'ko','MarkerFaceColor','w');
% plot(0,0,'k*','LineWidth',1)

% plot( [lob hib]  ,[(lnorm(2) - (lnorm(1) - lob)*eivec(1,2)/eivec(1,1)) ...
% (lnorm(2) - (lnorm(1) - hib)*eivec(1,2)/eivec(1,1))] ,'r--');
% plot( [0 lnorm(1)]  ,[0 lnorm(2)],'k--','LineWidth',1);


% grid on
% axis square
% set(gca,'YTickLabel',[]);
% axis([lob hib lob hib])

for jj = 1 : length(x_in)
        
    x0 = x_in(jj);
    y0 = y_in(jj);
    z0 = z_in(jj);
    
    m0 = [x0;y0;z0];
    epsl = epsl_in;
    
    r = wt_func(m0,0,epsl);   
    R = spdiags(r,0,3,3);

    if pp == 2

        sp = 1/(r'*r);
        sq = 1;%/max(abs(  Wx'*Wx * m0 ) );
        
        if isinf(sq)

            sq = 1;

        end
        
    else
    sp = 1;
    sq = 1;            

    end
    
    MOF = ( alpha_s * sp*(R'*R) + alpha_x /2 * sq*(Wx'*Wx));

    % Plot unit gradient at the starting model
%     uv = -sp * grad_l2(m0,R); %uv*scale(1)
% %             uv = uv/norm(uv);
%     quiver3(x0,y0,z0,uv(1),uv(2),uv(3),0.1,'LineWidth',1,'Color','g','MaxHeadSize',0.5);
% 
%     uv = -sq * grad_l2(m0,W); %uv*scale(1)
% %             uv = uv/norm(uv);
%     quiver3(x0,y0,z0,uv(1),uv(2),uv(3),0.1,'LineWidth',1,'Color','r','MaxHeadSize',0.5);
%        
%     uv = G'*d - G'*G*m0; %uv*scale(1)
%             uv = uv/norm(uv);
%     quiver3(x0,y0,z0,uv(1),uv(2),uv(3),0.1,'LineWidth',1,'Color','b','MaxHeadSize',0.5);
        

    beta_in = beta;

    m_iter = [x0;y0;z0];
    
    phi_out = l2_func(speye(3),m_iter); 
    
    phi = sum( ( G * m_iter - d ).^2 ) + ( m_iter' * beta_in * MOF * m_iter );
    
    solves = 1;       % Count number of GN solves
    rddm = 1;       % Mesure the relative change in rel|dm|
    ddm = [1 1];    % Mesure of change in model update |dm|
    
    %% IRLS ITERATIONS
    while rddm > 1e-3 && solves < 100

        if epsl > 1e-5
        epsl = epsl/5;
        end
%         if ii ~=1
%         beta_in = beta_in * (phid_l2 / norm(G*[x_iter;y_iter]-d+1e-4)+1e-4) ;
%         end
        
%         beta_in=beta_in/2;

r =wt_func(m0,0,epsl);
        R = spdiags(r,0,3,3);
        

        if pp == 2
            
            sp = 2/(r'*r);
            sq = 1;%/max(abs( Wx'*Wx * m_iter ) );
            
            if isinf(sq)
                
                sq = 1;
                
            end
            
        else
            sp = 1;
            sq = 1;            
            
        end


        
        phi_in = l2_func(sqrt(sp) * R,m_iter) + l2_func(sqrt(sq) /2 * Wx,m_iter);
     
        scale = phi_out/phi_in;
        
        MOF = scale*( alpha_s * sp*(R'*R) + alpha_x /2 * sq*(Wx'*Wx));

        A = (G'*G + beta_in*MOF);
        
%         RHS = G'*d;
%         
%         m_new = A\RHS;
%         
%         dm = m_new - m_iter;
        RHS = -(G'*(G*m_iter-d) + beta_in*MOF*m_iter);

        [dm,normr,count]=CGiter(zeros(3,1),A,RHS);
%         dm = A\RHS;
         gamma = 2;

        % Reduce step length in order to reduce phid
        phi_temp = 0;
        objfunc = @(x) sum( ( G * x - d ).^2 ) + ( x' * beta_in * MOF * x );
        m = m_iter;
        
        while (phi_temp > phi || gamma == 2) && gamma>1e-4

            gamma = 0.5 * gamma;

            gdm = gamma * dm;

            ddm(2) = norm(gdm);

            m_temp = m + gdm;

            phi_temp = objfunc(m_temp);

        end
          
        if solves == 1

            rddm = 1;
            ddm(1) = ddm(2);

        else

            rddm = ddm(2)/ddm(1);

        end
%         plot3([m_iter(1) m_iter(1)+ dm(1)],[m_iter(2) m_iter(2)+dm(2)],[m_iter(3) m_iter(3)+dm(3)],'k')
        
        m_iter = m_temp;

%         mis_func([x_iter y_iter])

        solves = solves + 1;

        phi_out = scale * phi_in;

    end

%      mis_func([x_iter y_iter]);
     
%     plot(x_iter,y_iter,'k^','MarkerSize',8);
%     text(x0+0.05,y0+(-1)^jj*0.1,['\beta=' num2str(beta_in)],...
%             'BackgroundColor','w','EdgeColor','k')

%% FOR PLOTTING ONLY

        if pp == 2
            
            sp = 1/max(abs( R'*R * m_iter ) );
            sq = 1/max(abs( Wx'*Wx * m_iter ) );
            
            if isinf(sq)
                
                sq = 1;
                
            end
            
        else
            sp = 1;
            sq = 1;            
            
        end
        % Plot unit gradient at the starting model
%        uv = -sp * grad_l2(m_iter,R); %uv*scale(1)
% %             uv = uv/norm(uv);
%     quiver3(m_iter(1),m_iter(2),m_iter(3),uv(1),uv(2),uv(3),0.1,'LineWidth',1,'Color','g','MaxHeadSize',0.5);
% 
%     uv = -sq * grad_l2(m_iter,W); %uv*scale(1)
% %             uv = uv/norm(uv);
%     quiver3(m_iter(1),m_iter(2),m_iter(3),uv(1),uv(2),uv(3),0.1,'LineWidth',1,'Color','r','MaxHeadSize',0.5);
%        
%     uv = G'*d - G'*G*m_iter; %uv*scale(1)
% %             uv = uv/norm(uv);
%     quiver3(m_iter(1),m_iter(2),m_iter(3),uv(1),uv(2),uv(3),0.1,'LineWidth',1,'Color','b','MaxHeadSize',0.5);
        
%     R = spdiags(wt_func([x_iter,y_iter],0,epsl),0,2,2);
% 
%     IRLS = l2_func(R,[X(:) Y(:)]) + l2_func(W,[X(:) Y(:)]);
%     IRLS = reshape(IRLS,size(X,1),size(X,2))*scale;

    if jj == 1
        
%         contour(xx,yy,misfit +  beta_in*IRLS,10,'LineColor','k','LineStyle','-'); hold on

    else

%         contour(xx,yy,misfit +  beta_in*IRLS,10,'LineColor','k','LineStyle','--'); hold on

    end
%             
%     if jj == 1
        aa=quiver3(x0,y0,z0,m_iter(1) - x0,m_iter(2) - y0,m_iter(3) - z0,0,'LineWidth',1,'MaxHeadSize',0.1);
        plot3(m_iter(1),m_iter(2),m_iter(3),'ko','MarkerFaceColor','k');
%          text(x_iter+0.05,y_iter+0.1,['[' num2str(round(100*x_iter)/100) ' ; ' num2str(round(100*y_iter)/100) ']'],...
%         'BackgroundColor','w','EdgeColor','k')
%     else
% %         quiver(x0,y0,[x_iter - x0],[y_iter - y0],0.9,'LineWidth',2,'Color','r','MaxHeadSize',0.25);
%         
% 
%         plot(x_iter,y_iter,'k^','MarkerSize',8);
% %          text(x_iter-0.2,y_iter-0.1,['[' num2str(num2str(round(100*x_iter)/100)) ' ; ' num2str(num2str(round(100*y_iter)/100)) ']'],...
% %         'BackgroundColor','w','EdgeColor','k')
%     
% 
%     end
    
end

if pp == 1

    xlabel('(d)')
    text(1.05,0.4,'\boldmath$l_1$', 'interpreter', 'latex','FontSize',14);
    
else

    xlabel('(e)')
    text(1.05,0.4,'\boldmath$l_0$', 'interpreter', 'latex','FontSize',14);

end

axis equal
axis square
% axis([lob hib lob hib lob hib])

end


