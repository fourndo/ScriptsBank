% Plotting the gradient of the objective function as a function of epsilon
clear all
close all

m=1e-6:1e-5:0.1;
eps = 10.^-(0:0.5:3);
Ekb = @(p,ep) (x.^2 + ep.^2).^(p/2);
r = @(x,p,ep) 1./(x.^(2-p) + ep);

minsup = @(x,p,ep) x./(x.^2 + ep.^2).^(1-p/2);

wght = @(x,p,ep) 1./(x.^2 + ep.^2).^(1-p/2);

gradfx = @(x,p,ep) ((p-1)*x.^(2-p) + ep)./(x.^(2-p) + ep).^2;
inflex = @(p,ep) ((3-p)*ep/abs(1-p)).^(1/(2-p));

func = @(x,p,ep)((x+ep).^(2-p) - x.*(2-p).*(x+ep).^(1-p))./(x+ep).^(4-2*p);

curv =  @(x,p,ep) (( (2-p)*(1-p) * x.^(3-2*p) + (2-p)*(p-3)*ep*x.^(1-p) ) ./...
    ( (x.^(2-p) + ep).^2 + ( ( (p-1)*x.^(2-p) + ep ) ./ (x.^(2-p) + ep) ).^2 ).^(3/2)).^-1;

curv2 =  @(x,p,ep) (x.^(1-p) * (2-p) .* ( x.^(2-p) - 3*ep - p*x.^(2-p) + p*ep)) ./...
    (x.^(2-p) + ep).^3 ./...
    (1 + (( (x.^(2-p) + ep) - (2-p)*x.^(2-p)) ./ (x.^(2-p) + ep).^2 ).^2 ).^(3/2);

% dkdt = @(p,ep)(ep.^(0.5) - ep ) .^(1/((2-p)));
dkdt = @(p,ep) ep .^(1/(2*(2-p)));
polyn = @(a,ep) [1 0 2*ep-a 0 (ep^2 + a *ep)];

% dkdt = @(p,ep)(ep).^(1/((2-p)));
d2fdx2 = @(x,p,ep) ( (2-p)*(1-p) * x.^(3-2*p) + (2-p)*(p-3)*ep*x.^(1-p) ) ./...
    (x.^(2-p) + ep).^3;

% d2fdx2_2 = @(x,p,ep) (x.^(1-p) * (2-p) .* ( x.^(2-p) - 3*ep - p*x.^(2-p) + p*ep)) ./...
%     (x.^(2-p) + ep).^3;

% approx = @(x,p,ep) ( ( (1-p)*x.^(-p).*(2-p).*( x.^(2-p) - 3*ep -p*x.^(2-p) +p*ep) +...
%     x.^(1-p) .* (2-p) .* ( (2-p)*x.^(1-p) - (2-p)*p*x.^(1-p))).*(x.^(2-p) +ep)-...
%      x.^(1-p).*(2-p) .* ( x.^(2-p) - 3*ep -p*x.^(2-p) +p*ep) .* 3 .* (2-p) .* x.^(1-p))./...
%      (x.^(2-p) + ep).^4;

p = 0.0;
figure;
for ee = 1 : length(eps)
    
%     y= Ekb(1,eps(ee));
%     norm(y)
figure

    dy1 = minsup(m,p,eps(ee));
%     r1 = r(m,p,eps(ee));
    scale_y1 = 1/max(dy1);
    
%     y2 = minsup(m,p+0.5,eps(ee));
% %     r2 = r(m,p+0.5,eps(ee));
%     scale_y2 = 1/max(y2);
    
    
    dy2 = minsup(m,p+1,eps(ee));
%     r3 = r(m,p+1,eps(ee));
    scale_y2 = 1/max(dy2);
    
    w1 = wght(m,p,eps(ee));
    w2 = wght(m,p+1,eps(ee));
    
%     dy1 = d2fdx2(m,p,eps(ee));
%     scale_dy1 = 1/max(y1);
    
%     dy2 = d2fdx2(m,p+0.5,eps(ee));
%     scale_dy2 = 1/max(y2);
    
%     dy3 = d2fdx2(m,p+1.0,eps(ee));
%     scale_dy3 = 1/max(y3);
%     grad_y1 = scale_y1*gradfx(m,p,eps(ee));
        
%     polynom = polyn(scale_y1,eps(ee));
%     r = roots(polynom);
    
% y1 = log(abs(curv(m,p,eps(ee))));
%     y2 = log(abs(curv2(m,p,eps(ee))));
% y2 = func(m,p,eps(ee));
%     plot(m,y,'b');hold on
%     y1 = approx(m,p,eps(ee));
%     y2 = d2fdx2(m,p,eps(ee));
    
    [ax,h1,h2] = plotyy(m,scale_y1*dy1,m,w1/max(w1));
%     plot(m,scale_y1*y1,'b',m,scale_y2*y2,'r--');hold on
    set(h2,'LineStyle','--','color','blue')
%     alpha = minsup(alpha,p,eps(ee));

%     figure;plot(grad_y1)
%     x_inflex = inflex(p,eps(ee));
%     x_curv_y1 = (alpha - 2*eps(ee) + sqrt( (2*eps(ee)-alpha)^2 - 4*(eps(ee)^2 + alpha*eps(ee))))/2;
%     x_curv_y1 = dkdt(p,eps(ee));
%     x_curv_y2 = dkdt(p+0.5,eps(ee)); 
%     x_curv_y3 = dkdt(p+1,eps(ee));   
%     plot(ax(1),max(r),scale_y1*minsup(max(r),p,eps(ee)),'ro');

    hold(ax(1));
%     axis(ax(1),'equal')
%     xlim(ax(2),[0 0.5])
    plot(ax(1),m,scale_y2*dy2,'r');
    
    
    hold(ax(2));
    plot(ax(2),m,w2/max(w2),'r--',[eps(ee) eps(ee)],[0 1],':','LineWidth',3);
    
%     xlim(ax(2),[0 eps(ee)]);
%     xlim(ax(1),[0 eps(ee)]);
%     plot(ax(1),x_curv_y2,scale_y2*minsup(x_curv_y2,p+0.5,eps(ee)),'ro');
%     
%     plot(ax(1),m,scale_y3*y3,'b');
%     plot(ax(1),x_curv_y3,scale_y3*minsup(x_curv_y3,p+1,eps(ee)),'go');
    
%     xlim(ax(2),[0 0.5])
%     axis(ax(2),'equal')
    
    
%     xlim([0 0.1])
%     axis equal
    
%     figure
%     func_A = d3fdx3_A(m,p,eps(ee)) - d3fdx3_B(m,p,eps(ee));
%     plot(m,func_A,'r');hold on
    grid on
%     func_B = ;
%     plot(m,func_B,'g');
    
    
%     y2 = minsup(p);
%     norm(y2) 
%     plot(x,y2,'g:');hold on
    
end
    