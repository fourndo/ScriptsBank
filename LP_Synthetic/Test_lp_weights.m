% Plotting the gradient of the objective function as a function of epsilon
clear all
close all

addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB

x=-1:1e-5:1;

% eps = 10.^-[0 0.5 0.1 0.01];
eps = [5 0.5 0.25 0.01];
p = [0 0.5 1 1.5 2];

pcode1 = { 'k-' 'k--' 'k:' 'k.-' };

set(figure, 'Position', [50 0 775 1000]);
% plot_h = tight_subplot(2, 1, [.1 .1],[.1 0.1]);
% plot_h = subplot(2, 1);
subplot(3, 1 , 1);
ylim([-1e-3 1]);
yy = get(gca,'Ylim');
plot([1e-2 1e-2],yy,'r-.');hold on
plot([-1e-2 -1e-2],yy,'r-.')
AX = legend('\epsilon = 10^{-2}','Location','SouthEast'); 
for pp = 1 : length(p)
    
    phim = x.^2 ./ (x.^(2) + (1e-2).^2).^(1 - p(pp)/2);

    
    plot(x, phim,'k','LineWidth',p(pp)+0.25);hold on
    text(0.4-0.05,0.4.^2 ./ (0.4.^(2) + (1e-2).^2).^(1 - p(pp)/2),['p =' num2str(p(pp))],...
        'BackgroundColor','w','EdgeColor','k')
    
    grid on
end


annotation('rectangle',...
    [0.482 0.71 0.07 0.035],'LineWidth',1);

ylim([-1e-3 1]);
ylabel('$\phi_m$', 'interpreter', 'latex','FontSize',16)
set(get(gca,'YLabel'),'Rotation',360);
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') - [0.001 0.0 0.00]);

axes('Position',[.32 .51 .4 .2])
xlim([-1e-2 1e+2]);
x=-2e-2:1e-5:2e-2;
yy = get(gca,'Ylim');
plot([1e-2 1e-2],yy,'r-.');hold on
plot([-1e-2 -1e-2],yy,'r-.')
for pp = 1 : length(p)
        
    phim = x.^2 ./ (x.^(2) + (1e-2).^2).^(1 - p(pp)/2);
    box on
    plot(x, phim,'k','LineWidth',p(pp)+0.25);hold on
    
    if p(pp) == 0
        
        text(-0.001-3e-4,0.001.^2 ./ (0.001.^(2) + (1e-2).^2).^(1 - p(pp)/2)-5e-4,['p =' num2str(p(pp))],...
        'BackgroundColor','w','EdgeColor','k')
    
    elseif p(pp) == 0.5
        
        text(-0.00325-5e-3,0.00325.^2 ./ (0.00325.^(2) + (1e-2).^2).^(1 - p(pp)/2)-5e-4,['p =' num2str(p(pp))],...
        'BackgroundColor','w','EdgeColor','k')
    
    elseif p(pp) == 1
        
        text(-0.0125-3e-3,0.0125.^2 ./ (0.0125.^(2) + (1e-2).^2).^(1 - p(pp)/2)-5e-4,['p =' num2str(p(pp))],...
        'BackgroundColor','w','EdgeColor','k')
        
    else
        
        text(-0.02+0.001,0.02.^2 ./ (0.02.^(2) + (1e-2).^2).^(1 - p(pp)/2),['p =' num2str(p(pp))],...
        'BackgroundColor','w','EdgeColor','k')

    end
    grid on
    
end
set(gca,'linewidth',2)
ylim([-1e-3 0.01]);
ylabel('$\phi_m$', 'interpreter', 'latex','FontSize',16)
set(get(gca,'YLabel'),'Rotation',360);
ylabh = get(gca,'YLabel');
text(0,-4e-3,'(a)', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center');
xlabel('$m$', 'interpreter', 'latex','FontSize',14)
set(ylabh,'Position',get(ylabh,'Position') - [0.001 0.0 0.00]);
% axis square



axes('Position',[.10 .10 .35 .3])

x=-1:1e-5:1;

epsl = 1e-2;

x=-1:1e-5:1;
for pp = 1 : length(p)
    
    r = 1./(x.^(2) + epsl.^2).^(1-p(pp)/2);
    
%     phim = phim./max(x.*phim);
    semilogy(x, r,'k','LineWidth',p(pp)+0.25);hold on
    text(0.05,1 ./ (0.05.^(2) + (1e-2).^2).^(1 - p(pp)/2),['p =' num2str(p(pp))],...
        'BackgroundColor','w','EdgeColor','k')
%     plot(eps(ee),eps(ee) ./ (2* eps(ee).^2),'b*');

end
xlim([-0.01 0.25])
ylim([0.5 1e+4])
grid on
yy = get(gca,'Ylim');
plot([1e-2 1e-2],yy,'r-.')
ylabel('$\mathbf{R}$', 'interpreter', 'latex','FontSize',16)
set(get(gca,'YLabel'),'Rotation',360);
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') - [0.01 0 0]);
xlabel('$m$', 'interpreter', 'latex','FontSize',14)
text(.108,0.075,'(c)', 'interpreter', 'latex','FontSize',14);


axes('Position',[.60 .10 .35 .3])

x=-1:1e-5:1;
for pp = 1 : length(p)
    
    r = x./(x.^(2) + epsl.^2).^(1-p(pp)/2);
    
%     phim = phim./max(x.*phim);
    semilogy(x, r,'k','LineWidth',p(pp)+0.25);hold on
    text(0.05,0.05 ./ (0.05.^(2) + (1e-2).^2).^(1 - p(pp)/2),['p =' num2str(p(pp))],...
        'BackgroundColor','w','EdgeColor','k')
%     plot(eps(ee),eps(ee) ./ (2* eps(ee).^2),'b*');

end

xlim([-0.01 0.25])
ylim([1e-3 1e+2])
grid on
yy = get(gca,'Ylim');
plot([10^(-2) 10^(-2)],yy,'r-.')
ylabel('$\mathbf{g}(m)$', 'interpreter', 'latex','FontSize',16)
set(get(gca,'YLabel'),'Rotation',360);
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') - [0.01 0 0]);
xlabel('$m$', 'interpreter', 'latex','FontSize',14)
text(.108,1e-4,'(b)', 'interpreter', 'latex','FontSize',14);

%% Second plot
set(figure, 'Position', [50 0 775 400]);

axes('Position',[.60 .20 .35 .75])

x=-10:1e-5:10;

epsl = 1;

for pp = 1 : length(p)
    
    r = 1./(x.^(2) + epsl.^2).^(1-p(pp)/2);
    
    dphi_dm = (x).*r;
%     scale = 1./(max(sqrt(r)));
    eta = 2*max(x)*epsl;
    scale = (eta)^((1-p(pp)/2));
    dphi_dm = dphi_dm * scale;
    
%     dphi_dm = (x*(1e-2).^(1-p(pp)/2))./(x.^(2) + 1e-2.^2).^(1-p(pp)/2);
    
%     phim = phim./ ( max(x.*phim));
%     phim = phim.^0.6;


     plot(x, dphi_dm,'k','LineWidth',p(pp)+0.25);hold on
     
     if p(pp) == 2 || p(pp)==1.5 || p(pp)==1
         
              text(0.7*max(x),(0.7*max(x)*(eta).^(1-p(pp)/2)) ./ ((0.7*max(x)).^(2) + (epsl).^2).^(1 - p(pp)/2),['p =' num2str(p(pp))],...
        'BackgroundColor','w','EdgeColor','k')
         
     elseif p(pp)==0
         
     text(0.1*max(x),(0.1*max(x)*(eta).^(1-p(pp)/2)) ./ ((0.1*max(x)).^(2) + (epsl).^2).^(1 - p(pp)/2),['p =' num2str(p(pp))],...
        'BackgroundColor','w','EdgeColor','k')
    
     else
    
      text(0.7*max(x),(0.7*max(x)*(eta).^(1-p(pp)/2)) ./ ((0.7*max(x)).^(2) + (epsl).^2).^(1 - p(pp)/2),['p =' num2str(p(pp))],...
        'BackgroundColor','w','EdgeColor','k')
     end
end
text(sqrt(eta),sqrt(eta),['$\sqrt{\eta}$'],'BackgroundColor','w','EdgeColor','k','interpreter', 'latex','FontSize',14,'HorizontalAlignment','left','VerticalAlignment','bottom')
plot(sqrt(eta),sqrt(eta),'ko')
% subplot(2,1,2)
% for ee = 1 : length(eps)
%     
%     r1 = 1 ./ (x.^(2) + eps(ee).^2);
%     
% %     phim = phim./max(x.*phim);
%     semilogy(x, r1,'k-','LineWidth',3-abs(log10(eps(ee)))/2);hold on
% %     plot(eps(ee),eps(ee) ./ (2* eps(ee).^2),'b*');
% 
% end
xlim([0 max(x)])
ylim([0 max(x)])
grid on
yy = get(gca,'Ylim');
plot([10^(-2) 10^(-2)],yy,'r-.')
ylabel('$\mathbf{\hat g}(m)$', 'interpreter', 'latex','FontSize',16)
set(get(gca,'YLabel'),'Rotation',360);
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') - [0.01 0 0]);
xlabel('$m$', 'interpreter', 'latex','FontSize',14)
text(.108,-0.25,'(b)', 'interpreter', 'latex','FontSize',14);

axes('Position',[.10 .20 .35 .75])


x=-1:1e-5:1;
for pp = 1 : length(p)
    
    r = x./(x.^(2) + 1e-2.^2).^(1-p(pp)/2);
    
%     phim = phim./max(x.*phim);
    semilogy(x, r,'k','LineWidth',p(pp)+0.25);hold on
    text(0.05,0.05 ./ (0.05.^(2) + (1e-2).^2).^(1 - p(pp)/2),['p =' num2str(p(pp))],...
        'BackgroundColor','w','EdgeColor','k')
%     plot(eps(ee),eps(ee) ./ (2* eps(ee).^2),'b*');

end

xlim([-0.01 0.25])
ylim([1e-3 1e+2])
grid on
yy = get(gca,'Ylim');
plot([10^(-2) 10^(-2)],yy,'r-.')
ylabel('$\mathbf{g}(m)$', 'interpreter', 'latex','FontSize',16)
set(get(gca,'YLabel'),'Rotation',360);
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') - [0.01 0 0]);
xlabel('$m$', 'interpreter', 'latex','FontSize',14)
text(.108,1e-4,'(a)', 'interpreter', 'latex','FontSize',14);

%% Third plot
set(figure, 'Position', [50 0 775 1000]);
x=0:1e-5:1;
% splot2 = tight_subplot(1,2,0.05 , 0.05,0.05);
% axes(splot2(1));

% subplot(3,1,1)
% for ee = 1 : length(eps)
%     
%     r1 = 1 ./ (x.^(2) + eps(ee).^2);
%     
% %     phim = phim./max(x.*phim);
%     semilogy(x, r1,'k-','LineWidth',3-abs(log10(eps(ee)))/2);hold on
% %     plot(eps(ee),eps(ee) ./ (2* eps(ee).^2),'b*');
% 
% end
% plot([-1 1],[0 0],'k-');
% plot([-1 1],[0 0],'k--');
% % semilogy(-1:0.05:1, 1,'bo');
% grid on
% xlim([-0.1 0.1])
% ylim([1e-1 1e+4])
% % title('$|\frac{\partial \phi(m,\epsilon)}{\partial m}| , p = 0$', 'interpreter', 'latex','FontSize',18)
% % aa=get(gca,);
% set(get(gca,'YLabel'),'Rotation',360);
% text(-0.01,5e-3,'(a)', 'interpreter', 'latex','FontSize',14);
% 
% subplot(3,1,1)
% for ee = 1 : length(eps)
%     
%     r2 = 1 ./ (x.^(2) + eps(ee).^2).^0.5;
%     
% %     phim = phim./max(x.*phim);
%     semilogy(x, r2,'k--','LineWidth',3-abs(log10(eps(ee)))/2);hold on
% %     plot(eps(ee),eps(ee) ./ (2* eps(ee).^2).^0.5,'r*');
% end
% 
% AX= legend('\epsilon = 10^{-0}','\epsilon = 10^{-1}','\epsilon = 10^{-3}','p = 0','p = 1','p = 2','Location','NorthWest');
% ylabel('$r_i(m,\epsilon)$', 'interpreter', 'latex','FontSize',14)
% % xlim([-0.25 0.25])
% ylim([1e-1 1e+7])
% % set(gca,'YTickLabel',[]);
% grid on
% xlabel('$m$', 'interpreter', 'latex','FontSize',14)
% set(figure, 'Position', [10 250 1000 500])

% splot2 = tight_subplot(1,2,0.01 , 0.1,0.1);
subplot(2,1,2)
for ee = 1 : length(eps)
    
    dphidm = x ./ (x.^(2) + eps(ee).^2);
    
%     phim = phim./max(x.*phim);
    plot(x, dphidm/max(dphidm),'k-','LineWidth',2);hold on
    plot(eps(ee),eps(ee)/(2*eps(ee)^2)/max(dphidm),'ko','MarkerFaceColor','w');

    text(eps(ee)+0.02,eps(ee)/(2*eps(ee)^2)/max(dphidm)+0.025,['\epsilon =' num2str(eps(ee))],...
        'BackgroundColor','w','EdgeColor','k')
end
plot(0.75,0.76,'ko','MarkerFaceColor','w');

text(0.75+0.02,0.76+0.025,['\epsilon =' num2str(eps(1))],...
    'BackgroundColor','w','EdgeColor','k')
    
xlim([0 1])
ylim([0 1.1])
% AX= legend('\epsilon = 10^{-0}','\epsilon = 10^{-1}','\epsilon = 10^{-2}','\epsilon = 10^{-3}','p = 0','p = 1','p = 2','Location','SouthEast');
ylabel('$\theta^2\mathbf{g_p}(m)$', 'interpreter', 'latex','FontSize',14)
set(get(gca,'YLabel'),'Rotation',360);
xlabel('$m$', 'interpreter', 'latex','FontSize',14)
text(0.48,-.18,'(b)', 'interpreter', 'latex','FontSize',14);
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') - [0.05 0 0]);
title('$p = 0$', 'interpreter', 'latex','FontSize',16)
% semilogy(0:0.01:1, 0:0.01:1,'bo');
grid on
% plot([0 1],[0.5 0],'b--')

% aa=get(gca,);
% set(get(gca,'YLabel'),'Rotation',360);

subplot(2,1,1)
for ee = 1 : length(eps)
    
    dphidm = x ./ (x.^(2) + eps(ee).^2).^0.5;
    
%     phim = phim./max(x.*phim);
    plot(x, dphidm/max(dphidm),'k','LineWidth',2);hold on
    plot(eps(ee),eps(ee)/(2*eps(ee)^2)^0.5/max(dphidm),'ko','MarkerFaceColor','w');
    text(eps(ee)+0.02,eps(ee)/(2*eps(ee)^2)^0.5/max(dphidm)+0.025,['\epsilon =' num2str(eps(ee))],...
        'BackgroundColor','w','EdgeColor','k')
end
plot(0.75,0.756,'ko','MarkerFaceColor','w');

text(0.75+0.02,0.756+0.025,['\epsilon =' num2str(eps(1))],...
    'BackgroundColor','w','EdgeColor','k')
% AX= legend('\epsilon = 10^{-0}',,'\epsilon = 10^{-1}','\epsilon = 10^{-2}','\epsilon = 10^{-3}','p = 0','p = 1','p = 2','Location','SouthEast');
ylabel('$\theta^2\mathbf{g_p}(m)$', 'interpreter', 'latex','FontSize',14)
set(get(gca,'YLabel'),'Rotation',360);
xlim([0 1])
ylim([0 1.1])
title('$p = 1$', 'interpreter', 'latex','FontSize',16)
% set(gca,'YTickLabel',[]);
grid on
xlabel('$m$', 'interpreter', 'latex','FontSize',14)
text(0.48,-0.18,'(a)', 'interpreter', 'latex','FontSize',14);
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') - [0.05 0 0]);
% plot([0 1],[0.5 0],'b--')

% for ee = 1 : length(eps)
%     
%     subplot(2,1,2)
%     dphidm = eps(ee) ./ (eps(ee).^(2) + eps(ee).^2);
%     plot(eps(ee),dphidm/max(dphidm),'x');
%     
%     subplot(2,1,1)
%     dphidm = eps(ee) ./ (eps(ee).^(2) + eps(ee).^2).^0.5;
%     plot(eps(ee),dphidm/max(dphidm),'x');
% %     plot(eps(ee),eps(ee) ./ (2* eps(ee).^2).^0.5,'r*');
% end