function [eps_p,eps_q] = get_eps(m,pct_i,Gx,Gy,Gz)
% Function [eps_p,eps_q] = get_eps(m,pct_i,Gx,Gy,Gz)
% Finds epsilon values for lp and lq norm based on the distribution of
% model parameters.
% Input:
% m : model from the l2-norm inversion
% pct_i: Base percentile increment used to compute curvature
% Gx: Gradient operator in the x-direction. 
% Gy (optional): Gradient operator in the y-direction. Assume 1D problem if
% empty
% Gz (optional): Gradient operator in the z-direction. Assume 2D problem if
% empty


%% SCRIPT STARTS HERE
% Check if 1D, 2D or 3D problem


if ~isempty(Gz) && ~isempty(Gy)
    
    gradm = sqrt( (Gx * m).^2 + (Gy * m).^2 + (Gz * m).^2 );

elseif isempty(Gz) && ~isempty(Gy)
    
    gradm = sqrt( (Gx * m).^2 + (Gy * m).^2 );
    
else
    
    gradm = abs(Gx * m);
    
end

% gradm = gradm(gradm~=0);
% m = m(m~=0);

mcell = length(m);

[msort,idx] = sort(abs(m));

figure(100); 
subplot(1,2,1)
plot(msort);
hold on
axis square

pct = 0:pct_i:100;
mpct = prctile(abs(m),pct);
pctm = floor(mcell/100)*pct_i;

plot(0:pctm:mcell,mpct,'rx')
%                             ylim([0 100])
%                             set(figure(100), 'Position', [50 200 750 750])
%                             plot(x,model,'LineWidth',2);axis([x(1) x(end) -0.1 0.5]);hold on
%                             axis square
grid on

%% Eps for model gradient term

mcell = length(gradm);

dm2d2x = (mpct(1:end-2) + mpct(3:end) - 2*mpct(2:end-1));
dmdx = ((mpct(2:end-1)+mpct(3:end))/2 - (mpct(1:end-2)+mpct(2:end-1))/2);
curv = abs(dm2d2x./(1 + dmdx.^2).^(3/2));
mpct = mpct(2:end-2);
% Find the point of maximum increase in curvature
[~,idx] = max(abs(curv(1:end-1) - curv(2:end)));

plot(idx*pctm,mpct(idx),'bo');
set(gca,'XTick',0:pctm:mcell)
set(gca,'XTickLabel',pct)
xlabel('Percentile')
ylabel('$m$','interpreter','latex','FontSize',12)

eps_p = mpct(idx);%prctile(abs(m),);
text(idx*pctm,mpct(idx),'$\epsilon_p$','interpreter','latex','FontSize',12,'HorizontalAlignment','left','VerticalAlignment','top')

[msort,idx] = sort(abs(gradm));
subplot(1,2,2)
plot(msort);
hold on
axis square
%                             pct = 0:10:100;
mpct = prctile(abs(gradm),pct);

plot(0:pctm:mcell,mpct,'rx')
%                             ylim([0 100])
%                             set(figure(100), 'Position', [50 200 750 750])
%                             plot(x,model,'LineWidth',2);axis([x(1) x(end) -0.1 0.5]);hold on
%                             axis square
grid on

dm2d2x = (mpct(1:end-2) + mpct(3:end) - 2*mpct(2:end-1));
dmdx = ((mpct(2:end-1)+mpct(3:end))/2 - (mpct(1:end-2)+mpct(2:end-1))/2);
curv = abs(dm2d2x./(1 + dmdx.^2).^(3/2));
mpct = mpct(2:end-2);
% Find the point of maximum increase in curvature
[~,idx] = max(abs(curv(1:end-1) - curv(2:end)));

plot(idx*pctm,mpct(idx),'bo');

eps_q = mpct(idx);%prctile(abs(m),);
ylabel('$\nabla m$','interpreter','latex','FontSize',12)
text(idx*pctm,mpct(idx),'$\epsilon_q$','interpreter','latex','FontSize',12,'HorizontalAlignment','left','VerticalAlignment','top')
%                     eps_s = std(m)/4;
set(gca,'XTick',0:pctm:mcell)
set(gca,'XTickLabel',pct)
xlabel('Percentile')