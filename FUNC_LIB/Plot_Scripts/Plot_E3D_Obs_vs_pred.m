% Plot E3D time channels for observed and predicted data
work_dir    = 'C:\LC\Private\dominiquef\Projects\4414_Minsim\Modeling\VTEM\Inv3_10pc_1em13';
addpath ..\..\FUNC_LIB;
dsep = '\';

close all
%% INPUT PARAMETERS
% obsE3D  = 'nanoTEM_new_topo.obs';
timefile = '..\Times.dat';
obsfile = '..\VTEM_FLIN_EM1D_flt20m.obs';
prefile = 'Inv_PRED_iter7.pre';

% Number of time channels per figure (tcpf)
tcpf = 8;

%% SCRIPT STARTS HERE
% Load E3DD obs file for time channels
% [~,d] = read_E3DTD_obs([work_dir dsep obsE3D]);
% 
% tc = d{1}{1}(:,4);
tc = load([work_dir dsep timefile]);
% Load EM1DTM files
dobs = convert_EM1D_2_E3D_pred([work_dir dsep obsfile],tc);

dpre = convert_EM1D_2_E3D_pred([work_dir dsep prefile],tc); 

%%
tc = unique(dobs(:,4));

% Find non-zero time channels
index = zeros(1,length(tc));
for ii = 1 : length(tc)
    
    % Only plot 
    if sum(~isnan(dobs(dobs(:,4)==tc(ii),end)))~=0
        
        index(ii) = 1;
        
    end
    
end

% Figure the limits of plots and create grids
ntimes = sum(index);
index = find(index);


xmin = min(dobs(:,1));
xmax = max(dobs(:,1));

ymin = min(dobs(:,2));
ymax = max(dobs(:,2));

dx = 10;
dy = 10;

x = xmin:dx:xmax;
y = ymin:dy:ymax;
[Y,X] = ndgrid(y,x);

Y = flipud(Y);

%% Scart ploting

count = 0;
figcount = 0;
for ii = 1 : ntimes
    
    
    
    % If more than 6 time channels, plot on a new figure
    if mod(count,tcpf) == 0
        
        set(figure(2*(figcount)+1), 'Position', [25 50 1800 900]);
%         set(figure(2*(figcount)+2), 'Position', [25 50 1800 900]);
        count = 0;
        figcount = figcount+1;
    end
    
    count = count + 1;
    % Extract data
    subdata = dobs( dobs(:,4) == tc(index(ii)),: );
    subpred = dpre( dpre(:,4) == tc(index(ii)),: );
    
    
%     
%     F_o = scatteredInterpolant(subdata(:,1),subdata(:,2),abs(log10(subdata(:,end))),'linear');
%     F_p = scatteredInterpolant(subpred(:,1),subpred(:,2),abs(log10(subpred(:,end))),'linear');
%     
%     data_2D = F_o(Y,X);
%     pred_2D = F_p(Y,X);
    
    data_2D = flipud(griddata(subdata(:,2),subdata(:,1),((subdata(:,end))),Y,X,'natural'));
%     pred_2D = flipud(griddata(subpred(:,2),subpred(:,1),((subpred(:,end))),Y,X,'natural'));
    
    figure(2*(figcount-1)+1)
    subplot(2,tcpf/2,count)
    h = imagesc(x,y,data_2D); hold on
    contour(x,y,data_2D,'k')
    contour(x,y,data_2D,[0],'r')
%     scatter(subpred(:,1),subpred(:,2),1,'k.')
    set(h,'alphadata',~isnan(data_2D));
    colormap(jet)
    colorbar
    title(['\bf T:' num2str(tc(index(ii)))])
    caxis([min(data_2D(:)) max(data_2D(:))])
    axis equal tight
    set(gca,'Ydir','normal')
    set(gca,'YTickLabel',[])
    
%     figure(2*(figcount-1)+2)
%     subplot(2,tcpf/2,count)
%     h = imagesc(x,y,(pred_2D));hold on
%     contour(x,y,pred_2D,'k')
%     contour(x,y,data_2D,[0],'r')
% %     scatter(subpred(:,1),subpred(:,2),1,'k.')
%     set(h,'alphadata',~isnan(pred_2D));
%     colormap(jet)
%     colorbar
%     title(['\bf Pred:' num2str(tc(index(ii)))])
%     caxis([min(data_2D(:)) max(data_2D(:))])
%     axis equal tight
%     set(gca,'YTickLabel',[])
%     set(gca,'Ydir','normal')
end