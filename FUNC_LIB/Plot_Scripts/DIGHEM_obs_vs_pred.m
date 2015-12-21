% Function: plot_obs_vs_pred_2D
% Plot 2D interpolated observed vs predicted data from EM1D inversion
% Line number is infered from the x,y location of each sounding, assuming
% that the data is ordered as it was surveyed.
% 
% Script assumes that each frequency has the same number of soundings. Will
% need more work to figure out a way to drop this assumption.
% 
% INPUTS:
% work_dir: directory for the files
% obsfile: Observed data matrix
% prefile: Predicted data matrix
% linefile(optional) : list of line names used for labeling
%
% Last update: August 23, 2015
% D Fournier
% fourndo@gmail.com

clear all
close all

'.\Functions';

%% INPUT FILES
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\DIGHEM\1D';

obsfile = 'Inv_Data_XYZ_ds15m.dat';

% predfile = 'FWR_pred_final.dat';
predfile = 'Inv_PRED_iter25.pre';
outline = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\EOSC 556b (2013) EM\TKC_Project\Modelling\Geology_Objects\DO27_Outline.dat';
hydro = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\EOSC 556b (2013) EM\TKC_Project\Modelling\Geology_Objects\Lake_trace.dat';

Out_grav = load([outline]);
Out_lake = load([hydro]);

%% Load data
set(figure, 'Position', [25 50 1800 900])
set(figure, 'Position', [25 50 1800 900])
data = load([work_dir '\' obsfile]);
pred = load([work_dir '\' predfile]);

freq = unique(data(:,4));

for ii = 1 : length(freq)
    
    indx = data(:,4) == freq(ii);
    sub_data = data(indx,:);
    
    indx = pred(:,4) == freq(ii);
    sub_pred = pred(indx,:);

    
    if ii == 1
    %% Set coordinates for plot
    xmin = min(sub_data(:,1));
    xmax = max(sub_data(:,1));

    ymin = min(sub_data(:,2));
    ymax = max(sub_data(:,2));

    dx = 10;
    dy = 10;

    x = xmin:dx:xmax;
    y = ymin:dy:ymax;
    [Y,X] = ndgrid(y,x);
    
%     Y = flipud(Y);
    end
    
%     F_I = TriScatteredInterp(sub_pred(:,2),sub_pred(:,1),sub_pred(:,6),'natural');
%     F_R = TriScatteredInterp(sub_pred(:,2),sub_pred(:,1),sub_pred(:,5),'natural');
    
    pred_I = griddata(sub_pred(:,2),sub_pred(:,1),sub_pred(:,6),Y,X);
    pred_R = griddata(sub_pred(:,2),sub_pred(:,1),sub_pred(:,5),Y,X);
    
%     F_I = TriScatteredInterp(sub_data(:,2),sub_data(:,1),sub_data(:,6),'natural');
%     F_R = TriScatteredInterp(sub_data(:,2),sub_data(:,1),sub_data(:,5),'natural');
    
    data_I = griddata(sub_data(:,2),sub_data(:,1),sub_data(:,6),Y,X);
    data_R = griddata(sub_data(:,2),sub_data(:,1),sub_data(:,5),Y,X);
    
%     U_I = TriScatteredInterp(sub_data(:,2),sub_data(:,1),sub_data(:,8),'natural');
%     U_R = TriScatteredInterp(sub_data(:,2),sub_data(:,1),sub_data(:,7),'natural');
    
    uncert_I = griddata(sub_data(:,2),sub_data(:,1),sub_data(:,8),Y,X);
    uncert_R = griddata(sub_data(:,2),sub_data(:,1),sub_data(:,7),Y,X);
    
    res_I = (data_I - pred_I) ./ uncert_I;
    res_R = (data_R - pred_R) ./ uncert_R;
%%    
    figure(1)
    ax1 = axes('Position',[-0.1+(ii-1)*0.325 0.55 0.4 0.4]);
    h = imagesc(x,y,data_R); hold on
    set(h,'alphadata',~isnan(data_R));
    title(['\bf Obs: ' num2str(freq(ii)) ' Real'])
    caxis([min(sub_data(:,5)) max(sub_data(:,5))])
    colormap(jet)
    colorbar
    set(ax1,'YTickLabel',[])
    set(ax1,'YDir','normal')
    set(ax1,'XTickLabel',[])
    plot(Out_grav(:,1),Out_grav(:,2),'k--','LineWidth',2) 
    plot(Out_lake(:,1),Out_lake(:,2),'r.','MarkerSize',2) 
    axis equal tight
    axis([xmin xmax ymin ymax])
    grid on
    
    ax1 = axes('Position',[0.05+(ii-1)*0.325 0.55 0.4 0.4]);
    h=imagesc(x,y,data_I); hold on 
    set(h,'alphadata',~isnan(data_I));
    title(['\bf Obs: ' num2str(freq(ii)) ' Imag'])
    caxis([min(sub_data(:,6)) max(sub_data(:,6))])
    colorbar
    set(ax1,'YTickLabel',[])
    set(ax1,'YDir','normal')
    set(ax1,'XTickLabel',[])
    plot(Out_grav(:,1),Out_grav(:,2),'k--','LineWidth',2) 
    plot(Out_lake(:,1),Out_lake(:,2),'r.','MarkerSize',2) 
    axis equal tight
    axis([xmin xmax ymin ymax])
    grid on
    
    ax1 = axes('Position',[-0.1+(ii-1)*0.325 0.05 0.4 0.4]);
    h = imagesc(x,y,pred_R); hold on
    set(h,'alphadata',~isnan(data_R));
    title(['\bf Pred: ' num2str(freq(ii)) ' Real'])
    caxis([min(sub_data(:,5)) max(sub_data(:,5))])
    colorbar
    set(ax1,'YTickLabel',[])
    set(ax1,'YDir','normal')
    set(ax1,'XTickLabel',[])
    plot(Out_grav(:,1),Out_grav(:,2),'k--','LineWidth',2) 
    plot(Out_lake(:,1),Out_lake(:,2),'r.','MarkerSize',2) 
    axis equal tight
    axis([xmin xmax ymin ymax])
    grid on
    
    ax1 = axes('Position',[0.05+(ii-1)*0.325 0.05 0.4 0.4]);
    h=imagesc(x,y,pred_I); hold on
    set(h,'alphadata',~isnan(data_I));
    title(['\bf Pred: ' num2str(freq(ii)) ' Imag'])
    caxis([min(sub_data(:,6)) max(sub_data(:,6))])
    colorbar
    set(ax1,'YTickLabel',[])
    set(ax1,'YDir','normal')
    set(ax1,'XTickLabel',[])
    plot(Out_grav(:,1),Out_grav(:,2),'k--','LineWidth',2) 
    plot(Out_lake(:,1),Out_lake(:,2),'r.','MarkerSize',2) 
    axis equal tight
    axis([xmin xmax ymin ymax])
    grid on


    
    figure(2)
    colormap(jet)
    ax1 = axes('Position',[-0.1+(ii-1)*0.325 0.55 0.4 0.4]);
    h = imagesc(x,y,data_R); hold on
    set(h,'alphadata',~isnan(data_R));
    title(['\bf Obs: ' num2str(freq(ii)) ' Real'])
    caxis([min(sub_data(:,5)) max(sub_data(:,5))])
    colorbar
    set(ax1,'YTickLabel',[])
    set(ax1,'YDir','normal')
    set(ax1,'XTickLabel',[])
    plot(Out_grav(:,1),Out_grav(:,2),'k--','LineWidth',2) 
    plot(Out_lake(:,1),Out_lake(:,2),'r.','MarkerSize',2) 
    axis equal tight
    axis([xmin xmax ymin ymax])
    grid on
    
    ax1 = axes('Position',[0.05+(ii-1)*0.325 0.55 0.4 0.4]);
    h=imagesc(x,y,data_I); hold on
    set(h,'alphadata',~isnan(data_I));
    title(['\bf Obs: ' num2str(freq(ii)) ' Imag'])
    caxis([min(sub_data(:,6)) max(sub_data(:,6))])
    colorbar
    set(ax1,'YTickLabel',[])
    set(ax1,'YDir','normal')
    set(ax1,'XTickLabel',[])
    plot(Out_grav(:,1),Out_grav(:,2),'k--','LineWidth',2) 
    plot(Out_lake(:,1),Out_lake(:,2),'r.','MarkerSize',2) 
    axis equal tight
    axis([xmin xmax ymin ymax])
    grid on
    
    ax1 = axes('Position',[-0.1+(ii-1)*0.325 0.05 0.4 0.4]);
    h = imagesc(x,y,res_R); hold on
    set(h,'alphadata',~isnan(data_R));
    title(['\bf Pred: ' num2str(freq(ii)) ' Real'])
    caxis([-5 5])
    colorbar
    set(ax1,'YTickLabel',[])
    set(ax1,'YDir','normal')
    set(ax1,'XTickLabel',[])
    plot(Out_grav(:,1),Out_grav(:,2),'k--','LineWidth',2) 
    plot(Out_lake(:,1),Out_lake(:,2),'r.','MarkerSize',2) 
    axis equal tight
    axis([xmin xmax ymin ymax])
    grid on
    
    ax1 = axes('Position',[0.05+(ii-1)*0.325 0.05 0.4 0.4]);
    h=imagesc(x,y,res_I); hold on
    set(h,'alphadata',~isnan(data_I));
    title(['\bf Pred: ' num2str(freq(ii)) ' Imag'])
    caxis([-5 5])
    colorbar
    set(ax1,'YTickLabel',[])
    set(ax1,'YDir','normal')
    set(ax1,'XTickLabel',[])
    plot(Out_grav(:,1),Out_grav(:,2),'k--','LineWidth',2) 
    plot(Out_lake(:,1),Out_lake(:,2),'r.','MarkerSize',2) 
    axis equal tight
    axis([xmin xmax ymin ymax])
    grid on


%%
end

