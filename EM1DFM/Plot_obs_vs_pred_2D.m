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
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Bookpurnong\Modelling\EM1DFM';

obsfile = 'Inv_Data_XYZ_ds15m.dat';

predfile = 'Inv_PRED_iter4.pre';
% predfile = 'Inv_PRED_iter6.pre';

%% Load data
set(figure, 'Position', [25 50 1800 900])
data = load([work_dir '\' obsfile]);
pred = load([work_dir '\' predfile]);

freq = unique(data(:,4));
nfreq = length(freq);
for ii = 1 : nfreq
    
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
    
    Y = flipud(Y);
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
    
    
    set(figure, 'Position', [25 50 1800 900])
    subplot(231)
    h = imagesc(x,y,data_R);
    set(h,'alphadata',~isnan(data_R));
    title(['\bf' num2str(freq(ii)) 'Real'])
    caxis([min(sub_data(:,5)) max(sub_data(:,5))])
    colorbar
    colormap(jet)
    axis equal
    grid on
    
    subplot(232)
    h=imagesc(x,y,pred_R);
    set(h,'alphadata',~isnan(data_R));
    title('\bfPredicted Real')
    caxis([min(sub_data(:,5)) max(sub_data(:,5))])
    colorbar
    axis equal
    grid on
    
    subplot(233)
    h=imagesc(x,y,( data_R-pred_R ) ./ uncert_R );
    set(h,'alphadata',~isnan(data_R));
    title('\bfNormalized Residual')
    caxis([-4 4])
    colorbar
    axis equal tight
    grid on

    subplot(234)
    h=imagesc(x,y,data_I);
    set(h,'alphadata',~isnan(data_I));
    title(['\bf' num2str(freq(ii)) ' Imag'])
    caxis([min(sub_data(:,6)) max(sub_data(:,6))])
    colorbar
    axis equal
    grid on
    
    subplot(235)
    h=imagesc(x,y,pred_I);
    set(h,'alphadata',~isnan(data_I));
    title('\bfPredicted Imag')
    caxis([min(sub_data(:,6)) max(sub_data(:,6))])
    colorbar
    axis equal
    grid on
    
    subplot(236)
    h=imagesc(x,y, ( data_I-pred_I ) ./ uncert_I);
    set(h,'alphadata',~isnan(data_I));
    caxis([-4 4])
    title('\bfNormalized Residual')
    colorbar
    axis equal
    grid on

    figure(1)
    subplot(2,nfreq,ii)
    h=imagesc(x,y,data_R);
    set(h,'alphadata',~isnan(data_R));
    title(['\bf' num2str(freq(ii)) ' Real'])
    caxis([min(sub_data(:,5)) max(sub_data(:,5))])
    colorbar
    colormap(jet)
    axis equal
    grid on
    
    subplot(2,nfreq,ii+nfreq)
    h=imagesc(x,y,data_I);
    set(h,'alphadata',~isnan(data_I));
    title(['\bf' num2str(freq(ii)) ' Imag'])
    caxis([min(sub_data(:,6)) max(sub_data(:,6))])
    colorbar
    axis equal
    grid on
end

