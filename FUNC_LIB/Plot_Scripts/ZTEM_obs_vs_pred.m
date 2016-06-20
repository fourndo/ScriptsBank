% Function: ZTEM_obs_vs_pred
% Plot 2D interpolated observed vs predicted data from ZTEM observation.
%
% INPUTS:
% work_dir: directory for the files
% obsfile: Observed data file
% prefile: Predicted data file
%
% Last update: March 13, 2016
% D Fournier
% fourndo@gmail.com

clear all
close all

addpath '..\';

%% INPUT FILES
work_dir = 'C:\LC\Private\dominiquef\Projects\4414_Minsim\Modeling\ZTEM\FLIN';

obsfile = 'FLIN_ZTEM_rot_flt50m_v2.obs';

% predfile = 'FWR_pred_final.dat';
predfile = 'Run7_Start1em2_SMOOTHMOD\dpred_08.txt';


%% Load data

[data,ndv,dtype] = load_ZTEM_data([work_dir '\' obsfile]);
pred = load([work_dir '\' predfile]);

freq = unique(data(:,1));

for ii = 1 : length(freq)
    
    set(figure, 'Position', [25 50 1800 900])
    
    indx = data(:,1) == freq(ii);
    sub_data = data(indx,:);
    sub_data = sub_data(2:end,:);
    
    sub_pred = pred(indx,:);
    sub_pred =sub_pred(2:end,:);
    
    if ii == 1
        %% Set coordinates for plot
        xmin = min(sub_data(:,2));
        xmax = max(sub_data(:,2));

        ymin = min(sub_data(:,3));
        ymax = max(sub_data(:,3));

        dx = 20;
        dy = 20;

        x = xmin:dx:xmax;
        y = ymin:dy:ymax;
        [Y,X] = ndgrid(y,x);
    
%     Y = flipud(Y);
    end
    
%     F_I = TriScatteredInterp(sub_pred(:,2),sub_pred(:,1),sub_pred(:,6),'natural');
%     F_R = TriScatteredInterp(sub_pred(:,2),sub_pred(:,1),sub_pred(:,5),'natural');
    
    % Cycle through each component
    
    for jj = 1 : 2
        
        
        % Only keep values
        indx = sub_pred(:,4+(jj-1)) ~= str2num(ndv);
        pred_I = griddata(sub_pred(indx,2),sub_pred(indx,1),sub_pred(indx,5+2*(jj-1)),Y,X);
        
        indx = sub_pred(:,5+(jj-1)) ~= str2num(ndv);
        pred_R = griddata(sub_pred(indx,2),sub_pred(indx,1),sub_pred(indx,4+2*(jj-1)),Y,X);

    %     F_I = TriScatteredInterp(sub_data(:,2),sub_data(:,1),sub_data(:,6),'natural');
    %     F_R = TriScatteredInterp(sub_data(:,2),sub_data(:,1),sub_data(:,5),'natural');
        indx = sub_data(:,7+4*(jj-1)) ~= str2num(ndv);
        data_I = griddata(sub_data(indx,3),sub_data(indx,2),sub_data(indx,7+4*(jj-1)),Y,X);
        
        indx = sub_data(:,5+4*(jj-1)) ~= str2num(ndv);
        data_R = griddata(sub_data(indx,3),sub_data(indx,2),sub_data(indx,5+4*(jj-1)),Y,X);

    %     U_I = TriScatteredInterp(sub_data(:,2),sub_data(:,1),sub_data(:,8),'natural');
    %     U_R = TriScatteredInterp(sub_data(:,2),sub_data(:,1),sub_data(:,7),'natural');

        uncert_I = griddata(sub_data(indx,3),sub_data(indx,2),sub_data(indx,8+4*(jj-1)),Y,X);
        uncert_R = griddata(sub_data(indx,3),sub_data(indx,2),sub_data(indx,6+4*(jj-1)),Y,X);

        res_I = (data_I - pred_I) ./ uncert_I;
        res_R = (data_R - pred_R) ./ uncert_R;

        %% Plot data
        ax1 = axes('Position',[0.0+(jj-1)*0.4 0.7 0.2 0.2]);
        h = imagesc(x,y,data_R); hold on
        set(h,'alphadata',~isnan(data_R));
        
        if jj == 1
            title(['\bf Obs X: ' num2str(freq(ii)) ' Real'])            
        else
            title(['\bf Obs Y: ' num2str(freq(ii)) ' Real'])
        end
        
%         caxis([-5 5])
        colorbar
        set(ax1,'YTickLabel',[])
        set(ax1,'YDir','normal')
        set(ax1,'XTickLabel',[])
        axis equal tight
        axis([xmin xmax ymin ymax])
        grid on
        
        ax1 = axes('Position',[0.2+(jj-1)*0.4 0.7 0.2 0.2]);
        h = imagesc(x,y,data_I); hold on
        set(h,'alphadata',~isnan(data_R));
        
        if jj == 1
            title(['\bf Obs X: ' num2str(freq(ii)) ' Imag'])
            
        else
            title(['\bf Obs Y: ' num2str(freq(ii)) ' Imag'])
        end
%         caxis([-5 5])
        colorbar
        set(ax1,'YTickLabel',[])
        set(ax1,'YDir','normal')
        set(ax1,'XTickLabel',[])
        axis equal tight
        axis([xmin xmax ymin ymax])
        grid on

        %% Plot predicted
        ax1 = axes('Position',[0.0+(jj-1)*0.4 0.4 0.2 0.2]);
        h = imagesc(x,y,pred_R); hold on
        set(h,'alphadata',~isnan(pred_R));
        
        if jj == 1
            title(['\bf Pred X: ' num2str(freq(ii)) ' Real'])            
        else
            title(['\bf Pred Y: ' num2str(freq(ii)) ' Real'])
        end
        
%         caxis([-5 5])
        colorbar
        set(ax1,'YTickLabel',[])
        set(ax1,'YDir','normal')
        set(ax1,'XTickLabel',[])
        axis equal tight
        axis([xmin xmax ymin ymax])
        grid on
        
        ax1 = axes('Position',[0.2+(jj-1)*0.4 0.4 0.2 0.2]);
        h = imagesc(x,y,pred_I); hold on
        set(h,'alphadata',~isnan(pred_R));
        
        if jj == 1
            title(['\bf Pred X: ' num2str(freq(ii)) ' Imag'])
            
        else
            title(['\bf Pred Y: ' num2str(freq(ii)) ' Imag'])
        end
%         caxis([-5 5])
        colorbar
        set(ax1,'YTickLabel',[])
        set(ax1,'YDir','normal')
        set(ax1,'XTickLabel',[])
        axis equal tight
        axis([xmin xmax ymin ymax])
        grid on
    end


%%
end

