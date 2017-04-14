% clear all
% close all

addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB

work_dir = 'C:\Users\DominiqueFournier\Downloads';
obsfile = 'TotalField.dat';
predfile = 'dpred_014.txt';
fsfile = 'predict.txt';

% outline = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\DO27_Outline.dat';
% hydro = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Lake_trace.dat';
% 
% Out_grav = load([outline]);
% Out_lake = load([hydro]);

% Load data
[trx,data] = load_E3D_obs([work_dir '\' obsfile]);

% Load pred
pred = importdata([ work_dir '\' predfile]);

% Load free space
[~,temp] = load_E3D_obs([work_dir '\' fsfile]);
freespace = temp(:,1:16);
% freespace = load_E3D_pred([ work_dir '\' fsfile]);

%% 2D plots
freq = unique(data(:,1));

set(figure(1), 'Position', [25 50 1800 900])
set(figure(2), 'Position', [25 50 1800 900])


for ii = 1 : length(freq)
    
    index = data(:,1) == freq(ii);
    subdata = data(index,:);
    
%     indx = pred(:,4) == freq(ii);
    subpred = pred(index,:);

    % Convert back real to ppm  
    subdata(:,25) = (subdata(:,25) - freespace(index,15) ) ./ freespace(index,15) * 1e+6;
    subpred(:,14) = (subpred(:,14) - freespace(index,15) ) ./ freespace(index,15) * 1e+6;
    
    % Convert back imag to ppm
    subdata(:,26) =subdata(:,26)./ freespace(index,15) * 1e+6; 
    subdata(:,27) =subdata(:,27)./ freespace(index,15) * 1e+6; 
    subdata(:,28) =subdata(:,28)./ freespace(index,15) * 1e+6; 
    subpred(:,15) =subpred(:,15)./ freespace(index,15) * 1e+6;
    
    if ii == 1
    %% Set coordinates for plot
    xmin = min(subdata(:,2));
    xmax = max(subdata(:,2));

    ymin = min(subdata(:,3));
    ymax = max(subdata(:,3));

    dx = 10;
    dy = 10;

    x = xmin:dx:xmax;
    y = ymin:dy:ymax;
    [Y,X] = ndgrid(y,x);
    
%     Y = flipud(Y);
    end
    
%     F_I = TriScatteredInterp(sub_pred(:,2),sub_pred(:,1),sub_pred(:,6),'natural');
%     F_R = TriScatteredInterp(sub_pred(:,2),sub_pred(:,1),sub_pred(:,5),'natural');
    
    pred_I = griddata(subpred(:,2),subpred(:,1),(abs(subpred(:,15))),Y,X,'natural');
    pred_R = griddata(subpred(:,2),subpred(:,1),(abs(subpred(:,14))),Y,X,'natural');
    
%     F_I = TriScatteredInterp(sub_data(:,2),sub_data(:,1),sub_data(:,6),'natural');
%     F_R = TriScatteredInterp(sub_data(:,2),sub_data(:,1),sub_data(:,5),'natural');
    
    data_I = griddata(subdata(:,3),subdata(:,2),(abs(subdata(:,27))),Y,X,'natural');
    data_R = griddata(subdata(:,3),subdata(:,2),(abs(subdata(:,25))),Y,X,'natural');
    
%     U_I = TriScatteredInterp(sub_data(:,2),sub_data(:,1),sub_data(:,8),'natural');
%     U_R = TriScatteredInterp(sub_data(:,2),sub_data(:,1),sub_data(:,7),'natural');
    
    uncert_I = griddata(subdata(:,3),subdata(:,2),subdata(:,28),Y,X,'natural');
    uncert_R = griddata(subdata(:,3),subdata(:,2),subdata(:,26),Y,X,'natural');
    
    res_I = ((data_I - pred_I)) ./ uncert_I;
    res_R = ((data_R - pred_R)) ./ uncert_R;
%%    
    figure(1)
    %ax1 = axes('Position',[-0.1+(ii-1)*0.325 0.55 0.35 0.35]);
    ax1 = subplot(2,6,(ii-1)*2 + 1);
    h = imagesc(x,y,data_R); hold on
    set(h,'alphadata',~isnan(data_R));
    scatter(data(index,2),data(index,3),'k.')
    title(['\bf Obs: ' num2str(freq(ii)) ' Real'])
    caxis([min(data_R(:)) max(data_R(:))])
    colormap(jet)
    colorbar('SouthOutside')
    set(ax1,'YTickLabel',[])
    set(ax1,'YDir','normal')
    set(ax1,'XTickLabel',[])
%     plot(Out_grav(:,1),Out_grav(:,2),'k--','LineWidth',2) 
%     plot(Out_lake(:,1),Out_lake(:,2),'r.','MarkerSize',2) 
    axis equal tight
    axis([xmin xmax ymin ymax])
    grid on
    
    %ax1 = axes('Position',[0.05+(ii-1)*0.325 0.55 0.35 0.35]);
    ax1 = subplot(2,6,(ii-1)*2 + 2);
    h=imagesc(x,y,data_I); hold on 
    set(h,'alphadata',~isnan(data_I));
    scatter(data(index,2),data(index,3),'k.')
    title(['\bf Obs: ' num2str(freq(ii)) ' Imag'])
    caxis([min(data_I(:)) max(data_I(:))])
    colorbar('SouthOutside')
    set(ax1,'YTickLabel',[])
    set(ax1,'YDir','normal')
    set(ax1,'XTickLabel',[])
%     plot(Out_grav(:,1),Out_grav(:,2),'k--','LineWidth',2) 
%     plot(Out_lake(:,1),Out_lake(:,2),'r.','MarkerSize',2) 
    axis equal tight
    axis([xmin xmax ymin ymax])
    grid on
    
    %ax1 = axes('Position',[-0.1+(ii-1)*0.325 0.05 0.35 0.35]);
    ax1 = subplot(2,6,(ii-1)*2 + 7);
    h = imagesc(x,y,pred_R); hold on
    set(h,'alphadata',~isnan(data_R));
    scatter(data(index,2),data(index,3),'k.')
    title(['\bf Pred: ' num2str(freq(ii)) ' Real'])
    caxis([min(data_R(:)) max(data_R(:))])
    colorbar('SouthOutside')
    set(ax1,'YTickLabel',[])
    set(ax1,'YDir','normal')
    set(ax1,'XTickLabel',[])
%     plot(Out_grav(:,1),Out_grav(:,2),'k--','LineWidth',2) 
%     plot(Out_lake(:,1),Out_lake(:,2),'r.','MarkerSize',2) 
    axis equal tight
    axis([xmin xmax ymin ymax])
    grid on
    
    %ax1 = axes('Position',[0.05+(ii-1)*0.325 0.05 0.35 0.35]);
    ax1 = subplot(2,6,(ii-1)*2 + 8);
    h=imagesc(x,y,pred_I); hold on
    set(h,'alphadata',~isnan(data_I));
    scatter(data(index,2),data(index,3),'k.')
    title(['\bf Pred: ' num2str(freq(ii)) ' Imag'])
    caxis([min(data_I(:)) max(data_I(:))])
    colorbar('SouthOutside')
    set(ax1,'YTickLabel',[])
    set(ax1,'YDir','normal')
    set(ax1,'XTickLabel',[])
%     plot(Out_grav(:,1),Out_grav(:,2),'k--','LineWidth',2) 
%     plot(Out_lake(:,1),Out_lake(:,2),'r.','MarkerSize',2) 
    axis equal tight
    axis([xmin xmax ymin ymax])
    grid on


    
    figure(2)
    colormap(jet)
    %ax1 = axes('Position',[-0.1+(ii-1)*0.325 0.55 0.35 0.35]);
    ax1 = subplot(2,6,(ii-1)*2 + 1);
    h = imagesc(x,y,data_R); hold on
    set(h,'alphadata',~isnan(data_R));
    scatter(data(index,2),data(index,3),'k.')
    title(['\bf Obs: ' num2str(freq(ii)) ' Real'])
    caxis([min(data_R(:)) max(data_R(:))])
    colorbar('SouthOutside')
    set(ax1,'YTickLabel',[])
    set(ax1,'YDir','normal')
    set(ax1,'XTickLabel',[])
%     plot(Out_grav(:,1),Out_grav(:,2),'k--','LineWidth',2) 
%     plot(Out_lake(:,1),Out_lake(:,2),'r.','MarkerSize',2) 
    axis equal tight
    axis([xmin xmax ymin ymax])
    grid on
    
    %ax1 = axes('Position',[0.05+(ii-1)*0.325 0.55 0.35 0.35]);
    ax1 = subplot(2,6,(ii-1)*2 + 2);
    h=imagesc(x,y,data_I); hold on
    set(h,'alphadata',~isnan(data_I));
    scatter(data(index,2),data(index,3),'k.')
    title(['\bf Obs: ' num2str(freq(ii)) ' Imag'])
    caxis([min(data_I(:)) max(data_I(:))])
    colorbar('SouthOutside')
    set(ax1,'YTickLabel',[])
    set(ax1,'YDir','normal')
    set(ax1,'XTickLabel',[])
%     plot(Out_grav(:,1),Out_grav(:,2),'k--','LineWidth',2) 
%     plot(Out_lake(:,1),Out_lake(:,2),'r.','MarkerSize',2) 
    axis equal tight
    axis([xmin xmax ymin ymax])
    grid on
    
    %ax1 = axes('Position',[-0.1+(ii-1)*0.325 0.05 0.35 0.35]);
    ax1 = subplot(2,6,(ii-1)*2 + 7);
    h = imagesc(x,y,res_R); hold on
    set(h,'alphadata',~isnan(data_R));
    scatter(data(index,2),data(index,3),'k.')
    title(['\bf Residual: ' num2str(freq(ii)) ' Real'])
    caxis([-5 5])
    colorbar('SouthOutside')
    set(ax1,'YTickLabel',[])
    set(ax1,'YDir','normal')
    set(ax1,'XTickLabel',[])
%     plot(Out_grav(:,1),Out_grav(:,2),'k--','LineWidth',2) 
%     plot(Out_lake(:,1),Out_lake(:,2),'r.','MarkerSize',2) 
    axis equal tight
    axis([xmin xmax ymin ymax])
    grid on
    
    %ax1 = axes('Position',[0.05+(ii-1)*0.325 0.05 0.35 0.35]);
    ax1 = subplot(2,6,(ii-1)*2 + 8);
    h=imagesc(x,y,res_I); hold on
    set(h,'alphadata',~isnan(data_I));
    scatter(data(index,2),data(index,3),'k.')
    title(['\bf Residual: ' num2str(freq(ii)) ' Imag'])
    caxis([-5 5])
    colorbar('SouthOutside')
    set(ax1,'YTickLabel',[])
    set(ax1,'YDir','normal')
    set(ax1,'XTickLabel',[])
%     plot(Out_grav(:,1),Out_grav(:,2),'k--','LineWidth',2) 
%     plot(Out_lake(:,1),Out_lake(:,2),'r.','MarkerSize',2) 
    axis equal tight
    axis([xmin xmax ymin ymax])
    grid on


%%
end

%% 1D plot

% %% Assign line number to obs location
% freq = unique(data(:,1));
% ccode = ['r' 'b' 'k'];
% 
% index = data(:,1)==freq(1);
% lineID = xy_2_lineID(data(index,2),data(index,3));
% line_num = unique(lineID);
% nlines = length(line_num);
%     
% for jj = 1 : nlines    
%     
%     set(figure, 'Position', [25 100 1800 900])
%     
%     for ii = 1 : length(freq)
%     
%     index = data(:,1)==freq(ii);
%     
%     subdata = data(index,:) ;
%     subpred = pred(index,:);
%     
%     % Convert back real to ppm  
%     subdata(:,25) = (subdata(:,25) - freespace(index,15) ) ./ freespace(index,15) * 1e+6;
%     subpred(:,14) = (subpred(:,14) - freespace(index,15) ) ./ freespace(index,15) * 1e+6;
%     
%     % Convert back imag to ppm
%     subdata(:,26) =subdata(:,26)./ freespace(index,15) * 1e+6; 
%     subdata(:,27) =subdata(:,27)./ freespace(index,15) * 1e+6; 
%     subdata(:,28) =subdata(:,28)./ freespace(index,15) * 1e+6; 
%     subpred(:,15) =subpred(:,15)./ freespace(index,15) * 1e+6;
%     
%     tempIN = subdata(:,25)./subpred(:,14);
%     tempQUAD = subdata(:,27)./subpred(:,15);
%     
%     lindex = lineID == jj;
%     
%     subplot(2,1,2);title('\bfQuadrature');
%     errorbar(subdata(lindex,2),subdata(lindex,27) , subdata(lindex,28) ,'Color',ccode(ii),'LineStyle','-','Marker','+'); hold on
%     plot(subdata(lindex,2),(subpred(lindex,15) ),'Color',ccode(ii),'LineStyle',':','Marker','o'); hold on
%     
%     subplot(2,1,1);title('\bfIn-Phase');
%     errorbar(subdata(lindex,2),( subdata(lindex,25) ), subdata(lindex,26) ,'Color',ccode(ii),'LineStyle','-','Marker','+'); hold on
%     plot(subdata(lindex,2),( subpred(lindex,14) ),'Color',ccode(ii),'LineStyle',':','Marker','o'); hold on
%     
% %     figure(100);
% %     subplot(2,1,2);(subdata(lindex,2),tempQUAD); hold on
% %     subplot(2,1,1);plot(subdata(lindex,2),tempIN); hold on
% %     set(h1,'LineStyle','-','color','red')
% %     set(h2,'LineStyle',':','color','red')
% %     set(h3,'LineStyle','-','color','black')
% %     set(h4,'LineStyle',':','color','black')
%         
% %         plot(,'o:','Color',ccode(ii,:)); hold on
% %         plot(,'*:','Color',ccode(ii,:)); hold on
%     end
%     
% end
% 
% legend('Obs 900','Pred 900','Obs 7200','Pred 7200','Obs 56k','Pred 56k','Location','NorthWest');
% figure(1);legend('Obs 900','Pred 900','Obs 7200','Pred 7200','Obs 56k','Pred 56k','Location','NorthWest');
% figure(2);legend('Obs 900','Pred 900','Obs 7200','Pred 7200','Obs 56k','Pred 56k','Location','NorthWest');