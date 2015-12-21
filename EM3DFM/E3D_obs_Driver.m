clear all
close all

addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB

work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\DIGHEM\3D\DO18\Inv10_wsqrt_10pct_1DReffix';
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\DIGHEM\3D\DO18\Inv8_wsqrt_10pct_1ppmfloor';
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\DIGHEM\3D\DO27\Inv41_10pct_1ppmfloor\newM_ref';
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\DIGHEM\3D\DO18\Inv9_wsqrt_10pct_1DRef';
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\DIGHEM\3D\DO27\Inv37_12stn_shallow_deep\REF1D_ax5_pct_flr';
obsfile = '..\Data_TotalField_10pct_1flr.dat';
predfile = 'dpred_07.txt';
fsfile = '..\e3d_data_FREESPACE.txt';

% Load data
[trx,data] = load_E3D_obs([work_dir '\' obsfile]);

% Load pred
pred = load_E3D_pred([ work_dir '\' predfile]);

freespace = load_E3D_pred([ work_dir '\' fsfile]);


%% Assign line number to obs location
freq = unique(data(:,1));
ccode = ['r' 'b' 'k'];

index = data(:,1)==freq(1);
lineID = xy_2_lineID(data(index,2),data(index,3));
line_num = unique(lineID);
nlines = length(line_num);
    
for jj = 1 : nlines    
    
    set(figure, 'Position', [25 100 1800 900])
    
    for ii = 1 : length(freq)
    
    index = data(:,1)==freq(ii);
    
    subdata = data(index,:) ;
    subpred = pred(index,:);
    
    % Convert back real to ppm  
    subdata(:,25) = (subdata(:,25) - freespace(index,15) ) ./ freespace(index,15) * 1e+6;
    subpred(:,15) = (subpred(:,15) - freespace(index,15) ) ./ freespace(index,15) * 1e+6;
    
    % Convert back imag to ppm
    subdata(:,26) =subdata(:,26)./ freespace(index,15) * 1e+6; 
    subdata(:,27) =subdata(:,27)./ freespace(index,15) * 1e+6; 
    subdata(:,28) =subdata(:,28)./ freespace(index,15) * 1e+6; 
    subpred(:,16) =subpred(:,16)./ freespace(index,15) * 1e+6;
    
    tempIN = subdata(:,25)./subpred(:,15);
    tempQUAD = subdata(:,27)./subpred(:,16);
    
    lindex = lineID == jj;
    
    subplot(2,1,2);title('\bfQuadrature');
    errorbar(subdata(lindex,2),subdata(lindex,27) , subdata(lindex,28) ,'Color',ccode(ii),'LineStyle','-','Marker','+'); hold on
    plot(subdata(lindex,2),(subpred(lindex,16) ),'Color',ccode(ii),'LineStyle',':','Marker','o'); hold on
    
    subplot(2,1,1);title('\bfIn-Phase');
    errorbar(subdata(lindex,2),( subdata(lindex,25) ), subdata(lindex,26) ,'Color',ccode(ii),'LineStyle','-','Marker','+'); hold on
    plot(subdata(lindex,2),( subpred(lindex,15) ),'Color',ccode(ii),'LineStyle',':','Marker','o'); hold on
    
%     figure(100);
%     subplot(2,1,2);(subdata(lindex,2),tempQUAD); hold on
%     subplot(2,1,1);plot(subdata(lindex,2),tempIN); hold on
%     set(h1,'LineStyle','-','color','red')
%     set(h2,'LineStyle',':','color','red')
%     set(h3,'LineStyle','-','color','black')
%     set(h4,'LineStyle',':','color','black')
        
%         plot(,'o:','Color',ccode(ii,:)); hold on
%         plot(,'*:','Color',ccode(ii,:)); hold on
    end
    
end

legend('Obs 900','Pred 900','Obs 7200','Pred 7200','Obs 56k','Pred 56k','Location','NorthWest');
figure(1);legend('Obs 900','Pred 900','Obs 7200','Pred 7200','Obs 56k','Pred 56k','Location','NorthWest');
figure(2);legend('Obs 900','Pred 900','Obs 7200','Pred 7200','Obs 56k','Pred 56k','Location','NorthWest');