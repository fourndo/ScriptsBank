
clear all
close all

addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB

% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\VTEM\3D\DO27\Inv1';
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\VTEM\3D\DO27\Inv1';
% filename = 'recv_h3dtd.txt';

obsfile = 'VTEM_AVG_h3d.obs';
prefile = 'dpred_03_02.txt';

dpre = load([work_dir '\' prefile]);

scale = 1;

[dobs,tx,index] = read_H3D_obs([work_dir '\' obsfile]);

% Find all unique stations, even if variable time channels
[stn_id,loc,~] = unique(index,'stable');

x = dobs(loc,1);
y = dobs(loc,2);

ndv = 99999;

lineID = xy_2_lineID(x, y);

line_num = unique(lineID);

nstn = length(stn_id);
nline = length(line_num);

% Find number of unique time channels
tc = unique(dobs(:,4));
ntc = length(tc);

obs = nan(nstn,ntc);
pre = nan(nstn,ntc);

% Re-structure data
for ii = 1 : nstn
    
    for jj = 1 : ntc
        
        grope = dobs(:,1) == x(ii) & dobs(:,2) == y(ii) & dobs(:,4)==tc(jj);
        
        if sum(grope)~=0
            
            obs(ii,jj) = dobs(grope,21);
            pre(ii,jj) = dpre(grope,13)*scale;
            
        end
        
    end
    
end

%% Read data and plot all time channels


subind = 1 ;
scale = 1;
clf;

pos_data = dobs; pos_data(dobs<0) = NaN;
neg_data = dobs; neg_data(dobs>0) = NaN;
set(figure(ii), 'Position', [0 0 2000 1000])

count= 1;
for ii = 1 : 2 : nline
    
    grope = lineID == line_num(ii);
    subplot(4,1,count)
    % Grab subset of data for the line
%     subdata = dobs(
    for jj = 1 : ntc;

    %     if subind~= ceil(ii/5)
    %         
    %         subind = ceil(ii/5);
    %         subplot(6,1,subind)
    %         
    %     end

           
           semilogy( x(grope) , abs(obs(grope,jj)),'b.:' ) ;hold on
           semilogy( x(grope) , abs(pre(grope,jj)),'r-o' ) ;hold on

%            semilogy(X_L150 , scale*(Hz_pred( ii : ntc : end)),'k'); hold on 
    %         semilogy( X_L150 , abs(Hz_pred( ii : nt : end)).*exp(-time(ii)) , 'r'); hold on

    end
    
    title('\bfObserved VTEM vs predicted using DIGHEM model (DO-27)');
    xlabel('\bfEasting (m)')
    ylabel('\bfdBz/dt');
    legend('OBS','PRED')
    count = count +1;

end

