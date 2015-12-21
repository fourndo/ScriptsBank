%% Script for hypercube 

clear all
close all

addpath ..\func_lib
addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB

work_dir = 'C:\LC\Private\dominiquef\Projects\3796_AGIC_Research\HyperCube\Projects\Foootprint';

% rule_file = 'Zones_Target\Coverage\Zone5_RAW_cpmx4_pure0p1_cov0p5.csv';
% rule_file = 'Zones_Target\Purity\Zone2_RAW_cpmx4_pure0p25_cov0p01.csv';

% rule_file = 'Uranium_Target\Purity\U_target_Zone1_cpmx4_pure0p25_cov0p01.csv';
% rule_file = 'Uranium_Target\Coverage\U_Target_Zone1_cpmx4_pure0p1_cov0p5.csv';

rule_file = 'Uranium_Target\Default_quantile10\U_target_Zone1_cpmx4_DEFAULT.csv';

data_file = 'Data\HypecubeDatasetU_nt\HypecubeDatasetU_tsf.csv';

% Directory seperator (Windows '\' || Unix '/')
dsep = '\';

% Specify column for output variable
out_var = 17;

% Specify columns used for clustering (omit output variable)
actv_var    = [1:16 18 19];
data_head   = [2:17 19:41];

% Define color palette for clusters
cpal = [0 0 0;52 170 220;76 217 100;255 204 0;255 59 48]/255;

%% Load data
% Load data
data    = importdata([work_dir dsep data_file]);

head    = data.textdata(1,data_head);
d       = data.data;
ndata   = size(d,1);

%% Load rules
[~,ruleID,~,varID,Lbound,Ubound] = get_rules_v2([work_dir dsep rule_file],head,7);

% Target variable
target  = d(:,out_var) == 1;
d_sub   = d(:,actv_var);

%% Compute purity, lift, coverage for all rules in file

nrules  = size(varID,2);

% Cycle through the rules and select data variable
for ii  = 1 :  nrules
    
   nvar = size(varID{ii},2);
   idx  = ones(ndata,1);
   
   % Select data within all bounds of the rule
   for  jj = 1 : nvar
       
      idx = idx .* (d_sub(:,varID{ii}(1,jj)) >= Lbound{ii}(jj) & d_sub(:,varID{ii}(1,jj)) <= Ubound{ii}(jj));
            
   end
   
      purity(ii)    = sum(idx & target) / sum(idx);
      lift(ii)      = purity(ii) / (sum(target) / size(d_sub,1));
      coverage(ii)  = sum(idx & target) / sum(target);
    
end



%% Divide Uraninum into 5 zones based on percentile distribution
pct = prctile(d_sub(:,1),[0 20 40 60 80 100]);

% Create targeting vector for each zone
target = zeros(size(d_sub,1),1);
for ii = 1 : 5
    
    target( d_sub(:,1) >= pct(ii) & d_sub(:,1) <= pct(ii+1)) = ii;
    
    fprintf('%i samples in zone %i\n',sum(target==ii),ii);
    
end

%% Compute k-means on data
d_sub_norm = d_sub(:,1:18);

% Normalize data before clustering (whitening)
for ii = 1 : size(d_sub_norm,2)
    
d_sub_norm(:,ii) = d_sub_norm(:,ii)/std(d_sub_norm(:,ii));

end

% Compute kmeans
[IDX,cmod] = kmeans(d_sub_norm,5);

% Plot all clusters in 2D with uranium on x-axis
figure;
ha = tight_subplot(5, 4, 0.05, 0.05, 0.05);

for ii = 1 : size(d_sub_norm,2)
    
    axes(ha(ii));
    scatter(log(d(:,1)),log(d(:,ii)),5,IDX); hold on
    title(head{ii})
    axis equal tight
    
end

colormap(cpal)


%% Filter all the rules for only complexity 3
% Group the rules with same variables

sub_rules = [];
for ii = 1 : nrules
    
    if size( varID{ii}, 2 ) == 3
        
        sub_rules = [sub_rules;ii sort(varID{ii}(1,:))];
        
    end
    
end

grp_rules = [];
indx = ones(size(sub_rules,1),1);
count = 0;
for ii = 1 : size(sub_rules,1)
    
    if indx(ii)==1
        
        count = count + 1 ;
        temp = sub_rules(:,2)==sub_rules(ii,2) &...
            sub_rules(:,3)==sub_rules(ii,3) &...
            sub_rules(:,4)==sub_rules(ii,4);
        
        indx(temp) = 0;
        
        grp_rules{count} =sub_rules(temp,1)';
        
    end
     
    
end

%% Plot rules as 3D boxes overlayed on scatter
for ii = 1 : size(grp_rules,2)
        
    % Plot cloud of point
    ax1 = figure(2);
    set(ax1,'Position',[50 50 800 800]);
    scatter3(log(d(:,varID{grp_rules{ii}(1)}(1,1))),log(d(:,varID{grp_rules{ii}(1)}(1,2))),log(d(:,varID{grp_rules{ii}(1)}(1,3))),15,IDX,'filled'); hold on
    xlabel(head{varID{grp_rules{ii}(1)}(1,1)});
    ylabel(head{varID{grp_rules{ii}(1)}(1,2)});
    zlabel(head{varID{grp_rules{ii}(1)}(1,3)});
    colormap(cpal)
    view([-30 30]);
    axis square
    
    % Plot Hypercube box in 3D, compensate for the 0's in log space by
    % replacing with smallest value for the lower bound
    for jj = 1 : length(grp_rules{ii})
        
        kk = grp_rules{ii}(jj);
        
        if Lbound{kk}(1)==0
            temp = d(:,varID{kk}(1,1)) > 0;

            xbox = [log(min(d(temp,varID{kk}(1,1)))) log(Ubound{kk}(1))];
        else

            xbox = [log(Lbound{kk}(1)) log(Ubound{kk}(1))];

        end

        if Lbound{kk}(2)==0
            temp = d(:,varID{kk}(1,2)) > 0;

            ybox = [log(min(d(temp,varID{kk}(1,2)))) log(Ubound{kk}(2))];
        else

            ybox = [log(Lbound{kk}(2)) log(Ubound{kk}(2))];

        end

        if Lbound{kk}(3)==0
            temp = d(:,varID{kk}(1,3)) > 0;

            zbox = [log(min(d(temp,varID{kk}(1,3)))) log(Ubound{kk}(3))];
        else

            zbox = [log(Lbound{kk}(3)) log(Ubound{kk}(3))];

        end

        ccode = randn(1,3);
        ccode = abs(ccode) / norm(ccode);
        
        patch([xbox(1) xbox(1) xbox(2) xbox(2)],[ybox(1) ybox(2) ybox(2) ybox(1)],[zbox(1) zbox(1) zbox(1) zbox(1)],ccode);
        patch([xbox(1) xbox(1) xbox(2) xbox(2)],[ybox(1) ybox(2) ybox(2) ybox(1)],[zbox(2) zbox(2) zbox(2) zbox(2)],ccode);

        patch([xbox(1) xbox(1) xbox(2) xbox(2)],[ybox(1) ybox(1) ybox(1) ybox(1)],[zbox(1) zbox(2) zbox(2) zbox(1)],ccode);
        patch([xbox(1) xbox(1) xbox(2) xbox(2)],[ybox(2) ybox(2) ybox(2) ybox(2)],[zbox(1) zbox(2) zbox(2) zbox(1)],ccode);

        patch([xbox(1) xbox(1) xbox(1) xbox(1)],[ybox(1) ybox(2) ybox(2) ybox(1)],[zbox(1) zbox(1) zbox(2) zbox(2)],ccode);
        patch([xbox(2) xbox(2) xbox(2) xbox(2)],[ybox(1) ybox(2) ybox(2) ybox(1)],[zbox(1) zbox(1) zbox(2) zbox(2)],ccode);
        
        text(xbox(1),ybox(1),num2str(kk),'FontSize',10);
        
        alpha(0.1);

    end

    title('Click on figure if you want to save | Press ENTER to exit')

    gin = ginput(1);

    if ~isempty(gin)

        copyobj(gcf, 0)
        temp = ['Nb. of rules: ' num2str(length(grp_rules{ii})) '--> Variables: ' head{varID{kk}(1,:)}];
        tlt = title(temp);

    end

    close(ax1);
        
end