%% Script for hypercube 

clear all
close all

addpath ..\func_lib
addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB

work_dir = 'C:\LC\Private\dominiquef\Projects\3796_AGIC_Research\HyperCube\Projects\Foootprint';

rule_file{1} = 'Zones_Target\Coverage\Zone5_RAW_cpmx4_pure0p1_cov0p5.csv';
rule_file{2} = 'Zones_Target\Purity\Zone5_RAW_cpmx4_pure0p25_cov0p01.csv';
rule_file{3} = 'Zones_Target\Default_quantile10\Zone5_RAW_cpmx4_DEFAULT.csv';

% rule_file = 'Uranium_Target\Purity\U_target_Zone1_cpmx4_pure0p25_cov0p01.csv';
% rule_file = 'Uranium_Target\Coverage\U_Target_Zone1_cpmx4_pure0p1_cov0p5.csv';

% rule_file{1} = 'Zones_Target\Default_quantile10\Zone5_RAW_cpmx4_DEFAULT.csv';
% rule_file{2} = 'Uranium_Target\Default_quantile10\U_target_Zone5_cpmx4_DEFAULT.csv';

data_file = 'Data\HypecubeDatasetU_nt\HypecubeDatasetU_tsf.csv';

% Directory seperator (Windows '\' || Unix '/')
dsep = '\';

% Specify column for output variable
out_var = 17;
targ_var = 5;

% Specify columns used for clustering (omit output variable)
actv_var    = [1:16 18 19];
data_head   = [2:17 19:41];

% Define color palette for clusters
cpal = [0 0 0;52 170 220;76 217 100;255 204 0;255 59 48;255 153 200]/255;

%% Load data
% Load data
data    = importdata([work_dir dsep data_file]);

head    = data.textdata(1,data_head);
d       = data.data;
ndata   = size(d,1);

% Target variable
target  = d(:,out_var) == targ_var;
d   = d(:,actv_var);
%% Divide Uraninum into 5 zones based on percentile distribution
% pct = prctile(d_sub(:,1),[0 20 40 60 80 100]);
% 
% % Create targeting vector for each zone
% target = zeros(size(d_sub,1),1);
% for ii = 1 : 5
%     
%     target( d_sub(:,1) >= pct(ii) & d_sub(:,1) <= pct(ii+1)) = ii;
%     
%     fprintf('%i samples in zone %i\n',sum(target==ii),ii);
%     
% end

%% Compute k-means on data
d_sub_norm = d(:,1:18);

% Normalize data before clustering (whitening)
for ii = 1 : size(d_sub_norm,2)
    
d_sub_norm(:,ii) = d_sub_norm(:,ii)/std(d_sub_norm(:,ii));

end

% Compute kmeans
[IDX,cmod] = kmeans(d_sub_norm,5);

% cpal = [255 100 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255;
cdot = interp1(0:5',cpal,IDX,'linear');

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

%% Load rules

for oo = 1 : size(rule_file,2)
    
    [~,rulesID,~,vartemp,Lb,Ub] = get_rules_v2([work_dir dsep rule_file{oo}],head,7);

    varID{oo} = vartemp;
    Lbound{oo} = Lb;
    Ubound{oo} = Ub;
    
    nrules = length(rulesID);

    %% Compute purity, lift, coverage for all rules in file

    % Cycle through the rules and select data variable
    for ii  = 1 :  nrules

       nvar = size(varID{oo}{ii},2);
       idx  = ones(ndata,1);

       % Select data within all bounds of the rule
       for  jj = 1 : nvar

          idx = idx .* (d(:,varID{oo}{ii}(1,jj)) >= Lbound{oo}{ii}(jj) & d(:,varID{oo}{ii}(1,jj)) <= Ubound{oo}{ii}(jj));

       end

          purity(oo,ii)    = sum(idx & target) / sum(idx);
          lift(oo,ii)      = purity(oo,ii) / (sum(target) / size(d,1));
          coverage(oo,ii)  = sum(idx & target);

    end

    %% Filter all the rules for only complexity 3
    % Group the rules with same variables

    sub_rules = [];
    for ii = 1 : nrules

        if size( varID{oo}{ii}, 2 ) == 3

            sub_rules = [sub_rules;ii sort(varID{oo}{ii}(1,:))];

        end

    end

    indx = ones(size(sub_rules,1),1);
    count = 0;
    for ii = 1 : size(sub_rules,1)

        if indx(ii)==1

            count = count + 1 ;
            temp = sub_rules(:,2)==sub_rules(ii,2) &...
                sub_rules(:,3)==sub_rules(ii,3) &...
                sub_rules(:,4)==sub_rules(ii,4);

            indx(temp) = 0;

            grp_rules{oo}{count} =sub_rules(temp,1)';

        end

    end
    
    ngroup(oo) = count;
    
end

%% Look for identical rules in various models
% Variable keeper holds the group id number associated with each
% model that use the same variables.
% keeper(1 , 1) = [model1_groupx, model2_groupy]


keeper = [];
count1 = 1;


for ii = 1 : ngroup(1)
    
    var_set = varID{1}{ grp_rules{1}{ii}(1) }(1,:);
    count2 = 1;
    keeper(count1,1) = ii;
     
    for kk = 2 : size(grp_rules,2)
        
        for jj = 1 : ngroup(kk)

            if sum(var_set == varID{kk}{ grp_rules{kk}{jj}(1) }(1,:)) == 3

                 count2 = count2 + 1;
                 keeper(count1,kk) = jj;

            end

        end
    
    end
    
    if count2 == size(grp_rules,2)
        
        count1 = count1 + 1;
        
    end
    
end


%% Plot rules as 3D boxes overlayed on scatter only for rules with same
% variables found in both models
% Strutured cell-array
% varid{ii} = model
% varid{ii}{jj} = rule within model
% varid{ii}{jj}(kk) = variable within the rule, within model




for oo = 1 : size(keeper,1)
        
    % Plot cloud of point using model 1 (only rules with same variables in
    % both models are targeted, so it does not matter which one we are 
    % using here.)
    
    ax1 = figure(2);
    set(ax1,'Position',[50 50 800 800]);
    
    % Plot the location of target variable
    scatter3(log10(d(target,varID{1}{grp_rules{1}{keeper(oo,1)}(1)}(1,1))),...
        log10(d(target,varID{1}{grp_rules{1}{keeper(oo,1)}(1)}(1,2))),...
        log10(d(target,varID{1}{grp_rules{1}{keeper(oo,1)}(1)}(1,3))),40,'k^','filled');
    
    hold on
    
    scatter3(log10(d(:,varID{1}{grp_rules{1}{keeper(oo,1)}(1)}(1,1))),...
        log10(d(:,varID{1}{grp_rules{1}{keeper(oo,1)}(1)}(1,2))),...
        log10(d(:,varID{1}{grp_rules{1}{keeper(oo,1)}(1)}(1,3))),15,cdot,'filled');

    

    
    xlabel(head{varID{1}{grp_rules{1}{keeper(oo,1)}(1)}(1,1)});
    ylabel(head{varID{1}{grp_rules{1}{keeper(oo,1)}(1)}(1,2)});
    zlabel(head{varID{1}{grp_rules{1}{keeper(oo,1)}(1)}(1,3)});
    colormap(cpal)
    view([-30 30]);
    axis square
    
    for ii = 1 : size(grp_rules,2)   
           
    
        % Plot Hypercube box in 3D, compensate for the 0's in log space by
        % replacing with smallest value for the lower bound
        for jj = 1 : length(grp_rules{ii}{keeper(oo,ii)})

            kk = grp_rules{ii}{keeper(oo,ii)}(jj);

            if Lbound{ii}{kk}(1)==0

                temp = d(:,varID{ii}{kk}(1,1)) > 0;

                xbox = [log10(min(d(temp,varID{ii}{kk}(1,1)))) log10(Ubound{ii}{kk}(1))];

            else

                xbox = [log10(Lbound{ii}{kk}(1)) log10(Ubound{ii}{kk}(1))];

            end

            if Lbound{ii}{kk}(2)==0

                temp = d(:,varID{ii}{kk}(1,2)) > 0;

                ybox = [log10(min(d(temp,varID{ii}{kk}(1,2)))) log10(Ubound{ii}{kk}(2))];

            else

                ybox = [log10(Lbound{ii}{kk}(2)) log10(Ubound{ii}{kk}(2))];

            end

            if Lbound{ii}{kk}(3)==0

                temp = d(:,varID{ii}{kk}(1,3)) > 0;

                zbox = [log10(min(d(temp,varID{ii}{kk}(1,3)))) log10(Ubound{ii}{kk}(3))];

            else

                zbox = [log10(Lbound{ii}{kk}(3)) log10(Ubound{ii}{kk}(3))];

            end

%             ccode = randn(1,3);
%             ccode = abs(ccode) / norm(ccode);

            ccode = zeros(1,3);
            ccode(ii) = 1;
            
            patch([xbox(1) xbox(1) xbox(2) xbox(2)],[ybox(1) ybox(2) ybox(2) ybox(1)],[zbox(1) zbox(1) zbox(1) zbox(1)],ccode);
            patch([xbox(1) xbox(1) xbox(2) xbox(2)],[ybox(1) ybox(2) ybox(2) ybox(1)],[zbox(2) zbox(2) zbox(2) zbox(2)],ccode);

            patch([xbox(1) xbox(1) xbox(2) xbox(2)],[ybox(1) ybox(1) ybox(1) ybox(1)],[zbox(1) zbox(2) zbox(2) zbox(1)],ccode);
            patch([xbox(1) xbox(1) xbox(2) xbox(2)],[ybox(2) ybox(2) ybox(2) ybox(2)],[zbox(1) zbox(2) zbox(2) zbox(1)],ccode);

            patch([xbox(1) xbox(1) xbox(1) xbox(1)],[ybox(1) ybox(2) ybox(2) ybox(1)],[zbox(1) zbox(1) zbox(2) zbox(2)],ccode);
            patch([xbox(2) xbox(2) xbox(2) xbox(2)],[ybox(1) ybox(2) ybox(2) ybox(1)],[zbox(1) zbox(1) zbox(2) zbox(2)],ccode);

            text(xbox(1),ybox(1),num2str(kk),'FontSize',10);

            alpha(0.1);

        end

        

    end

    title('Click on figure if you want to save | Press ENTER to exit')

        gin = ginput(1);

        if ~isempty(gin)

            copyobj(gcf, 0)
            temp = ['--> Variables: ' head{varID{ii}{kk}(1,:)}];
            tlt = title(temp);

        end

        close(ax1);
        
end