%% Script for hypercube 
% Load sets of rules and plot all rules of complexity 3.
% The rules are group by dominant variable pairs and plotted in the same
% window. 
%
% Top: displays the points in the database and visual
% representation of the hypercubes
%
% Middle: displays in 2-D the rules centroid, coverage(marker size) and
% lift(marker color) for the dominant variables
%
% Bottom: displays in 2-D the rules centroid, coverage(marker size) and
% lift(marker color) for the first and third variables


clear all
close all

addpath ..\func_lib
addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB

work_dir = 'C:\LC\Private\dominiquef\Projects\3796_AGIC_Research\HyperCube\Projects\Foootprint';

% rule_file{1} = 'Zones_Target\Coverage\Zone5_RAW_cpmx4_pure0p1_cov0p5.csv';
% rule_file{1} = 'Zones_Target\Purity\Zone5_RAW_cpmx4_pure0p25_cov0p01.csv';
% rule_file{1} = 'Zones_Target\Default_quantile10\Zone5_RAW_cpmx4_DEFAULT.csv';

% rule_file{1} = 'Zones_Target\Coverage\Zone2_RAW_cpmx4_pure0p1_cov0p5.csv';
% rule_file{2} = 'Zones_Target\Purity\Zone2_RAW_cpmx4_pure0p25_cov0p01.csv';

rule_file{1} = 'Uranium_Target\Purity\U_target_Zone5_cpmx4_pure0p25_cov0p01.csv';
% rule_file{2} = 'Uranium_Target\Coverage\U_Target_Zone1_cpmx4_pure0p1_cov0p5.csv';

% rule_file{1} = 'Zones_Target\Default_quantile10\Zone5_RAW_cpmx4_DEFAULT.csv';
% rule_file{2} = 'Uranium_Target\Default_quantile10\U_target_Zone5_cpmx4_DEFAULT.csv';

data_file = 'Data\HypecubeDatasetU_nt\HypecubeDatasetU_tsf_EDT.csv';

% Directory seperator (Windows '\' || Unix '/')
dsep = '\';

% [out_var,targ_var,
% Specify column for output variable
out_var = 17;
targ_var = 5;

% Specify columns used for clustering (omit output variable)
actv_var    = [1:16 18 19];
data_head   = [1:16 18 19];

% Define color palette for clusters
cpal = [0 0 0;52 170 220;76 217 100;255 204 0;255 59 48;255 153 200]/255;

% Define marker type for plot
mtype = ['o' '^' 's'];

% View angle
cang = [45 15];
    
%% Load data
% Load data
data    = importdata([work_dir dsep data_file]);

head    = data.textdata(1,data_head);
d       = data.data;
ndata   = size(d,1);

% Target variable
target  = d(:,out_var) == targ_var;
d   = d(:,actv_var);


%% Load rules and sort by groups within models

% For each model (set of rules)
for oo = 1 : size(rule_file,2)
    
    % Parce the csv output from Hypercube in matrix form
    [~,rulesID,~,vartemp,Lb,Ub] = get_rules_v2([work_dir dsep rule_file{oo}],head,7);

    % Variable index and bounds for each rules
    varID{oo} = vartemp; % Row1: Variable number Row2: Continuous (1) | discrete(0)
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
       
        base_pur = (sum(target) / size(d,1));
        N = sum(idx);
        coverage(oo,ii)  = sum(idx & target);
        purity(oo,ii)    = coverage(oo,ii) / N;
        lift(oo,ii)      = purity(oo,ii) / base_pur;
        z_sc(oo,ii)      = sqrt(N)*(purity(oo,ii) - base_pur) / sqrt(purity(oo,ii) + (1- base_pur)); 
        
        for  jj = 1 : nvar

            temp = d(idx & target,varID{oo}{ii}(1,jj));

            cnt_R(jj,ii) = mean(log10(temp(temp>0)));

        end
         

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

            R_grp{oo}{count} =sub_rules(temp,1)';

        end

    end
    
    ngroup(oo) = count;
    
end

% Rank rules by z-score
[~,R_rank] = sort(z_sc(sub_rules(:,1)),'descend');

%% Create color scheme
% lpal = [0 0 1;1 0 0];
clim = [min(z_sc(z_sc~=0)) mean([min(z_sc(z_sc~=0)) max(z_sc(:))]) max(z_sc(:))];
cmrk = interp1(clim,[0 0 255;0 255 0;255 0 0]/255,min(z_sc(z_sc~=0)):1e-2:max(z_sc(:)));
cdot = interp1(0:5',cpal,data.data(:,out_var),'linear');



%% PLOT TOP RULES
% For the top 10 rules, plot the rule and scatter plot
% Compute distance between centroid of rule vs target
for ii = 1 : size(R_grp,2)
    
    for jj = 1 : 10
        
        %% Compute centroid of target and centroid of rule
        indx = sub_rules(R_rank(jj),1);
        
        A = [d(target,varID{ii}{indx}(1,1)) d(target,varID{ii}{indx}(1,2)) d(target,varID{ii}{indx}(1,3))]; 
        
        cnt_tag = [mean(log10(A(A(:,1)>0,1))) mean(log10(A(A(:,2)>0,2))) mean(log10(A(A(:,3)>0,3)))];
        
        
        set(figure,'Position',[50 50 1900 900]);
        axb = axes('Position',[0.675 0.575 .4 .4]);

        scatter3(log10(d(target,varID{ii}{indx}(1,1))),...
            log10(d(target,varID{ii}{indx}(1,2))),...
            log10(d(target,varID{ii}{indx}(1,3))),...
            40,'r^','filled');

        hold on

        scatter3(log10(d(target==0,varID{ii}{indx}(1,1))),...
            log10(d(target==0,varID{ii}{indx}(1,2))),...
            log10(d(target==0,varID{ii}{indx}(1,3))),...
            15,'k','filled');

        scatter3((cnt_tag(1)),...
            (cnt_tag(2)),...
            (cnt_tag(3)),...
            100,'cs','filled','MarkerEdgeColor','k');
        
        scatter3((cnt_R(1,indx)),...
            (cnt_R(2,indx)),...
            (cnt_R(3,indx)),...
            100,'g^','filled','MarkerEdgeColor','k');
        
        
        xlabel(head{varID{ii}{indx}(1,1)});
        ylabel(head{varID{ii}{indx}(1,2)});
        zlabel(head{varID{ii}{indx}(1,3)});
        set(get(gca,'YLabel'),'Rotation',315);
        set(get(gca,'XLabel'),'Rotation',45);

        axis([min(log10(d(:,varID{ii}{indx}(1,1)))) ...
              max(log10(d(:,varID{ii}{indx}(1,1)))) ...
              min(log10(d(:,varID{ii}{indx}(1,2)))) ...
              max(log10(d(:,varID{ii}{indx}(1,2)))) ...
              min(log10(d(:,varID{ii}{indx}(1,3)))) ...
              max(log10(d(:,varID{ii}{indx}(1,3))))]);

        %             colormap(cdot)
        view(cang);
        axis square tight
        title(['Zone ' num2str(targ_var)])
        
        %% PLOT CUBE IN 3D
        if Lbound{ii}{indx}(1)==0

            temp = d(:,varID{ii}{indx}(1,1)) > 0;

            xbox = [log10(min(d(temp,varID{ii}{indx}(1,1))))...
                log10(Ubound{ii}{indx}(1))];

        else

            xbox = [log10(Lbound{ii}{indx}(1))...
                log10(Ubound{ii}{indx}(1))];

        end

        if Lbound{ii}{indx}(2)==0

            temp = d(:,varID{ii}{indx}(1,2)) > 0;

            ybox = [log10(min(d(temp,varID{ii}{indx}(1,2))))...
                log10(Ubound{ii}{indx}(2))];

        else

            ybox = [log10(Lbound{ii}{indx}(2))...
                log10(Ubound{ii}{indx}(2))];

        end

        if Lbound{ii}{indx}(3)==0

            temp = d(:,varID{ii}{indx}(1,3)) > 0;

            zbox = [log10(min(d(temp,varID{ii}{indx}(1,3))))...
                log10(Ubound{ii}{indx}(3))];

        else

            zbox = [log10(Lbound{ii}{indx}(3))...
                log10(Ubound{ii}{indx}(3))];

        end

        %             ccode = randn(1,3);
        %             ccode = abs(ccode) / norm(ccode);

        ccode = zeros(1,3);
        ccode(ii) = 1;

        patch([xbox(1) xbox(1) xbox(2) xbox(2)],[ybox(1) ybox(2) ybox(2) ybox(1)],[zbox(1) zbox(1) zbox(1) zbox(1)],ccode,'FaceAlpha',0.1);
        patch([xbox(1) xbox(1) xbox(2) xbox(2)],[ybox(1) ybox(2) ybox(2) ybox(1)],[zbox(2) zbox(2) zbox(2) zbox(2)],ccode,'FaceAlpha',0.1);

        patch([xbox(1) xbox(1) xbox(2) xbox(2)],[ybox(1) ybox(1) ybox(1) ybox(1)],[zbox(1) zbox(2) zbox(2) zbox(1)],ccode,'FaceAlpha',0.1);
        patch([xbox(1) xbox(1) xbox(2) xbox(2)],[ybox(2) ybox(2) ybox(2) ybox(2)],[zbox(1) zbox(2) zbox(2) zbox(1)],ccode,'FaceAlpha',0.1);

        patch([xbox(1) xbox(1) xbox(1) xbox(1)],[ybox(1) ybox(2) ybox(2) ybox(1)],[zbox(1) zbox(1) zbox(2) zbox(2)],ccode,'FaceAlpha',0.1);
        patch([xbox(2) xbox(2) xbox(2) xbox(2)],[ybox(1) ybox(2) ybox(2) ybox(1)],[zbox(1) zbox(1) zbox(2) zbox(2)],ccode,'FaceAlpha',0.1);
    
        %% Plot all varables in 2D scatter
        
        % Pick the variable with largest simple lift as x-axis
        tlift = 0;
        varid = 0;
        idx  = ones(ndata,1);
        
        for kk = 1 : 3
            
            varin = varID{ii}{indx}(1,kk);
                        
            idin = d(:,varin) >= Lbound{oo}{indx}(kk) & d(:,varin) <= Ubound{oo}{indx}(kk);
            
            idx = idx.*idin;
            
            temp = sum(idin & target)/sum(idin) / base_pur;
            
            if temp > tlift
                
                varid = kk;
                tlift = temp;
                
            end
                
            
        end
        
        % Plot all 2-D scatters
        countx = 0;
        county = 0;
        for kk = 1 : size(d,2)
            
            % Skip if the same variable as the x-axis
            if kk == varID{ii}{indx}(1,varid)
                
                continue
                
            end
            
            % Increment the window plots
            if countx > 5
                
                countx = 0;
                county = county+1;
                
            end
            
            axb = axes('Position',[-.025+countx*0.15 0.05+county*0.28 .2 .2]);

            % Plot all points that are target variables
            scatter(log10(d(target,varID{ii}{indx}(1,varid))),...
                log10(d(target,kk)),...
                20,'r^','filled'); hold on
        
            % Plot all other points 
            scatter(log10(d(target==0,varID{ii}{indx}(1,varid))),...
            log10(d(target==0,kk)),...
            20,'ko','filled');
        
           % Plot points included in the current rule
            scatter(log10(d(idx==1,varID{ii}{indx}(1,varid))),...
            log10(d(idx==1,kk)),...
            10,'go','filled');
        
            axis([min(log10(d(:,varID{ii}{indx}(1,varid)))) ...
              max(log10(d(:,varID{ii}{indx}(1,varid)))) ...
              min(log10(d(:,kk))) ...
              max(log10(d(:,kk)))]);
          
            xlabel(head{varID{ii}{indx}(1,varid)});
            ylabel(head{kk});
            set(get(gca,'YLabel'),'Rotation',90);

        

            axis square tight
        
            countx = countx + 1;
            
            % Add Hypercube in 2-D if variables
            ll = kk == varID{ii}{indx}(1,:);
            if sum( ll ) > 0
                
                if Lbound{ii}{indx}(varid)==0

                temp = d(:,varID{ii}{indx}(1,varid)) > 0;

                xbox = [log10(min(d(temp,varID{ii}{indx}(1,varid))))...
                    log10(Ubound{ii}{indx}(varid))];

                else

                xbox = [log10(Lbound{ii}{indx}(varid))...
                    log10(Ubound{ii}{indx}(varid))];

                end

                if Lbound{ii}{indx}(ll)==0

                temp = d(:,varID{ii}{indx}(1,ll)) > 0;

                ybox = [log10(min(d(temp,varID{ii}{indx}(1,ll))))...
                    log10(Ubound{ii}{indx}(ll))];

                else

                ybox = [log10(Lbound{ii}{indx}(ll))...
                    log10(Ubound{ii}{indx}(ll))];

                end
                
                patch([xbox(1) xbox(1) xbox(2) xbox(2)],[ybox(1) ybox(2) ybox(2) ybox(1)],ccode,'FaceAlpha',0.1);
                
            end
            
            % Compute the centroid of target vs rule, and distance vector
            A = [d(target,varID{ii}{indx}(1,varid)) d(target,kk)]; 

            cnt_tag = [mean(log10(A(A(:,1)>0,1))) mean(log10(A(A(:,2)>0,2)))];
            scatter(cnt_tag(1),cnt_tag(2),100,'cs','filled','MarkerEdgeColor','k');
           
            A = [d(idx==1,varID{ii}{indx}(1,varid)) d(idx==1,kk)]; 

            cnt_idx = [mean(log10(A(A(:,1)>0,1))) mean(log10(A(A(:,2)>0,2)))];
            scatter(cnt_idx(1),cnt_idx(2),100,'g^','filled','MarkerEdgeColor','k');
            
            plot([cnt_idx(1) cnt_tag(1)],[cnt_idx(2) cnt_tag(2)],'g','LineWidth',3) 
        end
    
    end
    
end

