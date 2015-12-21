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

rule_file{1} = 'Zones_Target\Coverage\Zone5_RAW_cpmx4_pure0p1_cov0p5.csv';
% rule_file{2} = 'Zones_Target\Purity\Zone5_RAW_cpmx4_pure0p25_cov0p01.csv';
% rule_file{1} = 'Zones_Target\Default_quantile10\Zone5_RAW_cpmx4_DEFAULT.csv';

% rule_file{1} = 'Zones_Target\Coverage\Zone2_RAW_cpmx4_pure0p1_cov0p5.csv';
% rule_file{2} = 'Zones_Target\Purity\Zone2_RAW_cpmx4_pure0p25_cov0p01.csv';

% rule_file{1} = 'Uranium_Target\Purity\U_target_Zone1_cpmx4_pure0p25_cov0p01.csv';
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
% d_sub_norm = d(:,1:18);
% 
% % Normalize data before clustering (whitening)
% for ii = 1 : size(d_sub_norm,2)
%     
% d_sub_norm(:,ii) = d_sub_norm(:,ii)/std(d_sub_norm(:,ii));
% 
% end
% 
% % Compute kmeans
% [IDX,cmod] = kmeans(d_sub_norm,5);
% 
% % cpal = [255 100 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255;



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
        coverage(oo,ii)  = sum(idx & target);
        purity(oo,ii)    = coverage(oo,ii) / sum(idx);
        lift(oo,ii)      = purity(oo,ii) / base_pur;
        z_sc(oo,ii)      = sqrt(size(d,1))*(purity(oo,ii) - base_pur) / sqrt(purity(oo,ii) + (1- base_pur)); 

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

%% Create color scheme for lift
% lpal = [0 0 1;1 0 0];
clim = [min(z_sc(z_sc~=0)) mean([min(z_sc(z_sc~=0)) max(z_sc(:))]) max(z_sc(:))];
cmrk = interp1(clim,[0 0 255;0 255 0;255 0 0]/255,min(z_sc(z_sc~=0)):1e-2:max(z_sc(:)));
cdot = interp1(0:5',cpal,data.data(:,out_var),'linear');
%% Rank variables on number of times they are used in all the models

var_count = zeros(size(actv_var));

% For each model
for ii = 1 : size(R_grp,2)
    
    % For each group of rules
    for jj = 1 : size(R_grp{ii},2)
        
        % For each rule
        for kk = 1 : length(R_grp{ii}{jj})
        
        var_count(varID{ii}{R_grp{ii}{jj}(kk)}(1,:)) = ...
            var_count(varID{ii}{R_grp{ii}{jj}(kk)}(1,:)) + 1;
       
        end
        
    end
    
end

% Rank the variables in order of most to least used
[~,var_rank] = sort(var_count,'descend');

%% Rank the group of rules based on the dominant variables
% For each model
for ii = 1 : size(R_grp,2)
    
    R_grp_Rank{ii} = zeros(1,size(R_grp{ii},2));
    % For each group of rules
    for jj = 1 : size(R_grp{ii},2)
        
        %Add points depending on the ranking of variables
        for kk = 1 : 3
            
            R_grp_Rank{ii}(jj) = R_grp_Rank{ii}(jj) +...
                (find(varID{ii}{R_grp{ii}{jj}(1)}(1,kk) == var_rank));
            
        end
        
    end
    % Get the ranking
    [~,R_grp_Rank{ii}] = sort(R_grp_Rank{ii});
end


%% SORTING PHASE
% For each variable (in rank order) plot the distribution of rules
% for all the models 

fig_indx = [];
kk = 0;
   
% For each model
for ii = 1 : size(R_grp,2)

    % For each group of rules
    for jj = 1 : size(R_grp{ii},2)

        FLAG = 1;

        % Check if the combination of variable has been used by another
        % model. If yes than plot in same graph, otherwise create a new
        % window

        var_temp = varID{ii}{R_grp{ii}{jj}(1)}(1,:);

        if ~isempty(fig_indx)

            nfig = length(unique(fig_indx(:,4)));

            for kk = 1 : nfig

                fig_temp = fig_indx(fig_indx(:,4)==kk,:);

                % Test if variables are re-used
                % If 1 variable: create a new figure
                % If 2 variable: create a new axis
                % If 3 variable: use same handle               
                var_test1 = (sum(var_temp(1) == fig_temp(:,1:3),2)) ;
                var_test2 = (sum(var_temp(2) == fig_temp(:,1:3),2)) ;
                var_test3 = (sum(var_temp(3) == fig_temp(:,1:3),2));

                test2 = [var_test1 var_test2 var_test3];

                if  sum(sum(test2,1)==size(fig_temp,1)) >=2

                    test3 = find(sum(test2,2)==3);

                    if ~isempty(test3)

                        fig_id = fig_temp(test3(1),4);
                        axs_id = fig_temp(test3(1),5);

                        fig_indx = [fig_indx;var_temp fig_id axs_id];
%                         shuffle = [find(var_temp==fig_indx(kk,1))...
%                             find(var_temp==fig_indx(kk,2))...
%                             find(var_temp==fig_indx(kk,3))];

%                             shuffle = 1:3;
                        FLAG = 0;

                        break

                    end

                    fig_id = fig_temp(1,4);
                    axs_id = max(fig_temp(:,5)) + 1;

                    fig_indx = [fig_indx;var_temp fig_id axs_id];

                    FLAG = 0;

                    break

                end

            end

        end
        % Have not found any rules like that -- Create new figure
        if FLAG == 1

            if isempty(kk)

                kk = 0;
            end

            % 
            fig_id = kk+1;
            axs_id = 1;
            
            fig_indx = [fig_indx;var_temp fig_id axs_id];
            set(figure(kk+1),'Position',[50 50 1900 800]);

            %% Add colorbar
            axb = axes('Position',[0.10 0.225 .15 .15]);
            cbar = colorbar(axb,'NorthOutside');
            colormap(axb,cmrk);

            set(cbar,'Ticks',[0 0.5 1])
            set(cbar,'TickLabels',round(clim*10)/10)
            
            set(gca,'Visible','off');
            text(0.50,2.1,'$Zscore$', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
            
            %% ADD LEGEND
            axb = axes('Position',[0.30 0.25 .15 .15]);
            temp = [0.1 0.2 0.4 0.6 0.8 1];
            scatter(temp,ones(6,1),(temp*100+5).^1.6);
            text(0.5,1.5,'Coverage')
            for pp = 1 : 6
                
                text(temp(pp),1,[num2str(temp(pp)*100) '%'],'HorizontalAlignment','center')
                
            end
            set(gca,'Visible','off');
            
        end
        
    end
    
end

%% PLOTTING SCRIPT
count = 0;
for ii = 1 : size(R_grp,2)

    % For each group of rules
    for jj = 1 : size(R_grp{ii},2)
        
        count = count + 1;
        
        fig_id = fig_indx(count,4);
        axs_id = fig_indx(count,5);
        
        % Find the two variables repeated in the figure
        fig_temp = fig_indx(fig_indx(:,4) == fig_id,:);
        
        fig_test = [sum(sum(fig_temp(1,1) == fig_temp(:,1:3))) ...
                    sum(sum(fig_temp(1,2) == fig_temp(:,1:3))) ...
                    sum(sum(fig_temp(1,3) == fig_temp(:,1:3)))];
           
        
        [~,var] = sort(fig_test,'descend');
        
        test_indx = [find(fig_indx(count,1:3)==fig_temp(1,var(1)))...
                            find(fig_indx(count,1:3)==fig_temp(1,var(2)))];
                        
        temp = zeros(1,3); temp(test_indx) = 1;
         
        
        shuffle = zeros(1,3);
        shuffle(1:2) = test_indx;
        shuffle(3) = find(temp==0); 
        
        % Find if it is the first time the window is used
        fig_temp = find(fig_indx(count,4)==fig_indx(:,4));
        [axs_temp,axs_indx] = unique(fig_indx(fig_temp,5));
        
        temp = fig_temp(axs_indx(axs_temp == fig_indx(count,5)));
        
        if temp == count
            
            figure(fig_id);
            ax{fig_id,axs_id,1} = axes('Position',[-0.025 + 0.135*(axs_id-1) .725 .25 .25]);
            
            scatter3(ax{fig_id,axs_id,1},...
                log10(d(target,varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(1)))),...
                log10(d(target,varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(2)))),...
                log10(d(target,varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(3)))),...
                40,'r^','filled');

            hold on

            scatter3(ax{fig_id,axs_id,1},...
                log10(d(target==0,varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(1)))),...
                log10(d(target==0,varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(2)))),...
                log10(d(target==0,varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(3)))),...
                15,'k','filled');

            xlabel(head{varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(1))});
            ylabel(head{varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(2))});
            zlabel(head{varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(3))});
            set(get(gca,'YLabel'),'Rotation',315);
            set(get(gca,'XLabel'),'Rotation',45);
            
            axis([min(log10(d(:,varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(1))))) ...
                  max(log10(d(:,varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(1))))) ...
                  min(log10(d(:,varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(2))))) ...
                  max(log10(d(:,varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(2))))) ...
                  min(log10(d(:,varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(3))))) ...
                  max(log10(d(:,varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(3)))))]);
            
%             colormap(cdot)
            view([-45 45]);
            axis square tight
            title(['Zone ' num2str(targ_var)])
            ax{fig_id,axs_id,2} = axes('Position',[-0.025 + 0.135*(axs_id-1) .375 .25 .25]);
            set(ax{fig_id,axs_id,2},'XTickLabel',[])
            ylabel(head{varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(2))});
            hold(ax{fig_id,axs_id,2},'on')  
%             axis square tight
            colormap(ax{fig_id,axs_id,2},jet);
            
            
            ax{fig_id,axs_id,3} = axes('Position',[-0.025 + 0.135*(axs_id-1) .05 .25 .25]);
%         axis(ax{fig_id,axs_id,3},[min(log10(d(:,varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(1))))) ...
%             max(log10(d(:,varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(1))))) ...
%             min(log10(d(:,varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(3))))) ...
%             max(log10(d(:,varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(3)))))],'square');
            xlabel(head{varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(1))});
            ylabel(head{varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(3))});
            hold(ax{fig_id,axs_id,3},'on')  
%             axis square tight
        end
        
        
        for kk = 1 : length(R_grp{ii}{jj})

            %% PLOT CUBE IN 3D
            indx = R_grp{ii}{jj}(kk);

            if Lbound{ii}{indx}(shuffle(1))==0

                temp = d(:,varID{ii}{indx}(1,shuffle(1))) > 0;

                xbox = [log10(min(d(temp,varID{ii}{indx}(1,shuffle(1)))))...
                    log10(Ubound{ii}{indx}(shuffle(1)))];

            else

                xbox = [log10(Lbound{ii}{indx}(shuffle(1)))...
                    log10(Ubound{ii}{indx}(shuffle(1)))];

            end

            if Lbound{ii}{indx}(shuffle(2))==0

                temp = d(:,varID{ii}{indx}(1,shuffle(2))) > 0;

                ybox = [log10(min(d(temp,varID{ii}{indx}(1,shuffle(2)))))...
                    log10(Ubound{ii}{indx}(shuffle(2)))];

            else

                ybox = [log10(Lbound{ii}{indx}(shuffle(2)))...
                    log10(Ubound{ii}{indx}(shuffle(2)))];

            end

            if Lbound{ii}{indx}(shuffle(3))==0

                temp = d(:,varID{ii}{indx}(1,shuffle(3))) > 0;

                zbox = [log10(min(d(temp,varID{ii}{indx}(1,shuffle(3)))))...
                    log10(Ubound{ii}{indx}(shuffle(3)))];

            else

                zbox = [log10(Lbound{ii}{indx}(shuffle(3)))...
                    log10(Ubound{ii}{indx}(shuffle(3)))];

            end

    %             ccode = randn(1,3);
    %             ccode = abs(ccode) / norm(ccode);

            ccode = zeros(1,3);
            ccode(ii) = 1;

            figure(fig_id);

            patch([xbox(1) xbox(1) xbox(2) xbox(2)],[ybox(1) ybox(2) ybox(2) ybox(1)],[zbox(1) zbox(1) zbox(1) zbox(1)],ccode,'FaceAlpha',0.1, 'Parent', ax{fig_id,axs_id,1});
            patch([xbox(1) xbox(1) xbox(2) xbox(2)],[ybox(1) ybox(2) ybox(2) ybox(1)],[zbox(2) zbox(2) zbox(2) zbox(2)],ccode,'FaceAlpha',0.1, 'Parent', ax{fig_id,axs_id,1});

            patch([xbox(1) xbox(1) xbox(2) xbox(2)],[ybox(1) ybox(1) ybox(1) ybox(1)],[zbox(1) zbox(2) zbox(2) zbox(1)],ccode,'FaceAlpha',0.1, 'Parent', ax{fig_id,axs_id,1});
            patch([xbox(1) xbox(1) xbox(2) xbox(2)],[ybox(2) ybox(2) ybox(2) ybox(2)],[zbox(1) zbox(2) zbox(2) zbox(1)],ccode,'FaceAlpha',0.1, 'Parent', ax{fig_id,axs_id,1});

            patch([xbox(1) xbox(1) xbox(1) xbox(1)],[ybox(1) ybox(2) ybox(2) ybox(1)],[zbox(1) zbox(1) zbox(2) zbox(2)],ccode,'FaceAlpha',0.1, 'Parent', ax{fig_id,axs_id,1});
            patch([xbox(2) xbox(2) xbox(2) xbox(2)],[ybox(1) ybox(2) ybox(2) ybox(1)],[zbox(1) zbox(1) zbox(2) zbox(2)],ccode,'FaceAlpha',0.1, 'Parent', ax{fig_id,axs_id,1});


            %% Plot the purity and coverage in 2D
            % Plot the center of rule scaled by the coverage
            % and scale the color by purity
            cmarker = interp1(clim,[0 0 255;0 255 0;255 0 0]/255,z_sc(ii,indx));

            plot([xbox(1) xbox(2)],[mean(ybox) mean(ybox)],'k','LineWidth',1,'Parent',ax{fig_id,axs_id,2})
            hold(ax{fig_id,axs_id,2},'on')   
            plot([mean(xbox) mean(xbox)],[ybox(1) ybox(2)],'k','LineWidth',1,'Parent',ax{fig_id,axs_id,2})              
            scatter(mean(xbox),mean(ybox),(coverage(ii,indx)/sum(target)*100+5)^1.6,cmarker,'Parent',ax{fig_id,axs_id,2},'Marker',mtype(ii),'LineWidth',1);             
%                 caxis(ax{fig_id,axs_id,2},[min(lift(lift~=0)) max(lift(:))])
            axis(ax{fig_id,axs_id,2},[min(log10(d(:,varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(1))))) ...
              max(log10(d(:,varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(1))))) ...
              min(log10(d(:,varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(2))))) ...
              max(log10(d(:,varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(2)))))],'square');

            plot([xbox(1) xbox(2)],[mean(zbox) mean(zbox)],'k','LineWidth',1,'Parent',ax{fig_id,axs_id,3})
            hold(ax{fig_id,axs_id,3},'on')
            plot([mean(xbox) mean(xbox)],[zbox(1) zbox(2)],'k','LineWidth',1,'Parent',ax{fig_id,axs_id,3})
            scatter(mean(xbox),mean(zbox),(coverage(ii,indx)/sum(target)*100+5)^1.6,cmarker,'Parent',ax{fig_id,axs_id,3},'Marker',mtype(ii),'LineWidth',1); 
%                 caxis(ax{fig_id,axs_id,3},[min(lift(lift~=0)) max(lift(:))])
            axis(ax{fig_id,axs_id,3},[min(log10(d(:,varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(1))))) ...
              max(log10(d(:,varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(1))))) ...
              min(log10(d(:,varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(3))))) ...
              max(log10(d(:,varID{ii}{R_grp{ii}{jj}(1)}(1,shuffle(3)))))],'square');

          colormap(ax{fig_id,axs_id,3},cmrk);
%           colorbar('North')
          
        end

    end

end
