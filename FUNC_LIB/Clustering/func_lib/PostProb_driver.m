% Weight of Evidence test
% Script loading the Ni-RIM data and replicating GOCAD's results

clear all
close all

addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB


index = '1';

% work_dir = ['C:\LC\Private\dominiquef\Projects\3796_AGIC_Research\HyperCube\Processing\MtDore_files_for_Hypercube\Targeting_2015\WofE\Training_Set_' index];
work_dir = ['C:\LC\Private\dominiquef\Projects\3796_AGIC_Research\HyperCube\Processing\MtDore_files_for_Hypercube\Targeting_2015\HyperCube\HighCoverage\Subset' index];
% work_dir = 'C:\LC\Private\dominiquef\Projects\3796_AGIC_Research\HyperCube\Processing\MtDore_files_for_Hypercube\Targeting_2015\HyperCube\LowFiltering\12K_rules';
% work_dir = 'C:\LC\Private\dominiquef\Projects\3796_AGIC_Research\HyperCube\Processing\MtDore_files_for_Hypercube\Targeting_2015\HyperCube\HighCoverage\Combined';
% work_dir = 'C:\LC\Private\dominiquef\Projects\3796_AGIC_Research\HyperCube\Processing\MtDore_files_for_Hypercube\Targeting_2015\HyperCube\LowFiltering\4K_rules';

% rulefile = ['weights_training_set' index '.csv'];
rulefile = ['Mira_' index '_CT0,01 - PI0,0001 - CI0,001 - C1-5.csv'];
% rulefile = 'Mira 1 - All - 12000.csv';
% rulefile = 'Mira 1 - P0.5 - 4000.csv';

datafile = '\..\..\..\Dataset_all_training.dat';
data_training = ['\..\..\..\buffer_zones_SUB' index '.dat'];

criteria_list = [work_dir '\..\..\..\Datamtrix_header.txt'];

nullcell = load([work_dir '\..\..\..\activecell.dat']);

% Choose WofE or HC
var = 'HC';
out_col = 19;
drange = 18;
flag = 0;
%% Reformat data entry
% Load data
database = load([work_dir '\' datafile]);
OUTvar = load([work_dir '\' data_training]); 
OUTvar = OUTvar(nullcell==1);
OUTvar(OUTvar==-99999) = flag;

ndata = size(database,1);

outcome = OUTvar;
nts = sum(outcome~=0);
nnts = ndata - nts;

% Rank vector for efficiency curve
rankit = (ndata:-1:1 ) / ndata * 100;

data = [database(:,1:drange) outcome]; 

%% Extract rule parameters from input table

switch var
    
    case 'WofE'
    [criteria_name,ruleID,rule_lift,rules,Lbound,Ubound] = get_WofE_rules([work_dir '\' rulefile],criteria_list);
    otherwise
    [criteria_name,ruleID,rule_lift,rules,Lbound,Ubound] = get_rules_v2([work_dir '\' rulefile],criteria_list,11);

end
nrules = length(ruleID);

%% Loading MIRA WofE result
% load([work_dir '\MIRA_WofE_rank'])
% figure; 
% plot( rankit , MIRA_rank, 'r');
% hold on
% grid on

%% Low risk (HyperCube)
% index = [6 15 27];
% condition = [1 1 1];
% Ubound = [49.54982 110.333 97.18526];
% Lbound = [44.91676 18.82986 95.6636];

%% Compute weights



switch var
    
    case 'WofE'
        
        % Pre-allocate space for weights
        Wp  = zeros(size(rules,2),1);
        Wm  = zeros(size(rules,2),1);
        C   = zeros(size(rules,2),1);
        stdC= zeros(size(rules,2),1);
        
        for ii = 1 : size(rules,2)



            if rules{ii}(1,2) == 1 

                % First compute number of data in each domain 
                %  (1:yes | training)  (0:no | non-training)

                N11 = sum( data(outcome~=0,rules{ii}(1,1)) >= Lbound{ii}(1) );
                N10 = sum( data(outcome==0,rules{ii}(1,1)) >= Lbound{ii}(1) );
                N01 = sum( data(outcome~=0,rules{ii}(1,1)) <= Lbound{ii}(1) );
                N00 = sum( data(outcome==0,rules{ii}(1,1)) <= Lbound{ii}(1) );

             elseif rules{ii}(1,2) == -1 

                N11 = sum( data(outcome~=0,rules{ii}(1,1)) <= Ubound{ii}(1) );
                N10 = sum( data(outcome==0,rules{ii}(1,1)) <= Ubound{ii}(1) );
                N01 = sum( data(outcome~=0,rules{ii}(1,1)) >= Ubound{ii}(1) );
                N00 = sum( data(outcome==0,rules{ii}(1,1)) >= Ubound{ii}(1) );

            elseif  rules{ii}(1,2) == 0    

                logic = zeros(size(data,1),1);
                for jj = 1 : length(Lbound{ii});
                    
                    logic = logic | data(:,rules{ii}(1,1)) == Lbound{ii}(jj);
%                     N10 = N10 + data(tsites==0,rules{ii}(1,1)) == Lbound{ii}(jj);
%                     N01 = N01 + data(tsites,rules{ii}(1,1)) ~= Lbound{ii}(jj);
%                     N00 = N00 + data(tsites==0,rules{ii}(1,1)) ~= Lbound{ii}(jj); 
                    
                end
                
                N11 = sum( logic(outcome~=0) == 1 );
                N10 = sum( logic(outcome==0) == 1 );
                N01 = sum( logic(outcome~=0) == 0 );
                N00 = sum( logic(outcome==0) == 0 );

            
            end

            N11 = sum(N11); if N11==0; N11=1; end
            N10 = sum(N10); if N10==0; N10=1; end
            N01 = sum(N01); if N01==0; N01=1; end
            N00 = sum(N00); if N00==0; N00=1; end
            
                sigma = sqrt( 1/N11 + 1 / N10 + 1/N01 + 1/N00 );

                Wp(ii)  = log( ( N11 / nts) / ( N10 / nnts) );

                Wm(ii)  = log( ( N01 / nts ) / ( N00 / nnts) );

                C(ii)   = Wp(ii)- Wm(ii);

                stdC(ii)= C(ii) / sigma;

        end
        
    otherwise
        
        
        Wp  = zeros(size(rules,2),1);
        Wm  = zeros(size(rules,2),1);
        C   = zeros(size(rules,2),1);
        stdC= zeros(size(rules,2),1);
        
        % Print progress
        progress=-1;
        tic
        for jj = 1 : size(rules,2)
            
            N11 = ones(nts,1)==1;
            N10 = ones(nnts,1)==1;
            N01 = ones(nts,1)==1;
            N00 = ones(nnts,1)==1;
            for ii = 1 : size(rules{jj},2)
                
                
                % Check if continuous variable
                if rules{jj}(2,ii) == 1

                    N11 = N11 & ( ( data(outcome~=0,rules{jj}(1,ii)) <= Ubound{jj}(ii) ) &...
                            ( data(outcome~=0,rules{jj}(1,ii)) >= Lbound{jj}(ii) ) );

                    N10 = N10 & (  ( data(outcome==0,rules{jj}(1,ii)) <= Ubound{jj}(ii) ) &...
                        ( data(outcome==0,rules{jj}(1,ii)) >= Lbound{jj}(ii) ) );

                    N01 = N01 & (  ( data(outcome~=0,rules{jj}(1,ii)) >= Ubound{jj}(ii) ) |...
                        ( data(outcome~=0,rules{jj}(1,ii)) <= Lbound{jj}(ii) ) );

                    N00 = N00 & (  ( data(outcome==0,rules{jj}(1,ii)) >= Ubound{jj}(ii) ) |...
                        ( data(outcome==0,rules{jj}(1,ii)) <= Lbound{jj}(ii) ) );
                    
                else
                    
                    N11 = N11 & ( data(outcome~=0,rules{jj}(1,ii)) == Lbound{jj}(ii) ) ;

                    N10 = N10 & ( data(outcome==0,rules{jj}(1,ii)) == Lbound{jj}(ii) );

                    N01 = N01 & ( data(outcome~=0,rules{jj}(1,ii)) ~= Lbound{jj}(ii) ) ;

                    N00 = N00 & ( data(outcome==0,rules{jj}(1,ii)) ~= Lbound{jj}(ii) );
                    
                    
                end

            end
        
        
            N11 = sum(N11); if N11==0; N11=1; end
            N10 = sum(N10); if N10==0; N10=1; end
            N01 = sum(N01); if N01==0; N01=1; end
            N00 = sum(N00); if N00==0; N00=1; end

            sigma = sqrt( 1/N11 + 1 / N10 + 1/N01 + 1/N00 );

            Wp(jj)   = log( ( N11 / nts) / ( N10 / nnts) );

            Wm(jj)   = log( ( N01 / nts ) / ( N00 / nnts) );

            C(jj)    = Wp(jj) - Wm(jj);

            stdC(jj) = C(jj) / sigma;
        
            % Print progress
            d_iter = floor(jj/size(rules,2)*100);
            if  d_iter > progress

                fprintf('Computed %i pct of data in %8.5f sec\n',d_iter,toc)
                progress = d_iter;

            end
    
        end        
            
end

%% Compute Post-Prob on learning 
[post_prob_L,top_rule,~,~] = comp_postprob(data,rules,Ubound,Lbound,Wp,Wm,out_col,var,'all');

[score_L,rank_PProb_L] = sort(post_prob_L,1,'descend'); 
efficiency_L = cumsum(outcome(rank_PProb_L)~=0) / nts *100;

% load([work_dir '\MIRA_WofE_rank'])
figure; 
plot( rankit , efficiency_L,'LineWidth',2);
xlabel('\bfRanked post-probability (%)');
ylabel('\bfTraining nodes (%)');
ylim([0 105])
set(gca,'XDir','Reverse')
title('\bfEfficiency of classification - TRAINING DATABASE')
% legend('WofE','Predictive Analytics','Location','SouthEast')
grid on

save([work_dir '\post_prob_L' index],'post_prob_L')
%% Compute Post-Prob on non-trained
outcome = (database(:,end)~=0 & OUTvar==0).*1;
data(:,end) = outcome;

[post_prob_U,~,~,~] = comp_postprob(data,rules,Ubound,Lbound,Wp,Wm,out_col,var,'all');

[score_U,rank_PProb_U] = sort(post_prob_U,1,'descend'); 
efficiency_U = cumsum(outcome(rank_PProb_U)) / sum(outcome) *100;
% load([work_dir '\MIRA_WofE_rank'])
figure; 
plot( rankit , efficiency_U,'LineWidth',2);
xlabel('\bfRanked post-probability (%)');
ylabel('\bfTraining nodes (%)');
ylim([0 105])
set(gca,'XDir','Reverse')
title('\bfEfficiency of classification - UN-TRAINED DATABASE')
% legend('WofE','Predictive Analytics','Location','SouthEast')
grid on

save([work_dir '\post_prob_U' index],'post_prob_U')
%% Compute Post-Prob on entire database
% data = database;
[post_prob_A,top_rule,OUTid,OUThit] = comp_postprob(database,rules,Ubound,Lbound,Wp,Wm,out_col,var,'all');
outcome = database(:,end)~=0;

post_prob_A = post_prob_A - min(post_prob_A);
post_prob_A = post_prob_A / max(post_prob_A);

[score_A,rank_PProb_A] = sort(post_prob_A,1,'descend'); 
efficiency_A = cumsum(outcome(rank_PProb_A)) / sum(outcome) *100;
% load([work_dir '\MIRA_WofE_rank'])
figure; 
plot( rankit , efficiency_A,'LineWidth',2);
xlabel('\bfRanked post-probability (%)');
ylabel('\bfTraining nodes (%)');
ylim([0 105])
set(gca,'XDir','Reverse')
title('\bfEfficiency of classification - FULL DATABASE')
% legend('WofE','Predictive Analytics','Location','SouthEast')
grid on

save([work_dir '\post_prob_A' index],'post_prob_A')

%% Count the number of time a new FOG is detected

OUTunique = zeros( length(rank_PProb_A) , 1 );
temp = database(rank_PProb_A,out_col);

GROUP_count = zeros( length(OUTid) , 1 );


for jj = 1 : length(temp)
    
    OUTunique(jj) = sum(GROUP_count);
    
    if GROUP_count(temp(jj)==OUTid)==0
        
        GROUP_count(temp(jj)==OUTid) = 1;
        
        
    end
        
    
end

% Normalize the unique detection on all detection
OUTunique = OUTunique / length(OUTid) * 100;


figure; 
plot( rankit , OUTunique,'LineWidth',2);
xlabel('\bfRanked post-probability (%)');
ylabel('\bfUnique Sites (%)');
set(gca,'XDir','Reverse')
ylim([0 105])
title('\bfEfficiency of classification - Detection of unique Sites')
legend('WofE','Predictive Analytics','Location','SouthEast')
grid on
%% Map back to 3D, each FOG individucally
% t_stamp = unique(data(:,31));
% load([work_dir '\XYZ']);
% for ii = 1 : length(t_stamp)
%    
%    index =  data(:,31) == t_stamp(ii);
%    
%    avg_data =  [xyzFOG(index,4:6)  post_prob(index)];
%    
%    save([work_dir '\postProb_FOG_' num2str(t_stamp(ii)) 'BStrap.dat'],'-ascii','avg_data');
%     
%     
% end

%% Test final rules against k-mean clustering
% load([work_dir '\XYZ']);
ncluster = 4;
training = database(OUTvar~=0,:);
training(training==-99999) = NaN;
idx = kmeans(training(:,1:drange),ncluster);
nnodes = zeros(ncluster,1);

kcluster = zeros(ncluster,length(OUTid));
% Extract the unique FOG numbers in each cluster
for ii = 1:ncluster
    
    temp = unique(training(idx==ii,out_col));
    kcluster(ii,1:length(temp)) = temp;
    nnodes(ii) = sum(idx==ii);
    
end

% Build kclusters with coordinates
% temp = xyzFOG(OUTvar~=0,4:6);

% Only output non-nans
% indx = isnan(idx)==0;

% idx_xyz = [temp(indx,:) idx(indx)];

% save([work_dir '\K_cluster_xyz.dat'],'-ascii','idx_xyz');

%% Write out final rules

% Load criteria names
% criteria_name = load( criteria_list);

% switch var
%     
%     case 'HC'
%     fid = fopen([work_dir '\FINAL_Rule_Selection.dat'],'w');
% 
%     for ll = 1 : length(top_rule)
% 
%                 jj = top_rule(ll);        
% 
%                 fprintf(fid,'Rule# %i:\n',jj);
%                 fprintf(fid,'Criteria \t Lower_Bound \t Upper_Bound\n');
% 
%                 for ii = 1 : size(rules{jj},2)
% 
%                     c_ID = rules{jj}(1,ii);
%     %                 if c_ID == 31
%     %                     
%     %                     continue
%     %                     
%     %                 end
%                     fprintf(fid,'%i %40s \t %f \t %f\n',c_ID,criteria_name{c_ID},Lbound{jj}(ii),Ubound{jj}(ii));
% 
% 
%                 end
% 
%                 vec_var = find(OUThit(top_rule(ll),:));
% 
%                 for rr = 1 : length(vec_var)
% 
%                     fprintf(fid,'%i ',OUTid(vec_var(rr)));
% 
%                 end
% 
% 
%                 fprintf(fid,'\n\n');
% 
% 
% 
%     end
% 
%     fclose(fid);
% end
%% 2D Shell models





 
idx_xyz = zeros(ndata,1);
for jj = 1 : size(kcluster,1)
    
    for ii = 1 : sum(kcluster(jj,:)~=0)
        
        
        idx_xyz(OUTvar==kcluster(jj,ii)) = jj;
        
    end
    
end

model = ones(length(nullcell),1) * -99999;
model(nullcell==1) = idx_xyz;



save([work_dir '\kcluster3D.dat'],'-ascii','model');

% save([work_dir '\kcluster2D.dat'],'-ascii','idx_xyz');

% idx_xyz = ones(ndata,1)*10;
% for jj = 1 : size(OUThit,1)
%     
%     temp = find(OUThit(jj,:));
%     
%     for ii = 1 : length(temp)
%         
%         
%         idx_xyz(OUTvar==temp(ii)) = min([unique(idx_xyz(OUTvar==temp(ii))) jj]);
%         
%     end
%     
% end
% idx_xyz(idx_xyz==10)=0;
% 
% model = ones(length(nullcell),1) * -99999;
% model(nullcell==1) = idx_xyz;
% 
% save([work_dir '\HC_rule_hit3D.dat'],'-ascii','model');

% save([work_dir '\HC_rule_hit2D.dat'],'-ascii','idx_xyz');

switch var
    
    case 'HC'
        
        model = ones(length(nullcell),1) * -99999;
        model(nullcell==1) = (post_prob_A);
        save([work_dir '\post_prob_' var index '.dat'],'-ascii','model');

    otherwise

        model = ones(length(nullcell),1) * -99999;
        model(nullcell==1) = (post_prob_A);
        save([work_dir '\post_prob_' var index '.dat'],'-ascii','model');
        
end

fclose all;