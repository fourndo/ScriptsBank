% Print a suite of effeciency curves

clear all
close all

% work_dir = 'C:\LC\Private\dominiquef\Projects\3796_AGIC_Research\HyperCube\Processing\MtDore_files_for_Hypercube\Targeting_2015\WofE';
work_dir = 'C:\LC\Private\dominiquef\Projects\3796_AGIC_Research\HyperCube\Processing\MtDore_files_for_Hypercube\Targeting_2015\HyperCube\HighCoverage';

nullcell = load([work_dir '\..\..\activecell.dat']);

OUTid = load([work_dir '\..\..\buffer_zones_ALL.dat']);
OUTid = OUTid(nullcell==1);
sites = unique(OUTid(OUTid~=0));
nsites = length(sites);

% Convert output variable to binary
outvar_all = OUTid;
outvar_all(outvar_all~=0) = 1;

idx = 1:5;


% Load training cells
outvar = zeros(sum(nullcell==1),length(idx));

for ii = 1 : length(idx)
    
    temp =  load([work_dir '\..\..\buffer_zones_SUB' num2str(idx(ii)) '.dat']);
    temp(temp~=0) = 1;
    outvar(:,ii) = temp(nullcell==1);
    
end

ndata = size(outvar,1);
rankit = (ndata:-1:1 ) / ndata * 100;

% Compute average efficiency
avg_PP_A = zeros(ndata,1);

for ii = 1 : length(idx)
    
    unuse_train = outvar(:,ii) == 0 & outvar_all~=0;
    
    % Load efficiency curves
    load([work_dir '\post_prob_L' num2str(idx(ii))]);
    load([work_dir '\post_prob_U' num2str(idx(ii))]);
    load([work_dir '\post_prob_A' num2str(idx(ii))]);
    
    avg_PP_A = avg_PP_A + post_prob_A;
    
    [score_L,rank_PProb_L] = sort(post_prob_L,1,'descend'); 
    efficiency_L = cumsum(outvar(rank_PProb_L,ii)) / sum(outvar(:,ii)) *100;
    
    [score_A,rank_PProb_A] = sort(post_prob_A,1,'descend'); 
    efficiency_A = cumsum(outvar_all(rank_PProb_A)) / sum(outvar_all) *100;
    
    [score_U,rank_PProb_U] = sort(post_prob_U,1,'descend'); 
    efficiency_U = cumsum(unuse_train(rank_PProb_U)) / sum(unuse_train) *100;

    figure(1); 
    plot( rankit , efficiency_L,'LineWidth',2);
    hold on
    
    figure(2); 
    plot( rankit , efficiency_A,'--','LineWidth',0.5);
    hold on
    
    figure(3); 
    plot( rankit , efficiency_U,'LineWidth',2);
    hold on
    
    
%     if ii == 1
%         
%         % Make 20 bines
%         bins = -0.01:0.1:1;
%         aa = histc(post_prob_A,bins);
%         
%         bb = histc(post_prob_U(non_train),bins);
%         
%         figure;
%         semilogy(bins,aa); hold on
%         semilogy(bins,bb,'r');
%         
%     end
    % legend('WofE','Predictive Analytics','Location','SouthEast')

    hold on
    
end

%% Compute efficency for average PostProb
avg_PP_A = avg_PP_A / length(idx);

[score_A,rank_PProb_A] = sort(avg_PP_A,1,'descend'); 
efficiency_A = cumsum(outvar_all(rank_PProb_A)) / sum(outvar_all) *100;
    
model = ones(length(nullcell),1) * -99999;
model(nullcell==1) = (avg_PP_A);
save([work_dir '\post_prob_AVG.dat'],'-ascii','model');

figure(2)
plot( rankit , efficiency_A,'r','LineWidth',3);
hold on

% Add percentile score as line
% pct_PP = prctile(score_A,[25 50 75]);
% pct_25th = find(pct_PP(1) == score_A);
% pct_50th = find(pct_PP(2) == score_A);
% pct_75th = find(pct_PP(3) == score_A);

% y1=get(gca,'ylim');
% plot([rankit(pct_25th(1)) rankit(pct_25th(1))],y1);
% plot([rankit(pct_50th(1)) rankit(pct_50th(1))],y1);
% plot([rankit(pct_75th(1)) rankit(pct_75th(1))],y1);
%% FIND LOCATION ON EFFICIENCY CURVE WHERE SITES ARE FOUND

knowledge = nan(nsites,1);
ranked_OUTid = OUTid(rank_PProb_A);
count =0 ;
for ii = 1 : length(ranked_OUTid)
    
    % Check if the current training is new knowledge
    % If yes, then flag it
    if sum(knowledge == ranked_OUTid(ii)) == 0 && ranked_OUTid(ii)~=0
        
        count = count +1;
        knowledge(count) = ranked_OUTid(ii);
        
        figure(2)
        plot(rankit(ii),efficiency_A(ii),'bo','MarkerFaceColor','r');
%         text(rankit(ii)-log(knowledge(count))^2,efficiency_A(ii),num2str(knowledge(count)));
    end
        
end

%% ADD Title and labels
figure(1)
xlabel('\bfPost-Probability Percentile');
ylabel('\bfTraining nodes (%)');
ylim([0 105])
set(gca,'XDir','Reverse')
title('\bfEfficiency of classification - LEARNING SET')
axis equal
axis square
    grid on
    
figure(2)
hold on
xlabel('\bfPost-Probability Percentile');
ylabel('\bfTraining nodes (%)');
ylim([0 105])
set(gca,'XDir','Reverse')
title('\bfEfficiency of classification - FULL SET')
axis equal
axis square
grid on
    
figure(3)
xlabel('\bfPost-Probability Percentile');
ylabel('\bfTraining nodes (%)');
ylim([0 105])
set(gca,'XDir','Reverse')
title('\bfEfficiency of classification - UNTRAINED SET')
axis equal
axis square
    grid on
