%% Script for hypercube 

clear all
close all

addpath ../func_lib
addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB

work_dir = 'C:\LC\Private\dominiquef\Projects\3796_AGIC_Research\HyperCube\Foootprint';
rule_file{1} = 'Zones_Target\Raw\Model\Zone5_Raw_Complx3_Quant5.csv';
rule_file{2} = 'Zones_Target\Logged\Model\Zone5_Logged_Complx3_Quant5.csv';
rule_file{3} = 'Zones_Target\Stdized\Model\Zone5_Std_Complx3_Quant5.csv';


% Load data
data = importdata([work_dir '\' 'Data\HypecubeDatasetU_nt\HypecubeDatasetU_tsf.csv']);

head = data.textdata(1,[2:17 19:end]);
d = data.data;

ndata = size(d,1);

% Target variable
target = d(:,17) == 1;
d = d(:,[1:16 18:end]);
    
%% Load rules
for jj = 1 : 3
    
    [~,rule_ID{jj},rule_lift{jj},varID{jj},Lbound{jj},Ubound{jj}] = get_rules_v2([work_dir '\' rule_file{jj}],head,7);

    % Compute purity, lift, coverage

    nrules = size(varID,2);

    for ii = 1 :  nrules

       nvar = size(varID{jj},2);
       vari = varID{jj}(1,ii);
       Lb = Lbound{jj}(ii);
       Ub = Ubound{jj}(ii);
       
       idx = ones(ndata,1);
       
       for  kk = 1 : nvar

          idx = idx .* (d(:,vari{1}(kk)) >= Lb{1}(kk) & d(:,vari{1}(kk)) <= Ub{1}(kk));

       end

          purity{jj}(ii) = sum(idx & target) / sum(idx);
          lift{jj}(ii) = purity(ii) / (sum(target) / size(d,1));
          coverage{jj}(ii) = sum(idx & target) / sum(target);

    end

end

