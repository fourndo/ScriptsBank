function [head,ruleID,rule_lift,aa,Lbound,Ubound] = get_rules_v2(rule_file,head,var_start)
    

% Load in Header file
% head = []; 
% count = 0;
% fid = fopen(header_file,'r');
% line = fgets(fid);
% 
% while line~=-1
%     
%     count = count+1;
%     head{count} = strtrim(line);
%     
%     line = fgets(fid);
%     
% end
% fclose(fid);

fid = fopen(rule_file,'r');
line = fgets(fid);

count = 0 ;
ruleID = [];
rule_lift = [];
rules = [];
Lbound = [];
Ubound{1} = [];
line = fgets(fid);
while line~=-1
    
    count = count + 1;
    
%     fprintf('%i\n',count)
    
    data = regexp(line,';','split');
    
    nvar = ( size(data,2)-(var_start)-1 ) / 3;
    
    
    ruleID(count) = count;
    
    Lb = [];
    Ub = [];
    indx = [];
    vartype = [];
    
    for ii = 1 : nvar

        if isempty(data{var_start + (ii-1)*3 +1}) == 0 && strcmp(data{var_start + (ii-1)*3 +1 },'Key_ID')==0
            
            % Look for variable number
            for jj = 1 : size(head,2)

                if strcmp(data{var_start + (ii-1)*3 + 1},head{jj})==1

                    indx = [indx jj];
                    break
                end

            end

            Lb = [Lb str2double(data{var_start +2 + (ii-1)*3})]; 

            % Check if variable is dicrete (0) or continuous (1)
            temp = str2double(data{var_start +3 + (ii-1)*3});

            if isnan(temp) == 1

                Ub = [Ub 0];
                vartype = [vartype 0];

            else

                Ub = [Ub temp];
                vartype = [vartype 1];
            end
            
        end

    end
    
    % Sort in ascending variable ID
    [~,oo] = sort(indx);
    
    rules{count} = [indx(oo);vartype(oo)];
    Ubound{count} = Ub(oo);
    Lbound{count} = Lb(oo);

    line = fgets(fid);
end
fclose(fid);            
aa = rules;


