function [head,ruleID,rule_lift,rules,Lbound,Ubound] = get_WofE_rules(rule_file,header_file)
    

% Conversion between GOCAD sorting and database
idx = [1 2 4 5 6 7 16 8 9 10 11 12 14 15 13 17 18 3];


fid = fopen(rule_file,'r');
line = fgets(fid);

count = 0 ;
ruleID = [];
rule_lift = [];
rules = [];
Lbound = [];
Ubound{1} = [];
head = [];

line = fgets(fid);
while line~=-1
    
    count = count + 1;
    
    data = regexp(line,',','split');
    
    ruleID(count) = count;

    % Check if variable is dicrete (0) or continuous (1)
    if strcmp(strtrim(data{5}), '+')==1

        Lbound{count} = str2double(data{4});
        Ubound{count} = str2double(data{3});
        rules{count}(1,1) = idx(count);
        rules{count}(1,2) = 1;

    elseif strcmp(strtrim(data{5}), '-')==1

        Lbound{count} = str2double(data{3});
        Ubound{count} = str2double(data{4});
        rules{count}(1,1) = idx(count);
        rules{count}(1,2) = -1;

    else
        
        temp = regexp(data{4},'\w*\d','match');
        
        Lbound{count} = zeros(size(temp,2),1);
        
        for jj = 1 : size(temp,2)
            
            Lbound{count}(jj) = str2double(temp{jj});
        
        end
        
        Ubound{count} = [];
        rules{count}(1,1) = idx(count);
        rules{count}(1,2) = 0;

    end
    
%     Ubound{count} = Ub;
%     Lbound{count} = Lb;

    line = fgets(fid);
end
fclose(fid);            


